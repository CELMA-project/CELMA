"""Decide whether an approver approved a numerics change, from facts GitHub sets.

Nothing is stored. GitHub itself records who applied a label and which head the
pull request had at that moment, and the pull request author cannot forge that:

* A pull request run of ``numerics-guard.yml`` may pass a numerics change only
  if its event is the label ``numerics-change-approved`` being applied
  (``labeled``) by a sender in the allowlist, the run is started (or re-run) by
  someone in the allowlist, and the pull request as it is NOW (the live API) is
  open, still has the label and still has the head and base of the event, since
  a re-run replays the old event. That run approves exactly that head and base. Every other event (opened, synchronize,
  reopened, edited, unlabeled, another label, a sender outside the allowlist)
  fails, so a push, a retarget or a reopen needs a new application of the label.
* After the merge, a push to ``main`` passes only if the LATEST
  ``pull_request_target`` run of the workflow for the merged pull request is a
  successful one for the label ``numerics-change-approved`` being applied, for
  the merged pull request's final head and for the base the push started from,
  triggered by someone in the allowlist (the original actor and the actor of any
  re-run), and the pull request still carried the label at the merge. A later
  event (an unlabel, a retarget, a push) makes a later run the latest, so it
  voids the approval. The run's name (``run-name`` in the workflow, from the
  base branch copy) carries event, action, pull request, head, base and label.

A copy of the workflow on a pull request branch cannot satisfy this. Its runs
have the event ``pull_request`` or ``push`` (or ``workflow_dispatch``), never
``pull_request_target``, which always runs the workflow of the base branch, and
the agent is their actor. An approval by another user, a status or a comment
that the agent posts is not read at all.
"""

import json
import os
import re
import sys
from dataclasses import dataclass
from typing import Any, Dict, List

LABEL = "numerics-change-approved"
TITLE = re.compile(
    r"numerics-guard event=(?P<event>\S+) action=(?P<action>\S*) pr=(?P<pr>\d*) "
    r"head=(?P<head>[0-9a-f]*) base=(?P<base>[0-9a-f]*) approval=(?P<approval>true|false)"
)


@dataclass(frozen=True)
class Push:
    """The push to a branch whose merge is being judged.

    Attributes
    ----------
    base : str
        The commit the push started from, the tip of the branch the merge was
        made on.
    after : str
        The commit the push ended on.
    branch : str
        The branch pushed to.
    """

    base: str
    after: str
    branch: str


def permitted(login: str, allowlist: List[str]) -> bool:
    """Tell whether a login may approve.

    Parameters
    ----------
    login : str
        GitHub login.
    allowlist : list
        Logins that may approve.

    Returns
    -------
    bool
        True if the login is in the list, ignoring case.
    """
    return bool(login) and login.lower() in [x.lower() for x in allowlist]


def pull_request_approved(
    payload: Dict[str, Any],
    triggering_actor: str,
    live: Dict[str, Any],
    allowlist: List[str],
) -> bool:
    """Decide a pull request run from its event, its actor and the live pull request.

    A re-run replays the original event payload, so the payload alone proves
    nothing about now: the run must also be started by a permitted user, and the
    pull request as it is now must still be open, carry the label and have the
    head and base the event approved.

    Parameters
    ----------
    payload : dict
        The event payload (``github.event``) of a ``pull_request_target`` run.
    triggering_actor : str
        Who started this run or re-run (``GITHUB_TRIGGERING_ACTOR``).
    live : dict
        The pull request now (``GET /repos/R/pulls/N``).
    allowlist : list
        Logins that may approve.

    Returns
    -------
    bool
        True only for the label being applied by a permitted sender, in a run
        started by a permitted user, while the pull request is still open with
        the label and with the head and base of the event.
    """
    event_pr = payload.get("pull_request") or {}
    live_labels = [x.get("name") for x in live.get("labels") or []]
    live_head = (live.get("head") or {}).get("sha")
    live_base = live.get("base") or {}
    event_base = event_pr.get("base") or {}
    return (
        payload.get("action") == "labeled"
        and (payload.get("label") or {}).get("name") == LABEL
        and permitted((payload.get("sender") or {}).get("login", ""), allowlist)
        and permitted(triggering_actor, allowlist)
        and live.get("state") == "open"
        and LABEL in live_labels
        and bool(live_head)
        and live_head == (event_pr.get("head") or {}).get("sha")
        and bool(live_base.get("sha"))
        and live_base.get("sha") == event_base.get("sha")
        and live_base.get("ref") == event_base.get("ref")
    )


def approving_run(
    run: Dict[str, Any], allowlist: List[str], pr_number: int, head: str, base: str
) -> bool:
    """Tell whether a workflow run is the approval of a pull request head.

    Parameters
    ----------
    run : dict
        An entry of the workflow-runs API.
    allowlist : list
        Logins that may approve.
    pr_number : int
        The merged pull request.
    head : str
        The final head of the merged pull request.
    base : str
        The tip of the base branch the merge started from.

    Returns
    -------
    bool
        True for a completed, successful ``pull_request_target`` run for this
        head and base, whose name says the approval label was applied, started
        and re-run only by permitted users.
    """
    match = TITLE.fullmatch(run.get("display_title") or "")
    if not match or run.get("event") != "pull_request_target":
        return False
    return (
        run.get("status", "completed") == "completed"
        and run.get("conclusion") == "success"
        and match["event"] == "pull_request_target"
        and match["action"] == "labeled"
        and match["approval"] == "true"
        and match["pr"] == str(pr_number)
        and match["head"] == head
        and match["base"] == base
        and permitted((run.get("actor") or {}).get("login", ""), allowlist)
        and permitted((run.get("triggering_actor") or {}).get("login", ""), allowlist)
    )


def latest_run_of(runs: List[Dict[str, Any]], pr_number: int) -> Dict[str, Any]:
    """Find the latest ``pull_request_target`` run for a pull request.

    Parameters
    ----------
    runs : list
        Runs of the workflow, any event, any state.
    pr_number : int
        The pull request.

    Returns
    -------
    dict
        The run with the highest run number among those named for the pull
        request, or an empty dict.
    """
    named = []
    for run in runs:
        match = TITLE.fullmatch(run.get("display_title") or "")
        if match and run.get("event") == "pull_request_target":
            if match["pr"] == str(pr_number):
                named.append(run)
    return max(named, key=lambda r: r.get("run_number", 0), default={})


def push_approved(
    runs: List[Dict[str, Any]],
    pulls: List[Dict[str, Any]],
    allowlist: List[str],
    push: Push,
) -> bool:
    """Decide a push to main from the merged pull request and the runs.

    Parameters
    ----------
    runs : list
        Runs of the workflow, any event, any state.
    pulls : list
        The merged pull requests associated with the pushed commit.
    allowlist : list
        Logins that may approve.
    push : Push
        The push: where it started, where it ended, which branch.

    Returns
    -------
    bool
        True if some pull request merged into this branch as exactly the pushed
        commit (``merge_commit_sha``, which GitHub sets to the merge, squash or
        last rebased commit) carried the label at the merge, and its latest run
        is the approval of its final head and base.
    """
    for pull in pulls:
        if (pull.get("base") or {}).get("ref") != push.branch:
            continue
        if not push.after or pull.get("merge_commit_sha") != push.after:
            continue
        labels = [x.get("name") for x in pull.get("labels") or []]
        latest = latest_run_of(runs, pull["number"])
        if LABEL in labels and approving_run(
            latest, allowlist, pull["number"], pull["head"]["sha"], push.base
        ):
            return True
    return False


def main(argv: List[str]) -> int:
    """Run a decision from the command line.

    Parameters
    ----------
    argv : list
        ``pr EVENT_JSON LIVE_PR_JSON`` (the triggering actor is in
        ``GITHUB_TRIGGERING_ACTOR``) or ``push RUNS_JSON PULLS_JSON BASE_SHA AFTER_SHA BRANCH``; the allowlist is
        in the environment variable ``ALLOWED_APPROVERS``, comma separated.

    Returns
    -------
    int
        0 if approved, 1 if not.
    """
    allowlist = [x for x in os.environ.get("ALLOWED_APPROVERS", "").split(",") if x]
    with open(argv[1], encoding="utf-8") as handle:
        first = json.load(handle)
    if argv[0] == "pr":
        with open(argv[2], encoding="utf-8") as handle:
            live = json.load(handle)
        actor = os.environ.get("GITHUB_TRIGGERING_ACTOR", "")
        approved = pull_request_approved(first, actor, live, allowlist)
    else:
        with open(argv[2], encoding="utf-8") as handle:
            approved = push_approved(
                first, json.load(handle), allowlist, Push(argv[3], argv[4], argv[5])
            )
    print("approved" if approved else "no approval for this change")
    return 0 if approved else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
