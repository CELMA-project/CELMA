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
  voids the approval. Labels other than the approval label never do: a run for
  one of them (``approval=false`` with the action ``labeled`` or ``unlabeled``)
  changes nothing about the head, the base or the approval, so it is not read
  as an event, and it only re-reads the live state (``pull_request_read``).
  The base of the approval is compared by the fork point with the head, not
  by its sha: ``pull_request.base.sha`` may lag behind the branch, and a merge
  of a branch that moved on is the same change. The run's name (``run-name`` in the workflow, from the
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
import subprocess  # nosec
import sys
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

LABEL = "numerics-change-approved"
TITLE = re.compile(
    r"numerics-guard event=(?P<event>\S+) action=(?P<action>\S*) pr=(?P<pr>\d*) "
    r"head=(?P<head>[0-9a-f]*) base=(?P<base>[0-9a-f]*) approval=(?P<approval>true|false)"
)
COVERS = Callable[[str, str, str], bool]


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


def git_covers(approved: str, now: str, head: str) -> bool:
    """Tell, from the git history of the working directory, if a base still covers.

    The change a pull request makes is its head against the fork point of the
    head and the base branch. An approval of the change against the base
    ``approved`` covers the change against the base ``now`` when the fork point
    with ``now`` is the fork point with ``approved`` or a descendant of it: the
    branch moved on, and the change did not grow.

    Parameters
    ----------
    approved : str
        The base sha the approval names.
    now : str
        The base sha to judge against.
    head : str
        The head of the pull request.

    Returns
    -------
    bool
        True if the change against ``now`` is inside the approved change; False
        if git cannot tell.
    """

    def run(*args: str) -> subprocess.CompletedProcess:  # type: ignore[type-arg]
        return subprocess.run(  # nosec
            ["git", *args],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )

    first = run("merge-base", approved, head)
    second = run("merge-base", now, head)
    if first.returncode != 0 or second.returncode != 0:
        return False
    old, new = first.stdout.decode().strip(), second.stdout.decode().strip()
    if not old or not new:
        return False
    return old == new or run("merge-base", "--is-ancestor", old, new).returncode == 0


def same_base(approved: str, now: str, head: str, covers: Optional[COVERS]) -> bool:
    """Tell whether the base an approval names is the base of the change now.

    Parameters
    ----------
    approved : str
        The base sha in the approving run.
    now : str
        The base sha to judge against.
    head : str
        The head of the pull request.
    covers : callable or None
        ``covers(approved, now, head)``, as ``git_covers``; ``None`` compares
        the shas only.

    Returns
    -------
    bool
        True if the shas are equal, or ``covers`` says the change did not grow.
    """
    if approved and approved == now:
        return True
    if covers is None or not approved or not now or not head:
        return False
    return covers(approved, now, head)


def approving_run(
    run: Dict[str, Any],
    allowlist: List[str],
    pr_number: int,
    change: Tuple[str, str],
    covers: Optional[COVERS] = None,
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
    change : tuple
        The final head of the merged pull request and the tip of the base
        branch the merge started from.
    covers : callable, optional
        Tells if an approved base covers the base now (``same_base``).

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
    head, base = change
    return (
        run.get("status", "completed") == "completed"
        and run.get("conclusion") == "success"
        and match["event"] == "pull_request_target"
        and match["action"] == "labeled"
        and match["approval"] == "true"
        and match["pr"] == str(pr_number)
        and match["head"] == head
        and same_base(match["base"], base, head, covers)
        and permitted((run.get("actor") or {}).get("login", ""), allowlist)
        and permitted((run.get("triggering_actor") or {}).get("login", ""), allowlist)
    )


def other_label_event(match: "re.Match[str]") -> bool:
    """Tell whether a run name is a label other than the approval label.

    Parameters
    ----------
    match : re.Match
        A match of ``TITLE``.

    Returns
    -------
    bool
        True for the action ``labeled`` or ``unlabeled`` of a label that is not
        the approval label.
    """
    return match["action"] in ("labeled", "unlabeled") and match["approval"] == "false"


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
        request, leaving out the runs for other labels, or an empty dict.
    """
    named = []
    for run in runs:
        match = TITLE.fullmatch(run.get("display_title") or "")
        if match and run.get("event") == "pull_request_target":
            if match["pr"] == str(pr_number) and not other_label_event(match):
                named.append(run)
    return max(named, key=lambda r: r.get("run_number", 0), default={})


def push_approved(
    runs: List[Dict[str, Any]],
    pulls: List[Dict[str, Any]],
    allowlist: List[str],
    push: Push,
    covers: Optional[COVERS] = None,
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
    covers : callable, optional
        Tells if an approved base covers the base now (``same_base``).

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
            latest,
            allowlist,
            pull["number"],
            (pull["head"]["sha"], push.base),
            covers,
        ):
            return True
    return False


def is_other_label_event(payload: Dict[str, Any]) -> bool:
    """Tell whether an event is another label being added or removed.

    Parameters
    ----------
    payload : dict
        The event payload (``github.event``).

    Returns
    -------
    bool
        True for ``labeled`` or ``unlabeled`` of a label that is not the
        approval label.
    """
    name = (payload.get("label") or {}).get("name")
    return payload.get("action") in ("labeled", "unlabeled") and name != LABEL


def pull_request_read(
    runs: List[Dict[str, Any]],
    live: Dict[str, Any],
    allowlist: List[str],
    covers: Optional[COVERS] = None,
) -> bool:
    """Judge a pull request as it is now, for an event that changes nothing.

    Another label added or removed neither approves nor voids: the pull request
    passes if it is open, still has the approval label and its latest event
    (leaving out other labels) is the approval of its current head and base.

    Parameters
    ----------
    runs : list
        Runs of the workflow, any event, any state.
    live : dict
        The pull request now (``GET /repos/R/pulls/N``).
    allowlist : list
        Logins that may approve.
    covers : callable, optional
        Tells if an approved base covers the base now (``same_base``).

    Returns
    -------
    bool
        True if the approval still holds.
    """
    labels = [x.get("name") for x in live.get("labels") or []]
    head = (live.get("head") or {}).get("sha") or ""
    base = (live.get("base") or {}).get("sha") or ""
    number = live.get("number")
    return bool(
        live.get("state") == "open"
        and LABEL in labels
        and head
        and base
        and number is not None
        and approving_run(
            latest_run_of(runs, number), allowlist, number, (head, base), covers
        )
    )


def main(argv: List[str]) -> int:
    """Run a decision from the command line.

    Parameters
    ----------
    argv : list
        ``pr EVENT_JSON LIVE_PR_JSON`` (the triggering actor is in
        ``GITHUB_TRIGGERING_ACTOR``), ``other-label EVENT_JSON`` (exit 0 if the
        event is another label), ``read RUNS_JSON LIVE_PR_JSON`` or
        ``push RUNS_JSON PULLS_JSON BASE_SHA AFTER_SHA BRANCH``; the allowlist
        is in the environment variable ``ALLOWED_APPROVERS``, comma separated.

    Returns
    -------
    int
        0 if approved, 1 if not.
    """
    allowlist = [x for x in os.environ.get("ALLOWED_APPROVERS", "").split(",") if x]
    with open(argv[1], encoding="utf-8") as handle:
        first = json.load(handle)
    if argv[0] == "other-label":
        return 0 if is_other_label_event(first) else 1
    if argv[0] == "read":
        with open(argv[2], encoding="utf-8") as handle:
            approved = pull_request_read(
                first, json.load(handle), allowlist, git_covers
            )
    elif argv[0] == "pr":
        with open(argv[2], encoding="utf-8") as handle:
            live = json.load(handle)
        actor = os.environ.get("GITHUB_TRIGGERING_ACTOR", "")
        approved = pull_request_approved(first, actor, live, allowlist)
    else:
        with open(argv[2], encoding="utf-8") as handle:
            approved = push_approved(
                first,
                json.load(handle),
                allowlist,
                Push(argv[3], argv[4], argv[5]),
                git_covers,
            )
    print("approved" if approved else "no approval for this change")
    return 0 if approved else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
