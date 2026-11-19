"""Decide whether Michael approved a numerics change, from facts GitHub sets.

Nothing is stored. GitHub itself records who applied a label and which head the
pull request had at that moment, and the pull request author cannot forge that:

* A pull request run of ``numerics-guard.yml`` may pass a numerics change only
  if its event is the label ``numerics-change-approved`` being applied
  (``labeled``) by a sender in the allowlist. That one run approves exactly the
  head and base its event carries. Every other event (opened, synchronize,
  reopened, edited, unlabeled, another label, a sender outside the allowlist)
  fails, so a push, a retarget or a reopen needs a new application of the label.
* After the merge, a push to ``main`` passes only if the workflow-runs API holds
  a successful ``pull_request_target`` run of the workflow for the merged pull
  request's final head, triggered by someone in the allowlist (the original
  actor and the actor of any re-run) for the label ``numerics-change-approved``.
  The run's name (``run-name`` in the workflow, from the base branch copy)
  carries event, action, pull request, head, base and label.

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
from typing import Any, Dict, List

LABEL = "numerics-change-approved"
TITLE = re.compile(
    r"numerics-guard event=(?P<event>\S+) action=(?P<action>\S*) pr=(?P<pr>\d*) "
    r"head=(?P<head>[0-9a-f]*) base=(?P<base>[0-9a-f]*) label=(?P<label>.*)"
)


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


def pull_request_approved(payload: Dict[str, Any], allowlist: List[str]) -> bool:
    """Decide a pull request run from its event payload.

    Parameters
    ----------
    payload : dict
        The event payload (``github.event``) of a ``pull_request_target`` run.
    allowlist : list
        Logins that may approve.

    Returns
    -------
    bool
        True only for the label being applied by a permitted sender. The
        approval is for the head and base of this very payload.
    """
    return (
        payload.get("action") == "labeled"
        and (payload.get("label") or {}).get("name") == LABEL
        and permitted((payload.get("sender") or {}).get("login", ""), allowlist)
    )


def approving_run(
    run: Dict[str, Any], allowlist: List[str], pr_number: int, head: str
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

    Returns
    -------
    bool
        True for a successful ``pull_request_target`` run for this head, whose
        name says the approval label was applied, started and re-run only by
        permitted users.
    """
    match = TITLE.fullmatch(run.get("display_title") or "")
    if not match or run.get("event") != "pull_request_target":
        return False
    return (
        run.get("conclusion") == "success"
        and match["event"] == "pull_request_target"
        and match["action"] == "labeled"
        and match["label"] == LABEL
        and match["pr"] == str(pr_number)
        and match["head"] == head
        and permitted((run.get("actor") or {}).get("login", ""), allowlist)
        and permitted((run.get("triggering_actor") or {}).get("login", ""), allowlist)
    )


def push_approved(
    runs: List[Dict[str, Any]], pulls: List[Dict[str, Any]], allowlist: List[str]
) -> bool:
    """Decide a push to main from the merged pull request and the runs.

    Parameters
    ----------
    runs : list
        Successful ``pull_request_target`` runs of the workflow.
    pulls : list
        The merged pull requests associated with the pushed commit.
    allowlist : list
        Logins that may approve.

    Returns
    -------
    bool
        True if some merged pull request's final head has an approving run.
    """
    return any(
        approving_run(run, allowlist, pull["number"], pull["head"]["sha"])
        for pull in pulls
        for run in runs
    )


def main(argv: List[str]) -> int:
    """Run a decision from the command line.

    Parameters
    ----------
    argv : list
        ``pr EVENT_JSON`` or ``push RUNS_JSON PULLS_JSON``; the allowlist is
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
        approved = pull_request_approved(first, allowlist)
    else:
        with open(argv[2], encoding="utf-8") as handle:
            approved = push_approved(first, json.load(handle), allowlist)
    print("approved" if approved else "no approval for this change")
    return 0 if approved else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
