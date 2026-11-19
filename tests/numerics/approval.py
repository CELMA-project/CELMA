"""Decide from the commit statuses whether a head was approved for a base.

Michael's approval is the label ``numerics-change-approved``. When he applies
it, the workflow records a commit status on the head commit that names the
merge base the head was compared with. A status belongs to one commit, so a
push gives a head without approval, and a status for another merge base (the
pull request was retargeted) does not count. The server sets the creation time
of a status and no commit date is involved.

Reads the JSON list of statuses of a commit (``GET /commits/SHA/statuses``,
newest first) on standard input; the exit code is 0 if approved.
"""

import json
import sys
from typing import Any, Dict, List, Optional

CONTEXT = "numerics-approval"


def description(merge_base: str) -> str:
    """Give the description the approval status carries.

    Parameters
    ----------
    merge_base : str
        Commit the head was compared with.

    Returns
    -------
    str
        The text of the status description.
    """
    return f"approved against {merge_base}"


def approved(statuses: List[Dict[str, Any]], merge_base: Optional[str]) -> bool:
    """Tell whether the newest approval status approves this head and base.

    Parameters
    ----------
    statuses : list
        Statuses of the head commit, newest first.
    merge_base : str or None
        Merge base of the pull request now, or ``None`` after the merge, when
        only the head itself has to be approved.

    Returns
    -------
    bool
        True if the newest status of the approval context is a success for
        this base.
    """
    for status in statuses:
        if status.get("context") != CONTEXT:
            continue
        if status.get("state") != "success":
            return False
        return merge_base is None or status.get("description") == description(
            merge_base
        )
    return False


def main(argv: List[str]) -> int:
    """Read statuses from standard input and report the verdict.

    Parameters
    ----------
    argv : list
        The merge base, or nothing.

    Returns
    -------
    int
        0 if approved, 1 if not.
    """
    ok = approved(json.load(sys.stdin), argv[0] if argv else None)
    print("approved" if ok else "no approval for this head and base")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
