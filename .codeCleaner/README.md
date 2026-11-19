# .codeCleaner

[Root](../README.md)

Old helper scripts for the code and the notebooks.
The continuous integration does not use them; formatting is now done by
`clang-format` and `ruff` through [pre-commit](../.pre-commit-config.yaml).

* [format-code.py](format-code.py) - Formats all `*.hxx` and `*.cxx` files
  with `clang-format` and `../.clang-format`. Run it from this directory.
* [convert-code.py](convert-code.py) - Runs `~/BOUT-dev/bin/bout_3to4.py`, the
  `BOUT++` 3 to 4 converter, on all `*.hxx` and `*.cxx` files. It does not
  upgrade to `BOUT++` 5.
* [cleanIpynb.py](cleanIpynb.py) - Strips the output from a notebook, for use
  as a git filter.
* [filterBranchCleanIpynb.sh](filterBranchCleanIpynb.sh) - Applies
  `cleanIpynb.py` to all notebooks in earlier commits with
  `git filter-branch`; the command is in the script.

## History

The repository was cleaned of accidentally committed binaries with
[bfg](https://rtyley.github.io/bfg-repo-cleaner/):

```bash
git clone --mirror https://github.com/CELMA-project/CELMA.git
wget -O bfg.jar http://repo1.maven.org/maven2/com/madgag/bfg/1.12.15/bfg-1.12.15.jar
java -jar bfg.jar --strip-blobs-bigger-than 20M CELMA.git
cd CELMA.git
git reflog expire --expire=now --all && git gc --prune=now --aggressive
git push
```

Pull requests are not cleaned, which leaves only the mirrored repository
large unless the pull request is removed, as described on
[this SO thread](http://stackoverflow.com/questions/34265266/remote-rejected-errors-after-mirroring-a-git-repository).
