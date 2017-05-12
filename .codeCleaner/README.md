# .codeCleaner

Contains scripts for cleaning up the code

* [format-code.py](format-code.py) - Formats all `*.hxx` and `*.cxx` files using
  `../.clang-format`
* [convert-code.py](convert-code.py) - Converts the code to the highest
  `BOUT++` version

**NOTE**:

The repo has been cleaned using these tools.

The old repo was still big as it contained accidentially commited binaries.
These were cleaned with

```
git clone --mirror https://github.com/CELMA-project/CELMA.git
wget -O bfg.jar http://repo1.maven.org/maven2/com/madgag/bfg/1.12.15/bfg-1.12.15.jar
java -jar bfg.jar --strip-blobs-bigger-than 20M CELMA.git
cd CELMA.git
git reflog expire --expire=now --all && git gc --prune=now --aggressive
git push
```

as described on [bfg](https://rtyley.github.io/bfg-repo-cleaner/)s homepage.

Note that pull requests are not cleaned (leaving only the mirrored repo large,
unless the PR is removed), as described on
[this SO thread](http://stackoverflow.com/questions/34265266/remote-rejected-errors-after-mirroring-a-git-repository).
