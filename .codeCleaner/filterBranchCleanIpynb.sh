#!/usr/bin/env bash

# Script meant to be used with git filter-branch to clean the output of *.ipynb
# files in earlier commits
# Call the following command from the repository's root folder
# export ROOT_REPO=${PWD} && git filter-branch -f --tree-filter 'bash ${ROOT_REPO}/.codeCleaner/filterBranchCleanIpynb.sh ${ROOT_REPO}' --tag-name-filter cat -- --all
# NOTE: For this to work, the path/to/repo cannot contain spaces

# Directory of calling dir
REPO_ROOT=$1

set -e

# Recall that this file will be called from the parent directory, so all paths
# are relative to the parent directory

echo "Processing commit: ${GIT_COMMIT}" >> ${REPO_ROOT}/.codeCleaner/filterBranch.log 2>&1
find . -iname "*.ipynb" | while read file; do
    echo "    Processing '${file}'" >> ${REPO_ROOT}/.codeCleaner/filterBranch.log 2>&1
    python ${REPO_ROOT}/.codeCleaner/cleanIpynb.py ${file} >> ${REPO_ROOT}/.codeCleaner/filterBranch.log 2>&1
done
