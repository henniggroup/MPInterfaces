#!/bin/bash

# run from master branch docs folder
# . ./update_gh-pages.sh

TMP_DOC_DIR="/tmp/mpint_doc"
make html
rm -rf $TMP_DOC_DIR
cp -r _build/html $TMP_DOC_DIR
git checkout gh-pages
cd ..
cp -r $TMP_DOC_DIR/* .
git status
echo "Commit and push, check http://henniggroup.github.io/MPInterfaces/"
