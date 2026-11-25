#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Amir H. Amanzadi
# Build the documentation site and publish it to the gh-pages branch.
# GitHub Pages serves that branch (Settings > Pages > Deploy from a branch > gh-pages, / (root)).
set -eu
ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT"
REMOTE="$(git remote get-url origin)"
SHA="$(git rev-parse --short HEAD)"
rm -rf site
uvx zensical==0.1.3 build --strict
touch site/.nojekyll   # the site is already built: tell Pages not to run Jekyll
cd site
git init -q -b gh-pages
git add -A
git commit -q -m "Publish docs from $SHA"
git push -q -f "$REMOTE" gh-pages   # one commit only: the branch never grows
rm -rf .git
echo "Published docs from $SHA to gh-pages."
