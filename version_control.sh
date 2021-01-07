#!/usr/bin/env bash

git add README.md
git add test.csv
git add Cargo.toml
git add Cargo.lock
git add src/
git add version_control.sh

git commit -m "$1"

git push