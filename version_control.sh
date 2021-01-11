#!/usr/bin/env bash

# don't add the local genome file I'm working on!
git add README.md
git add fw_out/
git add Cargo.toml
git add Cargo.lock
git add src/
git add version_control.sh

git commit -m "$1"

git push