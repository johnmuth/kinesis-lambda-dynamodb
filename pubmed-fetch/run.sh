#!/usr/bin/env bash

set -e

source /app/bin/activate

export LANG=en_GB.UTF-8

python pubmed_fetch.py $@
