#!/usr/bin/env bash

set -e

source /app/bin/activate

export LANG=en_GB.UTF-8

python test_pubmed_fetch.py $@
