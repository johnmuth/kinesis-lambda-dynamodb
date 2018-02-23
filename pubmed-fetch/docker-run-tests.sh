#!/usr/bin/env bash

set -e

docker-compose down
docker-compose run pubmed-fetch ./run_tests.sh $@
docker-compose down

