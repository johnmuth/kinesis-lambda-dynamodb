#!/usr/bin/env bash

set -e

docker-compose down
docker-compose run pubmed-search ./run.sh $@
docker-compose down

