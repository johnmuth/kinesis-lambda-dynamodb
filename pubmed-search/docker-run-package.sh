#!/usr/bin/env bash

set -e
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DOCKER_COMPOSE_FILE="$SCRIPT_DIR"/docker-compose.yml

docker-compose -f "$DOCKER_COMPOSE_FILE" down
docker-compose -f "$DOCKER_COMPOSE_FILE" build
docker-compose -f "$DOCKER_COMPOSE_FILE" run pubmed-search ./package.sh $@
docker-compose -f "$DOCKER_COMPOSE_FILE" down

