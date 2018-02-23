#!/usr/bin/env bash

set -e
set -u

./pubmed-search/docker-run-package.sh
./pubmed-fetch/docker-run-package.sh
aws s3 cp pubmed-fetch/dist/pubmed-fetch-lambda.zip s3://builds.evidentia.com/
aws s3 cp pubmed-search/dist/pubmed-search-lambda.zip s3://builds.evidentia.com/

