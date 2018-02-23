# kinesis-lambda-dynamodb

Demonstration using Kinesis, Lambda and DynamoDB together.

I want to do a search in PubMed, fetch the results and store them in DynamoDB.

I'm hoping Kinesis will be a good way to control the rate at which we write to DynamoDB.

## WIP!

## TODO

- don't hard-code pubmed-records s3 bucket
- create pubmed-records s3 bucket with cloudformation  

## Quick start

Clone this repo

Package the lambda functions and upload to S3

```bash
./pubmed-search/docker-run-package.sh
./pubmed-fetch/docker-run-package.sh
aws s3 cp pubmed-fetch/dist/pubmed-fetch-lambda.zip s3://my-lambda-zip-bucket/
aws s3 cp pubmed-search/dist/pubmed-search-lambda.zip s3://my-lambda-zip-bucket/
```

Create stack

```bash
./create-stack.sh
```

Invoke the pubmed-search lambda

```bash
./invoke-pubmed-search-lambda.sh
```

