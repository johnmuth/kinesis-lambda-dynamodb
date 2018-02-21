#!/usr/bin/env bash

stackName='kinesis-lambda-dynamodb'

echo "Calling create-stack"
aws cloudformation create-stack --capabilities CAPABILITY_NAMED_IAM --stack-name "$stackName" --template-body "file://cloudformation.$stackName.json"
echo "Waiting for create-stack to complete"
aws cloudformation wait stack-create-complete --stack-name "$stackName"
echo "create-stack complete"

