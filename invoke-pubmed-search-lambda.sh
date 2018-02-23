#!/usr/bin/env bash

search=`echo '("macular degeneration"[MeSH Terms] OR ("macular"[All Fields] AND "degeneration"[All Fields]) OR "macular degeneration"[All Fields]) AND ("loattrfree full text"[sb] AND "2013/02/25"[PDat] : "2018/02/23"[PDat])'|sed 's/"/\\\"/g'`
#search="macular degeneration"

payload='{"Records":[{"Sns":{"MessageAttributes":{"search":{"Value":"'"$search"'"}}}}]}'

echo "$payload" | jq '.'

echo aws lambda invoke --function-name PubmedSearchKLD2 --payload "$payload" --invocation-type Event outfile.txt
aws lambda invoke --function-name PubmedSearchKLD2 --payload "$payload"  --invocation-type Event outfile.txt


