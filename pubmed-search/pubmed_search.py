import boto3
import logging
from datetime import datetime
import json
from Bio import Entrez
import sys

Entrez.email = "support@evidentiapublishing.com"
Entrez.cache = "/tmp/"

logging.basicConfig()
logger = logging.getLogger('pubmed_search')
logger.setLevel(logging.DEBUG)

kinesis_client = boto3.client('kinesis', region_name='us-east-1')
stream_name = 'PubmedSearchResults'
partition_key = '1'
records_per_put = 500
results_per_record = 100

def pubmed_search(search):
    logger.info("Executing pubmed search")
    search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search, usehistory="y"), validate=False)
    count = int(search_results["Count"])
    web_env = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    logger.debug('Got {} results. web_env={} query_key={})'.format(str(count), web_env, query_key))
    record_num = 0
    put_num = 0
    for start_result_num in range(1, count, records_per_put * results_per_record):
        put_num += 1
        end_result_num = min(start_result_num + (records_per_put * results_per_record), count)
        records = []
        for first_record_num in range(start_result_num, end_result_num, results_per_record):
            record_num += 1
            last_result_num = min((first_record_num + results_per_record - 1), count)
            logger.debug('record_num={} put_num={} first_record_num={} last_record_num={}'.format(record_num, put_num, first_record_num, last_result_num))
            records.append({
                'Data': json.dumps({
                    'start': first_record_num,
                    'end': last_result_num,
                    'web_env': web_env,
                    'query_key': query_key,
                }),
                'PartitionKey': partition_key
            })

        kinesis_client.put_records(
            Records=records,
            StreamName=stream_name
        )


def get_message_attribute(event, attribute):
    if 'Records' not in event:
        raise Exception('No element Records in event: %s', json.dumps(event, default=json_serial))
    for record in event['Records']:
        if 'Sns' not in record or 'MessageAttributes' not in record['Sns'] or attribute not in record['Sns'][
            'MessageAttributes']:
            raise Exception('No element Records/Sns/MessageAttributes/{} in event so cannot process'.format(attribute))
        return record['Sns']['MessageAttributes'][attribute]['Value']


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    if isinstance(obj, datetime):
        serial = obj.isoformat()
        return serial
    raise TypeError("Type not serializable")


def handler(event, context):
    logger.debug('event: %s', json.dumps(event, default=json_serial))
    search = get_message_attribute(event, 'search')
    if search:
        pubmed_search(search)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        with open(sys.argv[1], 'r') as searchfile:
            search = searchfile.read().replace('\n', '')
            pubmed_search(search)
    else:
        print('Usage: ' + sys.argv[0] + ' <search file>')
