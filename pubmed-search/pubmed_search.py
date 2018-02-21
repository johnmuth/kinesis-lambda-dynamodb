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

def pubmed_search(search):
    logger.info("Executing pubmed search")
    search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search, usehistory="y"), validate=False)
    count = int(search_results["Count"])
    web_env = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    logger.debug('Got {} results. web_env={} query_key={})'.format(str(count), web_env, query_key))


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
