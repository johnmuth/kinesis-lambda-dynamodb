import boto3
import logging
from datetime import datetime
import json
from Bio import Entrez
import base64
import os
from urllib.error import HTTPError

Entrez.email = "support@evidentiapublishing.com"
Entrez.cache = "/tmp/"
output_bucket = 'pubmed.evidentia.com'
s3 = boto3.resource('s3', region_name='us-east-1')
batch_size = 100

logging.basicConfig()
logger = logging.getLogger('pubmed_fetch')
logger.setLevel(logging.DEBUG)

kinesis_client = boto3.client('kinesis', region_name='us-east-1')
stream_name = 'PubmedSearchResults'

def fetch_and_save_input_file(start, web_env, query_key, search_result_set_id):
    attempt = 1
    succeeded = False
    while attempt <= 3 and not succeeded:
        try:
            logger.info("Entrez.efetch(db=pubmed, retmode=xml, retstart={}, retmax={}, webenv={}, query_key={})".format(start, batch_size, web_env, query_key))
            fetch_handle = Entrez.efetch(db="pubmed",
                                         retmode="xml",
                                         retstart=start,
                                         retmax=batch_size,
                                         webenv=web_env,
                                         query_key=query_key)
            data = fetch_handle.read()
            fetch_handle.close()
            key = '{}.{}.xml'.format(search_result_set_id, start)
            local_file = '/tmp/{}'.format(key)
            logger.info("Writing to %s" % (local_file))
            out_handle = open(local_file, "w")
            out_handle.write(data)
            out_handle.close()
            s3uri = 's3://{}/{}'.format(output_bucket, key)
            bucket = s3.Bucket(output_bucket)
            logger.info("Uploading to %s" % (s3uri))
            bucket.upload_file(local_file, key)
            os.remove(local_file)
            succeeded = True

        except HTTPError as err:
            if 500 <= err.code <= 599:
                logger.warn("Received error from server %s" % err)
                logger.warn("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(3)
            else:
                logger.info("Received error from server %s" % err)
                raise

def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    if isinstance(obj, datetime):
        serial = obj.isoformat()
        return serial
    raise TypeError("Type not serializable")


def handler(event, context):
    logger.debug('event: %s', json.dumps(event, default=json_serial))
    if 'Records' not in event:
        raise Exception('No element Records in event: %s', json.dumps(event, default=json_serial))
    for record in event['Records']:
        if 'kinesis' not in record or 'data' not in record['kinesis']:
            raise Exception('No element Records/kinesis/data in event so cannot process')
        record_data = json.loads(base64.b64decode(record['kinesis']['data']))
        logger.debug('record_data: {}'.format(record_data))
        fetch_and_save_input_file(record_data['start'], record_data['web_env'], record_data['query_key'], record_data['search_result_set_id'])


if __name__ == '__main__':
    handler({}, {})
