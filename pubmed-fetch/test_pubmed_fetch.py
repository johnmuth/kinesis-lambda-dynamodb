import unittest
import pubmed_fetch
import logging
import json

logging.basicConfig()
logger = logging.getLogger('test_pubmed_fetch')
logger.setLevel(logging.DEBUG)

class PubmedFetchTestCase(unittest.TestCase):

    def test_handler(self):
        with open('kinesis.json') as json_data:
            event = json.load(json_data)
            pubmed_fetch.handler(event, {})

if __name__ == '__main__':
    unittest.main()