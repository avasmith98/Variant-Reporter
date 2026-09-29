import unittest
from unittest.mock import Mock, patch

import requests

import Variant_Reporter as reporter


class DomainTests(unittest.TestCase):
    @patch.object(reporter.time, 'sleep')
    @patch.object(reporter, 'send_request')
    def test_transient_failure_recovers(self, send, sleep):
        send.side_effect = [requests.HTTPError(response=Mock(status_code=500)),
                            Mock(json=Mock(return_value=[{
                                'type': 'Pfam', 'start': 1, 'end': 20,
                                'description': 'Kinase',
                            }]))]
        self.assertEqual(reporter.get_domains('ENSP00000275493'), [{
            'Source': 'Pfam', 'Start': '1', 'End': '20', 'Description': 'Kinase',
        }])
        self.assertEqual(send.call_count, 2)

    @patch.object(reporter.time, 'sleep')
    @patch.object(reporter, 'send_request')
    def test_persistent_outage_prevents_pdf_generation(self, send, sleep):
        send.side_effect = requests.HTTPError(response=Mock(status_code=500))
        with patch.object(reporter, 'get_results_dict',
                          side_effect=lambda _: reporter.get_domains('ENSP00000275493')), \
                patch.object(reporter, 'create_pdf_from_curated_data') as create_pdf:
            response = reporter.app.test_client().post('/', data={
                'hgvs_cdna_transcript_id': 'NM_005228.3(EGFR):c.2648T>C (p.Leu883Ser)',
            })
            self.assertEqual(response.status_code, 422)
            self.assertNotEqual(response.mimetype, 'application/pdf')
            self.assertIn(b'No report was generated', response.data)
            create_pdf.assert_not_called()
        self.assertEqual(send.call_count, 3)

    @patch.object(reporter.time, 'sleep')
    @patch.object(reporter, 'send_request')
    def test_invalid_payload_and_timeouts_are_unavailable(self, send, sleep):
        for response in (Mock(json=Mock(return_value={'error': 'unavailable'})),
                         Mock(json=Mock(side_effect=ValueError('not JSON'))),
                         Mock(json=Mock(return_value=[{'type': 'Pfam'}])),
                         Mock(json=Mock(return_value=[{'type': 'Pfam', 'start': 1, 'end': 2, 'description': None}]))):
            send.return_value = response
            with self.assertRaisesRegex(ValueError, 'No report was generated'):
                reporter.get_domains('ENSP00000275493')
        send.side_effect = requests.Timeout()
        with self.assertRaisesRegex(ValueError, 'No report was generated'):
            reporter.get_domains('ENSP00000275493')

    @patch.object(reporter, 'send_request')
    def test_empty_results_prevent_report(self, send):
        send.return_value = Mock(json=Mock(return_value=[]))
        with self.assertRaisesRegex(ValueError, 'No protein-domain data'):
            reporter.get_domains('ENSP00000275493')

    @patch.object(reporter, 'send_request')
    def test_permanent_http_error_is_not_retried(self, send):
        send.side_effect = requests.HTTPError(response=Mock(status_code=404))
        with self.assertRaisesRegex(ValueError, 'No report was generated'):
            reporter.get_domains('ENSP00000275493')
        self.assertEqual(send.call_count, 1)


if __name__ == '__main__':
    unittest.main()
