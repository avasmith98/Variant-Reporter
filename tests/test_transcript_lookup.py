import unittest
from unittest.mock import Mock, patch
import requests
import Variant_Reporter as reporter


class TranscriptLookupTests(unittest.TestCase):
    def transcript(self, accession='NM_005228.5', identifier='ENST00000275493'):
        return {'id': identifier, 'assembly_name': 'GRCh38',
                'MANE': [{'id': identifier, 'refseq_match': accession}]}

    @patch.object(reporter, 'send_request')
    def test_exact_mapping(self, send):
        send.return_value.json.return_value = {'Transcript': [self.transcript()]}
        self.assertEqual(reporter.get_ensembl_transcript_id('NM_005228.5', 'EGFR'), 'ENST00000275493')
        self.assertNotIn('/xrefs/', send.call_args.args[0])

    @patch.object(reporter, 'send_request')
    def test_no_default_when_missing_wrong_version_or_ambiguous(self, send):
        for transcripts in ([], [self.transcript('NM_005228.4')],
                            [self.transcript(), self.transcript(identifier='ENST_OTHER')],
                            [{'id': 'ENST00000275493'}]):
            with self.subTest(transcripts=transcripts):
                send.return_value.json.return_value = {'Transcript': transcripts}
                with self.assertRaisesRegex(ValueError, 'No unique exact'):
                    reporter.get_ensembl_transcript_id('NM_005228.5', 'EGFR')

    @patch.object(reporter, 'send_request')
    def test_wrong_assembly_or_inconsistent_mane_id_fails(self, send):
        transcript = self.transcript()
        transcript['assembly_name'] = 'GRCh37'
        send.return_value.json.return_value = {'Transcript': [transcript]}
        with self.assertRaises(ValueError):
            reporter.get_ensembl_transcript_id('NM_005228.5', 'EGFR')
        transcript['assembly_name'] = 'GRCh38'
        transcript['MANE'][0]['id'] = 'ENST_OTHER'
        with self.assertRaises(ValueError):
            reporter.get_ensembl_transcript_id('NM_005228.5', 'EGFR')

    @patch.object(reporter, 'send_request')
    def test_outage_fails(self, send):
        send.side_effect = requests.HTTPError('500')
        with self.assertRaises(requests.HTTPError):
            reporter.get_ensembl_transcript_id('NM_005228.5', 'EGFR')
