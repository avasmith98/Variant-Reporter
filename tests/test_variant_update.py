import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import Mock, patch
from xml.etree import ElementTree as ET
import requests
import Variant_Reporter as reporter

VARIANT = 'NM_005228.3(EGFR):c.2648T>C (p.Leu883Ser)'
FIXTURE = Path(__file__).parent / 'fixtures/egfr_clinvar.xml'


class VariantUpdateTests(unittest.TestCase):
    def setUp(self):
        self.root = ET.parse(FIXTURE).getroot()

    def hgvs(self, accession):
        return next(h for h in self.root.findall('.//HGVS')
                    if h.find('NucleotideExpression').get('sequenceAccessionVersion') == accession)

    def resolve(self, value=VARIANT, count='1', ids=None):
        search = Mock(json=Mock(return_value={'esearchresult': {'count': count, 'idlist': ['1050891'] if ids is None else ids}}))
        response = Mock(content=ET.tostring(self.root))
        with patch.object(reporter, 'get_current_version_hgvs_cdna_transcript_id', return_value='NM_005228.5'), \
             patch.object(reporter, 'send_request', side_effect=[search, response]) as send, \
             patch.object(reporter, 'verify_refseq_protein') as verify:
            result = reporter.resolve_variant(value)
            self.assertIn('NM_005228.3', send.call_args_list[0].args[0])
            return result

    def test_same_allele_update(self):
        result = self.resolve()
        self.assertEqual(result['transcript'], 'NM_005228.5')
        self.assertEqual(result['cdna'], 'c.2648T>C')
        self.assertEqual(result['protein'], '(p.Leu883Ser)')
        self.assertEqual(result['position'], 55192788)
        self.assertEqual(len(result['classifications']), 3)

    def test_changed_coordinates_and_protein_come_from_current_annotation(self):
        current = self.hgvs('NM_005228.5')
        current.find('NucleotideExpression').set('change', 'c.2651T>C')
        current.find('NucleotideExpression/Expression').text = 'NM_005228.5:c.2651T>C'
        current.find('ProteinExpression').set('change', 'p.Leu884Ser')
        result = self.resolve()
        self.assertEqual(result['cdna'], 'c.2651T>C')
        self.assertEqual(result['protein'], '(p.Leu884Ser)')

    def test_original_expression_must_be_present(self):
        self.hgvs('NM_005228.3').find('NucleotideExpression/Expression').text = 'NM_005228.3:c.2648T>A'
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            self.resolve()

    def test_current_version_must_be_present(self):
        self.hgvs('NM_005228.5').find('NucleotideExpression').set('sequenceAccessionVersion', 'NM_005228.4')
        with self.assertRaisesRegex(ValueError, 'No unique same-allele'):
            self.resolve()

    def test_ambiguous_current_mapping_fails(self):
        import copy
        self.root.find('.//HGVSlist').append(copy.deepcopy(self.hgvs('NM_005228.5')))
        with self.assertRaisesRegex(ValueError, 'No unique same-allele'):
            self.resolve()

    def test_wrong_gene_and_empty_search_fail(self):
        with self.assertRaises(ValueError):
            self.resolve(VARIANT.replace('EGFR', 'OTHER'))
        with self.assertRaises(ValueError):
            self.resolve(count='0', ids=[])
        with self.assertRaises(ValueError):
            self.resolve(count='101')

    def test_wrong_protein_fails(self):
        ET.SubElement(self.hgvs('NM_005228.3'), 'ProteinExpression', change='p.Leu883Val')
        with self.assertRaisesRegex(ValueError, 'protein change'):
            self.resolve()

    def test_missing_current_protein_fails(self):
        current = self.hgvs('NM_005228.5')
        current.remove(current.find('ProteinExpression'))
        with self.assertRaisesRegex(ValueError, 'annotation is missing'):
            self.resolve()

    def test_missing_classification_fails(self):
        record = self.root.find('ClassifiedRecord')
        record.remove(record.find('RCVList'))
        with self.assertRaisesRegex(ValueError, 'classification data is missing'):
            self.resolve()

    def test_original_protein_verified_from_codon(self):
        # Synthetic CDS with leucine CTG at codon 2, c.5T>C produces CCG (Pro).
        record = ET.fromstring('<GBSeq><GBSeq_sequence>ATGCTGTAA</GBSeq_sequence>'
                               '<GBSeq_feature-table><GBFeature><GBFeature_key>CDS</GBFeature_key>'
                               '<GBFeature_location>1..9</GBFeature_location></GBFeature>'
                               '</GBSeq_feature-table></GBSeq>')
        with patch.object(reporter, 'get_refseq_record', return_value=record):
            reporter.verify_refseq_protein('NM_1.1', 'c.5T>C', 'p.Leu2Pro')
            for cdna, protein in [('c.5T>C', 'p.Leu2Ser'), ('c.5A>C', 'p.Leu2Pro'),
                                  ('c.50T>C', 'p.Leu2Pro'), ('c.5del', 'p.Leu2Pro')]:
                with self.assertRaises(ValueError):
                    reporter.verify_refseq_protein('NM_1.1', cdna, protein)

    def test_error_response_never_creates_pdf(self):
        for failure, status in [(ValueError('Mapping not verified'), 422), (requests.Timeout(), 503)]:
            with patch.object(reporter, 'get_results_dict', side_effect=failure), \
                 patch.object(reporter, 'create_pdf_from_curated_data') as pdf:
                response = reporter.app.test_client().post('/', data={'hgvs_cdna_transcript_id': VARIANT})
                self.assertEqual(response.status_code, status)
                self.assertIn(b'No report was generated', response.data)
                pdf.assert_not_called()

    def test_invalid_input_fails_before_network(self):
        with patch.object(reporter, 'send_request') as send:
            with self.assertRaises(ValueError):
                reporter.resolve_variant('NM_005228(EGFR):c.2648T>C')
            send.assert_not_called()

    def test_report_uses_updated_annotations_and_keeps_original_input(self):
        variant = self.resolve()
        variant.update(cdna='c.2651T>C', protein='(p.Leu884Ser)')
        responses = {
            'resolve_variant': variant,
            'get_full_gene_name_and_ensembl_gene_id': ('epidermal growth factor receptor', 'ENSG00000146648'),
            'get_gene_start_end_chromosome': ('chr7', '55019017', '55211628'),
            'get_cytogenetic_band': '7p11.2',
            'get_high_protein_expression': 'skin',
            'get_ensembl_transcript_id': 'ENST00000275493',
            'get_transcript_details_and_ensembl_protein_id': ('ENSP00000275493', 10000, 1210, 28, 28, 22),
            'get_domains': [{'Source': 'Pfam', 'Description': 'Kinase', 'Start': '700', 'End': '1000'}],
        }
        with ExitStack() as stack:
            mocks = {name: stack.enter_context(patch.object(reporter, name, return_value=value))
                     for name, value in responses.items()}
            result = reporter.get_results_dict(VARIANT)
            self.assertEqual(result['Original variant'], VARIANT)
            self.assertEqual(result['Current HGVS ID'], 'NM_005228.5')
            self.assertEqual(result['cDNA change'], 'c.2651T>C')
            self.assertEqual(result['Protein change'], '(p.Leu884Ser)')
            mocks['get_ensembl_transcript_id'].assert_called_with('NM_005228.5', 'EGFR')
            response = reporter.app.test_client().post('/', data={'hgvs_cdna_transcript_id': VARIANT})
            self.assertEqual(response.status_code, 200)
            self.assertTrue(response.data.startswith(b'%PDF-'))

    def test_refseq_wrong_accession_or_version_fails(self):
        response = Mock(content=b'<GBSet><GBSeq><GBSeq_accession-version>NM_005228.4</GBSeq_accession-version></GBSeq></GBSet>')
        with patch.object(reporter, 'send_request', return_value=response):
            with self.assertRaises(ValueError):
                reporter.get_refseq_record('NM_005228.3')
            with self.assertRaises(ValueError):
                reporter.get_refseq_record('NM_999999')
