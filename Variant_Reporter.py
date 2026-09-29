#hgvs_cdna_transcript_id = "NM_005228.3(EGFR):c.2648T>C (p.Leu883Ser)"
import datetime
import requests
import re
import time
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, PageBreak, Spacer, Paragraph
from reportlab.lib import colors
from flask import Flask, request, render_template_string, send_file
from xml.etree import ElementTree as ET
from io import BytesIO
from urllib.parse import urlencode

def get_report_date():
    '''Gets the current date in the format 'Month Day, Year'.'''
    report_date = datetime.date.today().strftime('%B %d, %Y')
    return report_date
    
def send_request(url, headers=None):
    '''Sends a request to the specified URL and return the response.'''
    response = requests.get(url, headers=headers, timeout=45)
    response.raise_for_status()
    return response
    
def parse_variant(value):
    match = re.fullmatch(r"(NM_\d+\.\d+)\(([A-Za-z0-9_-]+)\):(c\.[^\s()]+)\s*\((p\.[^()]+)\)", value.strip())
    if not match:
        raise ValueError('Enter a versioned RefSeq transcript, gene, cDNA change and protein change.')
    return match.groups()


def get_refseq_record(accession):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
        'db': 'nucleotide', 'id': accession, 'rettype': 'gb', 'retmode': 'xml'})
    records = ET.fromstring(send_request(url).content).findall('.//GBSeq')
    if len(records) != 1:
        raise ValueError(f'Expected one RefSeq record for {accession}.')
    actual = records[0].findtext('GBSeq_accession-version')
    if not actual or (actual != accession if '.' in accession else actual.split('.')[0] != accession):
        raise ValueError(f'RefSeq did not return the requested transcript {accession}.')
    return records[0]


def get_current_version_hgvs_cdna_transcript_id(accession):
    return get_refseq_record(accession).findtext('GBSeq_accession-version')


def coding_sequence(record):
    features = [f for f in record.findall('./GBSeq_feature-table/GBFeature')
                if f.findtext('GBFeature_key') == 'CDS']
    if len(features) != 1:
        raise ValueError('Cannot verify the original protein change: missing or ambiguous CDS.')
    location = features[0].findtext('GBFeature_location', '')
    match = re.fullmatch(r'(\d+)\.\.(\d+)', location)
    sequence = record.findtext('GBSeq_sequence')
    if not match or not sequence:
        raise ValueError('Cannot verify the original protein change from this CDS.')
    start, end = map(int, match.groups())
    if not 1 <= start <= end <= len(sequence):
        raise ValueError('Invalid RefSeq CDS coordinates.')
    return sequence[start - 1:end].upper()


def verify_refseq_protein(accession, cdna, protein):
    """Verify a coding substitution directly when historical HGVS lacks protein data."""
    match = re.fullmatch(r'c\.(\d+)([ACGT])>([ACGT])', cdna)
    if not match:
        raise ValueError('No historical protein annotation; this variant cannot be verified automatically.')
    sequence = coding_sequence(get_refseq_record(accession))
    position = int(match[1]) - 1
    if not 0 <= position < len(sequence) or sequence[position] != match[2]:
        raise ValueError('The original cDNA reference base does not match RefSeq.')
    codon_start = position // 3 * 3
    codon = sequence[codon_start:codon_start + 3]
    changed = list(codon)
    changed[position % 3] = match[3]
    bases = 'TCAG'
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codons = [a + b + c for a in bases for b in bases for c in bases]
    table = dict(zip(codons, amino_acids))
    names = dict(zip('FLSY*CWPHQRIMTNKVADEG',
                     ['Phe', 'Leu', 'Ser', 'Tyr', 'Ter', 'Cys', 'Trp', 'Pro',
                      'His', 'Gln', 'Arg', 'Ile', 'Met', 'Thr', 'Asn', 'Lys',
                      'Val', 'Ala', 'Asp', 'Glu', 'Gly']))
    if codon not in table or ''.join(changed) not in table:
        raise ValueError('The original coding sequence contains an unverifiable codon.')
    before, after = table[codon], table[''.join(changed)]
    expected = f'p.{names[before]}{position // 3 + 1}{names[after]}'
    if before == after or protein != expected:
        raise ValueError('The supplied protein change does not match the original RefSeq codon.')


def resolve_variant(value):
    original, gene, cdna, protein = parse_variant(value)
    current = get_current_version_hgvs_cdna_transcript_id(original.split('.')[0])
    if (not re.fullmatch(r'NM_\d+\.\d+', current) or
            current.split('.')[0] != original.split('.')[0] or
            int(current.split('.')[1]) < int(original.split('.')[1])):
        raise ValueError('NCBI did not return a valid current version of the original transcript.')
    # Search the ORIGINAL expression, then require explicit same-allele evidence.
    expression = f'{original}:{cdna}'
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urlencode({
        'db': 'clinvar', 'term': f'"{expression}"', 'retmode': 'json', 'retmax': 100})
    search = send_request(url).json()['esearchresult']
    ids = search['idlist']
    if not ids or int(search['count']) != len(ids):
        raise ValueError('No complete ClinVar search result for the original variant.')
    matches = []
    for variation_id in ids:
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
            'db': 'clinvar', 'rettype': 'vcv', 'is_variationid': 'true', 'id': variation_id})
        root = ET.fromstring(send_request(url).content)
        for allele in root.findall('.//ClassifiedRecord/SimpleAllele'):
            if allele.get('VariationID') != variation_id:
                continue
            hgvs = allele.findall('./HGVSlist/HGVS')
            old = [h for h in hgvs if h.findtext('./NucleotideExpression/Expression') == expression]
            if old and gene in {g.get('Symbol') for g in allele.findall('./GeneList/Gene')}:
                matches.append((variation_id, allele, old, hgvs, root))
    if len(matches) != 1:
        raise ValueError('The original variant does not match exactly one ClinVar allele.')
    variation_id, allele, old, hgvs, root = matches[0]
    updated = [h for h in hgvs if h.find('NucleotideExpression') is not None and
               h.find('NucleotideExpression').get('sequenceAccessionVersion') == current]
    if len(updated) != 1:
        raise ValueError(f'No unique same-allele mapping to updated transcript {current}.')
    nucleotide = updated[0].find('NucleotideExpression')
    updated_cdna = nucleotide.get('change')
    updated_protein = updated[0].find('ProteinExpression')
    if not updated_cdna or not updated_cdna.startswith('c.') or updated_protein is None or not updated_protein.get('change'):
        raise ValueError('Updated cDNA or protein annotation is missing.')
    if nucleotide.findtext('Expression') != f'{current}:{updated_cdna}':
        raise ValueError('Inconsistent updated transcript annotation.')
    old_proteins = {h.find('ProteinExpression').get('change') for h in old
                    if h.find('ProteinExpression') is not None}
    if not old_proteins:
        verify_refseq_protein(original, cdna, protein)
        old_proteins = {protein}
    if old_proteins != {protein}:
        raise ValueError('The supplied protein change does not match the original variant.')
    locations = allele.findall('./Location/SequenceLocation')
    locations = [loc for loc in locations if loc.get('Assembly') == 'GRCh38']
    rsids = {x.get('ID') for x in allele.findall('./XRefList/XRef') if x.get('DB') == 'dbSNP' and x.get('Type') == 'rs'}
    if len(locations) != 1 or len(rsids) != 1:
        raise ValueError('Missing or ambiguous GRCh38 location or dbSNP identifier.')
    location = locations[0]
    if not location.get('start') or not location.get('stop') or not location.get('Chr'):
        raise ValueError('Incomplete GRCh38 location.')
    return {'transcript': current, 'gene': gene, 'cdna': updated_cdna,
            'protein': f"({updated_protein.get('change')})", 'clinvar_id': variation_id,
            'rsid': 'rs' + next(iter(rsids)), 'position': int(location.get('start')),
            'end': int(location.get('stop')), 'chromosome': location.get('Chr'),
            'classifications': get_clinvar(root)}


def get_full_gene_name_and_ensembl_gene_id(gene_symbol):
    '''Gets the full gene name and Ensembl gene ID for the specified gene symbol.'''
    url = f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers=headers)
    data = response.json()
    description = data.get('description')
    match = re.match(r"([^[]+)", description)
    full_gene_name = match.group(1).strip()
    ensembl_gene_id = data['id']
    return full_gene_name, ensembl_gene_id

def get_gene_start_end_chromosome(gene_symbol):
    '''Gets the genomic start/end positions and chromosome number of the specified gene symbol.'''
    url = f'https://api.genome.ucsc.edu/search?search={gene_symbol}&genome=hg38'
    response = send_request(url)
    data = response.json()
    for match in data['positionMatches'][0]['matches']:
        if gene_symbol in match['posName'] and 'ENST' in match['hgFindMatches']:
            position = match['position']
            chromosome, pos_range = position.split(':')
            start, end = pos_range.split('-')
            return chromosome, start, end
        
def get_cytogenetic_band(chromosome, start, end):
    '''Gets the cytogenetic band for the specified chromosome, start, and end positions.'''
    url = f'https://api.genome.ucsc.edu/getData/track?track=cytoBand;genome=hg38;chrom={chromosome};start={start};end={end}'
    response = send_request(url)
    data = response.json()
    cytogenetic_band = data['cytoBand']
    chromosome = cytogenetic_band[0]['chrom'][3:]
    cytoband = cytogenetic_band[0]['name']
    cytogenetic_band = chromosome + cytoband
    return cytogenetic_band

def get_high_protein_expression(ensembl_gene_id):
    '''Gets the tissues where the protein is highly expressed for the specified Ensembl gene ID.
    High protein expression is defined as having a level of 'high' for at least one tissue type in the Protein Atlas database.
    Example: Nasopharynx consits of basal cells, ciliated cells, etc. If 'basal cells' has high protein expression, then nasopharynx will be returned.'''
    url = f'https://www.proteinatlas.org/{ensembl_gene_id}.xml'
    response = send_request(url)
    root = ET.fromstring(response.content)
    high_protein_expression = []
    if not root.findall('.//data/level[@type="expression"]'):
        raise ValueError('Protein expression data is missing.')
    for data in root.findall('.//data'):
        tissue = data.find('tissue')
        levels = data.findall('level[@type="expression"]')
        if any(level.text.lower() == "high" for level in levels):
            high_protein_expression.append(tissue.text)
    if high_protein_expression:
        return(', '.join(high_protein_expression).lower())
    return('Protein not highly expressed.')

def get_ensembl_transcript_id(current_version_hgvs_cdna_transcript_id, gene_symbol):
    """Require one explicit, versioned RefSeq-to-Ensembl MANE mapping."""
    url = f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;mane=1'
    data = send_request(url, {'Content-Type': 'application/json'}).json()
    matches = {transcript['id'] for transcript in data.get('Transcript', [])
               for mapping in transcript.get('MANE', [])
               if mapping.get('refseq_match') == current_version_hgvs_cdna_transcript_id
               and mapping.get('id') == transcript.get('id')
               and transcript.get('assembly_name') == 'GRCh38'}
    if len(matches) != 1:
        raise ValueError(f'No unique exact Ensembl mapping for {current_version_hgvs_cdna_transcript_id}.')
    return matches.pop()

def get_transcript_details_and_ensembl_protein_id(ensembl_transcript_id, grch38_variant_position, variant_end, chromosome):
    '''Gets the Ensembl protein ID, transcript length, translation length, total exons, coding exons,
    exon number, and coding exons for the specified Ensembl transcript ID and variant position.'''
    
    url = f'https://rest.ensembl.org/lookup/id/{ensembl_transcript_id}?expand=1'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers)
    data = response.json()

    if (data.get('assembly_name') != 'GRCh38' or data.get('id') != ensembl_transcript_id
            or data.get('seq_region_name') != chromosome):
        raise ValueError('Ensembl returned a different transcript or assembly.')

    # Extract basic information
    if data['assembly_name'] == 'GRCh38':
        transcript_length = data['length']
        translation = data['Translation']
        translation_length = translation['length']
        ensembl_protein_id = translation['id']

    # Total number of GRCh38 exons
    total_exons = sum(1 for exon in data['Exon'] if exon.get('assembly_name') == 'GRCh38')

    # Number of coding exons in GRCh38
    coding_exons = sum(1 for exon in data['Exon']
                       if exon.get('assembly_name') == 'GRCh38' and 
                       int(exon['start']) <= int(translation['end']) and 
                       int(exon['end']) >= int(translation['start']))

    if data.get('strand') not in (1, -1):
        raise ValueError('Transcript strand is missing.')
    exons = sorted(data['Exon'], key=lambda exon: int(exon['start']), reverse=data['strand'] == -1)
    # Number exons in transcript order, including reverse-strand transcripts.
    exon_number = next((index + 1 for index, exon in enumerate(exons)
                        if exon.get('assembly_name') == 'GRCh38' and 
                        int(exon['start']) <= int(grch38_variant_position) <= variant_end <= int(exon['end'])), None)

    if exon_number is None:
        raise ValueError('Variant position is not in the verified transcript exons.')

    return ensembl_protein_id, transcript_length, translation_length, total_exons, coding_exons, exon_number


def get_domains(ensembl_protein_id):
    '''Gets protein domains, retrying temporary failures and failing if unavailable.'''
    url = f'https://rest.ensembl.org/overlap/translation/{ensembl_protein_id}'
    headers = {'Content-Type': 'application/json'}
    for attempt in range(3):
        try:
            response = send_request(url, headers)
            domain_data = response.json()
            if not isinstance(domain_data, list) or any(
                not isinstance(domain, dict) or
                not all(domain.get(key) is not None and domain.get(key) != ''
                        for key in ('type', 'description', 'start', 'end'))
                for domain in domain_data
            ):
                raise ValueError('Unexpected protein-domain response from Ensembl')
            break
        except (requests.RequestException, ValueError) as exc:
            status = getattr(getattr(exc, 'response', None), 'status_code', None)
            retryable = status is None or status in (429, 500, 502, 503, 504)
            if attempt == 2 or not retryable:
                raise ValueError(
                    f'Protein domains could not be retrieved for {ensembl_protein_id}. '
                    'No report was generated. Please try again later.'
                ) from exc
            time.sleep(0.5 * (2 ** attempt))
    if not domain_data:
        raise ValueError(
            f'No protein-domain data returned for {ensembl_protein_id}. '
            'No report was generated.'
        )
    domains = []
    for domain in domain_data:
        domains.append({
            'Source': domain['type'],
            'Description': domain['description'],
            'Start': str(domain['start']),
            'End': str(domain['end'])
        })
    return domains 

def get_clinvar(root):
    """Read classifications from the same variation record used for mapping."""
    rows = []
    for rcv in root.findall('.//ClassifiedRecord/RCVList/RCVAccession'):
        conditions = [node.text for node in rcv.findall('./ClassifiedConditionList/ClassifiedCondition')]
        classifications = rcv.find('RCVClassifications')
        if not conditions or not all(conditions) or classifications is None:
            raise ValueError('ClinVar condition or classification data is missing.')
        for classification in classifications:
            description = classification.findtext('Description')
            review = classification.findtext('ReviewStatus')
            if not description or not review:
                raise ValueError('ClinVar classification or review status is missing.')
            rows.append({'Variant classification': description,
                         'Condition': '; '.join(conditions),
                         'Variant more info': f'{classification.tag}: {review}'})
    if not rows:
        raise ValueError('ClinVar classification data is missing.')
    return rows


def get_results_dict(hgvs_cdna_transcript_id):
    '''Gets the results dictionary for the specified HGVS cDNA transcript ID. Also checks that the HGVS ID is valid
    and throws an error if it is not.'''
    report_date = get_report_date()
    variant = resolve_variant(hgvs_cdna_transcript_id)
    current_version_hgvs_cdna_transcript_id = variant['transcript']
    gene_symbol = variant['gene']
    rsID = variant['rsid']
    protein_change = variant['protein']
    cdna_change = variant['cdna']
    hgvs_cdna_transcript_id_formatted_2 = f'{cdna_change} {protein_change}'
    full_gene_name, ensembl_gene_id = get_full_gene_name_and_ensembl_gene_id(gene_symbol)
    chromosome, start, end = get_gene_start_end_chromosome(gene_symbol)
    cytogenetic_band = get_cytogenetic_band(chromosome, start, end)
    high_protein_expression = get_high_protein_expression(ensembl_gene_id)
    if chromosome != 'chr' + variant['chromosome']:
        raise ValueError('Variant and gene chromosomes do not match.')
    grch38_variant_position = variant['position']
    ensembl_transcript_id = get_ensembl_transcript_id(current_version_hgvs_cdna_transcript_id, gene_symbol)
    ensembl_protein_id, transcript_length, translation_length, total_exons, coding_exons, exon_number =  get_transcript_details_and_ensembl_protein_id(ensembl_transcript_id, grch38_variant_position, variant['end'], variant['chromosome'])
    domains = get_domains(ensembl_protein_id)
    cleaned_clinvar = variant['classifications']
        
    source_info = (
    f'https://www.proteinatlas.org/{ensembl_gene_id}-{gene_symbol}/tissue\n'
    f'https://www.ncbi.nlm.nih.gov/snp/?term={rsID}'
    )
    if len(cleaned_clinvar) > 0:
        source_info += f"\nhttps://www.ncbi.nlm.nih.gov/clinvar/variation/{variant['clinvar_id']}/"
    
    results_dict = {
        'Gene symbol': gene_symbol,
        'Report generated on': report_date,
        'Coding change': hgvs_cdna_transcript_id_formatted_2,
        'Full gene name': full_gene_name,
        'Cytogenetic band': cytogenetic_band,
        'High protein expression': high_protein_expression,
        'Original variant': hgvs_cdna_transcript_id,
        'Current HGVS ID': current_version_hgvs_cdna_transcript_id,
        'cDNA change' : cdna_change,
        'Protein change': protein_change,
        'Transcript length' : str(transcript_length) + ' base pairs',
        'Translation length': str(translation_length) + ' residues',
        'Total exons' : total_exons,
        'Coding exons': coding_exons,
        'Exon number' : exon_number,
        'Protein domains': domains,
        'Clinvar': cleaned_clinvar,
        'Sources': source_info
    }
    
    for field, value in results_dict.items():
        if value is None or value == '' or value == []:
            raise ValueError(f'Required report field is missing: {field}.')
    return results_dict

def curate_data_for_pdf(results_dict):
    """Curates the data for the PDF. Formats the data for the protein domains and Clinvar sections."""
    styles = getSampleStyleSheet()
    curated_data = []

    for key, value in results_dict.items():
        if key == 'Gene symbol':
            curated_data.append(('title', {'text': f'Variant Report: {results_dict["Full gene name"].title()} ({results_dict["Gene symbol"]})', 'style': 'Title'}))

        elif key == 'Protein domains':
            table_data = [['Protein Domains', '', '', '']]  # Add title row spanning all columns
            table_data.append(['Source', 'Description', 'Start', 'End'])
            for entry in value:
                table_data.append([
                    Paragraph(entry.get('Source', ''), styles['Normal']),
                    Paragraph(entry.get('Description', ''), styles['Normal']), 
                    Paragraph(entry.get('Start', ''), styles['Normal']),
                    Paragraph(entry.get('End', ''), styles['Normal']),
                ])
            curated_data.append(('Protein domains', {'data': table_data, 'columns': [100, 100, 100, 100]}))

        elif key == 'Clinvar':
            table_data = [['ClinVar', '', '']]  # Add title row spanning all columns
            table_data.append(['Classification', 'Condition', 'More info'])
            for entry in value:
                combined_info = Paragraph(
                    f"{entry.get('Condition', '')}<br/>{entry.get('Affected status', '')}<br/>{entry.get('Allele origin', '')}",
                    styles['Normal']
                )
                table_data.append([
                    Paragraph(entry.get('Variant classification', ''), styles['Normal']),
                    combined_info,
                    Paragraph(entry.get('Variant more info', ''), styles['Normal'])
                ])
            curated_data.append(('clinvar', {'data': table_data, 'columns': [125, 175, 200]}))

        else:
            curated_data.append((key, value))

    return curated_data

def append_story_elements(curated_data, order):
    '''Appends the story elements to the PDF.'''
    styles = getSampleStyleSheet()
    story = []

    domain_colors = [
        colors.HexColor('#60d5ca'),  # light sea green
        colors.HexColor('#c1d7e3'),  # light blue
        colors.HexColor('#fde3a3'),  # light orange
        colors.HexColor('#c7e9b3'),  # light green
        colors.HexColor('#fbc6c6'),  # light red
        colors.HexColor('#dfc7e5'),  # light purple
        colors.HexColor('#d9ecd9'),  # light sky blue
        colors.HexColor('#ffd6b3'),  # light peach
        colors.HexColor('#cce5ff'),  # light lavender blue
        colors.HexColor('#ffd699'),  # light sand
        colors.HexColor('#e5f2ff'),  # light powder blue
        colors.HexColor('#ffd6d6'),  # light pink
        colors.HexColor('#d6ffb3'),  # light lime
        colors.HexColor('#ffffe5'),  # light lemon
        colors.HexColor('#e5ffe5'),  # light mint
        colors.HexColor('#ffebd6'),  # light apricot
        colors.HexColor('#e5d6ff'),  # light mauve
        colors.HexColor('#d6ffff'),  # light aqua
        colors.HexColor('#ffebd6')   # light coral
    ]
    
    unique_descriptions = {}

    for key in order:
        for item in curated_data:
            if item[0] == key:
                if key == 'title':
                    title = Paragraph(item[1]['text'], styles[item[1]['style']])
                    story.append(title)
                    story.append(Spacer(1, 12))

                elif key == 'Protein domains':
                    table_data = item[1]['data']
                    
                    table_style = TableStyle([
                        ('SPAN', (0, 0), (-1, 0)),  
                        ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
                        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('FONTSIZE', (0, 0), (-1, 0), 14),
                        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                        ('BACKGROUND', (0, 1), (-1, 1), colors.grey),
                        ('TEXTCOLOR', (0, 1), (-1, 1), colors.whitesmoke),
                        ('FONTNAME', (0, 1), (-1, 1), 'Helvetica-Bold'),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black),
                        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                        ('WORDWRAP', (0, 0), (-1, -1)),
                        ('splitLongWords', (0, 0), (-1, -1), False)
                    ])

                    # Apply different colors to the domain rows based on unique descriptions
                    for i in range(2, len(table_data)):
                        description = table_data[i][1].getPlainText()

                        if description not in unique_descriptions:
                            color_index = len(unique_descriptions) % len(domain_colors)
                            unique_descriptions[description] = domain_colors[color_index]

                        color = unique_descriptions[description]
                        table_style.add('BACKGROUND', (0, i), (-1, i), color)

                    table = Table(table_data, colWidths=item[1]['columns'])
                    table.setStyle(table_style)
                    story.append(table)
                    story.append(Spacer(1, 12))

                elif key == 'clinvar':
                    table_data = item[1]['data']
                    
                    table_style = TableStyle([
                        ('SPAN', (0, 0), (-1, 0)),  
                        ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
                        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('FONTSIZE', (0, 0), (-1, 0), 14),
                        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                        ('BACKGROUND', (0, 1), (-1, 1), colors.grey),
                        ('TEXTCOLOR', (0, 1), (-1, 1), colors.whitesmoke),
                        ('ALIGN', (0, 1), (-1, -1), 'CENTER'),
                        ('FONTNAME', (0, 1), (-1, 1), 'Helvetica-Bold'),
                        ('BOTTOMPADDING', (0, 1), (-1, 1), 12),
                        ('BACKGROUND', (0, 2), (-1, -1), colors.beige),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black),
                        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                        ('WORDWRAP', (0, 0), (-1, -1)),
                        ('splitLongWords', (0, 0), (-1, -1), False)
                    ])
                    table = Table(table_data, colWidths=item[1]['columns'])
                    table.setStyle(table_style)
                    story.append(table)
                    story.append(Spacer(1, 12))

                else:
                    text = f'<b>{key}</b>: {item[1]}'
                    paragraph = Paragraph(text, styles['Normal'])
                    story.append(paragraph)
                    story.append(Spacer(1, 12))

    return story


def create_pdf_from_curated_data(story, filename):
    """Creates a PDF from the curated data and saves it to the specified filename. Adds page breaks if necessary."""
    doc = SimpleDocTemplate(filename, pagesize=letter)
    elements = []
    styles = getSampleStyleSheet()

    available_height = letter[1] - doc.topMargin - doc.bottomMargin
    current_height = 0

    for element in story:
        if isinstance(element, Table):
            table_width, table_height = element.wrap(doc.width, doc.height)
            
            if current_height + table_height > available_height:
                elements.append(PageBreak())
                current_height = 0

            elements.append(element)
            current_height += table_height
        
        elif isinstance(element, (Spacer, Paragraph)):
            element_width, element_height = element.wrap(doc.width, doc.height)
            
            if current_height + element_height > available_height:
                elements.append(PageBreak())
                current_height = 0

            elements.append(element)
            current_height += element_height

    if elements and isinstance(elements[-2], PageBreak):
        elements.pop(-2)
    
    doc.build(elements)


app = Flask(__name__)

html_template = """
<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
    <title>Variant Report Generator</title>
</head>
<body>
    <div class="container mt-5">
        <h1 class="text-center">Variant Report Generator</h1>
        {% if error %}<div class="alert alert-danger" role="alert">{{ error }} No report was generated.</div>{% endif %}
        <form method="post">
            <div class="form-group">
                <label for="hgvs_cdna_transcript_id">HGVS cDNA Transcript ID</label>
                <input type="text" class="form-control" id="hgvs_cdna_transcript_id" name="hgvs_cdna_transcript_id" placeholder="Enter HGVS cDNA Transcript ID" required>
            </div>
            <button type="submit" class="btn btn-primary">Generate Report</button>
        </form>
    </div>
</body>
</html>
"""

@app.route('/', methods=['GET', 'POST'])

def index():    
    if request.method == 'POST':
        hgvs_cdna_transcript_id = request.form['hgvs_cdna_transcript_id']
        try:
            results_dict = get_results_dict(hgvs_cdna_transcript_id)
        except (ValueError, KeyError, TypeError, ET.ParseError) as exc:
            return render_template_string(html_template, error=f'Variant validation failed: {exc}'), 422
        except requests.RequestException:
            return render_template_string(html_template, error='A required data service is unavailable. Please try again later.'), 503
        curated_data = curate_data_for_pdf(results_dict)
        ordered_keys = [
            'title',
            'Report generated on',
            'Original variant',
            'Current HGVS ID',
            'cDNA change',
            'Protein change',
            'Cytogenetic band',
            'High protein expression',
            'Transcript length',
            'Translation length',
            'Total exons',
            'Coding exons',
            'Exon number',
            'Protein domains',
            'clinvar',
            'Sources'
        ]
        story = append_story_elements(curated_data, ordered_keys)
        
        buffer = BytesIO()
        create_pdf_from_curated_data(story, buffer)
        buffer.seek(0)
        
        gene_symbol = results_dict['Gene symbol']
        return send_file(buffer, as_attachment=True, download_name=f'{gene_symbol}.pdf', mimetype='application/pdf')
    
    return render_template_string(html_template)

app.debug = True

def run_app():
    app.run()

if __name__ == '__main__':
    run_app()
