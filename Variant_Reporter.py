#NM_005228.3:c.2648T>C
import datetime
import requests
import re
from bs4 import BeautifulSoup
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, PageBreak, Spacer, Paragraph, Frame
from reportlab.lib import colors
from flask import Flask, request, render_template_string, send_file, abort
from xml.etree import ElementTree as ET
from io import BytesIO

def send_request(url, headers=None):
    '''Sends a request to the specified URL and return the response.'''
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response

def get_report_date():
    '''Gets the current date in the format 'Month Day, Year'.'''
    report_date = datetime.date.today().strftime('%B %d, %Y')
    return report_date
    
def format_hgvs_cdna_transcript_id_1(hgvs_cdna_transcript_id):
    '''Formats the HGVS cDNA transcript ID. Example: NM_005228.3:c.2648T>C -> NM_005228'''
    hgvs_cdna_transcript_id_formatted_1 = hgvs_cdna_transcript_id.split(':')[0].split('.')[0]
    return hgvs_cdna_transcript_id_formatted_1

def format_hgvs_cdna_transcript_id_2(hgvs_cdna_transcript_id):
    '''Formats the HGVS cDNA transcript ID. Example: NM_005228.3:c.2648T>C -> c.2648T>C'''
    hgvs_cdna_transcript_id_formatted_2 = hgvs_cdna_transcript_id.split(':')[1]
    return hgvs_cdna_transcript_id_formatted_2

def get_current_version_hgvs_cdna_transcript_id(hgvs_cdna_transcript_id_formatted_1):
    '''Gets the current version of the HGVS cDNA transcript ID. Example: NM_005228.3 -> NM_005228.5'''
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={hgvs_cdna_transcript_id_formatted_1}&rettype=gb&retmode=xml'
    response = send_request(url)
    root = ET.fromstring(response.content)
    for item in root.findall('.//GBSeq'):
        current_version_hgvs_cdna_transcript_id = item.find('GBSeq_accession-version').text
        return current_version_hgvs_cdna_transcript_id

def get_full_current_version_hgvs_cdna_transcript_id(current_version_hgvs_cdna_transcript_id, hgvs_cdna_transcript_id_formatted_2):
    '''Gets the full, current version of the HGVS cDNA transcript ID. Example: NM_005228.3:c.2648T>C -> NM_005228.5:c.2648T>C'''
    full_current_version_hgvs_cdna_transcript_id = current_version_hgvs_cdna_transcript_id + ':' + hgvs_cdna_transcript_id_formatted_2
    return full_current_version_hgvs_cdna_transcript_id

def validate_full_current_version_hgvs_cdna_transcript_id(full_current_version_hgvs_cdna_transcript_id):
    '''Checks if the specified HGVS ID is valid.'''
    url = f'https://rest.ensembl.org/vep/human/hgvs/{full_current_version_hgvs_cdna_transcript_id}?'
    headers = {"Content-Type": "application/json"}
    response = send_request(url, headers)
    
def get_ensembl_transcript_id(hgvs_cdna_transcript_id_formatted_1):
    '''Gets the Ensembl transcript ID for the specified HGVS cDNA transcript ID.'''
    url = f'https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{hgvs_cdna_transcript_id_formatted_1}?external_db=RefSeq_mRNA'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers)
    data = response.json()
    for entry in data:
        if entry['type'] == 'transcript':
            ensembl_transcript_id = entry['id']
            return ensembl_transcript_id

def get_gene_symbol(ensembl_transcript_id):
    '''Gets the gene symbol for the specified ensembl transcript ID.'''
    url = f'https://rest.ensembl.org/lookup/id/{ensembl_transcript_id}'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers)
    data = response.json()
    gene_symbol = data['display_name'].split('-')[0]
    return gene_symbol
    
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
    for data in root.findall('.//data'):
        tissue = data.find('tissue')
        levels = data.findall('level[@type="expression"]')
        if any(level.text.lower() == "high" for level in levels):
            high_protein_expression.append(tissue.text)
    if high_protein_expression:
        return(', '.join(high_protein_expression).lower())
    return('Protein not highly expressed.')

def get_amino_acid_change(current_version_hgvs_cdna_transcript_id, hgvs_cdna_transcript_id_split, ensembl_transcript_id):
    '''Gets the protein start, protein end, and amino acid change for the specified current version HGVS cDNA transcript ID,
    split HGVS cDNA transcript ID, and Ensembl transcript ID.'''
    url = f'https://rest.ensembl.org/vep/human/hgvs/{current_version_hgvs_cdna_transcript_id}:{hgvs_cdna_transcript_id_split}'
    headers = {"Content-Type": "application/json"}
    response = send_request(url, headers=headers)
    data = response.json()
    for consequence in data[0]['transcript_consequences']:
        if consequence['transcript_id'] == ensembl_transcript_id:
            variant_codon_start = consequence.get('protein_start', 'Not specified')
            variant_codon_end = consequence.get('protein_end', 'Not specified')  
            amino_acid_change = consequence.get('amino_acids', 'Not specified')
    return variant_codon_start, variant_codon_end, amino_acid_change

def get_clinvar_accession_id(full_current_version_hgvs_cdna_transcript_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={full_current_version_hgvs_cdna_transcript_id}&retmode=json"
    response = send_request(url)
    data = response.json()
    if 'esearchresult' in data and 'idlist' in data['esearchresult']:
        ids = data['esearchresult']['idlist']
        clinvar_accession_id = ids[0]
    return clinvar_accession_id

def get_rsID(clinvar_accession_id):
    '''Scrapes the rsID from the ClinVar Variation page.'''
    url = f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_accession_id}/'
    response = send_request(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    rsid_link = soup.find('a', href=lambda href: href and '/snp/rs' in href)
    rsID = rsid_link.text.strip()
    rsID = rsID.replace("dbSNP:", "")
    rsID = rsID.strip()
    return rsID
    
def get_grch38_variant_position(rsID):
    '''Gets the variant position (in GRCh38 currently) for the specified current version HGVS cDNA transcript ID's rsID.'''
    url = f"https://rest.ensembl.org/variation/human/{rsID}?content-type=application/json"
    response = send_request(url)
    data = response.json()
    # Check for mappings to GRCh38
    for mapping in data.get("mappings", []):
        if mapping.get("assembly_name") == "GRCh38":
            grch38_variant_position = mapping.get("start")
            return grch38_variant_position
    return None

def get_transcript_details_and_ensembl_protein_id(ensembl_transcript_id, grch38_variant_position):
    '''Gets the Ensembl protein ID, transcript length, translation length, total exons, coding exons,
    exon number, and coding exons for the specified Ensembl transcript ID and variant position.'''
    
    url = f'https://rest.ensembl.org/lookup/id/{ensembl_transcript_id}?expand=1'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers)
    data = response.json()
    
    # Extract basic information
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

    # Determine the exon number containing the variant position (GRCh38 only)
    exon_number = next((index + 1 for index, exon in enumerate(data['Exon'])
                        if exon.get('assembly_name') == 'GRCh38' and 
                        int(exon['start']) <= int(grch38_variant_position) <= int(exon['end'])), None)

    return ensembl_protein_id, transcript_length, translation_length, total_exons, coding_exons, exon_number


def get_pfam_smart_protein_domains(ensembl_protein_id):
    '''Gets the protein domains from the sources Pfam and Smart for the specified Ensembl protein ID.'''
    url = f'https://rest.ensembl.org/overlap/translation/{ensembl_protein_id}'
    headers = {'Content-Type': 'application/json'}
    response = send_request(url, headers)
    domain_data = response.json()
    domains = []
    for domain in domain_data:
        if domain.get('type') in ['Pfam', 'Smart']:
            domains.append({
                'Source': domain['type'],
                'Description': domain.get('description'),
                'Start': str(domain['start']),
                'End': str(domain['end'])
            })
    return domains 

def get_clinvar(rsID):
    '''Gets the variant classification, variant condition, and variant more info for the specified rsID. There may or may not be
    one or more Clinvar entries for the specified rsID. Also formats the output for the table in the pdf.'''
    url = f'https://www.ncbi.nlm.nih.gov/clinvar/?term={rsID}'
    response = send_request(url)
    clinvar = []
    soup = BeautifulSoup(response.text, 'html.parser')
    tables = soup.find_all('table')
    if len(tables) > 4:
        correct_table = tables[4]
        rows = correct_table.find_all('tr')
        for row in rows:
            cells = row.find_all('td')
            if len(cells) >= 5:
                classification_info = cells[0].text.strip()
                condition_info = cells[2].text.strip()
                more_info = cells[4].text.strip()
                if '(more)' in more_info:
                    more_info_parts = more_info.split("(more)")
                    if len(more_info_parts) > 1:
                        more_info = more_info_parts[1].strip()
                if '(less)' in more_info:
                    more_info = more_info.replace('(less)', '').strip()
                clinvar.append({
                    'Variant classification': classification_info,
                    'Variant condition': condition_info,
                    'Variant more info': more_info
                })
                    
    cleaned_clinvar = []
    for entry in clinvar:
        cleaned_entry = {}
        for key, value in entry.items():
            if key == 'Variant classification' or key == 'Variant more info':
                cleaned_value = ' '.join(value.split())
                cleaned_entry[key] = cleaned_value
            elif key == 'Variant condition':
                parts = [part.strip() for part in value.split('\n') if part.strip()]
                condition_dict = {}
                if len(parts) > 0:
                    condition_dict['Condition'] = 'Condition: ' + str(parts[0])
                if len(parts) > 1:
                    condition_dict['Affected status'] = parts[1]
                if len(parts) > 2:
                    condition_dict['Allele origin'] = parts[2]
                cleaned_entry.update(condition_dict)
            else:
                cleaned_entry[key] = value
        cleaned_clinvar.append(cleaned_entry)
    return cleaned_clinvar  
    
def get_results_dict(hgvs_cdna_transcript_id):
    '''Gets the results dictionary for the specified HGVS cDNA transcript ID. Also checks that the HGVS ID is valid
    and throws an error if it is not.'''
    report_date = get_report_date()
    hgvs_cdna_transcript_id_formatted_1 = format_hgvs_cdna_transcript_id_1(hgvs_cdna_transcript_id)
    hgvs_cdna_transcript_id_formatted_2 = format_hgvs_cdna_transcript_id_2(hgvs_cdna_transcript_id)
    current_version_hgvs_cdna_transcript_id = get_current_version_hgvs_cdna_transcript_id(hgvs_cdna_transcript_id_formatted_1)
    full_current_version_hgvs_cdna_transcript_id = get_full_current_version_hgvs_cdna_transcript_id(current_version_hgvs_cdna_transcript_id, hgvs_cdna_transcript_id_formatted_2)
    validate_full_current_version_hgvs_cdna_transcript_id(full_current_version_hgvs_cdna_transcript_id)
    ensembl_transcript_id = get_ensembl_transcript_id(hgvs_cdna_transcript_id_formatted_1)
    gene_symbol = get_gene_symbol(ensembl_transcript_id)
    full_gene_name, ensembl_gene_id = get_full_gene_name_and_ensembl_gene_id(gene_symbol)
    start, end, chromosome = get_gene_start_end_chromosome(gene_symbol)
    cytogenetic_band = get_cytogenetic_band(start, end, chromosome)
    high_protein_expression = get_high_protein_expression(ensembl_gene_id)
    variant_codon_start, variant_codon_end, amino_acid_change = get_amino_acid_change(current_version_hgvs_cdna_transcript_id, hgvs_cdna_transcript_id_formatted_2, ensembl_transcript_id)
    clinvar_accession_id = get_clinvar_accession_id(full_current_version_hgvs_cdna_transcript_id)
    rsID = get_rsID(clinvar_accession_id)
    grch38_variant_position = get_grch38_variant_position(rsID)
    ensembl_protein_id, transcript_length, translation_length, total_exons, coding_exons, exon_number = get_transcript_details_and_ensembl_protein_id(ensembl_transcript_id, grch38_variant_position)
    pfam_smart_protein_domains = get_pfam_smart_protein_domains(ensembl_protein_id)
    cleaned_clinvar = get_clinvar(rsID)
    source_info = (
    f'https://www.proteinatlas.org/{ensembl_gene_id}-{gene_symbol}/tissue\n'
    f'https://www.ncbi.nlm.nih.gov/snp/?term={rsID}'
    )
    if len(cleaned_clinvar) > 0:
        source_info += f'\nhttps://www.ncbi.nlm.nih.gov/clinvar/?term={rsID}'
    
    results_dict = {
        'Gene symbol': gene_symbol,
        'Report generated on': report_date,
        'Coding change': hgvs_cdna_transcript_id_formatted_2,
        'Full gene name': full_gene_name,
        'Cytogenetic band': cytogenetic_band,
        'High protein expression': high_protein_expression,
        'Current HGVS ID': current_version_hgvs_cdna_transcript_id,
        'Variant codon start': variant_codon_start,
        'Variant codon end' : variant_codon_end,
        'Amino acid change' : amino_acid_change,
        'Transcript length' : str(transcript_length) + ' base pairs',
        'Translation length': str(translation_length) + ' residues',
        'Total exons' : total_exons,
        'Coding exons': coding_exons,
        'Exon number' : exon_number,
        'Protein domains (Sources: Pfam and Smart)': pfam_smart_protein_domains,
        'Clinvar': cleaned_clinvar,
        'Sources': source_info
    }
    
    return results_dict

def curate_data_for_pdf(results_dict):
    """Curates the data for the PDF. Formats the data for the protein domains and Clinvar sections."""
    styles = getSampleStyleSheet()
    curated_data = []

    for key, value in results_dict.items():
        if key == 'Gene symbol':
            curated_data.append(('title', {'text': f'Variant Report: {results_dict["Full gene name"].title()} ({results_dict["Gene symbol"]})', 'style': 'Title'}))

        elif key == 'Protein domains (Sources: Pfam and Smart)':
            table_data = [['Protein Domains', '', '', '']]  # Add title row spanning all columns
            table_data.append(['Source', 'Description', 'Start', 'End'])
            for entry in value:
                table_data.append([
                    Paragraph(entry.get('Source', ''), styles['Normal']),
                    Paragraph(entry.get('Description', ''), styles['Normal']), 
                    Paragraph(entry.get('Start', ''), styles['Normal']),
                    Paragraph(entry.get('End', ''), styles['Normal']),
                ])
            curated_data.append(('protein_domains', {'data': table_data, 'columns': [100, 100, 100, 100]}))

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

                elif key == 'protein_domains':
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
        results_dict = get_results_dict(hgvs_cdna_transcript_id)
        curated_data = curate_data_for_pdf(results_dict)
        ordered_keys = [
            'title',
            'Report generated on',
            'Current HGVS ID',
            'Coding change',
            'Amino acid change',
            'Variant codon start',
            'Variant codon end',
            'Cytogenetic band',
            'High protein expression',
            'Transcript length',
            'Translation length',
            'Total exons',
            'Coding exons',
            'Exon number',
            'protein_domains',
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

run_app()
