import os
import logging
from flask import Flask, request, render_template, send_file
from html import escape
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from collections import Counter
import csv
import re

app = Flask(__name__)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def load_gene_data(file_path):
    df = pd.read_csv(file_path)
    df = df.dropna(subset=['Name', 'Protein Sequence'])
    return df

def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    score = aligner.score(seq1, seq2)
    similarity = (score / max(len(seq1), len(seq2))) * 100
    return similarity

def align_sequences_to_query(query, sequences):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligned_seqs = []
    max_len = 0

    for seq in sequences:
        alignment = aligner.align(query, seq)[0]
        aligned = alignment.aligned[1]
        aligned_seq = list("-" * len(query))
        for (start, end) in aligned:
            aligned_seq[start:end] = seq[start:end]
        aligned_str = "".join(aligned_seq)
        aligned_seqs.append(aligned_str)
        max_len = max(max_len, len(aligned_str))

    aligned_seqs = [s.ljust(max_len, '-') for s in aligned_seqs]
    return aligned_seqs

def compute_consensus_from_aligned(aligned_seqs):
    if not aligned_seqs:
        return ""
    max_len = max(len(seq) for seq in aligned_seqs)
    consensus = []
    for i in range(max_len):
        column = [seq[i] for seq in aligned_seqs if i < len(seq)]
        most_common = Counter(column).most_common(1)[0][0]
        consensus.append(most_common)
    return ''.join(consensus)

def analyze_sequence(protein_seq, gene_data, threshold=75):
    results = []
    pgpr_sequences = []

    for _, row in gene_data.iterrows():
        gene_name = row['Name']
        bacteria_name = row.get('Bacteria', 'Unknown')
        description = row.get('Description', 'No description available')
        sequence = row['Protein Sequence']

        if not isinstance(sequence, str) or not sequence.strip():
            continue

        if sequence.startswith(">"):
            sequence = ''.join(line.strip() for line in sequence.split('\n') if not line.startswith(">"))

        similarity = calculate_similarity(protein_seq, sequence)
        is_pgpr = similarity >= threshold

        if is_pgpr:
            pgpr_sequences.append(sequence)

        results.append({
            "Gene Name": gene_name,
            "Bacteria Name": bacteria_name,
            "Description": description,
            "Similarity (%)": similarity,
            "Is PGPR": is_pgpr
        })

    pgpr_aligned = align_sequences_to_query(protein_seq, pgpr_sequences) if pgpr_sequences else []
    consensus_sequence = compute_consensus_from_aligned(pgpr_aligned) if pgpr_aligned else None
    top_match = max(results, key=lambda x: x['Similarity (%)']) if results else None
    top_results = sorted(results, key=lambda x: x['Similarity (%)'], reverse=True)[:10]

    return top_results, top_match, consensus_sequence, pgpr_aligned

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    try:
        if request.method == 'POST':
            protein_seq = escape(request.form['protein_seq'].strip())
            if not protein_seq:
                return render_template('analyzer.html', error="Please provide a protein sequence.")

            gene_data = load_gene_data(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))
            if 'Name' not in gene_data.columns or 'Protein Sequence' not in gene_data.columns:
                return render_template('analyzer.html', error="CSV missing required columns.")

            results, max_gene, consensus_seq, pgpr_aligned = analyze_sequence(protein_seq, gene_data)

            chart_data = {
                'labels': [result['Gene Name'] for result in results],
                'data': [result['Similarity (%)'] for result in results]
            }

            return render_template(
                'analyzer.html',
                results=results,
                chart_data=chart_data,
                max_gene=max_gene,
                consensus_seq=consensus_seq,
                pgpr_aligned=pgpr_aligned
            )

        return render_template('analyzer.html')

    except Exception as e:
        logger.error(f"Error in analyzer: {str(e)}", exc_info=True)
        return render_template('analyzer.html', error=f"An error occurred: {str(e)}")

@app.route('/add_gene', methods=['GET', 'POST'])
def add_gene():
    try:
        if request.method == 'POST':
            gene_name = escape(request.form['gene_name'].strip())
            gene_id = escape(request.form['gene_id'].strip())
            bacteria_name = escape(request.form['bacteria_name'].strip())
            description = escape(request.form['description'].strip())
            protein_seq = escape(request.form['protein_seq'].strip())

            # Validate inputs
            if not all([gene_name, gene_id, bacteria_name, description, protein_seq]):
                return render_template('add_gene.html', error="All fields are required.")

            # Path to CSV file
            csv_path = os.path.join(os.path.dirname(__file__), 'PGPRgene.csv')

            # Read current CSV to determine next S. No.
            gene_data = load_gene_data(csv_path)
            next_sno = len(gene_data) + 1

            # Append new gene to CSV
            with open(csv_path, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([next_sno, gene_name, gene_id, description, bacteria_name, protein_seq])

            logger.info(f"Added gene: {gene_name} (ID: {gene_id}) to PGPRgene.csv")
            return render_template('add_gene.html', success="Gene added successfully!")

        return render_template('add_gene.html')

    except Exception as e:
        logger.error(f"Error in add_gene: {str(e)}", exc_info=True)
        return render_template('add_gene.html', error=f"Failed to add gene: {str(e)}")

@app.route('/download_csv')
def download_csv():
    try:
        csv_path = os.path.join(os.path.dirname(__file__), 'PGPRgene.csv')
        if not os.path.exists(csv_path):
            logger.error("CSV file not found")
            return render_template('index.html', error="CSV file not found.")
        
        logger.info("Serving PGPRgene.csv for download")
        return send_file(csv_path, as_attachment=True, download_name='PGPRgene.csv')
    
    except Exception as e:
        logger.error(f"Error in download_csv: {str(e)}", exc_info=True)
        return render_template('index.html', error=f"Failed to download CSV: {str(e)}")

@app.route('/convert_sequence', methods=['GET', 'POST'])
def convert_sequence():
    try:
        if request.method == 'POST':
            dna_sequence = escape(request.form['dna_sequence'].strip().upper())

            # Validate DNA sequence
            if not dna_sequence:
                return render_template('convert_sequence.html', error="Please provide a DNA sequence.")

            # Check for valid DNA characters (A, T, G, C)
            if not re.match(r'^[ATGC]+$', dna_sequence):
                return render_template('convert_sequence.html', error="Invalid DNA sequence. Use only A, T, G, C.")

            # Check if sequence length is divisible by 3
            if len(dna_sequence) % 3 != 0:
                return render_template('convert_sequence.html', error="DNA sequence length must be divisible by 3.")

            # Convert to protein sequence
            try:
                seq = Seq(dna_sequence)
                protein_seq = str(seq.translate(to_stop=True))
                if not protein_seq:
                    return render_template('convert_sequence.html', error="Translation resulted in an empty protein sequence.")
                
                logger.info(f"Converted DNA sequence to protein: {protein_seq[:50]}...")
                return render_template('convert_sequence.html', 
                                     dna_sequence=dna_sequence, 
                                     protein_sequence=protein_seq,
                                     success="DNA sequence successfully converted to protein sequence!")
            except Exception as e:
                logger.error(f"Translation error: {str(e)}")
                return render_template('convert_sequence.html', error=f"Failed to translate sequence: {str(e)}")

        return render_template('convert_sequence.html')

    except Exception as e:
        logger.error(f"Error in convert_sequence: {str(e)}", exc_info=True)
        return render_template('convert_sequence.html', error=f"An error occurred: {str(e)}")

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/project')
def project():
    return render_template('project.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)