import os
import logging
import numpy as np
import pandas as pd
from flask import Flask, request, render_template, send_file, flash, redirect, url_for
from html import escape
from Bio.Align import PairwiseAligner
import re

app = Flask(__name__)
app.secret_key = 'your_secret_key'

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

CSV_PATH = os.path.join(os.path.dirname(__file__), 'PGPRgene.csv')
# For Render deployment
# CSV_PATH = '/opt/render/project/src/PGPRgene.csv'

def load_gene_data(file_path):
    logger.debug(f"Loading gene data from {file_path}")
    try:
        df = pd.read_csv(file_path)
        logger.debug(f"CSV columns: {df.columns.tolist()}")
        df['Protein Sequence'] = df['Protein Sequence'].apply(
            lambda x: clean_sequence(x) if isinstance(x, str) else str(x).strip().upper()
        )
        for col in ['Name', 'Gene-ID', 'Description', 'Bacteria', 'Protein Sequence']:
            df[col] = df[col].fillna('')
        for col in ['Protein Sequence']:
            empty_rows = df[df[col] == '']
            if not empty_rows.empty:
                logger.warning(f"Rows with empty '{col}': {empty_rows.index.tolist()}")
        return df
    except Exception as e:
        logger.error(f"Error loading CSV: {str(e)}")
        raise

def clean_sequence(seq):
    """Remove FASTA headers, non-amino acid characters, and convert to uppercase."""
    if not seq or not isinstance(seq, str):
        return ""
    seq = seq.split('\n')[-1].strip() if '>' in seq else seq.strip()
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
    return seq.upper()

def calculate_similarity(seq1, seq2):
    try:
        seq1 = clean_sequence(seq1)
        seq2 = clean_sequence(seq2)
        if not seq1 or not seq2 or len(seq1) == 0 or len(seq2) == 0:
            logger.warning(f"Empty or invalid sequences after cleaning: seq1={seq1[:20]}..., seq2={seq2[:20]}...")
            return 0.0
        logger.debug(f"Comparing seq1={seq1[:20]}... (len={len(seq1)}), seq2={seq2[:20]}... (len={len(seq2)})")
        if seq1 == seq2:
            logger.debug(f"Exact match found: 100% similarity")
            return 100.0
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        score = aligner.score(seq1, seq2)
        norm_length = min(len(seq1), len(seq2))
        similarity = (score / norm_length) * 100 if norm_length > 0 else 0.0
        logger.debug(f"Score={score}, norm_length={norm_length}, similarity={similarity:.2f}%")
        if np.isnan(similarity):
            logger.warning(f"NaN similarity for seq1={seq1[:20]}..., seq2={seq2[:20]}...")
            return 0.0
        similarity = max(0, min(100, similarity))
        if similarity >= 99:
            logger.debug(f"High similarity: {similarity:.2f}%")
        return similarity
    except Exception as e:
        logger.error(f"Error in similarity calculation: {str(e)}")
        return 0.0

def analyze_sequence(protein_seq, gene_data, threshold=75):
    cleaned_input = clean_sequence(protein_seq)
    logger.debug(f"Cleaned input sequence: {cleaned_input[:50]}... (len={len(cleaned_input)})")
    if not cleaned_input:
        logger.error("Input sequence is empty after cleaning")
        return []
    results = []
    for idx, row in gene_data.iterrows():
        s_no = row['S. No.']
        gene_name = row['Name'] if row['Name'] and isinstance(row['Name'], str) and row['Name'].strip() else row['Gene-ID'] or f"Unknown_{s_no}"
        sequence = row['Protein Sequence']
        gene_id = row['Gene-ID']
        bacteria = row['Bacteria']
        logger.debug(f"Processing S. No. {s_no}, Name={gene_name}, Gene-ID={gene_id}, Seq={sequence[:20]}... (len={len(sequence)})")
        if not sequence or not isinstance(sequence, str) or sequence.strip() == '':
            logger.warning(f"Skipping row with invalid sequence for S. No. {s_no}")
            continue
        similarity = calculate_similarity(cleaned_input, sequence)
        results.append({
            "Gene Name": gene_name,
            "Gene-ID": gene_id,
            "Bacteria": bacteria,
            "Similarity (%)": similarity,
            "Is PGPR": similarity >= threshold
        })
    results = sorted(results, key=lambda x: x['Similarity (%)'], reverse=True)[:10]
    logger.debug(f"Top 10 results: {[r['Gene Name'] for r in results]}")
    return results

def is_duplicate_sequence(new_sequence, gene_data):
    cleaned_new_seq = clean_sequence(new_sequence)
    for seq in gene_data['Protein Sequence']:
        if cleaned_new_seq == seq:
            return True
    return False

@app.route('/')
def home():
    logger.debug("Accessing home page")
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    try:
        logger.debug("Accessing analyzer page")
        if request.method == 'POST':
            protein_seq = escape(request.form['protein_seq'].strip())
            logger.debug(f"Raw input sequence: {protein_seq[:50]}...")
            if not protein_seq:
                return render_template('analyzer.html', error="Please provide a protein sequence.")
            
            gene_data = load_gene_data(CSV_PATH)
            if not all(col in gene_data.columns for col in ['S. No.', 'Name', 'Gene-ID', 'Bacteria', 'Protein Sequence']):
                return render_template('analyzer.html', error="Invalid CSV format: Missing required columns")
            
            results = analyze_sequence(protein_seq, gene_data)
            if not results:
                return render_template('analyzer.html', error="No valid matches found. Check sequence or CSV data.")
            
            chart_data = {
                'labels': [result['Gene Name'] for result in results],
                'data': [result['Similarity (%)'] for result in results]
            }
            
            return render_template('analyzer.html', results=results, chart_data=chart_data)
        
        return render_template('analyzer.html')
    except Exception as e:
        logger.error(f"Error in analyzer: {str(e)}", exc_info=True)
        return render_template('analyzer.html', error=f"An error occurred: {str(e)}")

@app.route('/add_gene', methods=['GET', 'POST'])
def add_gene():
    try:
        if request.method == 'POST':
            name = escape(request.form.get('name', '').strip())
            gene_id = escape(request.form.get('gene_id', '').strip())
            description = escape(request.form.get('description', '').strip())
            bacteria = escape(request.form.get('bacteria', '').strip())
            protein_seq = escape(request.form.get('protein_seq', '').strip())
            
            if not all([name, gene_id, bacteria, protein_seq]):
                flash("All fields except Description are required.", "error")
                return render_template('add_gene.html')
            
            gene_data = load_gene_data(CSV_PATH)
            
            if is_duplicate_sequence(protein_seq, gene_data):
                flash("This protein sequence already exists in the database.", "error")
                return render_template('add_gene.html')
            
            new_row = {
                'S. No.': gene_data['S. No.'].max() + 1 if not gene_data['S. No.'].empty else 1,
                'Name': name,
                'Gene-ID': gene_id,
                'Description': description,
                'Bacteria': bacteria,
                'Protein Sequence': protein_seq
            }
            
            new_df = pd.DataFrame([new_row])
            new_df.to_csv(CSV_PATH, mode='a', header=False, index=False)
            logger.debug(f"Added new gene: {name}")
            flash("Gene added successfully!", "success")
            return redirect(url_for('add_gene'))
        
        return render_template('add_gene.html')
    except Exception as e:
        logger.error(f"Error in add_gene: {str(e)}", exc_info=True)
        flash(f"An error occurred: {str(e)}", "error")
        return render_template('add_gene.html')

@app.route('/download')
def download():
    try:
        logger.debug(f"Serving CSV download: {CSV_PATH}")
        return send_file(CSV_PATH, as_attachment=True, download_name='PGPRgene.csv')
    except Exception as e:
        logger.error(f"Error in download: {str(e)}")
        flash(f"Error downloading file: {str(e)}", "error")
        return redirect(url_for('home'))

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
    app.run(host='0.0.0.0', port=port, debug=True)