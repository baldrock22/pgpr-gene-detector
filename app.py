import os
import logging
from flask import Flask, request, render_template
from html import escape
import pandas as pd
from Bio.Align import PairwiseAligner

app = Flask(__name__)

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Load gene data from CSV
def load_gene_data(file_path):
    logger.debug(f"Loading gene data from {file_path}")
    df = pd.read_csv(file_path)
    logger.debug(f"CSV columns: {df.columns.tolist()}")

    # Check for duplicates based on Protein Sequence
    duplicates = df[df.duplicated(subset=['Protein Sequence'], keep=False)]
    if not duplicates.empty:
        logger.warning(f"Duplicate sequences found at S. No.: {duplicates['S. No.'].tolist()}")

    # Clean 'Protein Sequence' by removing FASTA headers and stripping whitespace
    df['Protein Sequence'] = df['Protein Sequence'].apply(
        lambda x: x.split('\n')[-1].strip() if isinstance(x, str) and '>' in x else x.strip()
    )
    return df

# Calculate similarity using PairwiseAligner (as in app.py)
def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    score = aligner.score(seq1, seq2)
    norm_length = max(len(seq1), len(seq2))
    similarity = (score / norm_length) * 100
    return score, norm_length, similarity

# Analyze the input sequence against the gene data
def analyze_sequence(protein_seq, gene_data):
    results = []
    for _, row in gene_data.iterrows():
        gene_name = row['Name']
        gene_id = row['Gene-ID']
        sequence = row['Protein Sequence']
        s_no = row['S. No.']

        if not sequence or not isinstance(sequence, str):
            logger.warning(f"Invalid sequence for gene {gene_name} (Gene-ID={gene_id}): {sequence}")
            continue

        logger.debug(f"Processing S. No. {s_no}, Name={gene_name}, Gene-ID={gene_id}, Seq={sequence[:20]}... (len={len(sequence)})")
        logger.debug(f"Comparing seq1={protein_seq[:20]}... (len={len(protein_seq)}), seq2={sequence[:20]}... (len={len(sequence)})")

        score, norm_length, similarity = calculate_similarity(protein_seq, sequence)
        logger.debug(f"Score={score}, norm_length={norm_length}, similarity={similarity:.2f}%")

        if similarity == 100:
            logger.debug("Exact match found: 100% similarity")
        else:
            logger.debug(f"High similarity but not identical: {similarity:.2f}%")

        results.append({
            "Gene Name": gene_name,
            "Similarity (%)": similarity,
            "Gene-ID": gene_id,
            "S. No.": s_no
        })

    # Sort by similarity in descending order and take top 10
    results = sorted(results, key=lambda x: x['Similarity (%)'], reverse=True)[:10]
    return results

@app.route('/')
def home():
    logger.debug("Accessing home page")
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    try:
        if request.method == 'POST':
            protein_seq = escape(request.form['protein_seq'].strip().upper())
            logger.debug(f"Received sequence: {protein_seq[:50]}...")
            if not protein_seq:
                return render_template('analyzer.html', error="Please provide a protein sequence.")

            logger.debug(f"Cleaned input sequence: {protein_seq} (len={len(protein_seq)})")

            gene_data = load_gene_data(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))
            if 'Name' not in gene_data.columns or 'Protein Sequence' not in gene_data.columns:
                return render_template('analyzer.html', error="Invalid CSV format: Missing 'Name' or 'Protein Sequence' columns")

            results = analyze_sequence(protein_seq, gene_data)

            # Prepare chart data for the frontend
            chart_data = {
                'labels': [result['Gene Name'] for result in results],
                'data': [result['Similarity (%)'] for result in results]
            }

            # Log top 10 results
            names = [result['Gene Name'] for result in results]
            similarities = [result['Similarity (%)'] for result in results]
            logger.debug(f"Top 10 results: {names} with similarities: {similarities}")

            return render_template('analyzer.html', results=results, chart_data=chart_data)

        return render_template('analyzer.html')
    except Exception as e:
        logger.error(f"Error in analyzer: {str(e)}", exc_info=True)
        return render_template('analyzer.html', error=f"An error occurred: {str(e)}")

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
