import os
import logging
from flask import Flask, request, render_template, jsonify
from html import escape
import pandas as pd
from Bio.Align import PairwiseAligner

app = Flask(__name__)

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def load_gene_data(file_path):
    logger.debug(f"Loading gene data from {file_path}")
    df = pd.read_csv(file_path)
    logger.debug(f"CSV columns: {df.columns.tolist()}")
    # Clean 'Protein Sequence' by removing FASTA headers
    df['Protein Sequence'] = df['Protein Sequence'].apply(
        lambda x: x.split('\n')[-1].strip() if isinstance(x, str) and '>' in x else x.strip()
    )
    return df

def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    score = aligner.score(seq1, seq2)
    similarity = (score / max(len(seq1), len(seq2))) * 100
    return similarity

def analyze_sequence(protein_seq, gene_data, threshold=75):
    results = []
    for _, row in gene_data.iterrows():
        gene_name = row['Name']  # Use 'Name' column from PGPRgene.csv
        sequence = row['Protein Sequence'].strip()
        if not sequence or not isinstance(sequence, str):
            logger.warning(f"Invalid sequence for gene {gene_name}: {sequence}")
            continue
        similarity = calculate_similarity(protein_seq, sequence)
        results.append({
            "Gene Name": gene_name,
            "Similarity (%)": similarity,
            "Is PGPR": similarity >= threshold
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
            protein_seq = escape(request.form['protein_seq'].strip())
            logger.debug(f"Received sequence: {protein_seq[:50]}...")
            if not protein_seq:
                return render_template('analyzer.html', error="Please provide a protein sequence.")
            
            gene_data = load_gene_data(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))
            if 'Name' not in gene_data.columns or 'Protein Sequence' not in gene_data.columns:
                return render_template('analyzer.html', error="Invalid CSV format: Missing 'Name' or 'Protein Sequence' columns")
            
            results = analyze_sequence(protein_seq, gene_data)
            
            chart_data = {
                'labels': [result['Gene Name'] for result in results],
                'data': [result['Similarity (%)'] for result in results]
            }
            
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