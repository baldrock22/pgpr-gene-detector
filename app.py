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

def load_gene_data(file_path):
    logger.debug(f"Loading gene data from {file_path}")
    df = pd.read_csv(file_path)
    logger.debug(f"CSV columns: {df.columns.tolist()}")
    duplicates = df[df.duplicated(subset=['Protein Sequence'], keep=False)]
    if not duplicates.empty:
        logger.warning(f"Duplicate sequences found at S. No.: {duplicates['S. No.'].tolist()}")
    df['Protein Sequence'] = df['Protein Sequence'].apply(
        lambda x: x.split('\n')[-1].strip() if isinstance(x, str) and '>' in x else x.strip()
    )
    return df

def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    score = aligner.score(seq1, seq2)
    norm_length = max(len(seq1), len(seq2))
    similarity = (score / norm_length) * 100
    return score, norm_length, similarity

@app.route('/')
def home():
    logger.debug("Accessing home page")
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    try:
        if request.method == 'POST':
            protein_seq = escape(request.form['protein_seq'].strip().upper())
            logger.debug(f"Raw input sequence: {protein_seq}")
            if not protein_seq:
                return render_template('analyzer.html', error="Please provide a protein sequence.")
            
            gene_data = load_gene_data(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))
            if 'Name' not in gene_data.columns or 'Protein Sequence' not in gene_data.columns:
                return render_template('analyzer.html', error="Invalid CSV format: Missing 'Name' or 'Protein Sequence' columns")
            
            logger.debug(f"Cleaned input sequence: {protein_seq} (len={len(protein_seq)})")
            results = []
            for _, row in gene_data.iterrows():
                gene_name = row['Name']
                gene_id = row['Gene-ID']
                sequence = row['Protein Sequence'].strip()
                if not sequence or not isinstance(sequence, str):
                    logger.warning(f"Invalid sequence for gene {gene_name}: {sequence}")
                    continue
                
                logger.debug(f"Processing S. No. {row['S. No.']}, Name={gene_name}, Gene-ID={gene_id}, Seq={sequence[:20]}... (len={len(sequence)})")
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
                    "Is PGPR": similarity >= 75
                })
            
            results = sorted(results, key=lambda x: x['Similarity (%)'], reverse=True)[:10]
            names = [result['Gene Name'] for result in results]
            similarities = [result['Similarity (%)'] for result in results]
            logger.debug(f"Top 10 results: {names} with similarities: {similarities}")
            
            chart_data = {
                'labels': names,
                'data': similarities
            }
            
            return render_template('analyzer.html', results=results, chart_data=chart_data)
        
        return render_template('analyzer.html')
    except Exception as e:
        logger.error(f"Error in analyzer: {str(e)}", exc_info=True)
        return render_template('analyzer.html', error=f"An error occurred: {str(e)}")

@app.route('/add_gene', methods=['GET', 'POST'])
def add_gene():
    if request.method == 'POST':
        gene_name = request.form.get('gene_name')
        protein_seq = request.form.get('protein_seq')
        new_entry = pd.DataFrame({
            'Name': [gene_name],
            'Protein Sequence': [protein_seq],
            'Gene-ID': ['TBD'],
            'S. No.': [len(pd.read_csv(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))) + 1]
        })
        new_entry.to_csv(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'), mode='a', header=False, index=False)
        return render_template('add_gene.html', success="Gene added successfully!")
    return render_template('add_gene.html')

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