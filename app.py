import os
from flask import Flask, request, render_template, jsonify
from html import escape
import pandas as pd
from Bio.Align import PairwiseAligner

app = Flask(__name__)

def load_gene_data(file_path):
    return pd.read_csv(file_path)

def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    score = aligner.score(seq1, seq2)
    similarity = (score / max(len(seq1), len(seq2))) * 100
    return similarity

def analyze_sequence(protein_seq, gene_data, threshold=75):
    results = []
    for _, row in gene_data.iterrows():
        gene_name = row['Gene Name']
        sequences = [row['Protein Sequence'].strip()]
        similarities = [calculate_similarity(protein_seq, seq) for seq in sequences]
        avg_similarity = sum(similarities) / len(similarities)
        results.append({
            "Gene Name": gene_name,
            "Similarity (%)": similarities[0],
            "Is PGPR": avg_similarity >= threshold
        })
    return results

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    if request.method == 'POST':
        protein_seq = escape(request.form['protein_seq'].strip())
        if not protein_seq:
            return render_template('analyzer.html', error="Please provide a protein sequence.")
        
        gene_data = load_gene_data(os.path.join(os.path.dirname(__file__), 'PGPRgene.csv'))
        results = analyze_sequence(protein_seq, gene_data)
        
        chart_data = {
            'labels': [result['Gene Name'] for result in results[:10]],
            'data': [result['Similarity (%)'] for result in results[:10]]
        }
        
        return render_template('analyzer.html', results=results, chart_data=chart_data)
    
    return render_template('analyzer.html')

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