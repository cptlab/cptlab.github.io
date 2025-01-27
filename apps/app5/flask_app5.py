from flask import Flask, render_template, request, redirect, url_for
from rdkit import Chem
import sqlite3
from markupsafe import Markup

app = Flask(__name__)

# Database connection function
def get_db_connection():
    conn = sqlite3.connect('smarts_patterns.db')
    conn.row_factory = sqlite3.Row  # To return rows as dictionaries
    return conn

# Display the main page with the search form
@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

# Route to handle the search and display results
@app.route('/results', methods=['POST'])
def results():
    smiles = request.form['smiles'].strip()
    error = None
    alerts = []

    if not smiles:
        error = "Please enter a SMILES string."
    else:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            error = "Invalid SMILES string. Please try again."
        else:
            # Fetch patterns from the database
            conn = get_db_connection()
            patterns = conn.execute('SELECT name, smarts FROM patterns').fetchall()
            conn.close()

            # Check for substructure matches
            for pattern in patterns:
                smarts = pattern['smarts']
                name = pattern['name']
                mol_pattern = Chem.MolFromSmarts(smarts)
                if mol.HasSubstructMatch(mol_pattern):
                    alerts.append(name)

    return render_template('results.html', smiles=smiles, error=error, alerts=alerts)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5005)

