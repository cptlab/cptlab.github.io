from flask import Flask, request, redirect, url_for, send_file, render_template
import pandas as pd
import os
from model import process_csv

app = Flask(__name__)

# Directory to save uploaded files
UPLOAD_FOLDER = 'uploads'
RESULT_FOLDER = 'results'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return 'No file part'

    file = request.files['file']
    if file.filename == '':
        return 'No selected file'

    if file and file.filename.endswith('.csv'):
        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)
        result_filepath = process_csv(filepath)
        return send_file(result_filepath, as_attachment=True)

    return 'Invalid file format. Please upload a CSV file.'


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5002)

