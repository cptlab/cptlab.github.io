from flask import Flask, render_template, request, redirect, url_for
from model import run_prediction_model

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    smiles = request.form['smiles']
    prediction = run_prediction_model(smiles)
    result = "Mutagen" if prediction > 0.5 else "Non-mutagen"
    return render_template('result.html', prediction=result)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001)

