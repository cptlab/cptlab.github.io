from flask import Flask, render_template, request, redirect, url_for
from model import run_prediction_model
from markupsafe import Markup


app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    smiles = request.form['smiles']
    prediction = run_prediction_model(smiles)
    prediction_value = round(prediction * 100, 3)

    result_text = "Mutagen" if prediction > 0.5 else "Non-mutagen"
    #result = f"<b>{result_text} ({prediction_percentage}%)</b><br>{ad_status} (Mahalanobis Distance: {round(distance, 2)})"
    result = f"<b>{result_text}, Mutagen probability: {prediction_value} % </b>"
    result_safe = Markup(result)
    return render_template('result.html', prediction=result_safe)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001)

