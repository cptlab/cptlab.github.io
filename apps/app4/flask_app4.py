from flask import Flask, render_template, request, redirect, url_for, send_from_directory
from model import run_prediction_model
import os
from io import BytesIO
from PIL import Image

app = Flask(__name__)

# Path to store temporary images
UPLOAD_FOLDER = 'static/temp'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    try:
        name = request.form['name']
        image = run_prediction_model(name)

        # Save the image to a temporary file
        image_filename = os.path.join(app.config['UPLOAD_FOLDER'], f"{name}.png")
        image.save(image_filename)

        # Generate the URL for the image
        image_url = url_for('send_image', filename=f"{name}.png", _external=True)

        return render_template('result.html', image_url=image_url)
    except ValueError as e:
        #Render error page.
        return render_template('error.html', error_message=str(e)), 500

@app.route('/Model4/static/temp/<filename>')
def send_image(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5003)

