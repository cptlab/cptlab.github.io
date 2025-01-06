from flask import Flask, request, render_template, redirect, url_for
import sqlite3

app = Flask(__name__)

# Database connection function
def get_db_connection():
    conn = sqlite3.connect('dvd_catalog.db')
    conn.row_factory = sqlite3.Row  # To return rows as dictionaries
    return conn

# Display menu for DVD DB
@app.route('/')
def index():
    return render_template('menu.html')

# Route to display all DVDs
@app.route('/all')
def show_all_dvds():
    conn = get_db_connection()
    dvds = conn.execute('SELECT * FROM dvds').fetchall()
    conn.close()
    return render_template('all.html', dvds=dvds)

# Route to add a new DVD
@app.route('/add', methods=['GET', 'POST'])
def add_dvd():
    if request.method == 'POST':
        title = request.form['title']
        genre = request.form['genre']
        year = request.form['year']
        director = request.form['director']
        rating = request.form['rating']

        conn = get_db_connection()
        conn.execute('INSERT INTO dvds (title, genre, year, director, rating) VALUES (?, ?, ?, ?, ?)',
                     (title, genre, year, director, rating))
        conn.commit()
        conn.close()

        return redirect(url_for('index'))

    return render_template('add_dvd.html')

# Route to search DVDs by title
@app.route('/search', methods=['GET', 'POST'])
def search():
    if request.method == 'POST':
        title = request.form['title']
        conn = get_db_connection()
        dvds = conn.execute('SELECT * FROM dvds WHERE title LIKE ?', ('%' + title + '%',)).fetchall()
        conn.close()
        return render_template('search_results.html', dvds=dvds)

    return render_template('search.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5004)
