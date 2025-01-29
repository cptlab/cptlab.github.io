import sqlite3
import pandas as pd
import argparse
import os

def add_smarts_to_db(input_csv, db_path='smarts_patterns.db'):
    # Check if the database exists
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"The database '{db_path}' does not exist in the current directory.")

    # Connect to the SQLite3 database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Read the input CSV file
    try:
        data = pd.read_csv(input_csv)
    except Exception as e:
        raise ValueError(f"Error reading the input CSV file: {e}")

    # Validate required columns
    required_columns = {'name', 'smarts', 'source'}
    if not required_columns.issubset(data.columns):
        raise ValueError(f"The input CSV must contain the following columns: {', '.join(required_columns)}")

    # Process each row in the CSV file
    for _, row in data.iterrows():
        source = row['source']
        name = row['name']
        smarts = row['smarts']

        # Add the source to the sources table if not already present
        cursor.execute('SELECT id FROM sources WHERE source = ?', (source,))
        result = cursor.fetchone()
        if result:
            source_id = result[0]
        else:
            cursor.execute('INSERT INTO sources (source) VALUES (?)', (source,))
            conn.commit()
            source_id = cursor.lastrowid

        # Add the SMARTS pattern to the patterns table
        cursor.execute(
            'INSERT INTO patterns (name, smarts, source_id) VALUES (?, ?, ?)',
            (name, smarts, source_id)
        )

    # Commit changes to the database
    conn.commit()
    print(f"Successfully added data from '{input_csv}' to '{db_path}'.")

    # Close the database connection
    conn.close()

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add SMARTS patterns from a CSV file to the smarts_patterns.db SQLite database. "
            "The input CSV must contain the columns: 'name', 'smarts', and 'source'."
        )
    )
    parser.add_argument('input_csv', help='Path to the input CSV file containing SMARTS patterns.')
    parser.add_argument('--db-path', default='smarts_patterns.db', help='Path to the SQLite database (default: smarts_patterns.db).')
    args = parser.parse_args()

    try:
        add_smarts_to_db(args.input_csv, args.db_path)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
