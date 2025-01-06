import sqlite3

# Step 1: Set up the SQLite database and table (runs once, only if dvd_catalog.db doesn't exist)
def setup_database():
    conn = sqlite3.connect('dvd_catalog.db')
    cursor = conn.cursor()

    cursor.execute('''
    CREATE TABLE IF NOT EXISTS dvds (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        title TEXT NOT NULL,
        genre TEXT NOT NULL,
        year INTEGER,
        director TEXT,
        rating REAL
    )
    ''')

    conn.commit()
    conn.close()

# Step 2: Define functions to interact with the database

def add_dvd(title, genre, year, director, rating):
    conn = sqlite3.connect('dvd_catalog.db')
    cursor = conn.cursor()

    cursor.execute('''
    INSERT INTO dvds (title, genre, year, director, rating)
    VALUES (?, ?, ?, ?, ?)
    ''', (title, genre, year, director, rating))

    conn.commit()
    conn.close()

def search_by_title(title):
    conn = sqlite3.connect('dvd_catalog.db')
    cursor = conn.cursor()

    cursor.execute('''
    SELECT * FROM dvds WHERE title LIKE ?
    ''', ('%' + title + '%',))

    results = cursor.fetchall()
    conn.close()

    return results

def filter_by_genre(genre):
    conn = sqlite3.connect('dvd_catalog.db')
    cursor = conn.cursor()

    cursor.execute('''
    SELECT * FROM dvds WHERE genre = ?
    ''', (genre,))

    results = cursor.fetchall()
    conn.close()

    return results

def display_all_dvds():
    conn = sqlite3.connect('dvd_catalog.db')
    cursor = conn.cursor()

    cursor.execute('SELECT * FROM dvds')
    results = cursor.fetchall()
    conn.close()

    return results

# Step 3: Create the CLI interface
def menu():
    setup_database()  # Ensure the database is set up when the script runs

    while True:
        print("\nDVD Catalog App")
        print("1. Add a new DVD")
        print("2. Search DVDs by Title")
        print("3. Filter DVDs by Genre")
        print("4. Display all DVDs")
        print("5. Exit")

        choice = input("Choose an option (1-5): ")

        if choice == '1':
            title = input("Enter DVD title: ")
            genre = input("Enter genre: ")
            year = int(input("Enter release year: "))
            director = input("Enter director: ")
            rating = float(input("Enter rating (0-10): "))
            add_dvd(title, genre, year, director, rating)
            print(f"{title} added to catalog.")

        elif choice == '2':
            title = input("Enter title to search: ")
            results = search_by_title(title)
            if results:
                print("\nSearch Results:")
                for dvd in results:
                    print(dvd)
            else:
                print("No DVDs found with that title.")

        elif choice == '3':
            genre = input("Enter genre to filter by: ")
            results = filter_by_genre(genre)
            if results:
                print("\nFiltered Results:")
                for dvd in results:
                    print(dvd)
            else:
                print(f"No DVDs found in the genre '{genre}'.")

        elif choice == '4':
            all_dvds = display_all_dvds()
            if all_dvds:
                print("\nAll DVDs in your catalog:")
                for dvd in all_dvds:
                    print(dvd)
            else:
                print("No DVDs in the catalog.")

        elif choice == '5':
            print("Exiting the catalog.")
            break
        else:
            print("Invalid option. Please choose again.")

# Run the program
if __name__ == "__main__":
    menu()
