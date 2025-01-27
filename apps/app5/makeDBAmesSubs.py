import sqlite3

# Connect to SQLite3 database (or create it if it doesn't exist)
conn = sqlite3.connect('smarts_patterns.db')
cursor = conn.cursor()

# Create a table to store the sources
cursor.execute('''
CREATE TABLE IF NOT EXISTS sources (
    id INTEGER PRIMARY KEY,
    source TEXT NOT NULL
)
''')

# Insert the source into the sources table
source = "Benigni R., C. Bossa, N. Jeliazkova, T. Netzeva, and A. Worth. (2008) The Benigni / Bossa rulebase for mutagenicity and carcinogenicity - a module of Toxtree. European Commission report EUR 23241"
cursor.execute('INSERT INTO sources (source) VALUES (?)', (source,))
conn.commit()

# Get the source_id of the inserted source
cursor.execute('SELECT id FROM sources WHERE source = ?', (source,))
source_id = cursor.fetchone()[0]

# Create a table to store the SMARTS patterns with a source_id column
cursor.execute('''
CREATE TABLE IF NOT EXISTS patterns (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    smarts TEXT NOT NULL,
    source_id INTEGER NOT NULL,
    FOREIGN KEY (source_id) REFERENCES sources(id)
)
''')

# Insert the SMARTS patterns into the table with the source_id
patterns = [
    ("acyl_halide", "C(=O)[Cl,Br,I]", source_id),
    ("alkyl_ester_sulphonic_phosphonic", "[C,c][S,P](O[C,c])", source_id),
    ("n_methylol_derivative", "[NX3][CH2][OH]", source_id),
    ("monohaloalkene", "C=C[Cl,Br,I]", source_id),
    ("s_n_mustard", "[NX3,SX2][CX4][CX4][Cl,Br,I]", source_id),
    ("propiolactone", "C1C(=O)O1", source_id),
    ("propiosultone", "C1C(=O)S1", source_id),
    ("epoxide", "C1OC1", source_id),
    ("aziridine", "C1NC1", source_id),
    ("aliphatic_halogen", "[CX4][Cl,Br,I]", source_id),
    ("alkyl_nitrite", "[CX4][NX2]O", source_id),
    ("alpha_beta_unsaturated_carbonyl", "C=CC(=O)", source_id),
    ("simple_aldehyde", "[CX3H1]=O", source_id),
    ("quinone", "C1=CC(=O)C=CC1=O", source_id),
    ("hydrazine", "NN", source_id),
    ("aliphatic_azo", "N=N", source_id),
    ("aliphatic_azoxy", "N=N[O]", source_id),
    ("isocyanate", "N=C=O", source_id),
    ("isothiocyanate", "N=C=S", source_id),
    ("alkyl_carbamate", "O=C(OC)[NX3]", source_id),
    ("alkyl_thiocarbamate", "O=C(SC)[NX3]", source_id),
    ("polycyclic_aromatic_hydrocarbon", "c1ccc2ccccc2c1", source_id),
    ("heterocyclic_polycyclic_aromatic_hydrocarbon", "c1ccc2[n,o,s]cccc2c1", source_id),
    ("alkyl_n_nitroso", "[NX3]N", source_id),
    ("aryl_n_nitroso", "[a][NX3]N", source_id),
    ("azide", "[NX2]N=[NX1]=[NX1]", source_id),
    ("triazene", "[NX3][NX2]=[NX2]", source_id),
    ("aliphatic_n_nitro", "[NX3][N+](=O)[O-]", source_id),
    ("alpha_beta_unsaturated_aliphatic_alkoxy", "C=CCOC", source_id),
    ("aromatic_nitroso", "[a]N", source_id),
    ("aromatic_ring_n_oxide", "[a][N+](=O)[O-]", source_id),
    ("nitro_aromatic", "[a][N+](=O)[O-]", source_id),
    ("primary_aromatic_amine", "[a][NH2]", source_id),
    ("hydroxyl_amine", "[NX3][OH]", source_id),
    ("hydroxyl_amine_ester", "[NX3][OX2][CX3]", source_id),
    ("aromatic_mono_dialkylamine", "[a][NX3][CX4]", source_id),
    ("aromatic_n_acyl_amine", "[a][NX3][CX3]=O", source_id),
    ("aromatic_diazo", "[a][N]=[N]", source_id),
    ("coumarin", "O=C1OCc2ccccc12", source_id),
    ("furocoumarin", "O=C1OCc2cc3ccccc3oc12", source_id)
]


cursor.executemany('INSERT INTO patterns (name, smarts, source_id) VALUES (?, ?, ?)', patterns)
conn.commit()

# Close the database connection
conn.close()
