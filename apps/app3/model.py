# model.py
from rdkit import Chem
from rdkit.Chem import Descriptors


def process_smiles(smiles):

    #Calculate the molecular weight using RDKit
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_weight = Descriptors.MolWt(mol)
            return mol_weight
        else:
            return 'Invalid SMILES'
    except Exception as e:
        return str(e)




def process_csv(filepath):
    import pandas as pd

    df = pd.read_csv(filepath)
    # Assuming the CSV has a column named 'smiles'
    df['result'] = df['SMILES'].apply(process_smiles)
    result_filepath = filepath.replace('uploads', 'results').replace('.csv', '_result.csv')
    df.to_csv(result_filepath, index=False)
    return result_filepath

