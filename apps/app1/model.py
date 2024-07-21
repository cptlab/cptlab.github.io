from rdkit import Chem

# Put required code for the model here.
# Then call the function inside run_prediction_model below.

def contains_aromatic_rings(mol):
    # Check for aromatic rings
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            return "The molecule contains aromatic rings."

    return "The molecule does not contain aromatic rings."


def run_prediction_model(smiles):
    # Flask will pass smiles txt to the function called run_prediction_model
    # and display what is returned.
    # First we check for valid smiles.
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError("Invalid smiles string")
    try:
        Chem.SanitizeMol(mol)
        prediction = contains_aromatic_rings(mol)
        return prediction
    except Exception as e:
        raise ValueError(f"Error processing smiles string: {str(e)}")


