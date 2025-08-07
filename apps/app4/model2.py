import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw


def extract_smiles(compound):
    # Defensive: check for _record and props key
    if not hasattr(compound, '_record'):
        return None
    props = compound._record.get('props', [])
    for prop in props:
        urn = prop.get('urn', {})
        # Look for the SMILES with label 'SMILES' and name 'Absolute' (canonical)
        if urn.get('label') == 'SMILES' and urn.get('name') == 'Absolute':
            value = prop.get('value', {})
            if isinstance(value, dict) and 'sval' in value:
                return value['sval']
    return None



# Function to draw a molecule by searching for its name on PubChem

def draw_mol_from_name(name):
    try:
        compound = pcp.get_compounds(name, 'name')
        if compound:
            smiles = extract_smiles(compound[0])
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"RDKit failed to parse SMILES: {smiles}")
                legend_text = f"User query: {name}"
                mol_img = Draw.MolToImage(mol, size=(300, 300), legend=legend_text)
                return mol_img
            else:
                raise ValueError(f"No SMILES found for compound '{name}'")
        else:
            raise ValueError(f"No compound found for '{name}'")
    except Exception as e:
        raise ValueError(f"Error retrieving molecule: {str(e)}")


# Make predictions
# Flask will pass the name text to this function

def run_prediction_model(name):
    try:
        output = draw_mol_from_name(name)  # Get the image from your model
        return output  # Return the PIL Image object
    except Exception as e:
        raise ValueError(f"Error processing molecule named: {str(e)}")

