import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

# Function to draw a molecule by searching for its name on PubChem
def draw_mol_from_name(name):
    try:
        compound = pcp.get_compounds(name, 'name')  # Search PubChem
        if compound:
            smiles = compound[0].canonical_smiles
            mol = Chem.MolFromSmiles(smiles)
            legend_text = f"User query: {name}"
            mol_img = Draw.MolToImage(mol, size=(300, 300), legend=legend_text)  # Create an image

            # No need to convert to PIL image; it already is one
            return mol_img
        else:
            raise ValueError(f"No compound found for the name: {name}")
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

