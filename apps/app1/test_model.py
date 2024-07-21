from model import run_prediction_model

test_smiles = "COc1ccc(C(=O)[O-])cc1"  # A valid SMILES string
invalid_smiles = "invalid"  # An invalid SMILES string

# Test with a valid SMILES string
try:
    result = run_prediction_model(test_smiles)
    print(f"Prediction for {test_smiles}: {result}")
except ValueError as e:
    print(f"Error: {e}")

# Test with an invalid SMILES string
try:
    result = run_prediction_model(invalid_smiles)
    print(f"Prediction for {invalid_smiles}: {result}")
except ValueError as e:
    print(f"Error: {e}")

