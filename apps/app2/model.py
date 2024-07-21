import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import AllChem

# Define your neural network model class
class Net(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.fc3 = nn.Linear(hidden_size, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        out = self.fc1(x)
        out = torch.relu(out)
        out = self.fc2(out)
        out = torch.relu(out)
        out = self.fc3(out)
        out = self.sigmoid(out)
        return out

# Load the model function
def load_model():
    input_size = 2048  # Size of Morgan fingerprint
    hidden_size = 128
    output_size = 1
    loaded_model = Net(input_size, hidden_size, output_size)
    checkpoint = torch.load("model.pt", map_location=torch.device('cpu'))
    loaded_model.load_state_dict(checkpoint)
    loaded_model.eval()
    return loaded_model

# Load the model
model = load_model()

# Function to convert SMILES string to Morgan fingerprints
def smiles_to_morgan_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    return fingerprint

# Define a function to convert a SMILES string to a tensor
def smiles_to_tensor(smiles):
    fingerprint = smiles_to_morgan_fingerprint(smiles)
    tensor = torch.Tensor(fingerprint).unsqueeze(0)
    return tensor

# Make predictions
# Flask will pass smiles txt to the function called run_prediction_model

def run_prediction_model(smiles):
    try:
        with torch.no_grad():
            input_tensor = smiles_to_tensor(smiles)
            output = model(input_tensor)
            prediction = output.item()
        return prediction
    except Exception as e:
        raise ValueError(f"Error processing smiles string: {str(e)}")


