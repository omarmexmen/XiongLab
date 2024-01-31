import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

def read_sdf_and_convert_to_csv(filename):
    # Path to the SDF file
    filepath = os.path.join(os.path.dirname(__file__), filename)
    
    # Check if the file exists
    if not os.path.isfile(filepath):
        print(f"File {filename} not found in the script's directory.")
        return
    
    # Read the SDF file
    sdf_data = Chem.SDMolSupplier(filepath)
    
    # Convert SDF data to DataFrame
    df = PandasTools.LoadSDF(filepath)
    
    # List to hold the data
    all_data = []
    
    # Extract SMILES, other information, and ID from each molecule
    for mol in sdf_data:
        if mol is not None:  # Check if the molecule is valid
            smiles = Chem.MolToSmiles(mol)
            mol_id = mol.GetProp('_Name') if mol.HasProp('_Name') else 'No ID'
            other_info = {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}  # Get other properties
            
            # Append to the list
            all_data.append({'SMILES': smiles, 'ID': mol_id, **other_info})
    
    # Create a DataFrame from the collected data
    all_data_df = pd.DataFrame(all_data)
    
    # Save the DataFrame to a CSV file
    output_file = os.path.join(os.path.dirname(__file__), 'convertedall.csv')
    all_data_df.to_csv(output_file, index=False)
    print(f"Data has been written to {output_file}")

# Call the function with 'acids.sdf', you can change the name if needed.
read_sdf_and_convert_to_csv('acids.sdf')
