from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def calculate_tanimoto_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return 0.0
    
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
    
    return round(DataStructs.TanimotoSimilarity(fp1, fp2), 6)

def find_substructure(scaffold, smiles):
    scaffold_mol = Chem.MolFromSmiles(scaffold)
    mol = Chem.MolFromSmiles(smiles)
    
    if scaffold_mol is None or mol is None:
        return False
    
    return mol.HasSubstructMatch(scaffold_mol)

def main():
    with open("analyzethis.txt", "r") as file1, open("PAINS.txt", "r") as file2:
        lines1 = file1.readlines()
        lines2 = file2.readlines()
        
    results_filename = "resultsscaffold.txt"
    with open(results_filename, "w") as output_file:
        pains_scaffolds = [line.strip() for line in lines2]
        
        for i, line1 in enumerate(lines1, start=1):
            line1 = line1.strip()
            compound_results = []  # To store the compound comparisons and smiles
            for j, line2 in enumerate(pains_scaffolds, start=1):
                line2 = line2.strip()
                score = calculate_tanimoto_similarity(line1, line2)
                found_substructure = find_substructure(line2, line1)
                
                if found_substructure or score > 0.3:
                    compound_results.append((j, score, line2, found_substructure))
            
            if compound_results:
                output_file.write(f"From compound {i} (SMILES: {line1}):\n")
                for result in compound_results:
                    pain_scaffold_msg = " (Exact scaffold found)" if result[3] else ""
                    output_file.write(f"  PAIN scaffold {result[0]} (SMILES: {result[2]}){pain_scaffold_msg}: {result[1]:.6f}\n")
                output_file.write("\n")
    
    # Find compounds from analyzethis.txt that are not mentioned in the result
    with open(results_filename, "r") as output_file:
        result_content = output_file.read()
    
    with open(results_filename, "a") as output_file:
        output_file.write("Unmentioned compounds:\n")
        for i, line1 in enumerate(lines1, start=1):
            line1 = line1.strip()
            mentioned = f"From compound {i} (SMILES: {line1})" in result_content
            if not mentioned:
                output_file.write(f"From compound {i} (SMILES: {line1})\n")
        output_file.write("\n")

if __name__ == "__main__":
    main()
