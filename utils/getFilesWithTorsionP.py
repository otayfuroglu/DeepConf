#
import os
import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
from rdkit.Chem import TorsionFingerprints
import pandas as pd
import tqdm 

def getTorsionPoints(rd_mol):
    torsion_points = []
    for torsions_list in TorsionFingerprints.CalculateTorsionLists(rd_mol):
        for torsions in torsions_list:
            if 180 in torsions:
                torsion_points.append(torsions[0][0])
    return(torsion_points)


#sdf_dir = "all_SDFs/"
#file_names = [file_name for file_name in os.listdir(sdf_dir) if ".sdf" in file_name]

csv_path = "sampl9.csv"
df = pd.read_csv(csv_path, header=None)

for i in tqdm.tqdm(range(len(df))):
    file_name = df[0][i]
    smile = df[1][i]
    rd_mol = Chem.MolFromSmiles(smile)
    torsions = getTorsionPoints(rd_mol)
    if len(torsions) <= 5:

        # for generate 3D structure
        rd_mol = Chem.AddHs(rd_mol)
        AllChem.EmbedMolecule(rd_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(rd_mol,1000)

        # write 3D structure as sdf
        with Chem.rdmolfiles.SDWriter(f"torPointUp5/{file_name}.sdf") as writer:
            writer.write(rd_mol)


