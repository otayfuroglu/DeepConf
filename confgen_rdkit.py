from rdkit import Chem
from rdkit.Chem import AllChem
import os, shutil



def getConfs():
    mol = next(Chem.SDMolSupplier(mol_path, removeHs=False))

    writer = Chem.SDWriter(f'{work_dir}/output_sorted_conformers.sdf')

    AllChem.EmbedMultipleConfs(mol,
                               numConfs=10,
                               maxAttempts=1000,
                               pruneRmsThresh=0.5,
                              )

    conformer_energies = []

    # Optimize each conformer and calculate energy
    for conf_id in range(mol.GetNumConformers()):
        # Option 1: UFF optimization and energy calculation
        #  ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)

         # Option 2: MMFF94 optimization
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)

        ff.Minimize()
        energy = ff.CalcEnergy()

        conformer_energies.append((conf_id, energy))

    conformer_energies.sort(key=lambda x: x[1])

    for conf_id, energy in conformer_energies:
        mol.SetProp('_Conformer_Energy', str(energy))
        mol.SetProp('_Name', str(f"conformer_{conf_id}"))
        writer.write(mol, confId=conf_id)
    writer.close()



struc_dir = "test_akocak"
for fl in os.listdir(struc_dir):
    if not fl.endswith(".sdf"):
        continue
    work_dir = fl.split('.')[0]
    print(work_dir)
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.mkdir(work_dir)
    mol_path = f"{struc_dir}/{fl}"
    getConfs()

