
import os
from collections import defaultdict
import numpy as np
import argparse
import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
from scipy.cluster.hierarchy import linkage, fcluster


def getClusterRMSDFromFiles(conf_dir, rmsd_thresh):
    suppl_list = {Chem.SDMolSupplier(f"{conf_dir}/{fl_name}"):
                  fl_name for fl_name in os.listdir(conf_dir)
                  if fl_name.endswith(".sdf")}
    mol_list = []
    for suppl, fl_name in suppl_list.items():
        mol = next(suppl)
        mol.SetProp("_Name", fl_name)
        mol_list.append(mol)

    n_mol=len(mol_list)
    if n_mol <= 1:
        print("Clustering do not applied.. There is just one conformer")
        return 0
    dist_matrix=np.empty(shape=(n_mol, n_mol))
    for i, mol1 in enumerate(mol_list):
        for j, mol2 in enumerate(mol_list):
            # first aling each conformer and then calculate rmsd
            #  Chem.rdMolAlign.AlignMol(mol1, mol2)
            dist_matrix[i, j] = Chem.rdMolAlign.GetBestRMS(mol1, mol2)
    linked = linkage(dist_matrix,'complete')
    cluster_conf = defaultdict(list)
    labelList = [mol.GetProp('_Name') for mol in mol_list]
    for key, fl_name in zip(fcluster(linked, rmsd_thresh, criterion='distance'), labelList):
        cluster_conf[key].append(fl_name)

        # save clusturedd files seperately
        directory = f"{conf_dir}/cluster_{key}"
        if not os.path.exists(directory):
            os.mkdir(directory)
        file_path = f"{directory}/{fl_name}"
        for mol in mol_list:
            if mol.GetProp('_Name') == fl_name:
                mol = mol
                break
        with Chem.rdmolfiles.SDWriter(file_path) as writer:
            writer.write(mol)

    #  return cluster_conf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Give something ...")
    #  parser.add_argument("-mof_num", "--mof_num",
    #                      type=int, required=True,
                        #  help="..")
    parser.add_argument("-confdir", type=str, required=True, help="..")
    parser.add_argument("-rmsd_thresh", type=str, required=True, help="..")
    args = parser.parse_args()

    getClusterRMSDFromFiles(conf_dir=args.confdir, rmsd_thresh=args.rmsd_thresh)
