#
from openbabel import openbabel, pybel
from ase import Atoms
from ase.io import write

import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
#  import os_util
from collections import defaultdict
import sys, os
import numpy as np
import pandas as pd

from multiprocessing import Pool, cpu_count
from itertools import product, repeat
from functools import wraps

from scipy.cluster.vq import kmeans, vq, whiten
from scipy.cluster.hierarchy import linkage, fcluster

from random import sample


NPROCS_ALL = int(cpu_count())
print("Number of total cpu core: ", NPROCS_ALL)


def calcFuncRunTime(func):
    import time

    @wraps(func)
    def wrapper(*args, **kwargs):
        s_time = time.time()
        func(*args, **kwargs)
        print(f"Function {func.__name__} executed in {(time.time()-s_time)/60:.5f} m")
    return wrapper


def calcRMSDsymm(pair_idx, mol_list):

    idx1 = pair_idx[0]
    idx2 = pair_idx[1]
    if idx1 < idx2:
        if len(mol_list) == 1:
            mol = mol_list[0]
            return (AllChem
                    .GetConformerRMS(mol,
                                     idx1,
                                     idx2,
                                     prealigned=False)
                   )
        else:
            return (Chem.rdMolAlign
                    .GetBestRMS(mol_list[idx1],
                                mol_list[idx2])
                   )


#  @calcFuncRunTime
def getDistMatrix(mol_list, conformerIds=None):

    n_mol=len(mol_list)
    if n_mol == 1 and conformerIds:
        n_mol = len(conformerIds)

    if n_mol <= 1:
        print("Clustering do not applied.. There is just one conformer")
        return None

    with Pool(NPROCS_ALL) as pool:
        results = pool.starmap(calcRMSDsymm,
                               zip(product(range(n_mol), repeat=2),
                                   repeat(mol_list)))

    ordered_all_rmsd = [result for result in results if result]
    return symmetricize(n_mol, ordered_all_rmsd)


def symmetricize(n: int, list1D: list) -> np.array:

    dist_matrix=np.zeros(shape=(n, n))
    i = 0
    for idx1 in range(n):
        for idx2 in range(n):
            if idx1 == idx2:
                dist_matrix[idx1, idx2] = 0.0
            elif idx1 < idx2:
                dist_matrix[idx1, idx2] = list1D[i]
                i += 1
    return dist_matrix + dist_matrix.T



    def run(self, fmax=0.05, steps=100):
        """Run the optimization until convergence or maximum steps."""
        for step in range(self.maxiter):
            positions, energy = self.step(self.atoms.get_forces())
            max_force = np.linalg.norm(self.atoms.get_forces(), axis=1).max()

            if max_force < fmax:
                print(f"Converged at step {step + 1} with max force {max_force:.6f} eV/A.")
                break
            else:
                print(f"Step {step + 1}: Energy = {energy:.6f} eV, Max force = {max_force:.6f} eV/A.")
        else:
            print("Reached maximum number of iterations without convergence.")


class confGen:
    """

    """

    def __init__(self, mol_path, addH, WORK_DIR):
        self.mol_path = mol_path
        self.WORK_DIR = WORK_DIR

        # for activete g16 optmization algorithm
        self.optG16 = False

        # set add missing H
        self.addH = addH

        # initialize RW mol
        self.rw_mol = None
        self._loadRWMol()

        # initialize calcultor
        self.calculator = None

        # initialize geom opt paprameters
        self.maxiter = None
        self.fmax = None

        # initialize optimization method
        self.opt_method = None

        # trial number
        self.n_trial = 1

    def getFileBase(self):
        return self.mol_path.split("/")[-1].split(".")[0]

    def increaseTrialNum(self):
        self.n_trial += 1

    def _getFileFormat(self, file_path=None):
        if file_path:
            return file_path.split(".")[-1]

        return self.mol_path.split(".")[-1]

    def _loadRWMol(self):
        self._loadMolWithOB()

    def _loadMolWithRW(self, mol_path, sanitize=True):
        rd_mol = next(Chem.SDMolSupplier(mol_path, sanitize=sanitize, removeHs=False))
        if sanitize is False:
            rd_mol.UpdatePropertyCache(strict=False)
        self.rw_mol = Chem.RWMol(rd_mol)

    def _rdKekuleizeError(self, rd_mol):
        # correction  kekuleize error (especially N in aromatic ring)
        print("\nWarning!: There is kekulize error, ingnored sanitize and kekulized for N atom which is in aromatic ring\n")
        for i, atom in enumerate(rd_mol.GetAtoms()):
            if atom.GetSymbol() == "N" and atom.GetIsAromatic():
                print("Aromatic N atom idex: ",i+1)
                atom.SetNumExplicitHs(1)
        return rd_mol

    def _loadMolWithOB(self):

        pb_mol = next(pybel.readfile(self._getFileFormat(), self.mol_path))
        tmp_file_name = f"{self.WORK_DIR}/tmp_ob_file.sdf"

        #add hydrogen with openbabel
        if self.addH:
            if self._getFileFormat() == "xyz":
                print("Error: Cannot add hydrogen atoms to XYZ file format!!!")
                exit(1)

            if self._getFileFormat() != "sdf":
                # coorection sfg for add true Hydrogen
                #  print(self._getFileFormat())
                pb_mol.write("sdf", tmp_file_name, overwrite=True)

                corr_tmp_file_name = "corr_tmp_ob_file.sdf"
                corr_part = "  0  0  0  0  0  0  0  0  0  0"
                with open(corr_tmp_file_name, "w") as corr_sdf:
                    with open(tmp_file_name) as lines:
                        for line in lines:
                            if len(line) == 70:
                                line = line[:40] + corr_part + "\n"
                            corr_sdf.write(line)

                pb_mol = next(pybel.readfile("sdf", corr_tmp_file_name))
                self._rmFileExist(tmp_file_name)
                self._rmFileExist(corr_tmp_file_name)

            pb_mol.addh()
            pb_mol.make3D()

        #  # openbabel file to rdkit mol2 file
        pb_mol.write("sdf", tmp_file_name, overwrite=True)

        # laod as RW file
        try:
            self._loadMolWithRW(tmp_file_name)
        except:
            self._loadMolWithRW(tmp_file_name, sanitize=False)
            self.rw_mol = self._rdKekuleizeError(self.rw_mol)

        # optmization for just added H
        if self.addH:
            self.setOptMethod(opt_method="LBFGS")
            self.setOptParams(fmax=0.05, maxiter=200)
            self.setANI2XCalculator()
            self.geomOptimization(fix_heavy_atoms=True)

    def addHwithRD(self):
        self.rw_mol = rdkit.Chem.rdmolops.AddHs(self.rw_mol, addCoords=True)

    def _rmFileExist(self, file_path):
        if os.path.exists(file_path):
            os.remove(file_path)

    def writeRWMol2File(self, file_path, **kwargs):

        # add missing H
        #  self.addHwithRD()

        file_format = self._getFileFormat(file_path)

        if file_format == "xyz":
            rdkit.Chem.rdmolfiles.MolToXYZFile(self.rw_mol, file_path)
        elif file_format == "pdb":
            rdkit.Chem.rdmolfiles.MolToPDBFile(self.rw_mol, file_path)
        elif file_format == "sdf":
            with Chem.rdmolfiles.SDWriter(file_path) as writer:
                for key, value in kwargs.items():
                    self.rw_mol.SetProp(key, str(value))
                    writer.write(self.rw_mol)

        else:
            print("Unknown file format")
            sys.exit(1)

    def _writeConf2File(self, mol, conformerId, file_path):
        with rdkit.Chem.SDWriter(file_path) as w:
            w.write(mol, conformerId)
            w.flush()
            w.close()

    def _getTorsionPoints(self):
        from rdkit.Chem import TorsionFingerprints
        torsion_points = []
        for torsions_list in TorsionFingerprints.CalculateTorsionLists(self.rw_mol):
            for torsions in torsions_list:
                if 180 in torsions:
                    torsion_points.append(torsions[0][0])
        return(torsion_points)

    def _getNumConfs(self, nfold, scaled=1):
        n_torsions = len(self._getTorsionPoints())
        if n_torsions >= 10:
            return 8096
        else:
            return int(scaled * nfold ** n_torsions)

    #  @calcFuncRunTime
    def _getClusterKmeansFromConfIds(self, conformerIds, dist_matrix, n_group):

        cluster_conf_id = defaultdict(list)
        whitened = whiten(dist_matrix)
        centroids, _ = kmeans(whitened, n_group)
        cluster, _ = vq(whitened,centroids)
        for key, value in zip(cluster, conformerIds):
            cluster_conf_id[key].append(value)

        return cluster_conf_id

    #  @calcFuncRunTime
    def _getClusterRMSDFromFiles(self, conf_dir, rmsd_thresh):

        mol_dict = {next(Chem.SDMolSupplier(f"{conf_dir}/{fl_name}", removeHs=False)):
                      fl_name for fl_name in os.listdir(conf_dir)
                      if fl_name.endswith(".sdf")}
        mol_list = []
        for mol, fl_name in mol_dict.items():
            mol.SetProp("_Name", fl_name)
            mol_list.append(mol)

        dist_matrix = getDistMatrix(mol_list, conformerIds=None)
        if dist_matrix is None:
            return 0

        linked = linkage(dist_matrix,'complete')
        labelList = [mol.GetProp('_Name') for mol in mol_list]
        cluster_conf = defaultdict(list)
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

        return cluster_conf

    def _getCluster_diffE(self, files_minE, diffE_thresh=0.001):
        n_files = len(files_minE)
        dist_matrix = np.zeros(shape=(n_files, n_files))
        for i , val1 in enumerate(files_minE.values()):
            for j, val2 in  enumerate(files_minE.values()):
                e_diff = val1 - val2
                dist_matrix[i, j] = abs(val1 - val2 )
        #  print(dist_matrix)
        linked = linkage(dist_matrix,'complete')
        label_list = list(files_minE.keys())
        cluster_conf = defaultdict(list)
        for key, fl_name in zip(fcluster(linked, diffE_thresh, criterion='distance'), label_list):
            cluster_conf[key].append(fl_name)
        return cluster_conf

    def _pruneOptConfs(self, cluster_conf, confs_energies, conf_dir, opt_prune_diffE_thresh):
        i = 0

        local_files_minE = {}
        # for rmsd filter
        print("Applied diff RMSD filter (Angstrom)")
        for fl_names in cluster_conf.values():
            for j, fl_name in enumerate(fl_names):
                #  e = float(confs_energies.loc[confs_energies["FileName"] == fl_name, " Energy(eV)"].item())
                e = float(confs_energies.loc[confs_energies["FileName"] == fl_name, "Energy(eV)"].item())
                if i == 0:
                    global_minE = e
                    global_minE_file = fl_name
                else:
                    if global_minE > e:
                        global_minE = e
                        global_minE_file = fl_name
                i += 1

                if j == 0:
                    minE = e
                    minE_file = fl_name
                else:
                    if minE > e:
                        minE = e
                        minE_file = fl_name

            fl_names.remove(minE_file)
            local_files_minE[minE_file] = minE

            if len (fl_names) != 0:
                for rm_file in fl_names:
                    print("Removed", rm_file)
                    os.remove(f"{conf_dir}/{rm_file}")

        # for the energy filter
        if len(local_files_minE) > 1:
            print("Applied diff Energy filter (eV/Atom)")
            cluster_conf = self._getCluster_diffE(local_files_minE, diffE_thresh=opt_prune_diffE_thresh)
            for fl_names in cluster_conf.values():
                if len(fl_names) > 1:
                    for fl_name in fl_names[1:]: # remove all file except first
                        if fl_name == f"pruned_{global_minE_file}": # if any candidate removed file is global min 
                            fl_name = fl_names[0] #  remove first file
                        print("Removed", fl_name)
                        os.remove(f"{conf_dir}/{fl_name}")
                        del local_files_minE[fl_name]

        local_files_minE_sorted = dict(sorted(local_files_minE.items(), key=lambda item: item[1]))

        with Chem.SDWriter(f"{conf_dir}/{self.getFileBase()}_output.sdf") as w:
            for fl_name, e in local_files_minE_sorted.items():
                mol = next(Chem.SDMolSupplier(f"{conf_dir}/{fl_name}", removeHs=False))
                mol.SetProp("Energy", str(e))
                mol.SetProp("_Name", fl_name)
                w.write(mol)
                os.remove(f"{conf_dir}/{fl_name}")

    def genGonformers(self, file_path,
                         numConfs=100,
                         ETKDG=False,
                         maxAttempts=10000,
                         pruneRmsThresh=0.1,
                         mmCalculator=False,
                         optimization_conf=False,
                         opt_prune_rms_thresh=0.2,
                         opt_prune_diffE_thresh=0.001,
                         saveConfs=True,
                         useExpTorsionAnglePrefs=True,
                         useBasicKnowledge=True,
                         enforceChirality=True,
                         nfold=2,
                         npick=2,
                         nscale=1,
                        ):

        import copy

        #  self.addHwithRD()
        print("Woking on conformer generation process")
        mol = copy.deepcopy(self.rw_mol)
        if numConfs == 0 or numConfs < self._getNumConfs(nfold, scaled=nscale):
            numConfs = self._getNumConfs(nfold, scaled=nscale)
            print(f"Maximum number of conformers setting to {numConfs}")

        if ETKDG:
            ps = rdkit.Chem.rdDistGeom.ETKDGv3()
            conformerIds = list(rdkit.Chem.rdDistGeom.EmbedMultipleConfs(
                mol,
                numConfs,
                ps
            ))
        else:
            conformerIds = list(rdkit.Chem.AllChem.EmbedMultipleConfs(
                mol,
                numConfs=numConfs,
                maxAttempts=maxAttempts,
                pruneRmsThresh=pruneRmsThresh,
                useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                useBasicKnowledge=useBasicKnowledge,
                enforceChirality=enforceChirality,
                numThreads=NPROCS_ALL,
            ))

        # file for saving energies
        file_csv = open("%s/all_confs_sp_energies.csv" %self.WORK_DIR, "w")
        print("FileName,Energy(eV)", file=file_csv)

        print("Number of generated conformation: %d" %len(conformerIds))

        #  for k-means clutering
        #  dist_matrix = self._getConfDistMatrix(mol, conformerIds)
        print("Obtaining pairwise distance distribution matrix")
        dist_matrix = getDistMatrix([mol], conformerIds)

        print("Processing k-means clustering")
        cluster_conf_id = self._getClusterKmeansFromConfIds(conformerIds, dist_matrix,
                                           n_group=self._getNumConfs(nfold, scaled=1)
                                          )
        print("Calculating SP energies")
        minEConformerIDs = []
        all_picked_confs = []

        for cluster, clustered_confIds in cluster_conf_id.items():

            if saveConfs:
                CONF_DIR = self.WORK_DIR + f"/confs_cluster_{cluster}"
                if not os.path.exists(CONF_DIR):
                    os.mkdir(CONF_DIR)

            for i, conformerId  in enumerate(clustered_confIds):
                #  print("%d. conformer processing..." %i)
                if saveConfs:
                    prefix = ""
                    conf_file_path = "%s/conf_%d.sdf"%(CONF_DIR, conformerId)
                    self._writeConf2File(mol, conformerId, conf_file_path)

                #create ase atoms
                ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
                if mmCalculator:
                    e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
                else:
                    e, _ = self._calcSPEnergy(mol, conformerId)

                if i == 0:
                    minE = e
                    minEConformerID = conformerId
                    minE_ase_atoms = ase_atoms
                else:
                    if minE > e:
                        minE = e
                        minEConformerID = conformerId
                        minE_ase_atoms = ase_atoms
                print("%sconf_%d.sdf,%s"%(prefix, conformerId, e), file=file_csv)

            minEConformerIDs.append(minEConformerID)

            # to pick conformer randomly
            picked = False
            while picked is False:
                rndConformerIDs = sample(clustered_confIds, npick)
                if minEConformerID not in rndConformerIDs:
                    picked = True

            picked_confs = [minEConformerID] + rndConformerIDs
            all_picked_confs += picked_confs

        # test
        assert len(minEConformerIDs) == len(cluster_conf_id.keys())
        # close to csv file
        file_csv.close()

        prefix = ""
        PICKED_CONF_DIR = self.WORK_DIR
        if optimization_conf:
            prefix = "opt_"

        PICKED_CONF_DIR = self.WORK_DIR + f"/{prefix}picked_confs"
        if not os.path.exists(PICKED_CONF_DIR):
            os.mkdir(PICKED_CONF_DIR)
        picked_file_csv = open(f"{PICKED_CONF_DIR}/{prefix}picked_confs_energies.csv", "w")
        print("FileName,Energy(eV),EnergyPerAtom(eV)", file=picked_file_csv)

        for i, conformerId  in enumerate(all_picked_confs):
            if optimization_conf:
                e, ase_atoms = self._geomOptimizationConf(mol, conformerId)
            else:
                e, ase_atoms = self._calcSPEnergy(mol, conformerId)
            conf_file_path = "%s/%sconf_%d.sdf"%(PICKED_CONF_DIR, prefix, conformerId)

            #  save optimized structure  with rdkit as sdf
            with Chem.rdmolfiles.SDWriter(conf_file_path) as writer:
                rwmol = self.aseAtoms2rwMol(ase_atoms)
                rwmol.SetProp("Energy", str(e))
                rwmol.SetProp("_Name", f"{prefix}conf_{conformerId}")
                writer.write(rwmol)

            print("%sconf_%d.sdf,%s,%s"%(prefix,
                                         conformerId,
                                         e,
                                         e/ase_atoms.get_number_of_atoms()),
                  file=picked_file_csv)
        picked_file_csv.close()

        # cluster and prune opitimzed confs by RMSD
        if optimization_conf:
            confs_energies = pd.read_csv(f"{PICKED_CONF_DIR}/{prefix}picked_confs_energies.csv")
            #  print(confs_energies)
            cluster_conf = self._getClusterRMSDFromFiles(PICKED_CONF_DIR, rmsd_thresh=opt_prune_rms_thresh)
            if cluster_conf != 0:
                self._pruneOptConfs(cluster_conf, confs_energies, PICKED_CONF_DIR, opt_prune_diffE_thresh)
            else:
                os.rename(f"{PICKED_CONF_DIR}/{confs_energies['FileName'][0]}", f"{PICKED_CONF_DIR}/{prefix}output.sdf")

    def _calcEnergyWithMM(self, mol, conformerId, minimizeIts):
        ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(
            mol,
            rdkit.Chem.AllChem.MMFFGetMoleculeProperties(mol),
            confId=conformerId)
        ff.Initialize()
        ff.CalcEnergy()
        results = {}
        if minimizeIts > 0:
            results["converged"] = ff.Minimize(maxIts=minimizeIts)
        results["energy_abs"] = ff.CalcEnergy()
        return results

    def setG16Calculator(self, label, chk, nprocs, xc, basis, scf, addsec=None, extra=None):
        from ase.calculators.gaussian import Gaussian
        self.optG16 = True

        self.calculator = Gaussian(
            label=label,
            #  chk=chk,
            nprocshared=nprocs,
            xc=xc,
            basis=basis,
            scf=scf,
            addsec=addsec,
            extra=extra,
        )

    def setANI2XCalculator(self):
        import torchani
        import torch
        print("Nuber of CUDA devices: ", torch.cuda.device_count())
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.calculator = torchani.models.ANI2x().to(device).ase()

    def setNequIPCalculator(self, model_path):
        from nequip.ase import NequIPCalculator
        print("Nuber of CUDA devices: ", torch.cuda.device_count())
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.calculator = NequIPCalculator.from_deployed_model(
            model_path=model_path, device=device)

    def _calcSPEnergy(self, mol, conformerId):
        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)

        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)
        ase_atoms.set_calculator(self.calculator)

        return ase_atoms.get_potential_energy(), ase_atoms 

    def calcSPEnergy(self):

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        ase_atoms= self.rwMol2AseAtoms()
        ase_atoms.set_calculator(self.calculator)
        return ase_atoms.get_potential_energy()

    def setOptParams(self, fmax, maxiter):
        self.maxiter = maxiter
        self.fmax = fmax

    def setOptMethod(self, opt_method):
        self.opt_method = opt_method.lower()

    def _getOptMethod(self, ase_atoms):
        if self.opt_method is None or self.opt_method=="lbfgs":
            from ase.optimize import LBFGS
            return LBFGS(ase_atoms)
        elif self.opt_method=="bfgs":
            from ase.optimize import BFGS
            return BFGS(ase_atoms)
        elif self.opt_method=="fire":
            from ase.optimize import FIRE
            return FIRE(ase_atoms)
        elif self.opt_method=="gpmin":
            from ase.optimize import GPMin
            return GPMin(ase_atoms)
        elif self.opt_method=="berny":
            from ase.optimize import Berny
            return Berny(ase_atoms)
        elif self.opt_method=="cg":
            from ase.optimize.sciopt import SciPyFminCG
            return SciPyFminCG(ase_atoms)
        elif self.opt_method=="newtonraphson":
            from ase_optmizer_newton_raphson import NewtonRaphson
            return NewtonRaphson(ase_atoms)

    def _geomOptimizationConf(self, mol, conformerId):
        from ase.calculators.gaussian import GaussianOptimizer, Gaussian

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)


        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)

        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)

        if self.optG16:
            dyn =  GaussianOptimizer(ase_atoms, self.calculator)
            dyn.run(steps=self.maxiter)
        else:
            ase_atoms.set_calculator(self.calculator)
            dyn = self._getOptMethod(ase_atoms)
            dyn.run(fmax=self.fmax, steps=self.maxiter)

        return ase_atoms.get_potential_energy(), ase_atoms

    def geomOptimization(self, fix_heavy_atoms=False):
        from ase.calculators.gaussian import GaussianOptimizer

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)

        ase_atoms = self.rwMol2AseAtoms()
        if fix_heavy_atoms:
            from ase.constraints import FixAtoms
            c = FixAtoms(indices=[atom.index for atom in ase_atoms if atom.symbol != 'H'])
            ase_atoms.set_constraint(c)

        if self.optG16:
            dyn =  GaussianOptimizer(ase_atoms, self.calculator)
            dyn.run(steps=self.maxiter)
        else:
            ase_atoms.set_calculator(self.calculator)
            #  self.dyn = LBFGS(ase_atoms)
            dyn = self._getOptMethod(ase_atoms)
            dyn.run(fmax=self.fmax, steps=self.maxiter)

        self.rw_mol = self.aseAtoms2rwMol(ase_atoms)
        return ase_atoms.get_potential_energy()

    def _rwConformer2AseAtoms(self, mol, conformerId):

        mol = mol.GetConformer(conformerId)

        atom_species = [atom.GetAtomicNum() for atom in mol.GetOwningMol().GetAtoms()]
        positions = mol.GetPositions()

        return Atoms(atom_species, positions)

    def rwMol2AseAtoms(self):

        atom_species = [atom.GetAtomicNum() for atom in self.rw_mol.GetAtoms()]

        conf = self.rw_mol.GetConformer()
        positions = [conf.GetAtomPosition(i) for i in range(len(atom_species))]

        return Atoms(atom_species, positions)

    def obMol2AseAtoms(self, ob_mol):
        from ase import Atom
        ase_atoms = Atoms()
        for i in range(ob_mol.NumAtoms()):
            obatom = ob_mol.GetAtom(i + 1)
            ase_atoms.append(Atom(obatom.GetAtomicNum(),
                              [obatom.GetX(),
                               obatom.GetY(),
                               obatom.GetZ()]
                             ))
        return ase_atoms

    def aseAtoms2rwMol(self, ase_atoms):

        write("tmp.pdb", ase_atoms)

        rd_mol = Chem.rdmolfiles.MolFromPDBFile("tmp.pdb", sanitize=True, removeHs=False)
        self._rmFileExist("tmp.pdb")

        try:
            return AllChem.AssignBondOrdersFromTemplate(self.rw_mol, rd_mol)
        except:
            print("Warnings: Can not assign bond borders!")
            return rd_mol


    def writeAseAtoms(self, file_path):
        ase_atoms = self.rwMol2AseAtoms()

        # write mol to xyz file by ase
        write(file_path, ase_atoms)


