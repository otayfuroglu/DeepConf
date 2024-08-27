
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np

from ase.io import read, write
from core_conf import confGen
import argparse
import os, sys, shutil
import multiprocessing
from itertools import product
import time

nprocs_all = int(multiprocessing.cpu_count())



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("structure_dir", type=str)
parser.add_argument("ignore_hydrogen", nargs="?", default="No") # args for bool
parser.add_argument("calculator_type", type=str)
parser.add_argument("optimization_method", nargs="?", default="No") # args for bool
parser.add_argument("optimization_conf", nargs="?", default="No") # args for bool
parser.add_argument("optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("pre_optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("genconformer", nargs="?", default="No") # args for bool
parser.add_argument("nprocs", type=int, default=nprocs_all)
parser.add_argument("thr_fmax", type=float, default=0.05)
parser.add_argument("maxiter", type=int, default=500)

parser.add_argument("ETKDG", nargs="?", default="No") # args for bool
parser.add_argument("num_conformers", type=int, default=50)
parser.add_argument("max_attempts", type=int, default=100)
parser.add_argument("prune_rms_thresh", type=float, default=0.2)
parser.add_argument("opt_prune_rms_thresh", type=float, default=0.2)
parser.add_argument("opt_prune_diffE_thresh", type=float, default=0.001)
parser.add_argument("nfold", type=int, default=2)
parser.add_argument("npick", type=int, default=2)


def calcFuncRunTime(func):
    import time
    def wrapper(*args, **kwargs):
        s_time = time.time()
        func(*args, **kwargs)
        print(f"Funtion {func.__name__} executed in {(time.time()-s_time)/60:.4f} m")
    return wrapper


def getBoolStr(string):
    string = string.lower()
    if "true" in string or "yes" in string:
        return True
    elif "false" in string or "no" in string:
        return False
    else:
        print("%s is bad input!!! Must be Yes/No or True/False" %string)
        sys.exit(1)


def setG16calculator(lig, file_base, label, WORK_DIR):
    lig.setG16Calculator(
            label="%s/g16_%s/%s"%(WORK_DIR, label, file_base),
            #  chk="%s.chk"%file_base,
            chk="",
            nprocs=nprocs,
            xc="HF",
            basis="3-21g",
            scf="maxcycle=100",
            extra="nosymm",
    )
    return lig


def setGenConformers(lig, out_file_path, mmCalculator):
    trial = lig.n_trial
    while trial <= 3:
        try:
            lig.genGonformers(
                file_path=out_file_path,
                numConfs=num_conformers,
                ETKDG=ETKDG,
                maxAttempts=max_attempts,
                pruneRmsThresh=prune_rms_thresh,
                mmCalculator=mmCalculator,
                optimization_conf=optimization_conf,
                opt_prune_rms_thresh=opt_prune_rms_thresh,
                opt_prune_diffE_thresh=opt_prune_diffE_thresh,
                nfold=nfold,
                npick=npick,)
        except:
            print(f"Trail {trial} failed, attempting new one... ")
            lig.increaseTrilNum()
            trial = lig.n_trial
            setGenConformers(lig, out_file_path, mmCalculator)
        finally:
            return lig
    else:
        print(f"{trial -1} attempts failed, Skipping...")
        return None


#  @calcFuncRunTime
def runConfGen(file_name):
    "Starting ligand preparetion process... "
    mol_path= "%s/%s"%(structure_dir, file_name)

    file_base = file_name.split(".")[0]
    #create destination directory
    WORK_DIR = file_base
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.mkdir(WORK_DIR)

    #Flags
    # default mm calculator set to False
    mmCalculator=False
    # default adding H is False
    addH = False

    # if desire adding H by openbabel
    prefix = ""
    if ignore_hydrogen:
        addH = True
        prefix += "addH_"
    if optimization_lig or optimization_conf:
        prefix += "opt_"

    # initialize confGen
    lig = confGen(mol_path, addH, WORK_DIR)
    lig.setOptMethod(optimization_method)
    #  lig.writeRWMol2File("test/test.xyz")

    if "ani2x" in calculator_type.lower():
        lig.setANI2XCalculator()
    elif "g16" in calculator_type.lower():
        lig = setG16calculator(lig, file_base, label="calculation", WORK_DIR=WORK_DIR)
    elif "uff" in calculator_type.lower():
        if optimization_conf:
            print("UFF calculator not support optimization")
            sys.exit(1)
        else:
            mmCalculator=True

    # set optimizetion parameters
    lig.setOptParams(fmax=thr_fmax, maxiter=args.maxiter)

    if pre_optimization_lig:
        print("Pre-Optimization process.. before generations")
        ase_atoms, _ = lig.geomOptimization()
        write(f"pre_opt_{file_base}.xyz", ase_atoms)


    if genconformer:
        out_file_path="%s/%sminE_conformer.sdf"%(WORK_DIR, prefix)
        lig = setGenConformers(lig, out_file_path, mmCalculator)
        if lig is None:
            return None

        print("Conformer generation process is done")
        if not optimization_conf and optimization_lig:
            print("Optimization for minumum energy conformer")
            _, e = lig.geomOptimization()
            _, e_file = open("%s/global_%s%s_energy.txt"%(WORK_DIR, prefix, file_base) , "w")
            #  lig.geomOptimization()

    else:
        out_file_path="%s/global_%s%s.sdf"%(WORK_DIR, prefix, file_base)
        # geometry optimizaton for ligand
        if  optimization_lig:
            _, e = lig.geomOptimization()
            e_file = open("%s/global_%s%s_energy.txt"%(WORK_DIR, prefix, file_base) , "w")
            print(e, " eV", file=e_file)
            #  ase_atoms = lig.rwMol2AseAtoms()

    # write minimun energy conformer to sdf file
    lig.writeRWMol2File(out_file_path)

    #  for the bug of reading sfd file which have charges in ase
    try:
        atoms = read(out_file_path)
    except:
        out_file_path="%s/%s%s.xyz"%(WORK_DIR, prefix, file_base)
        lig.writeRWMol2File(out_file_path)
        atoms = read(out_file_path)


if __name__ == "__main__":
    args = parser.parse_args()
    structure_dir = args.structure_dir
    calculator_type = args.calculator_type

    optimization_method = args.optimization_method

    optimization_conf = getBoolStr(args.optimization_conf)
    optimization_lig = getBoolStr(args.optimization_lig)
    pre_optimization_lig = getBoolStr(args.pre_optimization_lig)
    genconformer = getBoolStr(args.genconformer)
    ignore_hydrogen = getBoolStr(args.ignore_hydrogen)
    ETKDG = getBoolStr(args.ETKDG)

    nprocs = args.nprocs
    thr_fmax = args.thr_fmax
    maxiter = args.maxiter

    #get conformer generator parameters
    num_conformers = args.num_conformers
    max_attempts = args.max_attempts
    prune_rms_thresh = args.prune_rms_thresh
    opt_prune_rms_thresh = args.opt_prune_rms_thresh
    opt_prune_diffE_thresh = args.opt_prune_diffE_thresh
    nfold = args.nfold
    npick = args.npick

    file_names = [item for item in os.listdir(structure_dir) if not item.startswith(".")]
    failed_csv = open("failed_files.csv", "w")
    failed_csv.write("FileNames,\n")

    fl_timing = open("timings.csv", "w")
    print("FileName,Time(min.)", file=fl_timing)
    for file_name in file_names:
        file_base = file_name.split(".")[0]
        #  try:
        s_time = time.time()
        result = runConfGen(file_name)
        if result is None:
            continue

        print(file_name, ",", "%.4f"%((time.time()-s_time)/60), file=fl_timing)
        #  except:
        #      print("Error for %s file !!! Skipping...")
        #      failed_csv.write(file_name+",\n")
        #  break
    fl_timing.close()
    failed_csv.close()

