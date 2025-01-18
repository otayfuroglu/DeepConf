#! /usr/bin/env bash

ligPrep_DIR="/arf/scratch/otayfuroglu/ConfGen/"
PYTHON_DIR="$HOME/miniconda3/bin"

# struct_dir=./test_akocak
struct_dir=./geoms_several

# adding hydrogen if missing (yes/no) if yes, constraint heavy atoms and minimize hydrogens.
add_hydrogen=no

# set optizetions methods whichs availble in ASE (BFGS, LBFGS, GPMin, FIRE, Berny)
optimization_method=LBFGS
# optimization_method=newtonraphson
# optimization_method=gpmin

# optimization ligand if desired before conformer generation (yes/no)
pre_optimization_lig=no

# generate conformer if desired (yes/no)
genconformer=yes

sample_md=yes

#configuration for conformer generator parameters yes or no. No uses RDKit. if ETKG yes, max_attempts and prune_rms_threshold are redundant (not used).
ETKDG=no
num_conformers=1000
max_attempts=100000

# This is RMSD threshold for generation of the conformers at the beginning used in ETKG or torsion points by rdkit
prune_rms_thresh=0.0005

# this is RMSD threshold and  used for f-clustering after optimization of the picked conformers. Should be 0.1-0.5Angstrom usuually.
opt_prune_rms_thresh=0.5

# this is RMSD threshold and  used for f-clustering after optimization of the picked conformers. eV
opt_prune_diffE_thresh=0.01

# select caclulator type (ani2x/g16) for optimization conf
# caculator_type=g16
caculator_type="ani2x"

# perform geometry optimization for conformers if desired (yes/no)
optimization_conf=yes

# perform geometry optimization for orginal ligand if desired (yes/no)
optimization_lig=no

# set number of procssors for g16 calcultor (default=all cpu)
nprocs=1

# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.01

#maximum iteration for optimization
maxiter=50000

# number of fold for conformer generation
nfold=2
# to pick randomly extra conformer
npick=0

# to scale number of conformers
nscale=10

$PYTHON_DIR/python $ligPrep_DIR/runConfGen.py $struct_dir $add_hydrogen $caculator_type\
	$optimization_method $optimization_conf $optimization_lig $pre_optimization_lig $genconformer\
       	$nprocs $thr_fmax $maxiter $sample_md $ETKDG $num_conformers $max_attempts $prune_rms_thresh $opt_prune_rms_thresh $opt_prune_diffE_thresh $nfold $npick $nscale


