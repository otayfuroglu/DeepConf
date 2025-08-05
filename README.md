# Welcome to DeepConf!

DeepConf is a new tool leveraging ANI-ML potentials to explore low-energy conformational space of small molecules. It is particularly useful for generating bioactive conformers using machine learning-based geometry optimization. 
The tool has been tested with RDKit, TorchANI, and Gaussian16, and is suitable for preparing input geometries for docking or MD simulations.

DeepConf is an automated pipeline for generating and optimizing ligand conformers using ANI (ML), G16 (QM), or UFF (MM). It is designed to support high-throughput ligand processing and is inspired by workflows like DeepQM.
Please see the documentation below.

📚 Cite Us
DeepConf: Leveraging ANI-ML Potentials for Exploring Local Minima with Application to Bioactive Conformations
Omer Tayfuroglu, Irem Nur Zengin, Mehmet Serdar Koca, Abdulkadir Kocak
Journal of Chemical Information and Modeling, 2024, DOI: 10.1021/acs.jcim.4c02053

✅ Requirements
• Conda with Python ≥ 3.8
• RDKit
• ASE
• TorchANI
• PyTorch (≥ 0.4.1)
• Optional: Gaussian16 (G16) if QM optimization is required

💡 We recommend using a GPU for faster optimization with ANI-ML potentials.

⚙️ Installation

conda create --name DeepConf python=3.8
conda activate DeepConf
conda install -c conda-forge numpy tqdm rdkit ase pytorch torchani
conda install -c conda-forge numpy pandas tqdm ase pytorch torchani dftd3-python

# For G16 support:
Ensure Gaussian16 is installed and g16 command is available in your $PATH

Download the code:
git clone https://github.com/otayfuroglu/DeepConf.git

🚀 How To Use 
runConfGen.py is the main script for generating and optimizing conformers.
runConfGen.sh is for run runConfGen.py with parameters

1. Create a working directory:
mkdir <workdir>   #Make a folder called workdir.
cd <workdir>	 #Enter the <workdir> file.

mkdir <structures>  #Inside it, create a structures folder.

2. Copy ligand files:
cp /path_your_files/*.sdf /path/workdir/structures #Add all your .sdf or .pdb files into workdir/structures.

3. Copy the run script:
cp /path/DeepConf/runDeepConf.sh /path/workdir #Copy runDeepConf.sh from the DeepConf folder into workdir.

4. Edit the script:

In run_deepconf.sh, set the correct paths to your Python binary and DeepConf directory.
ligPrep_DIR="/home/path/DeepConf/" #Go to DeepConf directory and run "pwd" and write the path here
PYTHON_DIR="$HOME/path/bin" #Write the path when you run the code "which python" until the bin
struct_dir=structures  #Write your directory that consist structure files

5. Run the workflow:
Go to workdir, and run:
bash runDeepConf.sh

🔧 Parameters
structure_dir: Directory containing ligand structures (.sdf, .mol2, .xyz)
ignore_hydrogen: Whether to add hydrogen atoms (Yes/No)
calculator_type: Choose one of: ani2x, g16, uff
optimization_method: Method for optimization: BFGS, LBFGS, etc.
optimization_conf: Whether to optimize conformers after generation (Yes/No)
optimization_lig: Whether to optimize the ligand globally (Yes/No)
pre_optimization_lig: Whether to optimize before conformer generation (Yes/No)
genconformer: Enable conformer generation (Yes/No)
nprocs: Number of processors to use
thr_fmax: Force convergence threshold (default: 0.05)
maxiter: Maximum optimization steps
ETKDG: Use ETKDG method for RDKit conformer generation (Yes/No)
num_conformers: Number of conformers to generate
max_attempts: Maximum attempts for embedding conformers
prune_rms_thresh: RMSD threshold for pruning conformers
opt_prune_rms_thresh: Post-optimization RMSD pruning threshold
opt_prune_diffE_thresh: Post-optimization energy difference threshold (eV)
nfold: Factor to multiply conformers before pruning
npick: How many conformers to select after pruning
nscale: Energy scaling factor

▶️ Run Script
python runConfGen.py \
ligands yes ani2x BFGS yes yes yes yes \
8 0.05 500 yes 50 100 0.2 0.2 0.001 2 2 2

📂 Output Files
Each ligand gets its own directory:
ligand1/
├── pre_opt_ligand1.sdf
├── pre_opt_ligand1_energy.txt
├── minE_conformer.sdf
├── global_ligand1.sdf
├── global_ligand1_energy.txt

Additionally:
• timings.csv – Runtime per ligand
• failed_files.csv – Files that failed during execution

📌 Summary
DeepConf provides:
• ANI2x/ML optimization
• RDKit-based conformer generation
• Pruning with geometric and energetic thresholds
• UFF and G16 compatibility
• Easy integration into computational pipelines

🧾 License
MIT License
<img width="462" height="693" alt="image" src="https://github.com/user-attachments/assets/3d5481b2-65cd-49d0-a73d-15d8c01123ff" />
