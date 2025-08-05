# 🧠 DeepConf

**DeepConf** is a tool for exploring low-energy conformations of small molecules using machine learning potentials. It leverages ANI-ML, RDKit, TorchANI, and optionally Gaussian16 to generate bioactive conformers for downstream tasks like docking or MD simulations.

---

## 📘 Citation

> **DeepConf: Leveraging ANI-ML Potentials for Exploring Local Minima with Application to Bioactive Conformations**  
> Omer Tayfuroglu, Irem Nur Zengin, Mehmet Serdar Koca, Abdulkadir Kocak  
> *Journal of Chemical Information and Modeling*, 2024  
> DOI: [10.1021/acs.jcim.4c02053](https://doi.org/10.1021/acs.jcim.4c02053)

---

## ✅ Requirements

- Conda with Python ≥ 3.8
- [RDKit](https://www.rdkit.org/)
- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [TorchANI](https://aiqm.github.io/torchani/)
- PyTorch ≥ 0.4.1  
- *(Optional)*: [Gaussian16](https://gaussian.com/g16/) for QM-level optimization

💡 **Tip:** Use a GPU for faster ANI optimization.

---

## ⚙️ Installation

```bash
conda create --name DeepConf python=3.8
conda activate DeepConf
conda install -c conda-forge numpy pandas tqdm rdkit ase pytorch torchani dftd3-python
```

### For G16 support:
Ensure Gaussian16 is installed and g16 command is available in your $PATH

```bash
git clone https://github.com/otayfuroglu/DeepConf.git
```

## 🚀 How to Use

### Setup

```bash
mkdir workdir
cd workdir
mkdir structures
```

Copy your ligands:

```bash
cp /path/to/*.sdf structures/
```

Copy and configure the script:

```bash
cp /path/to/DeepConf/runDeepConf.sh .
```

Edit runDeepConf.sh:

```bash
ligPrep_DIR="/path/to/DeepConf"
PYTHON_DIR="/path/to/conda/envs/DeepConf/bin"
struct_dir="structures"
```

### Run

```bash
bash runDeepConf.sh
```

🔧 Parameters (for runConfGen.py)

| Parameter                | Description                             |
| ------------------------ | --------------------------------------- |
| `structure_dir`          | Input directory of ligands              |
| `ignore_hydrogen`        | Add hydrogens? (`yes`/`no`)             |
| `calculator_type`        | `ani2x`, `g16`, `uff`                   |
| `optimization_method`    | E.g. `BFGS`, `LBFGS`                    |
| `optimization_conf`      | Optimize conformers (`yes`/`no`)        |
| `optimization_lig`       | Global ligand optimization (`yes`/`no`) |
| `pre_optimization_lig`   | Pre-conformer ligand opt (`yes`/`no`)   |
| `genconformer`           | Generate conformers (`yes`/`no`)        |
| `nprocs`                 | Number of CPUs                          |
| `thr_fmax`               | Force threshold (default: 0.05)         |
| `maxiter`                | Max optimization steps                  |
| `ETKDG`                  | Use RDKit ETKDG (`yes`/`no`)            |
| `num_conformers`         | Number to generate                      |
| `max_attempts`           | Max embedding attempts                  |
| `prune_rms_thresh`       | RMSD pruning threshold                  |
| `opt_prune_rms_thresh`   | Post-opt RMSD threshold                 |
| `opt_prune_diffE_thresh` | Energy diff threshold (eV)              |
| `nfold`                  | Pre-prune conformer multiplier          |
| `npick`                  | Final conformer count                   |
| `nscale`                 | Energy scale factor                     |

### ▶️ Example Usage
```bash
python runConfGen.py \
ligands yes ani2x BFGS yes yes yes yes \
8 0.05 500 yes 50 100 0.2 0.2 0.001 2 2 2
```

### 📂 Output
Each ligand generates its own folder with outputs:

```bash
ligand1/
├── pre_opt_ligand1.sdf
├── pre_opt_ligand1_energy.txt
├── minE_conformer.sdf
├── global_ligand1.sdf
├── global_ligand1_energy.txt
```

Also:
- timings.csv – Processing time per ligand
- failed_files.csv – Ligands that failed during execution

# 📌 Summary
- ⚡ ANI-ML (ANI2x) conformer optimization
- 🧪 RDKit + ETKDG conformer generation
- 🧠 Pruning based on RMSD and energy thresholds
- 🧬 UFF and G16 support
- 🔗 Easy integration into pipelines
