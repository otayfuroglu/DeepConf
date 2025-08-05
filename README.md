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
### Run

```bash
bash runDeepConf.sh
```


