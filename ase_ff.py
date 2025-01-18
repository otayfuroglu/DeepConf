from ase import Atoms
from ase.calculators.calculator import Calculator, all_properties
from ase.units import kJ, mol
import openbabel
import numpy as np

class OpenBabelFFCalculator(Calculator):
    """
    A custom ASE calculator that uses Open Babel's force field library
    to calculate energy and forces.
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, forcefield="uff", **kwargs):
        super().__init__(**kwargs)
        self.forcefield = forcefield

    def calculate(self, atoms=None, properties=all_properties, system_changes=all_properties):
        super().calculate(atoms, properties, system_changes)

        # Convert ASE atoms to Open Babel molecule
        obmol = self._atoms_to_obmol(atoms)

        # Set up the force field
        obff = openbabel.OBForceField.FindForceField(self.forcefield)
        if obff is None:
            raise ValueError(f"Force field {self.forcefield} not found in Open Babel.")

        obff.Setup(obmol)

        # Compute energy
        energy = obff.Energy() * kJ / mol


        # Compute forces
        forces = np.zeros((obmol.NumAtoms(), 3))
        #  obff.GetGradient()  # Ensure gradients are updated
        for i in range(obmol.NumAtoms()):
            atom = obmol.GetAtom(i + 1)  # Open Babel uses 1-based indexing
            force = obff.GetGradient(atom)  # Ensure gradients are updated
            forces[i:] = np.array([force.GetX(), force.GetY(), force.GetZ()]) * kJ / mol

        self.results['energy'] = energy
        self.results['forces'] = forces

    #  @staticmethod
    #  def ase_to_obmol(atoms):
    #      """Convert ASE Atoms object to Open Babel molecule."""
    #      obmol = openbabel.OBMol()
    #      obconversion = openbabel.OBConversion()
    #      obconversion.SetInAndOutFormats("xyz", "xyz")
    #
    #      # Write ASE atoms to XYZ format
    #      xyz_data = atoms.get_positions()
    #      symbols = atoms.get_chemical_symbols()
    #      natoms = len(symbols)
    #      xyz_string = f"{natoms}\n\n"
    #      for symbol, coord in zip(symbols, xyz_data):
    #          xyz_string += f"{symbol} {coord[0]} {coord[1]} {coord[2]}\n"
    #
    #      # Read XYZ data into OBMol
    #      obconversion.ReadString(obmol, xyz_string)
    #      return obmol

    @staticmethod
    def _atoms_to_obmol(atoms):
        """Convert an Atoms object to an OBMol object.
        Parameters
        ==========
        Input
            atoms: Atoms
        Return
            obmol: OBMol
        """

        obmol = openbabel.OBMol()
        for atom in atoms:
            a = obmol.NewAtom()
            a.SetAtomicNum(int(atom.number))
            a.SetVector(atom.position[0], atom.position[1], atom.position[2])
        # Automatically add bonds to molecule
        obmol.ConnectTheDots()
        obmol.PerceiveBondOrders()
        return obmol
