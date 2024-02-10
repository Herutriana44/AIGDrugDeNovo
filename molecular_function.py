import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqUtils import IsoelectricPoint as IP
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
from scipy.constants import e, epsilon_0
from scipy.constants import Boltzmann
import datetime
import random
import string
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.constants import e
import pandas as pd
from rdkit.Chem import rdMolTransforms
from collections import Counter

class molecular_function:

    def calculate_distance_between_amino_acids(aa1, aa2):
        # Menghitung jarak antara dua pusat massa asam amino
        distance = abs(aa1 - aa2)
        return distance

    def generate_conformer(molecule_smiles):
        mol = Chem.MolFromSmiles(molecule_smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens for a more accurate 3D structure
        conformer = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)  # Generate a conformer
        return mol

    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms

    def calculate_molecular_center_of_mass(smiles):
        molecule = molecular_function.generate_conformer(smiles)
        try:
            if molecule is None:
                return None

            # Hitung pusat massa molekul
            center_of_mass = rdMolTransforms.ComputeCentroid(molecule.GetConformer())
            total_mass = Descriptors.MolWt(molecule)
            # print(center_of_mass)
            # print(total_mass)
            center_of_mass = sum([center_of_mass[i] * total_mass for i in range(len(center_of_mass))]) / total_mass
            # print(total_mass)
            
            return center_of_mass
        except Exception as e:
            print("Error:", str(e))
            return None

    def seq_to_smiles(seq):
        try:
            mol = Chem.MolFromSequence(seq)
            smiles = Chem.MolToSmiles(mol,kekuleSmiles=True)
            return str(smiles)
        except:
            return None
        
    def combine_epitope_with_adjuvant(epitope_sequence, adjuvant_smiles):
        # Konversi epitop ke molekul RDKit
        epitope_molecule = Chem.MolFromSequence(epitope_sequence)
        
        # Konversi SMILES adjuvant menjadi molekul RDKit
        adjuvant_molecule = Chem.MolFromSmiles(adjuvant_smiles)
        
        # Gabungkan epitop dan adjuvant
        combined_molecule = Chem.CombineMols(adjuvant_molecule, epitope_molecule)
        
        return combined_molecule
    
    def calculate_molecular_weight(molecule_smiles):
        try:
            #mol = generate_conformer(molecule_smiles)
            mol = Chem.MolFromSmiles(molecule_smiles)
            if mol is None:
                print("Gagal membaca molekul.")
                return None

            # Menghitung massa molekul
            molecular_weight = Descriptors.MolWt(mol)

            return molecular_weight
        except:
            return None
        
    def MolFromLongSequence(sequence, chunk_size=100):
        """Convert a long sequence into a Mol object by splitting into chunks and combining."""
        # Split the sequence into chunks of specified size
        chunks = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)]
        
        # Convert each chunk to a Mol object using MolFromSmiles
        mols = [Chem.MolFromSequence(chunk) for chunk in chunks if chunk]  # Exclude empty chunks
        
        # Combine the Mols into a single Mol
        combined_mol = Chem.Mol()
        for mol in mols:
            if mol:
                combined_mol = Chem.CombineMols(combined_mol, mol)
        
        return combined_mol
    
    def calculate_amino_acid_center_of_mass(sequence):
        try:
            amino_acid_masses = []
            for aa in sequence:
                try:
                    amino_acid_masses.append(Chem.Descriptors.MolWt(molecular_function.MolFromLongSequence(aa)))
                except:
                    return 0
                    break

            # Hitung pusat massa asam amino
            total_mass = sum(amino_acid_masses)
            center_of_mass = sum(i * mass for i, mass in enumerate(amino_acid_masses, start=1)) / total_mass

            return center_of_mass
        except:
            return 0
        
    def calculate_distance_between_amino_acids(aa1, aa2):
        # Menghitung jarak antara dua pusat massa asam amino
        distance = abs(aa1 - aa2)
        return distance
    
    def seq_to_smiles(seq):
        try:
            mol = molecular_function.MolFromLongSequence(seq)
            smiles = Chem.MolToSmiles(mol,kekuleSmiles=True)
            return str(smiles)
        except:
            return None
        
    def calculate_molecular_center_of_mass(smiles):
        molecule = molecular_function.generate_conformer(smiles)
        try:
            if molecule is None:
                return None

            # Hitung pusat massa molekul
            center_of_mass = rdMolTransforms.ComputeCentroid(molecule.GetConformer())
            total_mass = Descriptors.MolWt(molecule)
            # print(center_of_mass)
            # print(total_mass)
            center_of_mass = sum([center_of_mass[i] * total_mass for i in range(len(center_of_mass))]) / total_mass
            # print(total_mass)
            
            return center_of_mass
        except Exception as e:
            print("Error:", str(e))
            return None

    def calculate_molecular_weight(molecule_smiles):
        try:
            #mol = generate_conformer(molecule_smiles)
            mol = Chem.MolFromSmiles(molecule_smiles)
            if mol is None:
                print("Gagal membaca molekul.")
                return None

            # Menghitung massa molekul
            molecular_weight = Descriptors.MolWt(mol)

            return molecular_weight
        except:
            return None

class ms:
    def __init__(self):
        self.mass_of_argon = 39.948

    @staticmethod
    def attractive_energy(r, epsilon=0.0103, sigma=3.4):
        """
        Attractive component of the Lennard-Jones 
        interactionenergy.
        
        Parameters
        ----------
        r: float
            Distance between two particles (Å)
        epsilon: float 
            Negative of the potential energy at the 
            equilibrium bond length (eV)
        sigma: float 
            Distance at which the potential energy is 
            zero (Å)
        
        Returns
        -------
        float
            Energy of attractive component of 
            Lennard-Jones interaction (eV)
        """
        if r == 0:
            return 0
        return -4.0 * epsilon * np.power(sigma / r, 6)

    @staticmethod
    def repulsive_energy(r, epsilon=0.0103, sigma=3.4):
        """
        Repulsive component of the Lennard-Jones 
        interactionenergy.
        
        Parameters
        ----------
        r: float
            Distance between two particles (Å)
        epsilon: float 
            Negative of the potential energy at the 
            equilibrium bond length (eV)
        sigma: float 
            Distance at which the potential energy is 
            zero (Å)
        
        Returns
        -------
        float
            Energy of repulsive component of 
            Lennard-Jones interaction (eV)
        """
        if r == 0:
            return 0
        return 4 * epsilon * np.power(sigma / r, 12)
    
    @staticmethod
    def lj_energy(r, epsilon=0.0103, sigma=3.4):
        """
        Implementation of the Lennard-Jones potential 
        to calculate the energy of the interaction.
        
        Parameters
        ----------
        r: float
            Distance between two particles (Å)
        epsilon: float 
            Negative of the potential energy at the 
            equilibrium bond length (eV)
        sigma: float 
            Distance at which the potential energy is 
            zero (Å)
        
        Returns
        -------
        float
            Energy of the Lennard-Jones potential 
            model (eV)
        """
        if r == 0:
            return 0
        
        return ms.repulsive_energy(
            r, epsilon, sigma) + ms.attractive_energy(
            r, epsilon, sigma)
    
    @staticmethod
    def coulomb_energy(qi, qj, r):
        """
        Calculation of Coulomb's law.
        
        Parameters
        ----------
        qi: float
            Electronic charge on particle i
        qj: float
            Electronic charge on particle j
        r: float 
            Distance between particles i and j (Å)
            
        Returns
        -------
        float
            Energy of the Coulombic interaction (eV)
        """
        if r == 0:
            return 0
        
        energy_joules = (qi * qj * e ** 2) / (
            4 * np.pi * epsilon_0 * r * 1e-10)
        return energy_joules / 1.602e-19
    
    @staticmethod
    def bonded(kb, b0, b):
        """
        Calculation of the potential energy of a bond.
        
        Parameters
        ----------
        kb: float
            Bond force constant (units: eV/Å^2)
        b0: float 
            Equilibrium bond length (units: Å)
        b: float
            Bond length (units: Å)
        
        Returns
        float
            Energy of the bonded interaction
        """
        
        return kb / 2 * (b - b0) ** 2
    
    @staticmethod
    def lj_force(r, epsilon=0.0103, sigma=3.4):
        """
        Implementation of the Lennard-Jones potential 
        to calculate the force of the interaction.
        
        Parameters
        ----------
        r: float
            Distance between two particles (Å)
        epsilon: float 
            Potential energy at the equilibrium bond 
            length (eV)
        sigma: float 
            Distance at which the potential energy is 
            zero (Å)
        
        Returns
        -------
        float
            Force of the van der Waals interaction (eV/Å)
        """
        if r != 0:
            return 48 * epsilon * np.power(
                sigma, 12) / np.power(
                r, 13) - 24 * epsilon * np.power(
                sigma, 6) / np.power(r, 7)
        else:
            return 0
    @staticmethod
    def init_velocity(T, number_of_particles):
        """
        Initialise the velocities for a series of 
        particles.
        
        Parameters
        ----------
        T: float
            Temperature of the system at 
            initialisation (K)
        number_of_particles: int
            Number of particles in the system
        
        Returns
        -------
        ndarray of floats
            Initial velocities for a series of 
            particles (eVs/Åamu)
        """
        R = np.random.rand(number_of_particles) - 0.5
        return R * np.sqrt(Boltzmann * T / (
            ms.mass_of_argon * 1.602e-19))
    @staticmethod
    def get_accelerations(positions):
        """
        Calculate the acceleration on each particle
        as a  result of each other particle. 
        N.B. We use the Python convention of 
        numbering from 0.
        
        Parameters
        ----------
        positions: ndarray of floats
            The positions, in a single dimension, 
            for all of the particles
            
        Returns
        -------
        ndarray of floats
            The acceleration on each
            particle (eV/Åamu)
        """
        accel_x = np.zeros((positions.size, positions.size))
        for i in range(0, positions.size - 1):
            for j in range(i + 1, positions.size):
                r_x = positions[j] - positions[i]
                rmag = np.sqrt(r_x * r_x)
                force_scalar = ms.lj_force(rmag, 0.0103, 3.4)
                force_x = force_scalar * r_x / rmag
                accel_x[i, j] = force_x / ms.mass_of_argon
                accel_x[j, i] = - force_x / ms.mass_of_argon
        return np.sum(accel_x, axis=0)
    @staticmethod
    def update_pos(x, v, a, dt):
        """
        Update the particle positions.
        
        Parameters
        ----------
        x: ndarray of floats
            The positions of the particles in a 
            single dimension
        v: ndarray of floats
            The velocities of the particles in a 
            single dimension
        a: ndarray of floats
            The accelerations of the particles in a 
            single dimension
        dt: float
            The timestep length
        
        Returns
        -------
        ndarray of floats:
            New positions of the particles in a single 
            dimension
        """
        return x + v * dt + 0.5 * a * dt * dt
    
    @staticmethod
    def update_velo(v, a, a1, dt):
        """
        Update the particle velocities.
        
        Parameters
        ----------
        v: ndarray of floats
            The velocities of the particles in a 
            single dimension (eVs/Åamu)
        a: ndarray of floats
            The accelerations of the particles in a 
            single dimension at the previous 
            timestep (eV/Åamu)
        a1: ndarray of floats
            The accelerations of the particles in a
            single dimension at the current 
            timestep (eV/Åamu)
        dt: float
            The timestep length
        
        Returns
        -------
        ndarray of floats:
            New velocities of the particles in a
            single dimension (eVs/Åamu)
        """
        return v + 0.5 * (a + a1) * dt
    @staticmethod
    def run_md(dt, number_of_steps, initial_temp, x):
        """
        Run a MD simulation.
        
        Parameters
        ----------
        dt: float
            The timestep length (s)
        number_of_steps: int
            Number of iterations in the simulation
        initial_temp: float
            Temperature of the system at 
            initialisation (K)
        x: ndarray of floats
            The initial positions of the particles in a 
            single dimension (Å)
            
        Returns
        -------
        ndarray of floats
            The positions for all of the particles 
            throughout the simulation (Å)
        """
        positions = np.zeros((number_of_steps, 3))
        v = ms.init_velocity(initial_temp, 3)
        a = ms.get_accelerations(x)
        for i in range(number_of_steps):
            x = ms.update_pos(x, v, a, dt)
            a1 = ms.get_accelerations(x)
            v = ms.update_velo(v, a, a1, dt)
            a = np.array(a1)
            positions[i, :] = x
        return positions
    

    @staticmethod
    def lj_force_cutoff(r, epsilon, sigma):
        """
        Implementation of the Lennard-Jones potential 
        to calculate the force of the interaction which 
        is considerate of the cut-off.
        
        Parameters
        ----------
        r: float
            Distance between two particles (Å)
        epsilon: float 
            Potential energy at the equilibrium bond 
            length (eV)
        sigma: float 
            Distance at which the potential energy is 
            zero (Å)
        
        Returns
        -------
        float
            Force of the van der Waals interaction (eV/Å)
        """
        cutoff = 15 
        if r < cutoff:
            return 48 * epsilon * np.power(
                sigma / r, 13) - 24 * epsilon * np.power(
                sigma / r, 7)
        else:
            return 0