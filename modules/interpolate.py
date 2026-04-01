import numpy as np
from typing import List, Tuple, Optional, Union
import itertools


from geometry import GeometryOps, Quaternion
from molecule_class import Molecule
from psi4_interface import Psi4Calculator, Config


class Interpolator:
    """  
    A class for interpolation between two molecular minimum structures.
    """

    @staticmethod
    def slerp_interpolate(
            mol_a: "Molecule",
            mol_b: "Molecule",
            submolecule_indices: Optional[List[List[int]]] = None,
            lambdas: Optional[List[float]] = None,
            n_points: int = 3,
            calculator: Optional[Psi4Calculator] = None
        ) -> Union[List["Molecule"], Tuple[List["Molecule"], List[float]]]:
        """
        Generate interpolated structures between two molecules using rigid-body interpolation per submolecule.
        If a calculator is provided, it evaluates and returns the energies across the path.
        """
        if len(mol_a.coordinates) != len(mol_b.coordinates):
            raise ValueError("Molecules must have the same number of atoms for interpolation")
        
        if lambdas is None:
            # Angepasst, um Start (0.0) und Ende (1.0) mit einzuschließen
            lambdas = np.linspace(0, 1, n_points)
            
        if submolecule_indices is None:
            submolecule_indices = [list(range(len(mol_a.coordinates)))] 
            
        n_atoms = len(mol_a.coordinates)
        coords_a = mol_a.coordinates.copy()
        coords_b = mol_b.coordinates.copy()

        q_identity = Quaternion(1.0, 0.0, 0.0, 0.0)
        submol_data = []

        # Process each submolecule independently in absolute space
        for submol_idx in submolecule_indices:
            idx = np.array(submol_idx)
            submol_masses = mol_a.masses[idx]
            mass_col = submol_masses[:, np.newaxis] if submol_masses.ndim == 1 else submol_masses
            
            # 1. Compute COM in absolute frames A and B
            com_a = np.sum(mass_col * coords_a[idx], axis=0) / np.sum(submol_masses)
            com_b = np.sum(mass_col * coords_b[idx], axis=0) / np.sum(submol_masses)

            # 2. Extract local coordinates relative to respective COMs
            local_a = coords_a[idx] - com_a
            local_b = coords_b[idx] - com_b

            # 3. Calculate local rotation from A's orientation to B's orientation
            R_local = GeometryOps.kabsch_rotation(local_b, local_a)
            q_local = Quaternion.from_rotation_matrix(R_local)
            
            submol_data.append({
                "idx": idx,
                "com_a": com_a,
                "com_b": com_b,
                "local_a": local_a,
                "q_local": q_local
            })
        
        interpolated = []
        for lam in lambdas:
            coords_new = np.zeros((n_atoms, 3))

            for sd in submol_data:
                # Linear interpolation of Center of Mass
                com_t = (1 - lam) * sd["com_a"] + lam * sd["com_b"]

                # Slerp for smooth rotation
                q_t = Quaternion.slerp(q_identity, sd["q_local"], lam)
                R_t = q_t.to_rotation_matrix()

                # Rotate State A's rigid local geometry and move to new COM
                coords_new[sd["idx"]] = (R_t @ sd["local_a"].T).T + com_t

            new_mol = mol_a.copy()
            new_mol.coordinates = coords_new
            interpolated.append(new_mol)
            
        # Energieberechnung über das Interface
        if calculator is not None:
            print(f"\nEvaluating energies for {len(interpolated)} path structures...")
            results = calculator.batch_single_point_parallel(interpolated)
            energies = [res.energy if res.success else None for res in results]
            return interpolated, energies

        return interpolated

    @staticmethod
    def crossover_interpolate(
            mol_a: "Molecule",
            mol_b: "Molecule",
            submolecule_indices: Optional[List[List[int]]] = None,
            n_permutations: int = -1,
        ) -> Union[List["Molecule"], Tuple[List["Molecule"], List[float]]]:
        """
        Generates Crossovers between different Minimum Structures

        Algorithm:
        + Take two minima A and B with submolecules [[0,1,2], [3,4,5]]
        + Pick [[0,1,2]] from A and [[3,4,5]] from B combine this to a new structure
        + Do a similar permutation with [[3,4,5]] from A and [[0,1,2]] from B
        + Automate this

        Args:
            mol_a: First molecule (minimum structure)
            mol_b: Second molecule (minimum structure)
            submolecule_indices: List of lists, where each inner list contains the atom indices for a submolecule. If None, the entire molecule is treated as one submolecule.
            n_permutations: Number of random crossover permutations to generate. If -1, generates all
        """ 
        if len(mol_a.coordinates) != len(mol_b.coordinates):
            raise ValueError("Molecules must have the same number of atoms for interpolation")
        if submolecule_indices is None:
            submolecule_indices = [list(range(len(mol_a.coordinates)))] # default: whole molecule as one submolecule


        n_atoms = len(mol_a.coordinates)
        coords_a = mol_a.coordinates.copy()
        coords_b = mol_b.coordinates.copy()
        crossover_structures = []

        # 0 means to take from mol_a, 1 means to take from mol_b
        n_submols = len(submolecule_indices)
        all_combinations = list(itertools.product([0,1], repeat=n_submols))


        # Filter out all (0) or all (1)
        crossover_combinations = [comb for comb in all_combinations if 0 in comb and 1 in comb]

        # Deterministic Case
        if n_permutations >= len(crossover_combinations) or n_permutations == -1:
            for comb in crossover_combinations:
                new_coords = np.zeros((n_atoms, 3)) 
                # Map the submolecules based on current combinations
                for sub_idx, source_is_b in zip(submolecule_indices, comb):
                    idx_array = np.array(sub_idx)
                    if source_is_b:
                        new_coords[idx_array] = coords_b[idx_array]
                    else:
                        new_coords[idx_array] = coords_a[idx_array]
                
                # Create the new crossover molecule
                new_mol = mol_a.copy()
                new_mol.coordinates = new_coords
                crossover_structures.append(new_mol)
        else:
            # Make a Random choice of n_permutations from the possible combinations
            selected_combinations = np.random.choice(len(crossover_combinations), size=n_permutations, replace=False)
            for idx in selected_combinations:
                comb = crossover_combinations[idx]
                new_coords = np.zeros((n_atoms, 3)) 
                # Map the submolecules based on current combinations
                for sub_idx, source_is_b in zip(submolecule_indices, comb):
                    idx_array = np.array(sub_idx)
                    if source_is_b:
                        new_coords[idx_array] = coords_b[idx_array]
                    else:
                        new_coords[idx_array] = coords_a[idx_array]
                
                # Create the new crossover molecule
                new_mol = mol_a.copy()
                new_mol.coordinates = new_coords
                crossover_structures.append(new_mol)
    


        return crossover_structures  


            
    

if __name__ == "__main__":
    h2o_1 = """
    6
    Structure 10 | Energy: -152.473796 Hartree
    O -0.755372 -0.921122 -4.005642
    H -0.493420 -1.806856 -3.730554
    H -0.616336 -0.411307 -3.192303
    O -0.206197 0.784773 -1.718121
    H -0.965669 1.375017 -1.624367
    H 0.450896 1.356129 -2.137662 
    """
    h2o_2 = """
    6
    Structure 6 | Energy: -152.461845 Hartree
    O 0.933610 -1.962299 1.931146
    H 1.019776 -1.537261 2.792539
    H 1.182413 -2.871120 2.135874
    O -3.879426 -2.708703 -2.287702
    H -3.263914 -2.748691 -3.028802
    H -4.191267 -1.797892 -2.341029
    """
    submolecule_indices = [[0, 1, 2], [3, 4, 5]]
    mol1 = Molecule.from_xyz(h2o_1)
    mol2 = Molecule.from_xyz(h2o_2)
    
    # Richte den Psi4Calculator für schnelle Berechnungen ein
    config = Config(method="hf", basis="6-311g(d,p)", memory="2GB")
    calculator = Psi4Calculator(config=config, verbose=False)
    
    interpolator = Interpolator()
    n_interpolation_points = 20

    # Try out crossinterpolation
    
    crossover_mols = interpolator.crossover_interpolate(mol1, mol2, submolecule_indices=submolecule_indices, n_permutations=2)
    xyz_lines = []
    for idx, mol in enumerate(crossover_mols):
        mol.name = f"Crossover_Structure_{idx}"
        xyz_lines.append(mol.to_xyz_string())
    with open("crossover_structures.xyz", "w") as f:
        f.write("\n".join(xyz_lines))
    print("Saved -> crossover_structures.xyz")




    interpolated_mols, path_energies = interpolator.slerp_interpolate(
        mol1, mol2, 
        submolecule_indices=submolecule_indices, 
        n_points=n_interpolation_points,
        calculator=calculator
    )

    # 1. Speichere xyz Trajektorie und passe Datei-Header mit der Energie an
    xyz_lines = []
    for idx, (mol, energy) in enumerate(zip(interpolated_mols, path_energies)):
        if energy is not None:
            mol.name = f"Interpolated_Step_{idx} | Energy: {energy:.6f} Hartree"
        else:
            mol.name = f"Interpolated_Step_{idx} | Energy: FAILED"
        xyz_lines.append(mol.to_xyz_string())

    with open("interpolated_trajectory.xyz", "w") as f:
        f.write("\n".join(xyz_lines))
    print("Saved -> interpolated_trajectory.xyz")

    # 2. Generiere Plot des Energiepfads
    import matplotlib.pyplot as plt

    valid_indices = [i for i, e in enumerate(path_energies) if e is not None]
    valid_energies = [e for i, e in enumerate(path_energies) if e is not None]

    if valid_energies:    
       
        
        plt.figure(figsize=(8, 5))
        plt.plot(valid_indices, valid_energies, marker='o', linestyle='-', color='#1f77b4')
        plt.xlabel("Interpolation Step")
        plt.ylabel("Relative Energy (kcal / mol)")
        plt.title("SLERP Rigid-Body Interpolation Energy Profile")
        plt.grid(True, linestyle="--", alpha=0.7)
        
        plot_path = "interpolation_energy_profile.png"
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        print(f"Saved -> {plot_path}")
    else:
        print("Energy calculation failed for all structures. No plot generated.")
