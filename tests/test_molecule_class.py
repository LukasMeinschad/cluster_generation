import numpy as np
import pytest

from modules.molecule_class import Molecule

# Fixture is not a test itself but a setup function that can be used by multiple tests
@pytest.fixture
def water_labels():
    return ["O1", "H1", "H2"]

@pytest.fixture
def water_coords():
    return np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.9584, 0.0000, 0.0000],
            [-0.2396, 0.9271, 0.0000],
        ],
        dtype=np.float64,
    )

@pytest.fixture
def water_molecule(water_labels, water_coords):
    return Molecule.from_labels_and_coords(water_labels, water_coords, name="Water")

def test_from_labels_and_coords(water_molecule):
    assert water_molecule.get_number_of_atoms() == 3
    assert water_molecule.coordinates.shape == (3,3)

def test_periodic_lookups():
    assert Molecule.get_atomic_number("O1") == 8
    assert Molecule.get_atomic_number("H1") == 1
    assert Molecule.get_covalent_radius("O1") > 0
    assert Molecule.get_vdw_radius("H1") > 0

def test_fragmentation_returns_at_least_one_fragment():
    labels = [f"C{i+1}" for i in range(10)]
    coords = np.random.rand(10, 3) * 20.0
    mol = Molecule.from_labels_and_coords(labels, coords)
    frags = mol.fragment_by_connectivity()
    assert len(frags) >= 1


@pytest.mark.parametrize("n_points", [1000, 5000, 20000])
def test_mc_volume_positive(water_molecule, n_points):
    vol = water_molecule.compute_volume_mc(n_points=n_points)
    assert vol > 0.0


def test_ase_roundtrip(water_molecule):
    pytest.importorskip("ase")
    atoms = water_molecule.to_ase_atoms()
    mol_back = Molecule.from_ase_atoms(atoms, name="WaterBack")
    assert mol_back.get_number_of_atoms() == water_molecule.get_number_of_atoms()
    assert np.allclose(mol_back.coordinates, water_molecule.coordinates)