import numpy as np
import pytest
import random

from modules.transformations import Transformation 
from modules.molecule_class import Molecule
from modules.geometry import GeometryOps


@pytest.fixture(autouse=True)
def deterministic_seed():
    random.seed(0)
    np.random.seed(0)
    yield

@pytest.fixture
def water_coords():
    return np.array([
        [0.0000, 0.0000, 0.0000],
        [0.9584, 0.0000, 0.0000],
        [-0.2396, 0.9271, 0.0000],
    ], dtype=float)

@pytest.fixture
def water_labels():
    return ["O1", "H1", "H2"]

@pytest.fixture
def water_molecule(water_labels, water_coords):
    return Molecule.from_labels_and_coords(water_labels, water_coords, name="Water")

def test_get_center_centroid_and_com(water_molecule):
    T = Transformation()
    centroid = T.get_center(water_molecule, method="centroid")
    com = T.get_center(water_molecule, method="com")

    # Centroid is the mean of the coordinates
    expected_centroid = np.mean(water_molecule.coordinates, axis=0)
    assert np.allclose(centroid, expected_centroid), f"Expected centroid {expected_centroid}, got {centroid}"

    expected_com = GeometryOps.center_of_mass(water_molecule.coordinates, water_molecule.masses)
    assert np.allclose(com, expected_com), f"Expected COM {expected_com}, got {com}"

def test_get_center_coords_accepts_arrays():
    T = Transformation()
    coords = np.array([[0., 0., 0.],[1., 0., 0.]])
    centroid = T.get_center_coords(coords, method="centroid")
    expected_centroid = np.mean(coords, axis=0)
    assert np.allclose(centroid, expected_centroid), f"Expected centroid {expected_centroid}, got {centroid}"
    com = T.get_center_coords(coords, method="com", masses=np.array([1., 1.]))
    expected_com = GeometryOps.center_of_mass(coords, np.array([1., 1.]))
    assert np.allclose(com, expected_com), f"Expected COM {expected_com}, got {com}"

def test_translate_inplace_vs_copy(water_molecule):
    T = Transformation()
    v = np.array([1.0, 0.5, 0.0])
    orig = water_molecule.coordinates.copy()

    # In place 
    T.translate(water_molecule, v, in_place=True)
    assert np.allclose(water_molecule.coordinates, orig + v), "In-place translation failed"

    # Copy
    mol2 = Molecule.from_labels_and_coords(water_molecule.atom_labels, orig, name="WaterCopy")
    new_coords = T.translate(mol2, -v, in_place=False)
    assert np.allclose(new_coords, orig -v), "Copy translation failed"
    # Original molecule should be unchanged
    assert np.allclose(mol2.coordinates, orig), "Original molecule modified during copy translation"

def test_rotate_perserves_pairwise_distances(water_molecule):
    T = Transformation()
    coords0 = water_molecule.coordinates.copy()
    D0 = np.linalg.norm(coords0[None, :, :] - coords0[:, None, :], axis=-1)
    T.rotate(water_molecule, axis=np.array([0., 0., 1.]), angle_degree=30, in_place=True)
    coords1 = water_molecule.coordinates
    D1 = np.linalg.norm(coords1[None, :, :] - coords1[:, None, :], axis=-1)
    assert np.allclose(D0, D1), "Pairwise distances not preserved under rotation"

def test_rotate_quaternion_matches_rotate(water_molecule):
    T = Transformation()
    coords0 = water_molecule.coordinates.copy()
    mol_copy = Molecule.from_labels_and_coords(water_molecule.atom_labels, coords0, name="WaterCopy")
    axis = np.array([0., 0., 1.])
    angle = 45
    T.rotate(mol_copy, axis=axis, angle_degree=angle, in_place=True)
    mol_q = Molecule.from_labels_and_coords(water_molecule.atom_labels, coords0, name="WaterQ")
    T.rotate_quaternion(mol_q, axis=axis, angle_degree=angle, in_place=True)
    # Compare Rotation results
    assert np.allclose(mol_copy.coordinates, mol_q.coordinates), "Quaternion rotation does not match standard rotation"

def test_center_molecule_returns_and_inplace(water_molecule):
    T = Transformation()
    mol = Molecule.from_labels_and_coords(water_molecule.atom_labels, water_molecule.coordinates.copy(), name="WaterCopy")
    center, new_coords = T.center_molecule(mol, method="centroid", in_place=True)
    # Inplace should have moved coords to that the new centroid is close to 0 
    assert np.allclose(np.mean(mol.coordinates, axis=0), [0., 0., 0.]), f"Expected centroid at origin after centering, got {np.mean(mol.coordinates, axis=0)}"
    # Non-in-place should return center and new coords
    mol2 = Molecule.from_labels_and_coords(water_molecule.atom_labels, water_molecule.coordinates.copy(), name="WaterCopy2")
    center2, new_coords2 = T.center_molecule(mol2, method="centroid", in_place=False)
    assert np.allclose(center2, np.mean(mol2.coordinates, axis=0)), f"Expected returned center to match original centroid, got {center2}"
    assert np.allclose(new_coords2, mol2.coordinates - center2), "Expected new coords to be original coords minus center"

def test_compute_reference_frame(water_molecule):
    T = Transformation()
    com = GeometryOps.center_of_mass(water_molecule.coordinates, water_molecule.masses)
    ref = T.compute_reference_frame(water_molecule.coordinates, water_molecule.masses, center=com)
    # Axes should be orthonormal
    axes = ref.axes
    I = axes.T @ axes
    assert np.allclose(I, np.eye(3)), f"Expected orthonormal axes, got {I}"
    assert np.allclose(ref.origin, com), f"Expected origin to be at COM, got {ref.origin}"

def test_set_reference_frame(water_molecule):
    T = Transformation()
    mol = Molecule.from_labels_and_coords(water_molecule.atom_labels, water_molecule.coordinates.copy(), name="WaterCopy")
    ref = T.set_reference_frame(mol, method="centroid")
    # After set_reference frame the new centroid should be near zero
    centroid = np.mean(mol.coordinates, axis=0)
    assert np.allclose(centroid, [0., 0., 0.]), f"Expected centroid at origin after setting reference frame, got {centroid}"