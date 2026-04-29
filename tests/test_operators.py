import numpy as np
import random
import pytest

from modules.molecule_class import Molecule
from modules.operators import LocalOperators, NonLocalOperators
from modules.geometry import GeometryOps


@pytest.fixture(autouse=True)
def deterministic_seed():
    random.seed(0)
    np.random.seed(0)
    yield


@pytest.fixture
def water_dimer_labels():
    return ["O1", "H1", "H2", "O2", "H3", "H4"]


@pytest.fixture
def water_dimer_coords():
    return np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.9584, 0.0000, 0.0000],
            [-0.2396, 0.9271, 0.0000],
            [3.0000, 0.0000, 0.0000],
            [3.9584, 0.0000, 0.0000],
            [2.7604, 0.9271, 0.0000],
        ],
        dtype=np.float64,
    )


@pytest.fixture
def water_dimer(water_dimer_labels, water_dimer_coords):
    return Molecule.from_labels_and_coords(water_dimer_labels, water_dimer_coords, name="WaterDimer")


@pytest.fixture
def submolecule_indices():
    return [[0, 1, 2], [3, 4, 5]]


@pytest.fixture
def local_ops():
    return LocalOperators()


@pytest.fixture
def nonlocal_ops():
    return NonLocalOperators()


def _fragment_coms(mol, submolecule_indices):
    coms = []
    for idxs in submolecule_indices:
        coords = mol.coordinates[idxs]
        masses = mol.masses[idxs]
        coms.append(GeometryOps.center_of_mass(coords, masses))
    return np.vstack(coms)


def test_random_displacement_submolecule_changes_com(local_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = local_ops.random_displacement_submolecule(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert perturbed.coordinates.shape == water_dimer.coordinates.shape
    assert not np.allclose(before, after)


def test_random_rotation_submolecule_keeps_com(local_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = local_ops.random_rotation_submolecule(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert np.allclose(before, after)


def test_local_twist_submolecule_keeps_com(local_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = local_ops.local_twist_submolecule(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert np.allclose(before, after)


def test_correlated_displacement_changes_com(local_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = local_ops.correlated_displacement(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)


def test_small_principal_axis_rotation_keeps_com(local_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = local_ops.small_principal_axis_rotation(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert np.allclose(before, after)


def test_twist_operator_changes_com(nonlocal_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = nonlocal_ops.twist_operator(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)


def test_large_displacement_changes_com(nonlocal_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = nonlocal_ops.large_displacement(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)


def test_com_com_approach_operator_changes_com(nonlocal_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = nonlocal_ops.com_com_approach_operator(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)


def test_com_com_separation_operator_changes_com(nonlocal_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = nonlocal_ops.com_com_separation_operator(water_dimer, submolecule_indices, adaptive=False)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)


def test_exchange_operator_changes_com(nonlocal_ops, water_dimer, submolecule_indices):
    before = _fragment_coms(water_dimer, submolecule_indices)
    perturbed = nonlocal_ops.exchange_operator(water_dimer, submolecule_indices)
    assert isinstance(perturbed, Molecule)
    after = _fragment_coms(perturbed, submolecule_indices)
    assert not np.allclose(before, after)
