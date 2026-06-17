import numpy as np
import pytest
import os

from modules.molecule_class import Molecule
from modules.calculator import EnergyEvaluator


@pytest.fixture
def water_molecule():
    """Provides a simple water molecule for testing"""
    labels = ["O1", "H1", "H2"]
    coords = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.9584, 0.0000, 0.0000],
            [-0.2396, 0.9271, 0.0000],
        ],
        dtype=np.float64,
    )
    return Molecule.from_labels_and_coords(labels, coords, name="Water")


@pytest.fixture
def emt_evaluator():
    """Provides an EnergyEvaluator instance using the lightweight ASE EMT backend"""
    return EnergyEvaluator(backend="ase_emt")


# ====== Initialization tests ======

def test_initialization_ase_emt():
    """Test that EnergyEvaluator initializes correctly with the ase_emt backend"""
    evaluator = EnergyEvaluator(backend="ase_emt")
    assert evaluator.backend == "ase_emt"
    assert evaluator.calculator is not None


def test_initialization_invalid_backend():
    """Test that EnergyEvaluator raises ValueError for an unsupported backend"""
    with pytest.raises(ValueError, match="Unsupported backend"):
        EnergyEvaluator(backend="invalid_backend")


def test_conversion_factors():
    """Test that the unit conversion constants have the expected values"""
    assert np.isclose(EnergyEvaluator.EV_TO_HARTREE, 0.0367493)
    assert np.isclose(EnergyEvaluator.E_ANG_TO_DEBYE, 4.80320427)


def test_xtb_initialization():
    """Test EnergyEvaluator initialization with the xTB backend"""
    pytest.importorskip("xtb", reason="xTB is not installed")
    evaluator = EnergyEvaluator(backend="xtb", xtb_method="GFN2-xTB")
    assert evaluator.backend == "xtb"
    assert evaluator.calculator is not None


def test_psi4_initialization():
    """Test EnergyEvaluator initialization with the Psi4 backend"""
    pytest.importorskip("psi4", reason="Psi4 is not installed")
    evaluator = EnergyEvaluator(backend="psi4", qm_method="hf", qm_basis="6-31g")
    assert evaluator.backend == "psi4"
    assert evaluator.calculator is not None


def test_gpaw_initialization():
    """Test EnergyEvaluator initialization with the GPAW backend"""
    pytest.importorskip("gpaw", reason="GPAW is not installed")
    evaluator = EnergyEvaluator(backend="gpaw")
    assert evaluator.backend == "gpaw"
    assert evaluator.calculator is not None


# ====== PSI_SCRATCH tests ======

def test_ensure_psi_scratch_creates_default(monkeypatch, tmp_path):
    """Test that _ensure_psi_scratch sets a default scratch directory when PSI_SCRATCH is unset"""
    monkeypatch.delenv("PSI_SCRATCH", raising=False)
    monkeypatch.chdir(tmp_path)
    evaluator = EnergyEvaluator(backend="ase_emt")
    evaluator._ensure_psi_scratch()
    assert os.environ.get("PSI_SCRATCH") is not None
    assert os.path.isdir(os.environ["PSI_SCRATCH"])


def test_ensure_psi_scratch_accepts_valid_path(monkeypatch, tmp_path):
    """Test that _ensure_psi_scratch accepts a valid, writable PSI_SCRATCH path"""
    monkeypatch.setenv("PSI_SCRATCH", str(tmp_path))
    evaluator = EnergyEvaluator(backend="ase_emt")
    evaluator._ensure_psi_scratch()  # Should not raise


def test_ensure_psi_scratch_raises_on_invalid_path(monkeypatch):
    """Test that _ensure_psi_scratch raises PermissionError for a non-existent PSI_SCRATCH path"""
    monkeypatch.setenv("PSI_SCRATCH", "/nonexistent/path/that/does/not/exist")
    evaluator = EnergyEvaluator(backend="ase_emt")
    with pytest.raises(PermissionError):
        evaluator._ensure_psi_scratch()


# ====== Energy evaluation tests ======

def test_evaluate_returns_float(emt_evaluator, water_molecule):
    """Test that evaluate returns a scalar energy value"""
    energy = emt_evaluator.evaluate(water_molecule)
    assert isinstance(energy, float), "Energy should be a float"


def test_evaluate_energy_is_finite(emt_evaluator, water_molecule):
    """Test that the evaluated energy is a finite number"""
    energy = emt_evaluator.evaluate(water_molecule)
    assert np.isfinite(energy), "Energy should be a finite number"


def test_evaluate_energy_unit_conversion(emt_evaluator, water_molecule):
    """Test that the returned energy is scaled by EV_TO_HARTREE relative to a raw eV value"""
    from ase.calculators.emt import EMT
    atoms = water_molecule.to_ase_atoms()
    atoms.calc = EMT()
    energy_ev = atoms.get_potential_energy()
    energy_hartree = emt_evaluator.evaluate(water_molecule)
    assert np.isclose(energy_hartree, energy_ev * EnergyEvaluator.EV_TO_HARTREE)


# ====== Force evaluation tests ======

def test_evaluate_forces_shape(emt_evaluator, water_molecule):
    """Test that evaluate_forces returns an array of shape (n_atoms, 3)"""
    forces = emt_evaluator.evaluate_forces(water_molecule)
    n_atoms = water_molecule.get_number_of_atoms()
    assert forces.shape == (n_atoms, 3), f"Forces should have shape ({n_atoms}, 3), got {forces.shape}"


def test_evaluate_forces_are_finite(emt_evaluator, water_molecule):
    """Test that all force components are finite numbers"""
    forces = emt_evaluator.evaluate_forces(water_molecule)
    assert np.all(np.isfinite(forces)), "All force components should be finite"


def test_evaluate_forces_unit_conversion(emt_evaluator, water_molecule):
    """Test that returned forces are scaled by EV_TO_HARTREE relative to raw eV/Å values"""
    from ase.calculators.emt import EMT
    atoms = water_molecule.to_ase_atoms()
    atoms.calc = EMT()
    forces_ev = atoms.get_forces()
    forces_hartree = emt_evaluator.evaluate_forces(water_molecule)
    assert np.allclose(forces_hartree, forces_ev * EnergyEvaluator.EV_TO_HARTREE)


# ====== Geometry optimization tests ======

def test_optimize_geometry_bfgs_returns_molecule(emt_evaluator, water_molecule):
    """Test that geometry optimization with BFGS returns a Molecule with the same atom count"""
    optimized, energy = emt_evaluator.optimize_geometry(water_molecule, optimizer="BFGS")
    assert isinstance(optimized, Molecule), "Optimized result should be a Molecule instance"
    assert optimized.get_number_of_atoms() == water_molecule.get_number_of_atoms()
    assert isinstance(energy, float)


def test_optimize_geometry_lbfgs_returns_molecule(emt_evaluator, water_molecule):
    """Test that geometry optimization with LBFGS returns a Molecule with the same atom count"""
    optimized, energy = emt_evaluator.optimize_geometry(water_molecule, optimizer="LBFGS")
    assert isinstance(optimized, Molecule), "Optimized result should be a Molecule instance"
    assert optimized.get_number_of_atoms() == water_molecule.get_number_of_atoms()
    assert isinstance(energy, float)


def test_optimize_geometry_invalid_optimizer(emt_evaluator, water_molecule):
    """Test that an unsupported optimizer raises ValueError"""
    with pytest.raises(ValueError, match="Unsupported optimizer"):
        emt_evaluator.optimize_geometry(water_molecule, optimizer="gradient_descent")


def test_optimize_geometry_lowers_energy(emt_evaluator, water_molecule):
    """Test that geometry optimization produces a lower or equal energy than the initial structure"""
    energy_before = emt_evaluator.evaluate(water_molecule)
    _, energy_after = emt_evaluator.optimize_geometry(water_molecule, optimizer="BFGS")
    assert energy_after <= energy_before + 1e-6, (
        f"Optimized energy ({energy_after:.6f} Ha) should be <= initial energy ({energy_before:.6f} Ha)"
    )


# ====== Hessian tests ======

def test_compute_hessian_shape(emt_evaluator, water_molecule, tmp_path, monkeypatch):
    """Test that compute_hessian returns a 2D array of shape (3*n_atoms, 3*n_atoms)"""
    monkeypatch.chdir(tmp_path)  # Vibrations writes cache files to the working directory
    hessian = emt_evaluator.compute_hessian(water_molecule)
    n_atoms = water_molecule.get_number_of_atoms()
    assert hessian.ndim == 2, "Hessian should be a 2D array"
    assert hessian.shape == (3 * n_atoms, 3 * n_atoms), (
        f"Hessian should have shape ({3*n_atoms}, {3*n_atoms}), got {hessian.shape}"
    )


def test_compute_hessian_is_symmetric(emt_evaluator, water_molecule, tmp_path, monkeypatch):
    """Test that the Hessian matrix is symmetric"""
    monkeypatch.chdir(tmp_path)
    hessian = emt_evaluator.compute_hessian(water_molecule)
    assert np.allclose(hessian, hessian.T, atol=1e-6), "Hessian matrix should be symmetric"
