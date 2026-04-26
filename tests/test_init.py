import numpy as np
import pytest

from modules.box import SimulationBox
from modules.init import (
	ClusterInitializer,
	InitializationConfig,
	test_submolecule_partition_mapping,
)
from modules.molecule_class import Molecule


def _make_diatomic(name: str, offset: float = 0.0) -> Molecule:
	labels = ["H1", "H2"]
	coords = np.array(
		[
			[offset, 0.0, 0.0],
			[offset + 0.74, 0.0, 0.0],
		],
		dtype=np.float64,
	)
	return Molecule.from_labels_and_coords(labels, coords, name=name)


@pytest.fixture
def init_config() -> InitializationConfig:
	return InitializationConfig(
		optimize_submolecules=False,
		optimize_parallel=False,
		verbose=False,
	)


class TestConfigAndHelpers:
	def test_rank_candidates_selects_lowest_and_skips_failed(self, monkeypatch):
		energies = {"m1": -1.0, "m2": -5.0, "m3": -2.0}

		class FakeEvaluator:
			def __init__(self, **kwargs):
				pass

			def evaluate(self, molecule):
				if molecule.name == "m_fail":
					raise RuntimeError("failed evaluation")
				return energies[molecule.name]

		monkeypatch.setattr("modules.init.EnergyEvaluator", FakeEvaluator)

		cfg = InitializationConfig(verbose=False)
		molecules = [
			_make_diatomic("m1", 0.0),
			_make_diatomic("m_fail", 1.0),
			_make_diatomic("m2", 2.0),
			_make_diatomic("m3", 3.0),
		]

		ranked = cfg.rank_candidates(molecules, n_keep=2)

		assert len(ranked) == 2
		assert ranked[0][1].name == "m2"
		assert ranked[1][1].name == "m3"

	def test_partition_sphere_has_contiguous_full_coverage(self):
		parts = ClusterInitializer.partition_sphere(np.zeros(3), 2.0, 4)

		assert len(parts) == 4
		assert parts[0][0] == pytest.approx(0.0)
		assert parts[-1][1] == pytest.approx(2.0 * np.pi)

		for i in range(1, len(parts)):
			assert parts[i - 1][1] == pytest.approx(parts[i][0])

	def test_partition_sphere_raises_for_invalid_n(self):
		with pytest.raises(ValueError):
			ClusterInitializer.partition_sphere(np.zeros(3), 1.0, 0)


class TestPlacementAndPartitioning:
	def test_generate_partition_points_counts_and_shape(self):
		pts = ClusterInitializer.generate_partition_points(
			center=np.zeros(3),
			radius=2.0,
			n_partitions=2,
			n_theta=3,
			n_phi=4,
			n_r=5,
			spacing="sobol",
			sobol_scramble=True,
			sobol_seed=1,
		)

		assert len(pts) == 2
		assert pts[0].shape == (60, 3)
		assert pts[1].shape == (60, 3)

	def test_generate_random_configurations_validates_input(self, init_config):
		initializer = ClusterInitializer(init_config)
		box = SimulationBox(box_type="sphere", radius=3.0)

		with pytest.raises(ValueError):
			initializer._generate_random_configurations([], box, n_configurations=0)

		with pytest.raises(ValueError):
			initializer._generate_random_configurations([], box, n_configurations=1, placing_method="invalid")

	def test_generate_random_configurations_grid_defaults(self, monkeypatch, init_config):
		initializer = ClusterInitializer(init_config)
		box = SimulationBox(box_type="sphere", radius=3.0)
		submols = [_make_diatomic("a"), _make_diatomic("b", 2.0)]

		captured = []

		def fake_generate_configuration(**kwargs):
			captured.append(kwargs)
			return _make_diatomic("candidate", offset=10.0)

		monkeypatch.setattr(initializer, "_generate_configuration", fake_generate_configuration)

		out = initializer._generate_random_configurations(
			submolecules=submols,
			simulation_box=box,
			n_configurations=8,
			placing_method="grid",
			grid_spacing="sobol",
		)

		assert len(out) == 8
		assert len(captured) == 8
		assert all(call["placing_method"] == "grid" for call in captured)
		assert all(call["grid_spacing"] == "sobol" for call in captured)
		assert all(call["n_theta"] == 2 for call in captured)
		assert all(call["n_phi"] == 2 for call in captured)
		assert all(call["n_r"] == 2 for call in captured)

	def test_submolecule_partition_mapping_detects_valid_and_invalid_cases(self):
		box = SimulationBox(box_type="sphere", radius=10.0, center=np.zeros(3))
		submol_indices = [[0], [1]]

		good = Molecule.from_labels_and_coords(
			["H1", "H2"],
			np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64),
			name="good",
		)

		bad = Molecule.from_labels_and_coords(
			["H1", "H2"],
			np.array([[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=np.float64),
			name="bad",
		)

		assert test_submolecule_partition_mapping([good], submol_indices, box, verbose=False)
		assert not test_submolecule_partition_mapping([bad], submol_indices, box, verbose=False)


class TestInitializationPipelineWithMocks:
	def test_initialize_from_xyz_selects_lowest_energy(self, monkeypatch, init_config):
		initializer = ClusterInitializer(init_config)

		class DummySubmol:
			def __init__(self, idx):
				self._idx = idx

			def get_index_in_parent(self):
				return self._idx

		candidates = [
			_make_diatomic("c1", 0.0),
			_make_diatomic("c2", 1.0),
			_make_diatomic("c3", 2.0),
		]
		energies = {"c1": -1.0, "c2": -2.0, "c3": -0.5}

		class FakeEvaluator:
			def __init__(self, **kwargs):
				pass

			def evaluate(self, molecule):
				return energies[molecule.name]

		monkeypatch.setattr("modules.init.EnergyEvaluator", FakeEvaluator)
		monkeypatch.setattr(initializer, "_load_molecule", lambda _: _make_diatomic("root", 5.0))
		monkeypatch.setattr(initializer, "_fragment_molecule", lambda _: [DummySubmol([0]), DummySubmol([1])])
		monkeypatch.setattr(initializer, "_create_simulation_box", lambda _: SimulationBox(box_type="sphere", radius=4.0))
		monkeypatch.setattr(initializer, "_generate_random_configurations", lambda **kwargs: candidates)
		monkeypatch.setattr("modules.init.plt.savefig", lambda *args, **kwargs: None)

		selected, submol_indices, box = initializer.initialize_from_xyz(
			xyz_file="dummy.xyz",
			n_workers=2,
			n_configurations=3,
			placing_method="random",
			energy_backend="xtb",
		)

		assert [m.name for m in selected] == ["c2", "c1"]
		assert submol_indices == [[0], [1]]
		assert box.box_type == "sphere"

	def test_initialize_from_xyz_raises_when_all_energy_evals_fail(self, monkeypatch, init_config):
		initializer = ClusterInitializer(init_config)

		class DummySubmol:
			def get_index_in_parent(self):
				return [0]

		class AlwaysFailEvaluator:
			def __init__(self, **kwargs):
				pass

			def evaluate(self, molecule):
				raise RuntimeError("failed evaluation")

		monkeypatch.setattr("modules.init.EnergyEvaluator", AlwaysFailEvaluator)
		monkeypatch.setattr(initializer, "_load_molecule", lambda _: _make_diatomic("root", 0.0))
		monkeypatch.setattr(initializer, "_fragment_molecule", lambda _: [DummySubmol()])
		monkeypatch.setattr(initializer, "_create_simulation_box", lambda _: SimulationBox(box_type="sphere", radius=4.0))
		monkeypatch.setattr(
			initializer,
			"_generate_random_configurations",
			lambda **kwargs: [_make_diatomic("x", 0.0), _make_diatomic("y", 1.0)],
		)
		monkeypatch.setattr("modules.init.plt.savefig", lambda *args, **kwargs: None)

		with pytest.raises(RuntimeError, match="All energy evaluations failed"):
			initializer.initialize_from_xyz(
				xyz_file="dummy.xyz",
				n_workers=1,
				n_configurations=2,
				placing_method="random",
				energy_backend="xtb",
			)
