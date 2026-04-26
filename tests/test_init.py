import numpy as np 
import pytest

from modules.init import ClusterInitializer, InitializationConfig

@pytest.fixture
# Initialize Psi4 Backend
def psi4_initializer():
    config = InitializationConfig(
        backend="psi4",
        qm_method="hf",
        qm_basis="6-31g"
    )
    return ClusterInitializer(config)

@pytest.fixture
# Initialize XTB Backend
def xtb_initializer():
    config = InitializationConfig(
        backend="xtb",
        xtb_method="GFN2-xTB"
    )
    return ClusterInitializer(config)

