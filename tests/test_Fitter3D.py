

import pytest
import numpy as np

# from ROOT import Fitter3DGaussLcms

@pytest.fixture
def ROOT():
    return pytest.importorskip('ROOT')

@pytest.fixture
def gSystem(ROOT):
    from ROOT import gSystem
    return gSystem

@pytest.fixture
def Fitter3DGaussLcms(ROOT):
    from ROOT import Fitter3DGaussLcms
    return Fitter3DGaussLcms


@pytest.mark.parametrize("gamma, expected", [
    (1.0, 4.08248290463863),
])
def test_pseudoRinv(gamma, expected, Fitter3DGaussLcms):
    norm, lam, ro, rs, rl = 1.0, 1.0, 5.0, 5.0, 5.0
    params = Fitter3DGaussLcms.FitParams(np.array([norm, lam, ro, rs, rl]))
    assert np.isclose(params.PseudoRinv(gamma), expected)
