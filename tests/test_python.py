
from ROOT import gSystem


import pytest

@pytest.fixture
def ROOT():
    import ROOT
    return ROOT

@pytest.fixture
def FitterGaussOSL():
    from ROOT import FitterGaussOSL
    return FitterGaussOSL

def test_load_library():
    assert gSystem.Load("build/libFemtoFitter.so") == 0


def test_gauss(FitterGaussOSL):
    assert FitterGaussOSL.chi2()
