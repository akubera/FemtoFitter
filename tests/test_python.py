
import pytest


@pytest.fixture
def ROOT():
    return pytest.importorskip('ROOT')


@pytest.fixture
def gSystem(ROOT):
    from ROOT import gSystem
    return gSystem


@pytest.fixture
def FitterGaussOSL(ROOT):
    from ROOT import FitterGaussOSL
    return FitterGaussOSL


def test_load_library(gSystem):
    assert gSystem.Load("build/libFemtoFitter.so") >= 0


#def test_gauss(FitterGaussOSL):
#    assert FitterGaussOSL.calc_chi2()

