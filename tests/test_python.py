
import pytest


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


def test_load_library(gSystem):
    assert gSystem.Load("build/libFemtoFitter.so") >= 0


#def test_gauss(Fitter3DGaussLcms):
#    assert Fitter3DGaussLcms.calc_chi2()
