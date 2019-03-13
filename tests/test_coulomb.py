
import pytest

try:
    from ROOT import gSystem
except ImportError:
    pytest.skip("Could not load ROOT", allow_module_level=True)


@pytest.fixture
def kfile2():
    from ROOT import TFile
    tfile = TFile.Open("KFile2.root")
    if not tfile:
        raise FileNotFoundError(":-(")
    return tfile


@pytest.fixture
def FsiGamov():
    from ROOT import FsiGamov
    return FsiGamov


@pytest.fixture
def FsiKfile():
    from ROOT import FsiKFile
    return FsiKFile


@pytest.mark.parametrize("qinv, Rinv, ex", [
    (0.02, 7.3, 0.9758750214868225),
    (0.02, 7.3, 0.9758750214868225),
])
def test_kfile2(kfile2, qinv, Rinv, ex):
    k2ss = kfile2.Get("k2ss")
    assert k2ss.Interpolate(qinv, Rinv) == ex


@pytest.mark.parametrize("Rinv, qinv, ex", [
    (7.3, [0.02], 
          [0.9758750214868225]),
#    (4.2, [0.3, 0.003,],
#          [0.9946836731940754, .6,]),
])
def test_FSiKFile(FsiKfile, Rinv, qinv, ex):
    f = FsiKfile('KFile2.root')
    k = f.ForRadius(Rinv)

    for q, e in zip(qinv, ex):
        assert k(q) == e


@pytest.mark.parametrize("Rinv, qinv, ex", [
#   (4.2, [0.3, 0.003,],
#         [.98, .6,]),
])
def test_FSiKFile(FsiGamov, qinv, Rinv, ex):
    f = FsiGamov()
    k = f.ForRadius(Rinv)

    for q, e in zip(qinv, ex):
        assert k(q) == e


