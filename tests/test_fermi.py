from xray import fermi
from nose.tools import assert_almost_equal

def test_mu_T0():
  n = 1.0

  mu0 = fermi.fermi_energy(n)
  mu = fermi.mu_T(fermi.rhoe, n, 0)

  assert_almost_equal(mu, mu0)
