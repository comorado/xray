from xray import special
import numpy as np
from nose import tools as nt

def test_sph_j_of_0():
  nt.assert_almost_equal(special.sph_jn(0, 0), 1.0)
  nt.assert_almost_equal(special.sph_jn(0, [0]), 1.0)
  nt.assert_almost_equal(special.sph_jn(0, np.array([0])), 1.0)
  nt.assert_almost_equal(special.sph_jn(1, 0), 0.0)
  nt.assert_almost_equal(special.sph_jn(1, [0]), 0.0)
  nt.assert_almost_equal(special.sph_jn(1, np.array([0])), 0.0)

def test_sum_of_j_squared():
  x = np.linspace(0,5,11)
  y = np.sum([special.sph_jn(l, x)**2 * (2*l+1) for l in range(10)], 0)
  assert np.allclose(y, 1.0) 
