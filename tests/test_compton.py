import numpy as np
from xray import compton

def test_isotropic_profile():
  from xray.fermi import rhop, fermi_momentum

  N = 2.0
  V = 54.422 # (Be unit cell volume in a.u.)

  n = N/V
  pf = fermi_momentum(n)

  pq = np.linspace(0,2,101)


  J1 = compton.isotropic_profile(rhop, pq, pf) * V * (pf > pq)
  J2 = compton.fermi_profile(pq, N, V)

  print (J1-J2).min()
  assert np.allclose(J1, J2)
