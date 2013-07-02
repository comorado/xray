from xray import fermi, compton, const
import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser(description="Calculate Fermi gas Compton profiles")

parser.add_argument("-T", "--temperature", help="temperature (in eV)",
                    type=float, default=0.0)
parser.add_argument("-N", "--num_electrons", help="number of electrons",
                    type=float, default=1)
parser.add_argument("-V", "--volume", help="atomic volume (in a.u.)",
                    type=float, default=1.0)
parser.add_argument("-m", "--min", help="minimum p_q",
                    type=float, default=0.0)
parser.add_argument("-M", "--max", help="maximum p_q",
                    type=float, default=None)
parser.add_argument("-p", "--points", help="numer of p_q points",
                    type=int, default=1001)


args = parser.parse_args()
N = args.num_electrons
V = args.volume
n = N/V

T = args.temperature / const.HARTREE

#print N, V, args.temperature, T
pq_min = args.min
pq_max = args.max
npq = 1001 if args.points is None else args.points


pf = fermi.fermi_momentum(n)
if pq_max is None:
  pq_max = pf + 1.0 if T == 0 else pf + 3 * np.sqrt(2*T)
pq = np.linspace(pq_min, pq_max, npq)


if T == 0:
  jpq = compton.fermi_profile(pq, N, V)
else:
  mu = fermi.mu_T(N,V,T)
  jpq = compton.isotropic_profile(fermi.rhop, pq, (N,V,T,mu))
  #raise Exception("Unimplemented")

#header = """# Compton profile for:
## N = {N}
## V = {V}
## T = {T}
#""".format(dict(N=N,V=V,T=T))

#print(header)

np.savetxt(sys.stdout, np.c_[pq, jpq])
