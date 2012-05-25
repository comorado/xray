from xray import fermi, compton, const
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Calculate Fermi gas Compton profiles")

parser.add_argument("-T", "--temperature", help="temperature (in eV)",
                    type=float, default=0.0)
parser.add_argument("-N", "--num_electrons", help="number of electrons",
                    type=float, default=1)
parser.add_argument("-V", "--volume", help="number of electrons",
                    type=float, default=1.0)
parser.add_argument("-m" "--min", help="minimum p_q",
                    type=float, default=0.0)
parser.add_argument("-M" "--max", help="maximum p_q",
                    type=float, default=None)
parser.add_argument("-p" "--points", help="numer of p_q points",
                    type=int, default=1001)


args = parser.parse_args()

N = args.N
V = args.V
n = N/V 

T = args.T / const.HARTREE
pq_min = args.min
pq_max = args.max
npq = 1001 if args.points is None else args.points


pf = fermi.fermi_momentum(n)
if pq_max is None:
  pq_max = pf + 1.0 if T == 0 else pf + 3 * sqrt(2*T)
pq = np.linspace(pq_min, pq_max, npq)


mu = None
if T > 0:
  mu = 

if T == 0:
  jpq = compton.fermi_profile(args.N, args.V)
else:
  #jpq = compton.isotropic_profile(fermi.rhop, pq, ...)
  raise Exception("Unimplemented")
