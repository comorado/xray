from xray import nistdb, bragg
import numpy as np
from matplotlib import pyplot as pt, \
                       rcParams
from collections import defaultdict

transitions = ["Ka1", "Ka2", "Kb1", "La1", "La2", "Lb1", "Lg1", "Ma1", "Mb", "Mg"]
transitions = set([ nistdb.siegbahn_map[t] for t in transitions ]) # convert to IUPAC notation

lines = [l for l in nistdb.lines_in_range(3000,10000) if l.trans in transitions]
line_labels = [lambda l: "%s %s" % (l.elt_name, nistdb.iupac_map[l.trans])]
line_energies = [l.energy for l in lines]

xtals = ['Si', 'Ge', 'GaP']
hkl = bragg.generate_hkl_matrix(9)

# build list of reflections
reflections = defaultdict(lambda:defaultdict(lambda:[]))
for xtal in xtals:
  xtal_type, d = bragg.crystal_info[xtal]
  for m in bragg.filters[xtal_type](hkl):
    spacing = d/np.linalg.norm(m)
    mstr = '(%s)' % ''.join(str(mi) for mi in m)
    reflections[spacing][xtal].append(mstr)

spacings = np.array(sorted(reflections.keys()))
E1 = bragg.bragg_energy(90, spacings, degrees=True)
E2 = bragg.bragg_energy(70, spacings, degrees=True)

i = (E2 > 3000) & (E1 < 10000)
spacings = spacings[i]
E1 = E1[i]
E2 = E2[i]


# find intersections between emission lines and xtal reflections
intersections = [[],[],[]]
for E in line_energies:
  i = (E1 <= E) & (E <= E2)
  angles = bragg.bragg_angle(E, spacings[i], degrees=True)
  intersections[0].append(angles)
  intersections[1].append(np.ones_like(angles) * E)
  intersections[2].append(np.ones_like(angles)*40) # XXX cal dE/dtheta here
intersections = [np.hstack(x) for x in intersections]

styles = {
    'Si': dict(color='blue'),
    'Ge': dict(color='red'),
    'GaP': dict(color='green'),
    }
rcParams['lines.linewidth'] = 0.5

pt.figure(figsize=(20,30))
pt.hlines(line_energies, 70, 90)

angles = np.linspace(70,90,1001)
for s in spacings:
  ref = reflections[s]
  style_key = ' '.join(ref.keys())
  #style_key = ' '.join(x for x in ref if len(ref[x]) > 0)
  label = '; '.join([x + ' ' + ', '.join(ref[x]) for x in xtals if ref[x]])
  style = styles.get(style_key, dict(color='black'))

  ref_info = reflections[s]
  #label = build_reflection_label(ref_info)
  energy = bragg.bragg_energy(angles, s, degrees=True)
  pt.plot(angles, energy, **style)
  pt.text(90.2, energy[-1]-15, label, **style)

pt.scatter(intersections[0], intersections[1], intersections[2], facecolor='white', linewidth=0.5)

pt.xlim(70,90)
pt.ylim(3000,10000)
pt.savefig("poster.pdf")
