import freesasa, itertools, math
from collections import defaultdict

MAX = 236.0
CUT = 0.30
EX = {'44','49','200','421','234','215','272','113','88','33'}

# Genipin reach. The previous version used a 5.0 A floor, which counts pairs
# too close for the crosslinker to bridge: four of the six Lys377-Lys400
# cross-ring pairs sit at 8.4-8.5 A. Everything else in the project uses
# 10-16 A, so this does too.
LO, HI_GENIPIN = 5.0, 16.0

nz = {}
for l in open('4lnk_dodecamer.pdb'):
    if l.startswith('ATOM') and l[17:20] == 'LYS' and l[12:16].strip() == 'NZ':
        nz[(l[21], l[22:26].strip())] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))

r = freesasa.calc(freesasa.Structure('4lnk_dodecamer.pdb')).residueAreas()
rsa = {(c, n): a.total / MAX for c in r for n, a in r[c].items() if a.residueType == 'LYS'}
av = {k: v for k, v in nz.items() if rsa.get(k, 0) > CUT and k[1] not in EX}

ring = {c: 1 for c in 'ABCDEF'}
ring.update({c: 2 for c in 'GHIJKL'})
d = lambda a, b: math.sqrt(sum((i-j)**2 for i, j in zip(a, b)))

# cross-subunit pairs only; within-subunit coverage is multipoint2.py's job
allpairs = []
for a, b in itertools.combinations(sorted(av), 2):
    if a[0] == b[0]:
        continue
    allpairs.append((a, b, d(av[a], av[b])))

# ---- headline: the genipin window ----
same = [p for p in allpairs if LO <= p[2] <= HI_GENIPIN and ring[p[0][0]] == ring[p[1][0]]]
cross = [p for p in allpairs if LO <= p[2] <= HI_GENIPIN and ring[p[0][0]] != ring[p[1][0]]]
print(f"=== genipin span {LO:.0f}-{HI_GENIPIN:.0f} A (cross-subunit pairs) ===")
print(f"  within-ring  {len(same):3d}")
print(f"  CROSS-RING   {len(cross):3d}")
print()
print("  cross-ring pairs at genipin span:")
for a, b, x in sorted(cross, key=lambda p: p[2]):
    print(f"    {a[0]}/Lys{a[1]} - {b[0]}/Lys{b[1]}  {x:.1f} A")

# ---- sweep, both floors, so the window dependence is visible ----
print()
print("=== sweep: pair counts against upper bound, both floors ===")
print(f"{'upper':>6} | {'floor 10 A':>20} | {'floor 5 A':>20}")
print(f"{'(A)':>6} | {'within':>9}{'cross':>11} | {'within':>9}{'cross':>11}")
for UP in (16, 18, 20, 22, 25, 30):
    row = []
    for FLOOR in (10.0, 5.0):
        s = sum(1 for p in allpairs
                if FLOOR <= p[2] <= UP and ring[p[0][0]] == ring[p[1][0]])
        c = sum(1 for p in allpairs
                if FLOOR <= p[2] <= UP and ring[p[0][0]] != ring[p[1][0]])
        row += [s, c]
    print(f"{UP:>6} | {row[0]:>9}{row[1]:>11} | {row[2]:>9}{row[3]:>11}")

# ---- residue pairs stitching the rings, up to 30 A ----
print()
print("=== residue pairs stitching the rings, up to 30 A, by frequency ===")
res = defaultdict(int)
for a, b, x in sorted([p for p in allpairs
                       if ring[p[0][0]] != ring[p[1][0]] and p[2] <= 30],
                      key=lambda p: p[2]):
    res[tuple(sorted((a[1], b[1])))] += 1
for k, v in sorted(res.items(), key=lambda x: -x[1]):
    print(f"  Lys{k[0]}-Lys{k[1]}: {v}")
