import freesasa, itertools, math
from collections import defaultdict

MAXACC = {'LYS': 236.0}
# Genipin reach. The LOWER bound is steric, not chemical: a bifunctional
# reagent is limited by how far it can span, not by how close two amines
# may be, and genipin self-polymerises into short dimers and trimers. A
# 10 A floor here once excluded four of six real cross-ring contacts at
# 8.4-8.5 A. Report the window with every count.
SPAN_MIN, SPAN_MAX = 5.0, 16.0
RSA_CUT = 0.30

def load_lys_nz(pdb):
    out = {}
    for line in open(pdb):
        if line.startswith('ATOM') and line[17:20] == 'LYS' and line[12:16].strip() == 'NZ':
            out[(line[21], line[22:26].strip())] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out

def rsa_map(pdb):
    res = freesasa.calc(freesasa.Structure(pdb)).residueAreas()
    out = {}
    for ch in res:
        for num, a in res[ch].items():
            if a.residueType == 'LYS':
                out[(ch, num)] = a.total / MAXACC['LYS']
    return out

def dist(a, b):
    return math.sqrt(sum((i-j)**2 for i, j in zip(a, b)))

def analyse(pdb, exclude=()):
    nz  = load_lys_nz(pdb)
    rsa = rsa_map(pdb)
    avail = {k: v for k, v in nz.items()
             if rsa.get(k, 0) > RSA_CUT and k[1] not in exclude}

    print(f"\n=== {pdb} ===")
    print(f"available: {len(avail)} of {len(nz)}")

    nbrs = defaultdict(list)
    cross = []
    for a, b in itertools.combinations(sorted(avail), 2):
        d = dist(avail[a], avail[b])
        if SPAN_MIN <= d <= SPAN_MAX:
            nbrs[a].append((b, d))
            nbrs[b].append((a, d))
            if a[0] != b[0]:
                cross.append((a, b, d))

    hist = defaultdict(int)
    per_chain = defaultdict(int)
    for k in sorted(avail):
        n = len(nbrs[k])
        hist[n] += 1
        if n >= 1:
            per_chain[k[0]] += 1

    print("\nneighbours within reach -> number of lysines")
    for n in sorted(hist):
        print(f"  {n}: {hist[n]}")

    print("\nlysines with >=1 partner, per chain:")
    for ch in sorted(per_chain):
        print(f"  {ch}: {per_chain[ch]}")

    first = sorted(set(k[0] for k in avail))[0]
    print(f"\nbest anchor sites (chain {first}):")
    sel = [(k, len(nbrs[k])) for k in nbrs if k[0] == first]
    for k, n in sorted(sel, key=lambda x: -x[1])[:8]:
        print(f"  Lys{k[1]}: {n} partners")

    print(f"\ncross-chain reachable pairs: {len(cross)}")
    seen = set()
    for a, b, d in sorted(cross, key=lambda x: x[2])[:12]:
        key = (a[1], b[1])
        if key in seen: continue
        seen.add(key)
        print(f"  {a[0]}/Lys{a[1]} - {b[0]}/Lys{b[1]}  {d:.1f} A")

analyse('4lnk_dodecamer.pdb', exclude={'44','49','200','421','234','215','272','113','88','33'})
