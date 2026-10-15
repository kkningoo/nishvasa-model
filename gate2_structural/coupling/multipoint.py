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
    """Return {(chain, resi): (x,y,z)} for every Lys NZ atom."""
    out = {}
    for line in open(pdb):
        if line.startswith('ATOM') and line[17:20] == 'LYS' and line[12:16].strip() == 'NZ':
            ch, resi = line[21], line[22:26].strip()
            out[(ch, resi)] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
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
    print(f"total Lys: {len(nz)}   available (RSA>{RSA_CUT}, not excluded): {len(avail)}")

    # pairs within crosslinker reach
    pairs = []
    for a, b in itertools.combinations(sorted(avail), 2):
        d = dist(avail[a], avail[b])
        if SPAN_MIN <= d <= SPAN_MAX:
            pairs.append((a, b, d))

    print(f"pairs within {SPAN_MIN}-{SPAN_MAX} A: {len(pairs)}")

    # connected clusters
    adj = defaultdict(set)
    for a, b, _ in pairs:
        adj[a].add(b); adj[b].add(a)

    seen, clusters = set(), []
    for k in avail:
        if k in seen: continue
        stack, comp = [k], []
        while stack:
            n = stack.pop()
            if n in seen: continue
            seen.add(n); comp.append(n)
            stack.extend(adj[n] - seen)
        clusters.append(sorted(comp))

    multi = [c for c in clusters if len(c) > 1]
    print(f"clusters of 2+: {len(multi)}   singletons: {len(clusters)-len(multi)}")

    # per-chain coverage
    per_chain = defaultdict(int)
    for c in multi:
        for ch, _ in c:
            per_chain[ch] += 1
    print("\nmultipoint-capable lysines per chain:")
    for ch in sorted(per_chain):
        print(f"  {ch}: {per_chain[ch]}")

    chains = sorted(set(ch for ch, _ in avail))
    missing = [c for c in chains if per_chain.get(c, 0) == 0]
    print(f"\nchains with NO multipoint cluster: {missing if missing else 'none'}")

    print("\nlargest clusters:")
    for c in sorted(multi, key=len, reverse=True)[:6]:
        print(f"  n={len(c)}: " + ", ".join(f"{ch}/Lys{r}" for ch, r in c))

analyse('ppk2_minimised.pdb', exclude={'43', '164', '138'})
