import freesasa, itertools

MAXACC = {'LYS': 236.0}

def sasa(pdb, sel=None):
    s = freesasa.Structure(pdb)
    return freesasa.calc(s)

# ---- 1. lysine RSA on the assembled tetramer ----
s = freesasa.Structure('ppk2_minimised.pdb')
res = freesasa.calc(s).residueAreas()

print("=== Lysine RSA (assembled tetramer) ===")
counts = {}
for ch in sorted(res):
    exposed = 0
    for num, a in sorted(res[ch].items(), key=lambda x: int(x[0])):
        if a.residueType == 'LYS':
            rsa = a.total / MAXACC['LYS']
            if rsa > 0.30:
                exposed += 1
            if ch == 'A':
                flag = ' <<<' if num in ('43','73','164') else ''
            print(f"{ch} Lys{num:>4}  {rsa:.2f}")
    counts[ch] = exposed
print("\nexposed (RSA>0.30) per chain:", counts)
