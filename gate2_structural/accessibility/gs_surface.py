"""GS surface analysis on the dodecamer, with the isolated subunit as reference.

Produces three things the Model page cites and no current script generates:
  1. Per-lysine RSA on the assembly, and which lysines assembly buries
  2. Terminal accessibility, N- and C-, monomer vs assembly  (the 82% -> 9% figure)
  3. The availability count at RSA > 0.30

Run from the directory holding 4lnk_dodecamer.pdb.
Writes gs_surface.tsv alongside the printed summary.

    python3 gs_surface.py
"""

import freesasa
import os

# Tien et al. (2013) THEORETICAL maxima, Gly-X-Gly context.
# Note: hotspots.py currently uses MET 203, SER 130, THR 142 - a different
# table. Use these everywhere so RSA values are comparable across scripts.
TIEN = {
    'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167,
    'GLU': 223, 'GLN': 225, 'GLY': 104, 'HIS': 224, 'ILE': 197,
    'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240, 'PRO': 159,
    'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174,
}

ASSEMBLY = '4lnk_dodecamer.pdb'
MONOMER = 'gs_chainA_only.pdb'
AVAIL = 0.30

# Residues to report as termini. GS is 2-444 in 4LNK.
N_TERM = ['2', '3', '4', '5', '6']
C_TERM = ['440', '441', '442', '443', '444']


def extract_chain(src, dst, chain='A'):
    """Chain A alone, polymer only, so FreeSASA sees an isolated subunit."""
    with open(dst, 'w') as out:
        for line in open(src):
            if line.startswith('ATOM') and line[21] == chain:
                out.write(line)
        out.write('END\n')


def rsa(pdb, chain='A'):
    """{resnum: (resname, rsa)} for one chain."""
    areas = freesasa.calc(freesasa.Structure(pdb)).residueAreas()
    out = {}
    for num, a in areas[chain].items():
        m = TIEN.get(a.residueType)
        if m:
            out[num] = (a.residueType, a.total / m)
    return out


def main():
    if not os.path.exists(MONOMER):
        extract_chain(ASSEMBLY, MONOMER)
        print(f'wrote {MONOMER}')

    asm = rsa(ASSEMBLY, 'A')
    mono = rsa(MONOMER, 'A')
    shared = sorted(set(asm) & set(mono), key=int)

    rows = []
    for n in shared:
        rt, a = asm[n]
        _, m = mono[n]
        rows.append((n, rt, m, a, a - m))

    # ---- 1. lysines ----
    lys = [r for r in rows if r[1] == 'LYS']
    exposed_asm = [r for r in lys if r[3] > AVAIL]
    exposed_mono = [r for r in lys if r[2] > AVAIL]
    lost = [r for r in lys if r[2] > AVAIL and r[3] <= AVAIL]

    print(f'\n=== GS lysines, chain A ===')
    print(f'total {len(lys)} | exposed on monomer {len(exposed_mono)} '
          f'| exposed on assembly {len(exposed_asm)} | lost to assembly {len(lost)}')
    print(f'\n{"Lys":>6}{"monomer":>9}{"assembly":>10}{"delta":>8}')
    for n, rt, m, a, d in sorted(lys, key=lambda r: -r[3]):
        flag = '  <- buried by assembly' if (m > AVAIL and a <= AVAIL) else ''
        print(f'{n:>6}{m:9.2f}{a:10.2f}{d:+8.2f}{flag}')

    # ---- 2. termini ----
    print(f'\n=== Termini, chain A ===')
    print(f'{"res":>6}{"name":>6}{"monomer":>9}{"assembly":>10}{"delta":>8}')

    def block(label, nums):
        vals_m, vals_a = [], []
        for n in nums:
            if n in asm and n in mono:
                rt, a = asm[n]
                _, m = mono[n]
                vals_m.append(m); vals_a.append(a)
                print(f'{n:>6}{rt:>6}{m:9.2f}{a:10.2f}{a-m:+8.2f}')
        if vals_m:
            print(f'{label} mean: monomer {sum(vals_m)/len(vals_m):.2f}  '
                  f'assembly {sum(vals_a)/len(vals_a):.2f}')
            print(f'{label} range: monomer {min(vals_m):.2f}-{max(vals_m):.2f}  '
                  f'assembly {min(vals_a):.2f}-{max(vals_a):.2f}\n')

    block('N-terminus', N_TERM)
    block('C-terminus', C_TERM)

    print('The C-terminal mean is the 82% -> 9% figure cited for both the fusion\n'
          'rejection and N-terminal tag placement. Values above 1.0 at a free\n'
          'N-terminus are expected: Tien maxima assume a Gly-X-Gly context.')

    # ---- 3. full table ----
    with open('gs_surface.tsv', 'w') as f:
        f.write('resnum\tresname\trsa_monomer\trsa_assembly\tdelta\n')
        for n, rt, m, a, d in rows:
            f.write(f'{n}\t{rt}\t{m:.4f}\t{a:.4f}\t{d:+.4f}\n')
    print(f'\nwrote gs_surface.tsv  ({len(rows)} residues)')


if __name__ == '__main__':
    main()
