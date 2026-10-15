"""Chemical liability map. CORRECTED maxima.

Change from the previous version: MET, SER and THR now use Tien et al.
theoretical maxima (224, 155, 172) rather than 203, 130, 142. The old values
inflated methionine RSA by ~10% and threonine by ~21%, and did not match the
table the Model page cites. Every other residue was already correct.

Re-run and replace the liability RSA values on the page.
"""

import freesasa

MAXACC = {'ALA':129,'ARG':274,'ASN':195,'ASP':193,'CYS':167,'GLU':223,'GLN':225,
          'GLY':104,'HIS':224,'ILE':197,'LEU':201,'LYS':236,'MET':224,'PHE':240,
          'PRO':159,'SER':155,'THR':172,'TRP':285,'TYR':263,'VAL':174}

THREE2ONE = {v: k for k, v in {
 'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY',
 'H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER',
 'T':'THR','W':'TRP','Y':'TYR','V':'VAL'}.items()}

RSA_CUT = 0.20


def chain_seq(pdb, chain):
    seq, nums, seen = [], [], set()
    for l in open(pdb):
        if l.startswith('ATOM') and l[21] == chain and l[12:16].strip() == 'CA':
            r = l[22:26].strip()
            if r in seen:
                continue
            seen.add(r)
            seq.append(THREE2ONE.get(l[17:20], 'X'))
            nums.append(r)
    return ''.join(seq), nums


def rsa_all(pdb, chain):
    res = freesasa.calc(freesasa.Structure(pdb)).residueAreas()
    out = {}
    for num, a in res[chain].items():
        m = MAXACC.get(a.residueType)
        if m:
            out[num] = (a.residueType, a.total / m)
    return out


def scan(pdb, chain, critical, label):
    seq, nums = chain_seq(pdb, chain)
    rsa = rsa_all(pdb, chain)
    hits = []

    for i in range(len(seq) - 1):
        if seq[i] == 'N' and seq[i+1] in 'GSTNH':
            hits.append(('deamidation', f'N{nums[i]}{seq[i+1]}{nums[i+1]}', nums[i]))
        if seq[i] == 'D' and seq[i+1] in 'GST':
            hits.append(('isomerisation', f'D{nums[i]}{seq[i+1]}{nums[i+1]}', nums[i]))

    for i, aa in enumerate(seq):
        if aa in 'MCW':
            hits.append(('oxidation', f'{aa}{nums[i]}', nums[i]))

    print(f"\n=== {label} (chain {chain}) ===")
    rows = []
    for kind, motif, num in hits:
        rt, r = rsa.get(num, ('?', 0.0))
        if r >= RSA_CUT:
            rows.append((r, kind, motif, 'CRITICAL' if num in critical else ''))

    rows.sort(reverse=True)
    print(f"{'RSA':>5}  {'type':<14} {'site':<12} flag")
    for r, kind, motif, flag in rows:
        print(f"{r:5.2f}  {kind:<14} {motif:<12} {flag}")
    print(f"\nexposed liabilities: {len(rows)}   critical: {sum(1 for x in rows if x[3])}")


ppk2_crit = {'38','39','40','41','42','43','44','94','95','155','164','166','168','169'}
scan('ppk2_minimised.pdb', 'A', ppk2_crit, 'PPK2')

gs_crit = {'44','49','132','189','195','196','200','245','333'}
scan('gs_prot.pdb', 'A', gs_crit, 'GS')
