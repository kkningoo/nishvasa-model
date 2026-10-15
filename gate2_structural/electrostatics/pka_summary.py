import sys, re

MODEL = {'ASP':3.80,'GLU':4.50,'HIS':6.50,'CYS':9.00,
         'TYR':10.00,'LYS':10.50,'ARG':12.50,
         'N+':8.00,'C-':3.20}

def parse(f):
    rows, curve = [], []
    mode = None
    for line in open(f):
        s = line.strip()
        if s.startswith('SUMMARY OF THIS PREDICTION'):
            mode = 'sum'; continue
        if s.startswith('Free energy of'):
            mode = 'curve'; continue
        if s.startswith('The pH of optimum stability') or s.startswith('The pI'):
            print(s); mode = None; continue
        if mode == 'sum':
            p = s.split()
            if len(p) >= 4 and p[0] in MODEL:
                try: rows.append((p[0], p[1], p[2], float(p[3])))
                except ValueError: pass
        if mode == 'curve':
            p = s.split()
            if len(p) == 2:
                try: curve.append((float(p[0]), float(p[1])))
                except ValueError: pass
    return rows, curve

f = sys.argv[1]
rows, curve = parse(f)
print(f"\n=== {f} ===")
print(f"titratable residues: {len(rows)}\n")

shifted = [(r, r[3]-MODEL[r[0]]) for r in rows if abs(r[3]-MODEL[r[0]]) >= 1.5 and 0 <= r[3] <= 14]
shifted.sort(key=lambda x: -abs(x[1]))
print("strongly shifted (>=1.5 pH units from model):")
for (res, num, ch, pka), d in shifted[:25]:
    print(f"  {ch}/{res}{num:>4}  pKa {pka:5.2f}  (model {MODEL[res]:.2f}, shift {d:+.2f})")

if curve:
    print("\nfolding free energy vs pH (kcal/mol):")
    for ph, dg in curve:
        if abs(ph*2 - round(ph*2)) < 1e-6 and ph % 0.5 == 0:
            print(f"  pH {ph:4.1f}  {dg:7.2f}")
