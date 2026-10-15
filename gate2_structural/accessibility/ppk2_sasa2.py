import freesasa, itertools

def total(pdb):
    return freesasa.calc(freesasa.Structure(pdb)).totalArea()

single = {c: total(f'chain_{c}.pdb') for c in 'ABCD'}

print("\n=== Interface areas (buried per side, Å²) ===")
for c1, c2 in itertools.combinations('ABCD', 2):
    # build the pair file
    with open(f'pair_{c1}{c2}.pdb','w') as f:
        for c in (c1,c2):
            f.write(open(f'chain_{c}.pdb').read())
    pair = total(f'pair_{c1}{c2}.pdb')
    bsa = (single[c1] + single[c2] - pair) / 2
    print(f"{c1}-{c2}  {bsa:7.0f}")
