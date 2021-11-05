from math import log10

def random_genome(dna, gc_content):
    GC = 0
    n = len(dna)
    for c in dna:
        if c == 'G' or c == 'C':
            GC += 1
    AT = n - GC
    dna_calc = []


    for val in gc_content:
        p = GC * log10(val / 2) + AT * log10((1 - val) / 2)
        dna_calc.append(round(p, 3))
    return dna_calc

def random_genome1(dna, gc_content):
    dna = dna.upper()
    cg = len(dna.replace('A', '').replace('T', ''))
    at = len(dna.replace('G', '').replace('C', ''))
    print(cg, at)
    result = []
    for i in range(len(gc_content)):
        prob = cg * log10(float(gc_content[i]) / 2) + at * log10((1 - float(gc_content[i])) / 2)
        result.append(round(prob, 3))
    return result

print(random_genome("ACGATACAA", [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]))
print(random_genome1("ACGATACAA", [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]))
