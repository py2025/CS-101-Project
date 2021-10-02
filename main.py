genetic_code = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

#@param dna A string representing a DNA sequence
#returns a dictionary with the count of each nucleotide in the dna string
def s(dna):
    nucleotides = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for l in dna.upper():
        try:
            nucleotides[l] = nucleotides.get(l) + 1
        except (KeyError, TypeError):
            continue

    return nucleotides

#@param dna A string representing a DNA sequence
#returns the transcribed RNA string
def dna2rna(dna):
    rna = ''
    for n in dna.upper():
        if n == 'T':
            rna += 'U'
        else:
            rna += n

    return rna

#@param dna A string representing a DNA sequence
#returns the reverse complement of the DNA sequence
def reverse_complement(dna):
    r_dna = dna.upper()[::-1]
    complement = ''
    for n in r_dna:
        if n == 'T':
            complement += 'A'
        elif n == 'A':
            complement += 'T'
        elif n == 'C':
            complement += 'G'
        elif n == 'G':
            complement += 'C'
        else:
            complement += n

    return complement

#@param hom The number of organisms that are homozygous dominant
#@param het The number of organisms that are heterozygous
#@param rec The number of organisms that are homozygous recessive
#returns the probability for an offspring of two randomly selected organisms to possess a dominant allele
def mendels_law(hom, het, rec):
    totalal = hom + het + rec

    #array of discrete probabilities for every couple
    couples = [
        (hom * (hom - 1)) / (total - 1), #2 homozygous dom organisms
        (2 * hom * het) / (total - 1), #1 homozygous dom and 1 heterozygous organism
        (2 * hom * rec) / (total - 1), #1 homozygous dom and 1 homozygous rec organism
        (het * (het - 1) * 0.75) / (total - 1), #2 heterozygous organisms
        (het * rec) / (total - 1) #1 heterozygous and 1 homozygous rec organism
    ]

    return sum(couples) / total

#@param n Number of months
#@param k Number of rabbit pairs produced in the next month
#returns number of rabbit pairs present after n months
def fibonacci_rabbits(n, k):
    if n == 1 or n == 2:
        return 1
    else:
        return fibonacci_rabbits(n - 1, k) + fibonacci_rabbits(n - 2, k) * k

#@param dna_list List of DNA strings
#returns index of DNA string with highest GC-content and its GC-content percentage as a tuple
def gc_content(dna_list):
    GC_count = []
    total = []
    index = 0
    for dna in dna_list:
        GC = 0
        tot = 0
        for c in dna:
            if c == 'C' or c == 'G':
                GC += 1
                tot += 1
            else:
                tot += 1

        GC_count.append(GC)
        total.append(tot)

    max_GC = GC_count[0]
    index = 0
    for i in range(1, len(GC_count)):
        if GC_count[i] >= max_GC:
            max_GC = GC_count[i]
            index = i

    return (index, GC_count[index] / total[index] * 100)

#@param rna A string representing an RNA sequence
#returns corresponding amino acid string
def rna2codon(rna):
    amino = ''

    for i in range(0, int(len(rna) / 3)):
        try:
            print((rna[3 * i:3 * i + 3]))
            amino += genetic_code[(rna[3 * i:3 * i + 3])]
        except KeyError:
            continue

    return amino

#@param dna_snippet A substring of dna
#@param dna A string representing a DNA sequence
#returns all locations of the substring as a list of integers
def locate_substring(dna_snippet, dna):
    length = len(dna_snippet)
    indices = []
    for i in range(0, len(dna) - length + 1):
        if dna[i:i + length] == dna_snippet:
            indices.append(i)
    return indices

#@param dna1 A string representing a DNA sequence
#@param dna2 A string representing a DNA sequence
#returns the Hamming distance between the two strings
def hamming_dist(dna1, dna2):
    total = 0
    for i in range(0, len(dna1)):
        if dna1[i] != dna2[i]:
            total += 1

    return total

#@param genotypes A list with the number of couples in a population with a specific genotype pairing
#returns expected number of offspring displaing the dominant phenotype, assuming each couple has two offspring
def count_dom_phenotype(genotypes):
    e_val = 0
    couples = [
        2.0, #2 homozygous dom organisms
        2.0, #1 homozygous dom and 1 heterozygous organism
        2.0, #1 homozygous dom and 1 homozygous rec organism
        1.5, #2 heterozygous organisms
        1.0, #1 heterozygous and 1 homozygous rec organism
        0.0 #2 homozygous rec organisms
    ]

    for i in range(len(genotypes)):
        e_val += couples[i] * genotypes[i]

    return e_val

#@param protein A string representing a sequence of proteins
#returns total number of different source RNA strings
def source_rna(protein):
    total_c = []
    total = 3
    for c in protein:
        count = 0
        for k, v in genetic_code.items():
            if(v == c):
                count += 1
        total_c.append(count)
        count = 0

    for v in total_c:
        total *= v

    return total

#@param dna A string representing a DNA sequence
#@param intron_list A list of strings representing introns
#returns a protein string
def splice_rna(dna, intron_list):
    rna = dna2rna(dna)

    for intron in intron_list:
        rna = rna.replace(dna2rna(intron), '')

    protein_string = rna2codon(rna)

    return protein_string

