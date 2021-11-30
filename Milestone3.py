from math import factorial, log10

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
    total = hom + het + rec

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
def GC_content(dna_list):
    GC_count = {}
    index = 0
    max_GC = 0

    for dna in dna_list:
        GC = 0
        tot = 0
        for c in dna:
            if c == 'C' or c == 'G':
                GC += 1
            tot += 1
        GC_count[dna] = GC / tot * 100

    for key in GC_count:
        if GC_count[key] > max_GC:
            index = dna_list.index(key)
            max_GC = GC_count[key]

    return (index, max_GC)

#@param rna A string representing an RNA sequence
#returns corresponding amino acid string
def rna2codon(rna):
    amino = ''

    for i in range(0, int(len(rna) / 3)):
        try:
            if(genetic_code[(rna[3 * i: 3 * i + 3])] == '*'):
                break
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

#@param dna_motif A subsequence of DNA
#@param dna A DNA sequence
#returns the index positions of the motif characters for the first occurence of it within dna
def find_splice(dna_motif, dna):
    char_list = []
    [char_list.append(c) for c in dna_motif]
    pos_list = []
    current_index = -1 #index of last found character in DNA
    for char in char_list:
        indices = [i for i, c in enumerate(dna) if c == char] #finds all occurrences of c in dna
        for j in indices:
            if indices[-1] <= current_index:
                return []
            if j > current_index:
                pos_list.append(j)
                current_index = j
                break
    return pos_list

#@param dna_list A list of dna strings
#returns a string of dna bases common to the strings
def shared_motif(dna_list):
    dna1 = dna_list[0]
    common_motifs = []

    for base1 in range(len(dna1)): #start with one letter in the fragment
        for base2 in range(base1 + 1, len(dna1) + 1): #move on to next letter of fragment
            tempstring = dna1[base1:base2] #create a variable for the motif; not final, just to search through other fragments
            Available = True
            for motif1 in range(1, len(dna_list)): #look through other fragments for motif
                if tempstring not in dna_list[motif1]:
                    Available = False
                    break
            if Available:
                common_motifs.append(tempstring)
    final_motif = ""
    for substring in common_motifs:
        if len(substring) > len(final_motif): #make sure the returned motif is the longest one
            final_motif = substring
    return final_motif

#@param dna_dict A dictionary that maps DNA strings to their ROSALIND identifiers
#returns an adjacency list
def get_edges(dna_dict):
    adj_list = []
    #O(n^2)
    for key in dna_dict.keys():
        #this line is absolutely gross, but it gets the last three characters of the DNA string
        curr_end = dna_dict[key][-1:-4:-1][::-1]
        for key1 in dna_dict.keys():
            if key1 == key:
                continue
            curr_start = dna_dict[key1][0:3:1]
            if curr_end == curr_start:
                adj_list.append((key, key1))
    return adj_list

#helper method
#@param s1 A string
#@param s2 A string
#@param out_str A string to store resulting string after overlap
#returns maximum overlap of s1 and s2 and the resulting string
def find_overlaps(s1, s2, out_str):
    max_overlap = 0
    n = min([len(s1), len(s2)])

    for i in range(1, n + 1):
        if s1[-i:] == s2[:i]:
            if max_overlap < i:
                max_overlap = i
                out_str = s1 + s2[i:]

    for i in range(1, n + 1):
        if s2[-i:] == s1[:i]:
            if max_overlap < i:
                max_overlap = i
                out_str = s2 + s1[i:]
    return (max_overlap, out_str)

#@param dna_list A list of dna
#returns the shortest superstring of dna in dna_list
def assemble_genome(dna_list):
    n = len(dna_list)

    while n != 1:
        max_overlap = 0
        p = -1 #p and q are indices involved in max overlap
        q = -1
        res_str = ''
        temp_str = '' #keeps track of resulting string after overlap

        for i in range(n):
            for j in range(i + 1, n):
                overlap, temp_str = find_overlaps(dna_list[i], dna_list[j], temp_str)

                if max_overlap < overlap:
                    max_overlap = overlap
                    res_str = temp_str
                    p = i
                    q = j
        n -= 1

        if max_overlap == 0:
            dna_list[0] += dna_list[n]
        else:
            dna_list[p] = res_str
            dna_list[q] = dna_list[n]

    return dna_list[0]

#@param rna An RNA string
#returns the total possible number of perfect matchings of nucleotide bases
def perfect_match(rna):
    char_count = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    for c in rna.upper():
        try:
            char_count[c] = char_count.get(c) + 1
        except (KeyError, TypeError):
            continue
    if char_count['A'] != char_count['U'] and char_count['C'] != char_count['G']:
        return 0
    return factorial(char_count['A']) * factorial(char_count['C'])

#@param dna A DNA string
#@param gc_content An array of GC-content values
#returns the common logarithm of probability that a random string will match dna
def random_genome(dna, gc_content):
    GC = 0
    for c in dna:
        if c == 'G' or c == 'C':
            GC += 1
    AT = len(dna) - GC
    dna_calc = []

    for val in gc_content:
        p = GC * log10(val / 2) + AT * log10((1 - val) / 2)
        dna_calc.append(round(p, 3))
    return dna_calc

#@param dna A DNA string
#returns a list of (position, length) tuples
def rev_palindrome(dna):
    pos_len_list = []
    for i in range(len(dna)):
        for j in range(4, 13):
            substring = dna[i:i + j]
            if i + j > len(dna):
                continue
            if reverse_complement(substring) == substring:
                pos_len_list.append((i, j))
    return pos_len_list

#@param file_name A file to be opened
#returns a list of dna strings in the file
def load_file(file_name):
    dna_list = []
    with open(file_name) as dna_file:
        for dna in dna_file:
            dna_list.append(dna.strip('\n'))
    dna_file.close()
    return dna_list

#@param dna_list A list of dna with 8 characters of overlap between strings
#returns the shortest superstring of dna in dna_list with 8 character overlap
def assemble_genome2(dna_list):
    temp_dict2 = {}
    superstring = ''
    for dna in dna_list:
        matching_pre = False #flag for if matching prefix is found for each string
        for comp_dna in dna_list:
            if dna == comp_dna:
                continue
            else: #checks for suffix and prefix match of 8 characters via brute force
                if comp_dna[-1:-9:-1][::-1] == dna[:8]:
                    matching_pre = True
                    temp_dict2[comp_dna] = dna
        if matching_pre == False:
            superstring = dna
            prev_key = dna
    for i in range(len(temp_dict2)): #iterate through dictionary with every value becoming the next key
        superstring += temp_dict2[prev_key][8:]
        prev_key = temp_dict2[prev_key]
    return superstring
