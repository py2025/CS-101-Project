from math import factorial, log10

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

#helper method
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

print(rev_palindrome("TCAATGCATGCGGGTCTATATGCAT"))
