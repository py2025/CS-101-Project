from math import factorial

#@param dna_motif A subsequence of DNA
#@param dna A DNA sequence
#returns the index positions of the motif characters for the first occurence of it within dna
'''
def find_splice(dna_motif, dna):
    char_list = []
    [char_list.append(c) for c in dna_motif]
    ind_1 = 0 #index to start iterating through char_list
    for char in char_list:
        if char not in dna:
            return []
        pos = dna.find(char)
        pos_list.append(pos)
'''

def find_splice(dna_motif, dna):
    char_list = []
    [char_list.append(c) for c in dna_motif]
    pos_list = []
    deleted_dna = ''
    for char in char_list:
        if char not in dna:
            return []
        pos = dna[1:].find(char)
        pos_list.append(pos + len(deleted_dna) + 1) #accounts for indices lost from deletion
        deleted_dna += dna[:pos]
        dna = dna[pos:]
        print(dna)
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
#@param str_list A list of strings to concatenate with the substring they share
#returns the concatenated strings
def cat_strings(str_list):
    new_str_2 = ''
    if str_list[0].find(str_list[2]) == 0:
        new_str_2 += str_list[0].replace(str_list[2], '')
        return str_list[1] + new_str_2
    new_str_2 += str_list[1].replace(str_list[2], '')
    return str_list[0] + new_str_2

#helper method
#@param list1 a list of strings to check for matches with prefixes and suffixes
#returns concatenated strings with matching prefixes and suffixes
def check_for_match(list1):
    test_arr = []
    #O(n^4) Who doesn't love quadric time complexity?
    for dna in list1:
        for i in range(len(dna)):
            sub = dna[:i]
            for dna1 in list1:
                largest_common = 0
                for j in range(-2, -len(dna1) - 2, -1):
                    if dna1 == dna:
                        continue
                    if dna1[-1:j:-1][::-1] == sub and len(dna1[-1:j:-1][::-1]) > largest_common:
                        test_arr.append((dna, dna1, dna1[-1:j:-1][::-1]))
                        largest_common = len(dna1[-1:j:-1][::-1])
    list_of_cats = []
    for l in test_arr:
        list_of_cats.append(cat_strings(l))
    return list_of_cats

#@param dna_list A list of DNA strings
#returns the shortest superstring containing all given DNA strings
def assemble_genome(dna_list):
    no_dupes = []
    #second iteration concatenates strings again to ensure every permutation of strings is accounted for
    second_it = check_for_match(check_for_match(dna_list))
    [no_dupes.append(x) for x in second_it if x not in no_dupes]
    for i in range(len(no_dupes)):
        for s in dna_list:
            try:
                if no_dupes[i].find(s) == -1:
                    del no_dupes[i]
            except IndexError:
                continue
    min_len = len(no_dupes[0])
    index = 0
    for i in range(len(no_dupes)):
        if len(no_dupes[i]) < min_len:
            index = i
            min_len = len(no_dupes[i])
    return no_dupes[index]

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
