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

print(GC_content(["CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG","CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC","CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"]))
