#!/usr/bin/env python
# coding: utf-8

# In[ ]:



mass_table = {"G": 57.021464,
        "A": 71.037114,
        "S": 87.032028,
        "P": 97.052764,
        "V": 99.068414,
        "T": 101.047678,
        "C": 103.009184,
        "I": 113.084064,
        'L': 113.084064,
        'N': 114.042927,
        'D': 115.026943,
        'Q': 128.058578,
        'K': 128.094963,
        'E': 129.042593,
        'M': 131.040485,
        'H': 137.058912,
        'F': 147.068414,
        'R': 156.101111,
        'Y': 163.063329,
        'W': 186.079313}

gencode = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGC': 'C', 'UGU': 'C', 'UGA': 'Stop', 'UGG': 'W'}



def read_fasta(file_path):
    from itertools import groupby
    with open("file_path) as f:
        groups = groupby(f, key=lambda x: not x.startswith(">"))
        d = {}
        for k,v in groups:
        if not k:
            key, val = list(v)[0].rstrip(), "".join(map(str.rstrip,next(groups)[1],""))
            d[key] = val
    f.close()
    return (d)   

def calc_mass(prot_seq):
     
    result = 0
    for i in s:
        result += table[i]
    
    return(result)


def orf(rna_seq):

    
    result = []
    
    lst_start_codons = []
    
    decoded = []
    
    
    #finding start codons
    rna1 = rna_seq
      
    n_starts = rna_seq.count("AUG")
    
    while n_starts > 0:
        x = rna1.find("AUG")
        rna1 = rna1[x+3:]
        y = len(rna_seq) - len(rna1) - 3 
        lst_start_codons.append(y)
        n_starts -= 1
            
    
    #finding stop codons
    rna2 = rna_seq
    n_stops = rna_seq.count("UAA") + rna_seq.count("UAG") + rna_seq.count("UGA")
    lst_stop_codons = []
    
    while n_stops > 0:
        n_stops -= 1
        nearest_stops = []
        uaa = rna2.find("UAA")
        uag = rna2.find("UAG")
        uga = rna2.find("UGA")
        if uaa != -1:
            nearest_stops.append(uaa)
        if uag != -1:
            nearest_stops.append(uag)
        if uga != -1:
            nearest_stops.append(uga)
        if not nearest_stops:  #empty collections are false, check if empty
            break
        
        closest = min(nearest_stops)
        rna2 = rna2[closest + 3:]
        stop_codon = len(rna_seq) - len(rna2)
        lst_stop_codons.append(stop_codon)
        
    
    #finding start codons with adequate stop codons
    lst_start_stop = []
    lst_external_start = lst_start_codons.copy()
    
    for start in lst_start_codons:
        if start > max(lst_stop_codons):
            lst_external_start.remove(start)
            continue
        for stop in lst_stop_codons:
            if start < stop and (stop - start) % 3 == 0:
                lst_start_stop.append([start,stop])
                break  #we stop the process after finding the nearest stop codon
    
    #finding internal start codons
    internal_start = []
    
    for i in lst_start_stop:
        rna3 = rna_seq
        rna3 = rna3[i[0]:i[1]]
        
        for j in range(3,len(rna3),3): #start from 3 to exclude starting codon
            if gencode[rna3[j:j+3]] == "M":
                internal_start.append(j+i[0])
                
                
    
    #only start and stop with external start codons
    
    lst_start_stop_final = lst_start_stop.copy()
    
    for i in lst_start_stop:
        if i[0] in internal_start:
            lst_start_stop_final.remove(i)
    
    
    #decoding strands      
            
    for i in lst_start_stop_final:
        rna3 = rna_seq
        rna3 = rna3[i[0]:i[1]]
        
        for i in range(0, len(rna3),3):
            if gencode[rna3[i:i+3]] == "Stop":
                break
                
            decoded.append(gencode[rna3[i:i+3]])
           
        protein = "".join(decoded)
        result.append(protein)
        decoded.clear()
        
        

        
        
   
    result.sort()
    return (result)

def translate(rna_seq):

    start_codon = 0
    decoded = []
    for i in range(len(rna_seq)):
        if gencode[rna_seq[i:i+3]] == "M":
            start_codon = int(i)
            break
    for i in range(start_codon + 3, len(rna_seq), 3):
        if gencode[rna_seq[i:i+3]] == "Stop":
            break
        else:
            decoded.append(gencode[rna_seq[i:i+3]])

    
    protein = "".join(decoded)
    return (protein)
    

