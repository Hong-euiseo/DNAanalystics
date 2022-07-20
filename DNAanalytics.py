def revSeq(seq):
    data = seq[::-1]
    compDic = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    result = ''
    for i in data:
        result += compDic[i]
    return(result)

def DNAtoRNA(revComp_result):
    seq_RNA = revComp_result.replace('T', 'U')
    return seq_RNA

codonTable = {
"UUU" : "F",
"CUU" : "L",
"AUU" : "I",
"GUU" : "V",
"UUC" : "F",
"CUC" : "L",
"AUC" : "I",
"GUC" : "V",
"UUA" : "L",
"CUA" : "L",
"AUA" : "I",     
"GUA" : "V",
"UUG" : "L",
"CUG" : "L",
"AUG" : "M",
"GUG" : "V",
"UCU" : "S",
"CCU" : "P",
"ACU" : "T",
"GCU" : "A",
"UCC" : "S",
"CCC" : "P",
"ACC" : "T",
"GCC" : "A",
"UCA" : "S",
"CCA" : "P",
"ACA" : "T",
"GCA" : "A",
"UCG" : "S",
"CCG" : "P",
"ACG" : "T",
"GCG" : "A",
"UAU" : "Y",
"CAU" : "H",
"AAU" : "N",
"GAU" : "D",
"UAC" : "Y",
"CAC" : "H",
"AAC" : "N",
"GAC" : "D",
"UAA" : "Stop",
"CAA" : "Q",
"AAA" : "K",
"GAA" : "E",
"UAG" : "Stop",
"CAG" : "Q",
"AAG" : "K",
"GAG" : "E",
"UGU" : "C",
"CGU" : "R",
"AGU" : "S",
"GGU" : "G",
"UGC" : "C",
"CGC" : "R",
"AGC" : "S",
"GGC" : "G",
"UGA" : "Stop",
"CGA" : "R",
"AGA" : "R",
"GGA" : "G",
"UGG" : "W",
"CGG" : "R",
"AGG" : "R",
"GGG" : "G"
}

def RNAtoProtein(RNA_result):
    protein = ''
    if len(RNA_result)%3 == 0:
        for i in range(0, len(RNA_result), 3):
            codon = RNA_result[i:i+3]
            protein += codonTable[codon]
    return protein

def  main():
    # seq = "AAATTTGGGCGC"
    # 1. input sequence data
    seq = input("Enter the seqence information! : ")
    print('1. Input sequence information : ', seq)

    # 2. Complementary seqence data
    revComp_result = revSeq(seq)
    print('2. Complementary sequence information : ', revComp_result)
    
    # 3. DNA to RNA seq
    RNA_result = DNAtoRNA(revComp_result)
    print('3. RNA sequence information : ', RNA_result)

    # 4. RNA to Protein 
    Protein_result = RNAtoProtein(RNA_result)
    print('4. RNA to Protein sequence information: ', Protein_result)

    # 5. Save all results
    file_name = input('Please Enter the save file name here : ')
    file_name2 = str(file_name) + '.txt'
    file_result = open(file_name2, 'w')
    file_result.write('1. Input sequence information : ' + seq + '\n')
    file_result.write('2. Complementary sequence information : ' + revComp_result + '\n')
    file_result.write('3. RNA sequence information : ' + RNA_result + '\n')
    file_result.write('4. Protein sequence information : ' + Protein_result + '\n')
    file_result.close()
    print("Results saved on file to : ", file_name2)

main()

