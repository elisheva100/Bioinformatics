#Elisheva Heilbrun 208833012
#This program translates dna sequence into proteins.
from string import *

#This function generates the complemntary strand.
def Complementary_strand(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #dictionary for the complementary nucleotides.
    comp = "" #initializes the new sequence.
    #Build the complement sequence.
    for i in range(len(seq)):
        for k,v in complement.items():
            if (k == seq[i]):
                comp += v
  
    comp = comp[::-1] #Reverses the sequence
    
    return comp

#-------------------------------------------------------
#This function finds all the open reading frames on the sequence.
def Orf(seq):
    minimum = 300
    start_codon = 0 #initializes the start codon to poition 0. 
    stop_codon = 0 #initializes the stop codon to poition 0.
    stops = ["TAA","TGA","TAG"]
    orfs = [[],[],[]] #List for all the orfs.
    
    #Reads the equence on all the optional frames.
    for i in range(3):
        #Scans for all the codons
        for j in range(i,len(seq),3):
            codon = seq[j:j+3]#The codon is 3 nucleotides.
            if(start_codon == 0 and codon == "ATG"): #If the codon is a start codon and doesn't "contained" in other start codon
               start_codon = j+1 #The position of the start codon.
              
            if(start_codon > 0 and codon in stops): #if the codon is a stop  and codon it's start codon matches to it.
               stop_codon = j+3 #The position of the end of the stop codon.(according to the ORF finder)
                                  
            if(stop_codon > 0 and  start_codon > 0 ): #Matches the start and stop codons
                if(stop_codon - start_codon >= minimum): #If the is at least 300 nucleotides.
                    pair = start_codon,stop_codon #insert the matching codons into tuple.
                    orfs[i].append(pair) #Add to the list the matching pair.
                start_codon = 0 #initializes the start codon back to -0.
                stop_codon = 0 #initializes the stop codon back to -0.

    return orfs

#-----------------------------------------------------------------------------"
#This function prints  the orf sequence and it's details.
def Printing(orf,seq,index,sign):
    #ORF_seq = [] #List for the sequence of all the nucleotides.
    orf_seq = ""
    start = int(orf[0]) #The position of the start codon
    stop = int(orf[1]) #The position of the stop codon
    orf_seq += seq[start-1:stop] #Add to the string all the nucleotides on the orf
                    
    print "strand:" + str(sign) + ", frame:" + str(index) + ", start: " + str(start) + ", stop: " + str(stop) + ", length: " + str(stop-start +1) + " " + seq[start-1:start-1+6] + " " + seq[stop-6:stop]
    print "\nORF:\n"
    print orf_seq
    return orf_seq

#--------------------------------------------------------
#The function translates the dna sequence into protein according to the ORFs.
def Codons2Protein(orf):
    #Dictionary that contained all the amino acids encoding.
    code ={"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    protein = "" #Initializes the new protein sequence.
    #Reads 3 codons for encoding to amino acid.
    for i in range(0,len(orf),3):
        codon = orf[i:i+3]
        protein += code[codon] #The matchin amino acid.
    
    return protein
        
            

#Main
name = raw_input("Please enter the name of the sequence: ")
print "\n"
s = file(name,"r").read()
if(s.find(">") > -1):
    s = s[s.find("\n") + 1:] #Deletes the information line.
s = s.replace("\n","") #Removes all the end of lines.
s = s.upper() #Converts all the letters into capital letters.
orf = Orf(s) #All the ORFs on the leading strand.
index_l = 1 #Index for the leading strand.

#Print all the leding strand's orfs,they details and the traslation into proteins.
for i in range(len(orf)):

    if(not (orf[i])): #There is no orf/s on this specific frame.
        print "There is no orf/s on frame" + str(i+1) + " on the leading strand\n"
 
    else:
         for j in range(len(orf[i])):
            print"---------------------------------------------------------------"
            print "ORF" + str(index_l) + "\n"
            orfs_s = Printing(orf[i][j],s,i+1,"+")
            print "\nProtein:\n"
            print Codons2Protein(orfs_s) +"\n"
            index_l += 1 #Promotes the index
            
c_s = Complementary_strand(s) #gets the complementary strand.
c_orf = Orf(c_s) #All the ORFs on the complementary strand.
index_c = index_l #Index for the complementary strand.

#print all the complementary strand's orfs,they details and the traslation into proteins.
for i in range(len(orf)):

    if(not (c_orf[i])): #There is no orf/s on this specific frame.
        print "There is no orf/s on frame" + str(i+1) + " on the complementary strand\n"
 
    else:
         for j in range(len(c_orf[i])):
            print"---------------------------------------------------------------"
            print "ORF" + str(index_c) + "\n"
            orfs_s = Printing(c_orf[i][j],c_s,i+1,"-")
            print "\nProtein:\n"
            print Codons2Protein(orfs_s) +"\n"
            index_c += 1 #Promotes the index


