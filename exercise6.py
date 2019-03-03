#Elisheva Heilbrun 208833012
#This program executes the NW algorithm for global sequence alignment.
import numpy as np
from string import *
import random
    
#---------------------------------------------------------------
#This function returns the values for case of match or mismatch.
#---------------------------------------------------------------
def Diagonal(n1,n2,pt):
    if(n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']

#-----------------------------------------------------------------------------------------------------------------   
#This function gets the optional elements of the aligment matrix and returns the elements for the pointers matrix.
#-----------------------------------------------------------------------------------------------------------------
def Pointers(di,ho,ve):
    
    pointer = max(di,ho,ve) #based on python default maximum(return the first element).

    if(di == pointer):
        return 'D'
    elif(ve == pointer):
        return 'V'
    else:
        return 'H'    
    
#------------------------------------------------------------------------------------
#This function creates the alignment and pointers matrices - executes the NW algorithm
#------------------------------------------------------------------------------------
def NW(s1,s2,match = 1,mismatch = -1, gap = -2):
    penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap} #A dictionary for all the penalty valuse.
    n = len(s1) + 1 #The dimension of the matrix columns.
    m = len(s2) + 1 #The dimension of the matrix rows.
    al_mat = np.zeros((m,n),dtype = int) #Initializes the alighment matrix with zeros.
    p_mat = np.empty((m,n),dtype = str) #Initializes the alighment matrix with zeros.
    #Scans all the first rows element in the matrix and fill it with "gap penalty"
    for i in range(m):
        al_mat[i][0] = penalty['GAP'] * i
        p_mat[i][0] = 'V'
    #Scans all the first columns element in the matrix and fill it with "gap penalty"
    for j in range (n):
        al_mat[0][j] = penalty['GAP'] * j
        p_mat [0][j] = 'H'
    #Fill the matrix with the correct values.

    p_mat [0][0] = 0 #Return the first element of the pointer matrix back to 0.
    for i in range(1,m):
        for j in range(1,n):
            di = al_mat[i-1][j-1] + Diagonal(s1[j-1],s2[i-1],penalty) #The value for match/mismatch -  diagonal.
            ho = al_mat[i][j-1] + penalty['GAP'] #The value for gap - horizontal.(from the left cell)
            ve = al_mat[i-1][j] + penalty['GAP'] #The value for gap - vertical.(from the upper cell)
            al_mat[i][j] = max(di,ho,ve) #Fill the matrix with the maximal value.(based on the python default maximum)
            p_mat[i][j] = Pointers(di,ho,ve)
    
    return al_mat, p_mat,i,j,al_mat[i][j]

#------------------------------------------
#This function prints the alignment matrix.
#------------------------------------------
def Print_Al_matrix(mat, s1, s2):
    s1 = " " + s1 #Add space to the header.
    s1 = list(s1) #Turns the header into list.
    s2 = " " + s2 
    s2 = list(s2) 
    s = mat.__repr__() #Gets the matrix values so the interpreter will be ablr to read it.
    s = s.split("array(")[1] #Deletes the word "array"
    s = s.replace("      ", "")#Removes all the unnecessary spaces.
    s = s.replace("[[", " [") #Deletes the double "["
    s = s.replace("]])", "]") #Deletes the double "]"
    pos = [i for i, ltr in enumerate(s.splitlines()[0]) if ltr == ","] #Gets the positions of the end of the lines.
    pos[-1] = pos[-1]-1 #Adjusts the last position
    empty = " " * len(s.splitlines()[0])
    s = s.replace("],", "]")
    s = s.replace(",", "") #Deletes the ",".
    lines = [] #Initializes list for the matrix lines.
    #Adds the vertical header to the matrix lines.
    for i, l in enumerate(s.splitlines()):
        lines.append(s2[i] + l)
    s  ="\n".join(lines) #Separates the lines.
    empty = list(empty)
    #Adds the horizontal header
    for i, p in enumerate(pos):
        empty[p-i] = s1[i]
    s = "".join(empty) + "\n" + s
    print s
    
#----------------------------------------
#This function prints the pointers matrix
#----------------------------------------
def print_P_matrix(mat,cols,rows):
    cols = " " + cols #Add space to the header.
    rows = " " + rows
    cols =  " " + cols.replace("","   ")
    print cols
    i=0
    
    for line in mat:
        print rows[i],str(line)
        i += 1

#----------------------------------------------------------------
#This function returns the path according to the pointers matrix.
#----------------------------------------------------------------
def Path(matrix,cur1,cur2,path):

    if(matrix[cur1][cur2] == '0'):#Stop condition.
        return path
    elif(matrix[cur1][cur2] == 'D'):
        path += 'D'
        return Path(matrix,cur1-1,cur2-1,path)
    elif(matrix [cur1][cur2] == 'V'):
        path += 'V'
        return  Path(matrix,cur1-1,cur2,path)
    elif(matrix[cur1][cur2] == 'H'):
        path += 'H'
        return Path(matrix,cur1,cur2-1,path)

#--------------------------------------------
#This function prints the pairwise alignment
#--------------------------------------------
def Pairwise(s1,s2,p):
    p1 = "" #initializes an empty string for seq1
    p2 = "" #initializes an empty string for seq2
    j = 0 #An index for gaps
    k = 0 #An index for gaps
    length = len(p)
    #Scanning the path and creates the alignment according to the path.
    for i in range (length):
        if(p[i] == 'D'):#Match/Mismatch.
            p1 += s1[j]
            p2 += s2[k]
            j += 1
            k += 1
        elif(p[i] == 'H'):#Gap - Horizontal.
            p1 += s1[j]
            p2 += "-"
            j += 1
        elif(p[i] == 'V'):#Gap - Vertex.
            p1 += '-'
            p2 += s2[k]
            k += 1
    print "Pairwise Alingment:\n"
    print p1
    print p2
    
#--------------------------------------------
#This function generates random dna sequence.
#--------------------------------------------
def Random_DNA(length):
    dnalist = ['A','T','C','G']
    DNA = ""
    for i in range(length):
        DNA += random.choice(dnalist)
    return DNA

#--------------------------------------------------
#This function generates random protein sequence.
#--------------------------------------------------
def Random_amino(length):
    aalist = ['M','T','Q','A','P','F','L','S','V','E','G','H','I','W','R','D','K','N','Y','C']
    Amino = ""
    for i in range(length):
        Amino += random.choice(aalist)
    return Amino

#-------------------------------------------------------
#This function return the score of 2 sequences identity.
#-------------------------------------------------------
def Match(seq1,seq2):
    match = 0
    #Scans the sequences and add 1 if there is a match.
    for i in range (len(seq1)): #(Both sequences have the same length).
        if(seq1[i] == seq2[i]):
            match += 1
    return match

#------------------------------------------------------------------------------
#This function prints the average results of many dna sequences score with gaps.
#(By default parameters)
#------------------------------------------------------------------------------
def average_DNA(ammount,size,match = 1,mismatch = 0,gap = 0):
    dna1, dna2 = "",""
    avg = 0 #The score of all the matches.
    #Generate pairs of random sequences:
    for i in range (ammount):
      dna1 =  Random_DNA(size)
      dna2 =  Random_DNA(size)
      mat1,mat2,x,y,score = NW(dna1,dna2,match,mismatch,gap) #sends to the NW algorithm.
      score = int(score)
      avg += score
    print "Averrage of ",ammount," DNA's sequences sequence alignment identity:"
    print float(avg) / ammount
    print"---------------------------------------------------------------------"
    
#-----------------------------------------------------------------------------------
#This function prints the average results of many protein sequences score with gaps.
#(By default parameters)
#-----------------------------------------------------------------------------------
def average_amino(ammount,size,match = 1,mismatch = 0,gap = 0):
    amino1, amino2 = "",""
    avg = 0 #The score of all the matches.
    #Generate pairs of random sequences:
    for i in range (ammount):
      amino1 =  Random_amino(size)
      amino2 =  Random_amino(size)
      mat1,mat2,x,y,score = NW(amino1,amino2,match,mismatch,gap) #sends to the NW algorithm.
      score = int(score)
      avg += score
    print "Averrage of ",ammount," protein's sequences sequence alignment identity:"
    print float(avg) / ammount
    print"---------------------------------------------------------------------"

#---------------------------------------------------------------------------------
#This function prints the average results of many dna sequences score without gaps.
#---------------------------------------------------------------------------------
def Match_DNA(ammount,size):
    dna1,dna2 = "","" #Initializes the sequences.
    avg = 0 #The score of all the matches.
    #Generate pairs of random sequences:
    for i in range (ammount):
      dna1 =  Random_DNA(size)
      dna2 =  Random_DNA(size)
      score = Match(dna1,dna2)#Sends to the "MATCH" function.
      avg += score
    print "Averrage of ",ammount," DNA's sequences matching identity:(without gaps)"
    print float(avg) / ammount
    print"---------------------------------------------------------------------"
    
#--------------------------------------------------------------------------------------
#This function prints the average results of many protein sequences score without gaps.
#--------------------------------------------------------------------------------------
def Match_amino(ammount,size):
    amino1,amino2 = "",""
    avg = 0 #The score of all the matches.
    #Generate pairs of random sequences:
    for i in range (ammount):
      amino1 =  Random_amino(size)
      amino2 =  Random_amino(size)
      score = Match(amino1,amino2)#Sends to the "MATCH" function.
      avg += score
    print "Averrage of ",ammount," protein's sequences matching identity:(without gaps)"
    print float(avg) / ammount
    print"---------------------------------------------------------------------"
#------------------------------------------------------------------------------
#Main
#------------------------------------------------------------------------------
choice = raw_input("If you want to analyse sequences from file,press 1\nIf you want to compare randomal sequences, press 2: ")
if(choice == '1'): #You want to analyse from file.
    name = raw_input("Please Enter the name of the file where your sequences are written: ")
    s = open(name,"r")
    info = s.readline() #Removes the information line from the sequence.
    seq1,seq2 = "",""
    for line in s:
        if(line.startswith(">")): #The second information line.
            break
        else:
            seq1 += line
    for line in s:
        seq2 += line
    seq1 = seq1.replace("\n","") #Removes all the end of lines.
    seq2 = seq2.replace("\n","") #Removes all the end of lines.

    default = raw_input("If you want to insert paramters for sequence alignment,press 1.\nIf you want to use default parameters, press 2: ")
    if(default == '1'): #You want to insert parameters.
        print "Please enter parameters for the global sequence alighment:"
        match = int( raw_input("Enter value for match: "))
        mismatch = int(raw_input("Enter value for mis-match penalty: "))
        gap = int(raw_input("Enter value for gap penalty: "))
        matrix1,matrix2,x,y,score = NW(seq1,seq2,match,mismatch,gap)
    elif(default == '2'): #You want to use the default paramters.
        matrix1,matrix2,x,y,score = NW(seq1,seq2)
    else:
        print "ERROR\n"

    print "\nAlingment Matrix:\n"
    Print_Al_matrix(matrix1,seq1,seq2)
    print "\nPointers Matrix:\n"
    print_P_matrix(matrix2,seq1,seq2)
    path = "" #Initializes a string parameter for the path.
    p =  Path(matrix2,x,y,path)
    p = p[::-1] #Reverse the accepted path - that will start from the beginning.
    print "\nPath:\n"
    print p + "\n"
    Pairwise(seq1,seq2,p)

elif(choice == '2'): #You want to compare matches with gaps and without gaps.
    ch = raw_input("If you want to analyse dna sequences, press d\nIf you want to analyse protein sequences press p: ")
    if(ch == 'd'):
        ammount = int(raw_input("Enter how many pairs of sequences you want to raffle:  "))
        length = int(raw_input("Enter the length of the sequences: "))
        print "\n"
        average_DNA(ammount,length)
        Match_DNA(ammount,length)
    elif(ch == 'p'):
        ammount = int(raw_input("Enter how many pairs of sequences you want to raffle:  "))
        length = int(raw_input("Enter the length of the sequences: "))
        print "\n"
        average_amino(ammount,length)
        Match_amino(ammount,length)
    else:
        print "ERROR\n"    
else:
    print "ERROR\n"
