#This program executes some Statistical Functions on sequences.

from string import *
import random

#This function Generate Random Dna Sequence Data With equal base frequencies.
def dna(length):

    dna_list = [n for n in ''.join([ 'atgc' for i in range(length/4)])] #There are only 4 nucleotides. 
    random.shuffle(dna_list) #Shuffle the list.
    DNA = ''.join(dna_list) #Turn the list into string.
    #Writes the sequence into new file.
    f = open("Random_DNA.txt",'w')
    f.write(DNA)
    f.close()
    
    return "Random_DNA.txt"

#--------------------------------------------
#This function returns the frequency of every nucleotid in the sequence.
def bases_freq(dna_seq_file):

    frequ = {} #Creates empty dictionary.
    nucs = ['a','t','c','g'] #Initializes the list with all the nucleotides.
    freq = [0,0,0,0] #Initializes the list for finding the frequency in each part
    length = 0 #Sum of all the parts length.

    while(True):
        s = dna_seq_file.read(10000) #divieds the file into blocks.
        if(not s): #The file is over
            break
        if(s.find(">")> -1):
            s = s[s.find("\n")+1:] #Delete the information.
        s = s.replace("\n","") #Delete the end of lines.
        s = s.replace("N","") #Delete all the 'N' of chromosome 22.
        s = s.lower() #Turns all the letters into small letters.
        length += len(s)
        #Add to the list the frequency of every nucleotide in each block.
        for i in range(len(freq)):
            freq[i] += s.count(nucs[i])
    dna_seq_file.seek(0) #Initializes the file's obgect to the beggining.
    #print length
    #Add to the dictionary the frequency of each nucleotide in the whole sequence
    for i in range(len(nucs)):
        frequ[nucs[i]] = float(freq[i])/length
    frequ['gc'] = frequ['g'] + frequ['c'] #Add the "gc" content
    return frequ

#-------------------------------------------------------
#This function prints the nucleotides frequency as percentags.
def calculate (dic):

    nuc = ['a','t','c','g','gc'] #List of all the nucleotides names.
    pe = [] #Creates empty list for the frequences.
    #Append to the list the values by order:
    for i in range(len(nuc)):
        for k,v in dic.items():
            if(nuc[i] == k):
                pe.append("%.1f" %((v*100))) #Convert to percentages.              

    print "  a   t    c    g    gc"
    print " ".join(str(element) for element in pe)
    return
#-------------------------------------------
#This function creates new file with only the nucleotides of the sequence of the given file.
def new_file(dna_seq_file):

    f = open("dna_file.txt",'w')#Creat new file
    while(True):
        s = dna_seq_file.read(10000) #divieds the file into blocks
        if(not s): #The file is over
            break
        if(s.find(">")> -1):
            s = s[s.find("\n")+1:] #Delete the information.
        s = s.replace("\n","") #Delete the end of lines.
        s = s.replace("N","") #Delete all the 'N' of chromosome 22.
        s = s.lower() #Turns all the letters into small letters.
        f.write(s)
    dna_seq_file.seek(0) #Initializes the file's obgect to the beggining.
    f.close()
    return "dna_file.txt"

#------------------------------------------------------
#This fuction returns values for drawing a graph of "GC" content in sliding window.
def slidingwindowplot(dna_file,window_length):

    x = [] #X axis coordinates.
    y = [] #Y axis coordinates.
    name = new_file(dna_file)#The same file without information and end of lines.
    f = open(name,'r')
    #Reads the file acording to the window's size.
    while(True):
        s = f.read(int(window))
        if(not s):
            break
        content = float(s.count('g'))/len(s) + float(s.count('c'))/len(s)#"gc" conten in each window
        y.append(content)
        x.append(str(f.tell()).strip("L"))#The location of the file.
        
    return x,y

#--------------------------------------------------
#This function returns the ratio between the oberved and the expected frequencies of a word.
def representation(word, dna_seq_file):
    expect = 1.0 #paramter for calculating the expected frequency
    lst = [] #Creates empty list 
    freq = bases_freq(dna_seq_file) #Returns the frequencies of all the nucleotide
    observ = 0 #Paramter to find the observed frequency.
    length = 0 #Sum of all the parts length .
    while(True):
        s = dna_seq_file.read(10000) #divieds the file into blocks.
        if(not s): #The file is over.
            break
        if(s.find(">")> -1):
            s = s[s.find("\n")+1:] #Delete the information.
        s = s.replace("\n","") #Delete the end of lines.
        s = s.replace("N","") #Delete all the 'N' of chromosome 22.
        s = s.lower() #Turns all the letters into small letters.
        length += len(s) 
        observ += s.count(word)
    print(observ)
    observ = float(observ)/length
    dna_seq_file.seek(0)#Initializes the file's obgect to the beggining.
    #Add the frequency of every character in the word to the list:
    for i in range(len(word)):
        lst.append(freq[word[i]])
    #print lst   
    #Add the whole expected frequency:
    for i in range(len(lst)):
        expect *= lst[i]

    return("%.3f"%(observ / expect))
    
        
#Main:
ch = int(raw_input("If you want to analyze sequence from file press 1,\nIf you want to generate a random sequence press 2: "))

if(ch == 1):
    name = raw_input("Enter the name of the sequence that you want to search for: ")
elif(ch == 2):
    length = raw_input("Enter the length of the randomal DNA: ")
    name = dna(int(length))
    
f = open(name ,'r')
di = bases_freq(f)
print name.strip(".txt")
calculate (di)
window = raw_input("Enter a size for the sliding window: ")
x,y = slidingwindowplot(f,window)
c = open("plot.txt", "w") #Creat a file for what the function returns
c.write("X axis coordinates:\n")
c.write("\n".join(str(element) for element in x)) #Convert 'x' values into string
c.write("\n")
c.write("Y axis coordinates:\n")
c.write("\n".join(str(element) for element in y)) #Convert 'y' values into string
c.close()
print "gc representation:",representation('gc', f)
f.close()

