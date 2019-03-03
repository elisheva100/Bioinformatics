from string import *

def frequent (seq, k):
    kmers = {} #Create empty dictionary.
    #All the optionals k-mers.
    for i in range(len(seq)-k+1):  
        kmer = seq[i:i+k]
        if(kmer in kmers): #The k-mer allready exist:
            kmers[kmer] += 1 #Add to it's count the number of its appearance.
        else:
            kmers[kmer] = 1 #Add new element to the dictionary.

    freq = [] #Create empty list.
    max_value = max(kmers.values())#The value of the highest frequency.
    #Split the dictionary to keys and valuse.
    for k,v in kmers.items(): 
        if(v == max_value): #If the element has the maximum frequency:
            freq.append(k) #Add it to the list
    print " ".join(str(element) for element in freq) #Print the list as a continuous string
    #print max_value
#main
    
name = raw_input("Enter the name of the sequence that you want to search for: ")
number = input("Enter the k-mer's size: ")
s = file(name,"r").read()
#if(s.find(">") > - 1): # If the sequence has an information:
    #s = s[s.find("\n"):] #Delete the information.
    #s = s.replace("\n","") #Delete the end of lines.
#while (number != 0):
frequent(s,number)
    #number = input("Enter the k-mer's size: ")

#All the code lines that marked as comments are for section 3.
#(see at the output)
    
            
