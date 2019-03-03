#This program assembles genomic sequence.

#-----------------------------------------------------------------------
#This function generates a dictionary with the all the reads by indexes.
#-----------------------------------------------------------------------
def readDataFromFile(fileName):
    diction = {} #A dictionary for all the sequences.
    s = file(fileName).readlines()
    #Scans all the lines in the file and insert the sequences into a dictionary.
    for line in s:
        line = line.split() #Splits the line into index and sequence(a list).
        diction[line[0]] = line[1] #Insert the index into the "key's" dictionary and the sequence into the "value's" dictionary.

    return diction
    
#---------------------------------------------------------
#This function calculates the mean length of the sequence.
#---------------------------------------------------------
def meanLength(fileName):
    ammount,length = 0.0,0 #Initializes paramters for the reult.
    diction =  readDataFromFile(fileName) #Gets the dictionary with all the sequences.
    #Scans all the dictionary and returns the length of each value.
    for k,v in diction.items():
        length += len(v)
        ammount += 1

    return length/ammount #Returns the average.

#-----------------------------------------------------------------------
#This function gets 2 sequences and calculates the overlap between them.
#-----------------------------------------------------------------------
def getOverlap(left, right):
    
    overlap = "" #Initializes an empty string.
    #Scans the sequences to find overlaps.
    for i in range(len(left)):
        if (left [i:] == right[:len(left)-i]): #If there is an overlap.
            overlap =  left[i:] #Insert the overlap into the string.
            break

    return overlap

#---------------------------------------------------------------------
#This function gets a dictionary with all the reads
#and return a dictionary with all the overlaps for each sub sequence.
#---------------------------------------------------------------------
def getAllOverlaps(reads):
    
    d = {} #Initializes an empty dictionary.
    lap = 0 #Initializes a parameter for the length of the lap.
    #Scans all the dictionary which the sequences contained. 
    for k1,v1 in reads.items():
        #For each sequence scans all the other sequences.
        for k2,v2 in reads.items():
            if(k1 == k2): #If the index is the same index - skip over it.
                continue
            elif (k1 not in d): #If the index not in the new dictionary:               
                d[k1] = {} #Generate new sub-dictionary.
            lap = getOverlap(v1,v2)#Gets the overlap where v1 is the left and v2 is the right.
            d[k1][k2] = len(lap) #Add to the dictionary in the proper place the length of the overlap.

    return d

#----------------------------------------------------------------------
#This function prints a matrix with all the overlaps for each sequence.
#----------------------------------------------------------------------
def prettyPrint(overlaps):
    
    length = len(overlaps)
    print "  ",
    #Prints the horizontal header of the right reads.
    for i in range(length):
        print "% 3d" %(i+1),
    print
    #Prints the vertical header of left reads and the cells values.
    for i in range(length):
        print " " + str(i+1),
        for j in range(length):
            if(i==j): #If it's the same sequence.
                print "  -",
            else:
                print "% 3d" % overlaps[str(i+1)][str(j+1)], #Prints the number of overlapping bases for a left-right read pair
        print 

#-----------------------------------------------------------------
#This function return the name of the first read of the sequence.
#-----------------------------------------------------------------
def findFirstRead(overlaps):
    
    #Scans for all the reads of the "main" dictionary.
    for key1 in overlaps:
        flag = True #Initializes bolean flag.
        #Scans for all the reads on the "sub" dictionary.
        for key2 in overlaps[key1]:
            if(overlaps[key2][key1] > 3): #If it has significant overlaps: 
                flag = False     
        if(flag):
            return key1
        
#------------------------------------------------------
#This function gets the inner dictionary
#and returns the key associated with the largest value
#------------------------------------------------------
def findKeyForLargestValue(d):
    
    maxi = max(d.values()) #Gets the max value of the "sub" dictionary.
    #Scans the dictionary to find the matching key.
    for k,v in d.items():
        if(v == maxi):
           maxKey = k

    return maxKey

#---------------------------------------------------------------------------------          
#This function is a recursive function
#that returns a list of the reads by the order of their apperance in the sequence.
#---------------------------------------------------------------------------------
def findOrder(name, overlaps):
    
    if(max(overlaps[name].values()) < 3): #Stop condition(The first read).
        return [name]
    else:
        nextName = findKeyForLargestValue(overlaps[name]) #Finds the next read.
        return [name] + findOrder(nextName, overlaps) #Add to the list the next read.

#-----------------------------------------------    
#This function returns the whole genom sequence.
#-----------------------------------------------    
def assembleGenome(readOrder, reads, overlaps):
    
    genomic_seq = "" #Initializes a string for the genom sequence.
    max_overlap=[] #Initializes a list for finding the overlaps.
    #Scans all the sub sequences and removes the overlaps.
    for index in readOrder[:-1]: #Except from the last.
        print readOrder[:-1]
        #Scans all the overlaps in the dictionary
        for over in overlaps[index].values():
            if (over > 2): #If it's a significant overlap: 
                max_overlap.append(over) #Append to the list the overlaps
        overlap = max(max_overlap)
        max_overlap=[] #Reset the list.
        genomic_seq += reads[index][:-overlap] #Add the sub sequence without the overlap
    genomic_seq += reads[readOrder[-1]]#Add the last sub sequence.
    return genomic_seq
    
        
#Main
fileName = raw_input("Enter the name of the file where the sequences are: ")
fragments =  readDataFromFile(fileName)
print "\nThe raw fragments:\n"
for i in range(len(fragments)):
    print str(i+1) + ": " + fragments[str(i+1)]
print "\nThe mean length:\n"    
print "%.2f" % meanLength(fileName)
diction = readDataFromFile(fileName)
d = getAllOverlaps(diction)
print "\nThe overlaps matrix:\n"
prettyPrint(d)
name = findFirstRead(d)
findKeyForLargestValue(d)
order =  findOrder(name, d)
print "\nThe order of the fragments:\n\n" + str(order)
print"\nThe final joined sequence:\n"
print assembleGenome(order,diction,d)





