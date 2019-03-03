#Elisheva Heilbrun 208833012
from string import *

print ">gi number  Length  motive"
print "--------------------------"
amount = 0# count the sequences
s = file("seq.txt", "r").read()
parts = s.split(2*"\n")#Separates between the sequences.
for seq  in parts:
    info = seq[0:seq.find("|",seq.find("|")+1)]#till it's get the second "|".
    length = len(seq[seq.find("\n"):])- seq.count("\n")
    # The length of the sequence without information and without end of lines.
    subseq = seq.replace("\n","")#delete all the end of lines ("\n")
    last = subseq[-51:]#the last 50 characters
    motive = last.find("GGGCACCCG")#Checks if the asked motive exist
    if(motive > 0):#Found
        ans = '+'
    else:
        ans = '-'
    amount += seq.count(">")
    print info + "|" + "  " + str(length) + "      " + ans
print "Total of " + str(amount)+" sequences"

