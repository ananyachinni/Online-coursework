#!/usr/bin/env python3

"""
@author : Ananya Chinni Krishnan
Last edit on: 14-07-2023
"""
import re, sys, os
from collections import Counter
lengths = {}
ids = []
values = []
sequence = []

"""
A small program to try and open the input file, warn the user when the file does not open/exist.
"""
try:
   file = open(sys.argv[1])
except IOError:
       print("The file does not exist!! Please check again.")
	
length = 0
for lin in file:	
        lin = lin.rstrip()
        if lin.startswith('>'):
           		length = length + 1
                
                
print("Number of FASTA sequences in the input file are %i." % length)

file = open(sys.argv[1])

## Reference:
## https://bioinformatics.stackexchange.com/questions/16184/count-of-each-sequence-length-from-a-fasta-file-with-header-using-len-function-o
"""
Small program to parse a multi-FASTA file and write to a dict.
After which, splitting the dict items into lists for further usage.
"""
for line in file:
    line = line.strip()
    if line[0] == ">":
        header = line
        lengths[header] = ""
    else:
        sequ = line
        lengths[header] +=sequ


for lock, key in lengths.items():
    ids.append(lock)
    values.append(len(key))
    sequence.append(key) # Adding only the sequences to a list for ORF analysis


## Reference:
## https://stackoverflow.com/questions/73176254/how-to-remove-the-first-character-of-a-string-thats-in-a-list-python
ids = [go[1:] for go in ids] #List comprehension to remove the first character of every element of the list
##print(ids)


## Reference:
## https://pythonhow.com/how/create-a-dictionary-with-list-comprehension/
"""
Writing the above lists into a dictionary using dict comprehension
"""
final = {ids[i] : values[i] for i in range(len(values))}
##print(final)
matches = {ids[i] : sequence[i] for i in range(len(sequence))} # Writing the dict for only IDs and sequences (ORF analysis)
##print(matches)

## Reference:
## https://datagy.io/python-get-dictionary-key-with-max-value/
"""
Getting the longest and shortest sequence lengths and assign them back to their
identifiers
"""
max_measures = [{valu, key} for key, valu in final.items() if valu == max(final.values())]
print(max_measures,"\n")
min_measures = [{valu,key} for key, valu in final.items() if valu == min(final.values())]
print(min_measures,"\n")

"""
A function to get the value when the key is provided
"""
def find(key, dictionary):
    for k,v in dictionary.items():
        if k == key:
             yield v
## Reference:
## https://www.biostars.org/p/384922/
## https://stackoverflow.com/questions/67178286/you-cant-use-a-list-as-a-default-argument-but-can-you-use-a-tuple-instead
"""
Function for finding possible ORFs in a sequence
"""
def startstop_codon(dna,frame=0):
    dna = dna.upper()
    try:
        if frame not in (0,1,2):
            raise
        else:
            for i in range(frame, len(dna), 3):
                codon1 = dna[i:i+3]
                if codon1 == 'ATG':
                    position1 = i
                    for j in range(position1, len(dna), 3):
                        codon2 = dna[j:j+3]
                        if codon2 in ['TAA', 'TAG', 'TGA']:
                            position2 = j
                            yield (position1, dna[position1:position2+3], position2-position1+3)
                            break
    except:
        pass
        print("Wrong reading frame! Please provide frames with respect to the forward strand.")


def finding_all_ORFs(dictionary, frame):
    """
Function to find all possible ORFs for any frame in the forward strand in every sequence of a multi-FASTA file and
get the length & sequence ID of the longest ORF.
    """
    dump = {}
    for ids, sequences in dictionary.items():
        a = startstop_codon(sequences, frame)
        try:
            if a == None:
                a = 0
                raise
                continue
            else:
                result = list(a)
                dump.update({ids:result})
##            print(dump)
        except:
            pass
    max_orflength = {k: [z[2] for z in sorted(v,key=lambda x:x[-1])] for k,v in dump.items()} # To get the lengths alone from the list of tuples and sorting
                                                                                               # them in the reverse order, writing the values as a dict with the key
##    print(max_orflength)
    length = {k: [z for z in sorted(v,reverse = True)] for k,v in max_orflength.items()} #To arrange the lengths values list in descending order
##    print(length)
    answer = {k: [l for l in v[:1]] for k,v in length.items()} #To get only the first element of the lengths values list
##    print(answer)
    yoohoo = max(answer, key=answer.get) #To get the key of the max length value from the previous dict
##    print(yoohoo)
    print("The length of the longest ORF for this particular frame %i in this multi-FASTA file is %d and its id is %s.\n" %(frame,answer[yoohoo][0], yoohoo))
    ## Alternatively:
    start_pos = list(find(yoohoo,dump))
    max_orf_pos = [(excel[0],excel[2]) for excel in start_pos[0]]
    final = max_orf_pos[3]
    final_pos = final[0] + 1
    print("The length of the longest ORF for this particular frame %i in this multi-FASTA file is %d - position %d and its id is %s.\n" %(frame, answer[yoohoo][0], final_pos, yoohoo))


finding_all_ORFs(matches,2)



def finding_ORFs_specific(dictionary, frame):
    """
Function to obtain the position, length & sequence of the longest ORF in the user-defined sequence. If there are no ORFs present, this function can also
handle the exception and inform the user that there is no value.
    """
    sequences = [j for j in dictionary.values()]
    ids = [a for a in dictionary.keys()]
    if ident in ids:
        seqstr = dictionary.get(ident)
        p = startstop_codon(seqstr,frame)
        try:
           if p == None:
                p = 0
                raise
        except:
           pass
           print("No ORF found!\nIf statement exited.")
        else:
           result = list(p)
##            print(result)
           max_length=max(result, key = lambda item:item[2])
           print("> %s\n%s - length %i\nstarts at position %i.\n" %(ident,max_length[1],max_length[2],max_length[0])) # As max_length is a tuple, second index is added to access the elements
        


ident = input("Please enter the identifier:")
finding_ORFs_specific(matches, 2)


## Reference:
## http://binf.gmu.edu/swang36/NGS/python_genomics_project.html
## https://stackoverflow.com/questions/30692655/counting-number-of-strings-in-a-list-with-python
def find_repeats(zoo, repeat_len):
    """
    Function to find repeats in every sequence of a multi-FASTA file and choose the 5 most common repeats
    """
    repeats_list = []
    common = []
    for i in zoo:
        for j in range(len(i)-repeat_len):
            repeats_list.append(i[j:(j+repeat_len)])
##            print(repeats_list)
    trial = Counter(repeats_list) # Converting the repeats list into a Counter class object
    common = trial.most_common(5) # Finding the 5 most common repeats from the Counter class object
    print(common)

find_repeats(sequence,7)
			
file.close()
