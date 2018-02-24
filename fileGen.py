#!/usr/bin/python3

import sys
import random

# number of sequences, sequence length
if len(sys.argv) != 3:
  print("""Wrong number of agruments
    \tUseage: fileGen.py <number of sequences> <sequence length>""")
  sys.exit()
  
#open file
fileName = "{}by{}.FASTA".format(sys.argv[1], sys.argv[2])
fileObj = open(fileName, "w+")

#gen strings
seqLabel = ""
seq = ""
nucleotides = ['a' , 'c', 'g', 't']

for i in range(int(sys.argv[1])):
  seqLabel = ">seq {}\n".format(i+1)
  seq = ""
  
  #gen sequence
  for j in range(int(sys.argv[2])):
    seq += random.choice(nucleotides)
  seq += '\n'
  
  #write to file
  fileObj.write(seqLabel)
  fileObj.write(seq)
  
#close file
fileObj.close()
