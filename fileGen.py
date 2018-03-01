+#!/usr/bin/python3
+
+import sys
+import random
+
+def mutate(string):
+  nucleotides = ['a' , 'c', 'g', 't']
+  basesToChange = len(string) // 3
+  for i in range(basesToChange):
+    index = random.randrange(len(string))
+    string = string[:index] + random.choice(nucleotides) + string[index + 1:]
+  return string
+  
+
+# number of sequences, sequence length
+if len(sys.argv) != 4:
+  print("""Wrong number of agruments
+    \tUseage: fileGen.py <number of sequences> <sequence length> <motif to implant>""")
+  sys.exit()
+
+numOfSequences = int(sys.argv[1])
+sequenceLength = int(sys.argv[2])
+motif = sys.argv[3]
+
+#open file
+fileName = "{}by{}_{}.FASTA".format(numOfSequences, sequenceLength, motif)
+fileObj = open(fileName, "w+")
+
+#gen strings
+seqLabel = ""
+seq = ""
+nucleotides = ['a' , 'c', 'g', 't']
+
+for i in range(numOfSequences):
+  implant = mutate(motif)
+  seqLabel = ">seq {} | Implanted String: {} | Motif for Entire Set: {}\n".format(i+1, implant,motif)
+  seq = ""
+  
+  #gen sequence
+  for j in range(sequenceLength):
+    seq += random.choice(nucleotides)
+    
+  #implant the mutated motif into the sequence
+  implantStart = random.randint(0, sequenceLength - len(implant))
+  implantEnd = implantStart + len(implant)
+  seq = seq[:implantStart] + implant + seq[implantEnd:]
+  
+  #add newline to keep correct formating
+  seq += '\n'
+  
+  #write to file
+  fileObj.write(seqLabel)
+  fileObj.write(seq)
+  
+#close file
+fileObj.close()