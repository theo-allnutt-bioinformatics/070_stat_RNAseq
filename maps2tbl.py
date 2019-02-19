#!/usr/bin/python

import sys
import os
import re
import glob

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))
				  


folder = sys.argv[1] #working folder
outfile=sys.argv[2] #output file

g=open(outfile,'w')

filelist=glob.glob(folder)

filelist.sort(key=tokenize)

print filelist

data={}
allspecies=[]
filenames=[]
n1=0
print 'Collecting genes'
for i in filelist:
	print i
	file1 = open(i,'r')
	filename = i.split("/")[-1].split(".")[0]
	filenames.append(filename)
	data[filename]={}
	
	for k in file1:
		if k[0]<>"#":
			
			#print k
			
			mrna = k.split("\t")[0]
			
			
			ambig=int(k.split("\t")[-1].rstrip("\n"))
			unambig=int(k.split("\t")[-2])
			
			
			freq= ambig+unambig

			if mrna not in allspecies:
				allspecies.append(mrna)
				n1=n1+1
				print i, mrna.split(" ")[0],n1
			data[filename][mrna]=freq
		
	file1.close()
			
allspecies.sort()

#print "\n".join(str(x) for x in allspecies)	

g.write("\t"+"\t".join(str(x) for x in filenames)+"\n")

print 'writing file'
for i in allspecies:
	out = i+"\t"
	for name in filenames:
		
		if i in data[name].keys():
			out = out + str(data[name][i])+"\t"
			
		else:
			out = out + "0" +"\t"
	
	g.write(out+"\n")

