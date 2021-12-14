import sys

fileName = sys.argv[1]
infile = open(fileName, 'r')

toRemove = []

line = infile.readline()
line = line.strip().split('\t')
nElements = len(line)

toKeep = [ i for i in range(nElements) ]

for line in infile:
	line = line.strip().split('\t')
	toKeep = [ i for i in toKeep if line[i]!='NA' ]
infile.close()

infile = open(fileName, 'r')
for line in infile:
	line = line.strip().split('\t')
	res = "\t".join([ line[i] for i in toKeep ])
	print(res)
infile.close()
