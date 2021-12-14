import sys
filename = sys.argv[1]

infile = open(filename, "r")

line = infile.readline()
print(line.strip())

for line in infile:
	res = "\t".join( [ i if i!='2' else '0' for i in line.strip().split('\t') ] )
	print(res)
infile.close()


