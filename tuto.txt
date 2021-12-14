###################
## for hierfstat ##
###################
# produces a table in a dos format.
# dos format: AA=0; Aa=1; aa=2; one line per individual; one column per SNP
fasta2dos time python3 fasta2dos.py mytilus_renamed.fas > mytilus.dos 
fasta2dos time python3 fasta2dos.py test.fas > test.dos 

# the table .dos can be used by the library hierfstat in R.
# keep in mind that the fs.dosage function requires at least 2 populations (not the case for test.fas)
library(hierfstat)
library(data.table)
x=fread('mytilus.dos', h=T)
res=fs.dosage(x[,-c(1,2)], pop=t(t(x[,1])))
plot(res)

# if you want to produce a lighter .dos table by removing all SNPs containing at least one missing data (coded by NA)
python3 removeNA.py mytilus.dos > mytilus_cleaned.dos
python3 removeNA.py test.dos > test_cleaned.dos

##################
## for inbreedR ##
##################
# if you want to convert a .dos table (0/1/2) into a .inbreedR table (0 for AA or aa; 1 for Aa)
python3 hierf2inbreed.py mytilus.dos > mytilus.inbreedR
python3 hierf2inbreed.py test.dos > test.inbreedR

# the table .inbreedR can be used by the library inbreedR in R.
library(inbreedR)
library(data.table)
x=fread('mytilus.inbreedR', h=T)
g2_snps(x[,-c(1,2)], nperm = 100, nboot = 100, CI = 0.95)

##################################################################
## producing all figures + tables from .dos and .inbreedR files ##
##################################################################
Rscript inbred_stats.R input_hierfstat=mytilus.dos input_inbreedR=mytilus.inbreedR

