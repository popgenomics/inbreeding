from Bio.SeqIO import parse
import sys

fastaFileName = sys.argv[1]
#fastaFileName = "/home/croux/Programmes/DILS_sfs/RNAseqFGT_src/test.fas"
missing_data = ['N', 'n', '-', '.', ' ', '_', '?']

def get_SNPs(seq, populations, all_loci, loci_length, missing_data):
	# seq = dictionnary with all the parsed fasta
	# populations = vector of positions to work on
	
	# res: returned object
	# {'Contig69': [57, 108, 332, 374, 701, 1014, 1163, 1166],
	# 'Contig74': [1100, 1220],
	# 'Contig79': [503, 545, 716, 1484, 1511, 1563, 1592, 1647]}

	res = {}	
	
	# loop over loci
	for locus_tmp in range(len(all_loci)):
		locus_i = all_loci[locus_tmp]
		L_i = loci_length[locus_tmp]
	
		# loop over genomic positions
		list_of_SNPs = []
		alleles_1 = []
		alleles_2 = []
		for pos_tmp in range(L_i):
			list_of_alleles = []
			
			# loop over populations
			for pop_tmp in populations:
				# loop over individuals
				for ind_tmp in seq[locus_i][pop_tmp]:
					# loop over alleles
					for allele_tmp in seq[locus_i][pop_tmp][ind_tmp]:

#						seq['Contig69']['Emys_orbicularis']['GA03C']['Allele_1'][60]
						allele = seq[locus_i][pop_tmp][ind_tmp][allele_tmp][pos_tmp]
						if allele not in missing_data:
							list_of_alleles.append(allele)
				
			alleles = list(set(list_of_alleles))
			nAlleles = len(alleles)
			# only works on biallelic positions
			if nAlleles==2:
				list_of_SNPs.append(pos_tmp)
				alleles_1.append(alleles[0])
				alleles_2.append(alleles[1])
		if len(list_of_SNPs)>0:
			res[locus_i]={}
			res[locus_i]['SNPs']=list_of_SNPs
			res[locus_i]['Allele_1']=alleles_1
			res[locus_i]['Allele_2']=alleles_2
	return(res)


all_populations = []
all_individuals = {}
all_loci = []
loci_length = []

infile = parse(fastaFileName, 'fasta')
seq = {}
for seq_tmp in infile:
	locus = seq_tmp.id.split('|')[0]
	population = seq_tmp.id.split('|')[1]
	ind = seq_tmp.id.split('|')[2]
	allele = seq_tmp.id.split('|')[3]
	
	if locus not in seq:
		seq[locus] = {}
		if locus not in all_loci:
			all_loci.append(locus)
			loci_length.append( len(seq_tmp.seq) )
	
	if population not in seq[locus]:
		seq[locus][population] = {}
		if population not in all_populations:
			all_populations.append(population)
			all_individuals[population] = []
	
	if ind not in seq[locus][population]:
		seq[locus][population][ind] = {}
	
	if ind not in all_individuals[population]:
		all_individuals[population].append(ind)
	seq[locus][population][ind][allele] = seq_tmp.seq

SNPs = get_SNPs(seq=seq, populations=all_populations, all_loci=all_loci, loci_length=loci_length, missing_data=missing_data)


#### final output
res = "population\tindividual"
for locus_tmp in SNPs.keys():
	for pos_tmp in SNPs[locus_tmp]['SNPs']:
		res += "\t{locus}.{SNP}".format(locus=locus_tmp, SNP=pos_tmp)
print(res)

for pop_tmp in all_populations:
	for ind_tmp in all_individuals[pop_tmp]:
		res = "{pop}\t{ind}".format(pop=pop_tmp, ind=ind_tmp)
		for locus_tmp in SNPs.keys():
			for pos_tmp in range(len(SNPs[locus_tmp]['SNPs'])):
				pos_i = SNPs[locus_tmp]['SNPs'][pos_tmp]
				allele1 = seq[locus_tmp][pop_tmp][ind_tmp]['Allele_1'][pos_i]
				allele2 = seq[locus_tmp][pop_tmp][ind_tmp]['Allele_2'][pos_i]
#				print("{0}/{1}".format(allele1, allele2))	
				# if no N
				if allele1 in missing_data or allele2 in missing_data:
					genotype = "NA"
				else:
					# if homozygote
					if allele1==allele2:
						if allele1==SNPs[locus_tmp]['Allele_1'][pos_tmp]:
							genotype = 0
						else:
							genotype = 2
					# if heterozygote
					else:
						genotype = 1
				res += "\t{genotype}".format(genotype=genotype)
		print(res)

