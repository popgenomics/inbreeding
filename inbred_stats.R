library(data.table)
library(hierfstat)
library(inbreedR)

fileName_hierfstat = NULL
fileName_inbreedR = NULL
for(tmp in commandArgs()){
	tmp = strsplit(tmp, "=")
	if(tmp[[1]][1] == "input_hierfstat"){
		fileName_hierfstat=tmp[[1]][2]
	}
	if(tmp[[1]][1] == "input_inbreedR"){
		fileName_inbreedR=tmp[[1]][2]
	}
}

if(is.null(fileName_hierfstat)){
	print("no input file indicated by user (ex: ./inbred_stats.R input_hierfstat=input.dos input_inbreedR=input.inbreedR)")
	quit(status=1)
}

############################
######    hierfstat    #####
############################
# read the input file
x=fread(fileName_hierfstat, h=T)

# measures statistics
res=fs.dosage(x[,-c(1,2)], pop=t(t(x[,1])))
write.table(x = round(res$Fs,5), file = "results_hierfstat.txt", col.names=T, row.names=T, sep="\t", quote=F)

# Fst : pairwise Fst
Fst_table=res$Fst2x2

popA = c()
popB = c()
Fst = c()
for(popA_tmp in 1:ncol(Fst_table)){
	for(popB_tmp in 1:nrow(Fst_table)){
		popA = c(popA, colnames(Fst_table)[popA_tmp])
		popB = c(popB, rownames(Fst_table)[popB_tmp])
		Fst = c(Fst, Fst_table[popB_tmp, popA_tmp])
	}
}

Fst_table = tibble(popA=popA, popB=popB, Fst=round(Fst, 5))
write.table(x=Fst_table, file="results_Fst.txt", col.names=T, row.names=F, sep="\t", quote=F)
Fst_plot = Fst_table %>% dplyr::mutate(Fst=ifelse(Fst<0, 0, Fst)) %>% ggplot(aes(x=popA, y=popB, fill=Fst)) + geom_tile() + scale_fill_viridis_c(limits=c(0, 1))
ggsave(filename = "plot_Fst.pdf", plot = Fst_plot, bg="white")

# Fi
## Fi list of individual inbreeding coefficients, estimated with the reference being the population to which the individual belongs.
Fi_table=res$Fi
pop = c()
Fi = c()

for(pop_name in names(Fi_table)){
	for(ind in 1:length(Fi_table[[pop_name]])){
		pop = c(pop, pop_name)
		Fi = c(Fi, Fi_table[[pop_name]][ind])
	}
}
	
Fi_table = tibble(pop=pop, Fi=round(Fi, 5))
write.table(x=Fi_table, file="results_Fi.txt", col.names=T, row.names=F, sep="\t", quote=F)

Fi_table_2 = Fi_table %>% group_by(pop)  %>% summarise(Fi_median=round(median(Fi),5), Fi_mean=round(mean(Fi),5), Fi_sd=round(sd(Fi),5), Fi_min=round(min(Fi),5), Fi_max=round(max(Fi),5))
write.table(x=Fi_table_2, file="results_Fi_summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

if(min(Fi_table$Fi)<0){
	ymin = min(Fi_table$Fi)
}else{
	ymin = 0
}
if(max(Fi_table$Fi)<0){
	ymax = max(Fi_table$Fi)
}else{
	ymax = 1
}

Fi_plot = Fi_table %>% ggplot(aes(x=pop, y=Fi)) + geom_boxplot() + theme_bw() + ggtitle("Ind. Inb. coeff.") + ylim(ymin, ymax)
ggsave(filename = "plot_Fi.pdf", plot = Fi_plot, bg="white")

###########################
######    inbreedR    #####
###########################
filterSNP = function(x){
	nNA = sum(is.na(x))
	if(nNA>0){
		return(FALSE)
	}else{
		if(sum(x)==0){
			return(FALSE)
		}else{
			return(TRUE)
		}
	}
}

x_inbreedR=fread(fileName_inbreedR, h=T)
populations = sort(unique(as.character(t(t(x_inbreedR[,1])))))
pop = c()
g2 = c()
g2_se = c()
nInd = c()
nLocus = c()

pdf("plot_G2.pdf", bg="white")
for(pop_tmp in populations){
	print(pop_tmp)
	y = as.matrix(x_inbreedR[which(x_inbreedR[,1]==pop_tmp), -c(1,2)])
	y2 = apply(y, MARGIN=2, FUN="filterSNP")
	y3 = y[,y2]
	res_g2 = g2_snps(y3, nperm = 1000, nboot = 1000, CI = 0.95, parallel=T, ncores=4)
	
	pop = c(pop, pop_tmp)
	g2 = c(g2, res_g2$g2)
	g2_se = c(g2_se, res_g2$g2_se)
	nInd = c(nInd, res_g2$nobs)
	nLocus = c(nLocus, res_g2$nloc)
	
	plot(res_g2, main=pop_tmp)
}
dev.off()

res_tot = tibble(populations=pop, nIndividuals=nInd, nSNPs=nLocus, g2=round(g2, 5), g2_se=round(g2_se, 5))

write.table(x=res_tot, file="results_g2.txt", col.names=T, row.names=F, sep="\t", quote=F)

