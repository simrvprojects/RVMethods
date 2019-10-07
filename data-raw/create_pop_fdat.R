library(SimRVSequences)
#create mutation and genome set from slim output
sout = read_slim(file_path = "C:/Users/cnieuwoudt/Google Drive/sfuthesis-2.2.1/SlimData/slim_chr8_50Kfound.txt",
                 keep_maf = 0.001,
                 recomb_map = create_slimMap(hg_exons[hg_exons$chrom == 8, ]),
                 pathway_df = hg_apopPath[hg_apopPath$chrom == 8, ])

# sout contains 2 items
# the first is the haplotypes
# the second is the mutations dataset, which describes the haplotypes.
str(sout)

#tablulate the number of segments in each gene of interest in chromosome 8
table(hg_apopPath$gene[hg_apopPath$chrom == 8])
table(sout$Mutations$chrom[sout$Mutations$pathwaySNV])

#create is_CRV variable
table(sout[[2]]$afreq[sout[[2]]$pathwaySNV])

#the number of haplotypes that carry an SNV with the specified allele frequency
afreq_vals <- as.numeric(names(table(sout[[2]]$afreq[sout[[2]]$pathwaySNV])[1:15]))
#number of SNVs with specified allele frequency
numSnv <- as.numeric(table(sout[[2]]$afreq[sout[[2]]$pathwaySNV])[1:15])

data.frame(afreq = afreq_vals,
           totSNV = numSnv,
           tot_afreq = afreq_vals*numSnv,
           tot_haps = afreq_vals*numSnv*100000)

#option 2: use all SNVs in path with afreq in (0.00011,0.00010)
#cumulative allele frequency
sum(numSnv[c(11, 10)]*afreq_vals[c(11, 10)])
#the number of distinct haplotypes in the pool
sum(numSnv[c(11, 10)]*afreq_vals[c(11, 10)]*100000)
#the number of distinct SNVs
sum(numSnv[c(11, 10)])

sout$Mutations$is_CRV = FALSE
sout$Mutations$is_CRV[which(sout$Mutations$pathwaySNV & sout$Mutations$afreq == 0.00011)] = TRUE
sout$Mutations$is_CRV[which(sout$Mutations$pathwaySNV & sout$Mutations$afreq == 0.00010)] = TRUE
sum(sout$Mutations$afreq[sout$Mutations$is_CRV])

pop_fdat = sout
keep_cols = pop_fdat$Mutations$colID[pop_fdat$Mutations$is_CRV]
keep_rows = which(rowSums(pop_fdat$Haplotypes[, pop_fdat$Mutations$colID[pop_fdat$Mutations$is_CRV]])>0)
pop_fdat$Haplotypes[keep_rows, keep_cols]
rowSums(pop_fdat$Haplotypes[keep_rows, keep_cols])


class(pop_fdat)

save(pop_fdat, file="data/pop_fdat.rda", compress='xz')
