#######  check matches in WES data set where MP used

################################    rare variant in WES only

## step 1: read WES rare variant data set and do missingness filtering 

m.rare0 <- readRDS('/Volumes/Kexin/ASP/previous_doc/process_data_m_rare0.rds')
sum_missing <- function(x){
  sum(x == './.')
}

m.rare0_missingN <- apply(m.rare0[,c(11:52)],1,sum_missing)
m.rare0_all <- m.rare0
m.rare0 <- m.rare0[m.rare0_missingN<=4,]
dup_indi <- c(1:nrow(m.rare0))[duplicated(m.rare0[,c('CHROM','POS')])]
m.rare0_unique <- m.rare0[-dup_indi,]

##  step 2: read WGS data set and find shared 

shared <- data.frame()
for (chr in 1:22) {
  chr_data <- readRDS(paste0("/Volumes/Kexin/ASP/jointly_called_WGS/annotation/with_gene_symble/chr",chr,"_variants_ensembl.rds"))
  bi_snp <- chr_data$bi_snp
  multi_snp <- chr_data$multi_snp
  indel <- chr_data$indel
  pos <- unique(m.rare0$POS[m.rare0$CHROM==chr])
  
  bi_snp_share <- bi_snp[bi_snp$POS %in% pos,]
  multi_snp_shared <- multi_snp[multi_snp$POS %in% pos,]
  indel_shared <- indel[indel$POS %in% pos,]
  shared <- rbind(shared,bi_snp_share,multi_snp_shared,indel_shared)
}


saveRDS(shared,'/Volumes/Kexin/ASP/WES_WGS_compare/data/wgs_wes_rare_var_overlap_unfiltered.rds')


################################    all variant in WES 
m0=read.csv("/Volumes/Kexin/ASP/previous_doc/SisterPairs_recalibratedSNPs_INDELs_shellScript_lifted-hg19_noMeta.vcf_txt2csv.csv_vcfAnnotated.csv_colSelect_dep10_PASS_genotypeON_UKB.csv", stringsAsFactors=FALSE)


####### modify maf 

m0$X1000G_ALL[m0$X1000G_ALL=='.'] = NA
m0$X1000G_ALL[m0$X1000G_ALL==''] = NA
m0$ESP6500si_ALL[m0$ESP6500si_ALL=='.'] = NA
m0$ESP6500si_ALL[m0$ESP6500si_ALL==""] = NA

maf.1kg = as.numeric(m0$X1000G_ALL)
maf.esp = as.numeric(m0$ESP6500si_ALL)
maf.ukbb = as.numeric(m0$UKBB)

maf.1kg = pmin( maf.1kg, 1-maf.1kg)
maf.esp = pmin( maf.esp, 1-maf.esp)
maf.ukbb = pmin( maf.ukbb, 1-maf.ukbb)

max_maf = pmax(maf.1kg, maf.esp, maf.ukbb, na.rm=T)

m0$max_maf <- max_maf

####### find share

wgs_shared <- data.frame()
wgs_shared_unique <- data.frame()
wes_shared <- data.frame()
nshared <- c()
wes_n <- c()
wgs_n <- c()
for (chr in 1:22) {
  chr_data <- readRDS(paste0("/Volumes/Kexin/ASP/jointly_called_WGS/annotation/with_gene_symble/chr",chr,"_variants_ensembl.rds"))
  bi_snp <- chr_data$bi_snp
  multi_snp <- chr_data$multi_snp
  indel <- chr_data$indel
  pos <- unique(m0$POS[m0$CHROM==paste0('chr',chr)])
  wes_n <- c(wes_n,length(unique(pos)))
  wgs_n <- c(wgs_n, length(unique(c(bi_snp$POS,multi_snp$POS,indel$POS))))
  
  bi_snp_share <- bi_snp[bi_snp$POS %in% pos,]
  multi_snp_shared <- multi_snp[multi_snp$POS %in% pos,]
  indel_shared <- indel[indel$POS %in% pos,]
  nshared <- c(nshared,length(unique(c(bi_snp_share$POS,multi_snp_shared$POS,indel_shared$POS))))
  wgs_shared <- rbind(wgs_shared,bi_snp_share,multi_snp_shared,indel_shared)
  wgs_shared_unique <- rbind(wgs_shared_unique,bi_snp_share,multi_snp_shared[-duplicated(multi_snp_shared$POS),],indel_shared[-duplicated(indel_shared$POS),])
  wes_shared <- rbind(wes_shared,m0[m0$CHROM==paste0('chr',chr) & m0$POS %in% unique(c(bi_snp_share$POS,multi_snp_shared$POS,indel_shared$POS)),])
}
sum(wes_n)
sum(wgs_n)
sum(nshared)

sum(wgs_n)-sum(nshared)
sum(wes_n)-sum(nshared)

colnames(wes_shared)

wgs_shared$CHROM <- paste0('chr',wgs_shared$CHROM)
# sum((wes_shared$max_maf==0 | is.na(wes_shared$max_maf)) & wgs_shared$maf ==0)
WGS_WES_shared <- merge(x = wes_shared[,c('CHROM','POS','Func.refgene','Gene.refgene','max_maf')],
                             y = wgs_shared[,c('CHROM','POS','maf','gene_type','gene_symbol','effect_priority','typeseq_priority')],
                             by = c('CHROM','POS'), all.x= T)

sum((WGS_WES_shared$max_maf == 0 | is.na(WGS_WES_shared$max_maf)) & WGS_WES_shared$maf==0)
sum((WGS_WES_shared$max_maf > 0.05 & WGS_WES_shared$max_maf<=0.5 & WGS_WES_shared$maf>0.05 & WGS_WES_shared$maf<=0.5),na.rm=T)

plot(WGS_WES_shared$max_maf ~ WGS_WES_shared$maf,xlab='WGS_maf', ylab='WES_maf',main='Paired MAF in shared variants')
rare_indi <- which(WGS_WES_shared$max_maf>0 & WGS_WES_shared$max_maf<0.005)
plot(WGS_WES_shared[rare_indi,'max_maf'] ~ WGS_WES_shared[rare_indi,'maf'],
     xlim=c(0,0.005),ylim=c(0,0.005),
     xlab='WGS_maf', ylab='WES_maf',main='Paired MAF in rare variants MP used')



########   function agreement

table(WGS_WES_shared$Func.refgene)

index = c(1:nrow(WGS_WES_shared))[sapply(WGS_WES_shared$Func.refgene, grepl,pattern='upstream',fixed = TRUE) | 
                                    sapply(WGS_WES_shared$Func.refgene, grepl,pattern='downstream',fixed = TRUE)]

index = c(1:nrow(WGS_WES_shared))[sapply(WGS_WES_shared$Func.refgene, grepl,pattern='intronic',fixed = TRUE)]

wgs_type = (WGS_WES_shared$typeseq_priority[index])

sum(grepl('exonic',wgs_type,fixed=T))
sum(grepl('intergenic',wgs_type,fixed=T))
sum(grepl('intronic',wgs_type,fixed=T))
sum(grepl('downstream',wgs_type,fixed=T))+sum(grepl('upstream',wgs_type,fixed=T))
sum(grepl('UTR5',wgs_type,fixed=T))+sum(grepl('UTR3',wgs_type,fixed=T))
sum(grepl('slicing',wgs_type,fixed=T))



