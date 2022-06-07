##  perform LZ test on WGS and WES shared data 
############  read WGS data set
############

library(logistf)

EM <- function(para, cn, cs) {
  factor <- rep(NA, ncol(para))
  for(i in 1:ncol(para)) {#i <- 1 iterate over the positions
    u <- para[1,i] #number of unknown configuration (Double hets in IBD 1)
    if(u==0) {
      factor[i] <- NA
      next
    }
    
    #initialization
    kn <- para[2,i] #known non-shared variants (On ibd 0 or single variants on ibd 1)
    ks <- para[3,i] #known shared variants  (On ibd 2 or more than two variants on ibd 1)
    cn <- cn #total number of non-shared chromosomes
    cs <- cs # total number of shared chromosomes
    
    pn.init <- kn/(cn-u*2) #probability of rare variant on non-shared chromosome
    pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
    ps.init <- ks/(cs-u) #probability of rare variant on shared chromosome
    ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
    delta <- Inf
    iter <- 1
    
    while(delta > 10^-6) {
      #E step
      #us <- u*ps.cur/(pn.cur+ps.cur)
      #un <- u*pn.cur/(pn.cur+ps.cur)
      us <- u* ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
      un <- u* (1-ps.cur)*pn.cur^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
      #M-step
      pn.new <- (kn + 2*un)/cn
      ps.new <- (ks+us)/cs
      #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
      
      #check convergence
      delta <- max(abs(pn.cur - pn.new), abs(ps.cur - ps.new))
      pn.cur <- pn.new
      ps.cur <- ps.new
      
      #print(c(pn.cur, ps.cur, iter))
      #iter <- iter + 1
    }
    #c(pn.init, ps.init)
    factor[i] <- result <- c(ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2))
  }
  #output a correction factor for each position
  factor
}

wgs_shared_data <- readRDS("/Volumes/Kexin/ASP/previous_doc/wgs_shared_filtered.rds")

############   filter on quality and missingness
############

wgs_shared_data$het <- 1*(wgs_shared_data$HET==wgs_shared_data$NCALLED)
wgs_shared_data$homo <- 1*(wgs_shared_data$HOM.REF==wgs_shared_data$NCALLED | wgs_shared_data$HOM.VAR==wgs_shared_data$NCALLED)
wgs_shared_data <- wgs_shared_data[wgs_shared_data$het!=1,]
wgs_shared_data <- wgs_shared_data[wgs_shared_data$homo!=1,]
wgs_shared_data <- wgs_shared_data[wgs_shared_data$missing10!=1,]
wgs_shared_data <- wgs_shared_data[wgs_shared_data$QD_f!=1,]

wgs_shared_data$minor_allele_count <- NA
wgs_shared_data$minor_allele_count[wgs_shared_data$ALT==wgs_shared_data$minor_allele] <- wgs_shared_data$HET[wgs_shared_data$ALT==wgs_shared_data$minor_allele]+2*wgs_shared_data$'HOM.VAR'[wgs_shared_data$ALT==wgs_shared_data$minor_allele]  
wgs_shared_data$minor_allele_count[wgs_shared_data$REF==wgs_shared_data$minor_allele] <- wgs_shared_data$HET[wgs_shared_data$REF==wgs_shared_data$minor_allele]+2*wgs_shared_data$'HOM.REF'[wgs_shared_data$REF==wgs_shared_data$minor_allele]  

############   prs 
############
prs_all <-  readRDS("/Volumes/Kexin/ASP/PRS/Mavaddat/data/WGS_prs.rds")

chr = 22
wgs_shared_data_chr = wgs_shared_data[wgs_shared_data$CHROM==chr,]
############   IBD info
############  
ibd <-  read.table("/Volumes/Kexin/ASP/KING/IBD/no_maf_all_chr.segments",header = T)
ibd_f2_sis <- ibd [ibd$Chr==chr, ]
rm(ibd)
if (chr>=10){
  d<- 4
} else{
  d<-3
}
ibd_f2_sis$Startbp <- substring(ibd_f2_sis$StartSNP,d)
ibd_f2_sis$Startbp <- substr(ibd_f2_sis$Startbp,1,(nchar(ibd_f2_sis$Startbp)-4))

ibd_f2_sis$Stopbp <- substring(ibd_f2_sis$StopSNP,d)
ibd_f2_sis$Stopbp <- substr(ibd_f2_sis$Stopbp,1,(nchar(ibd_f2_sis$Stopbp)-4))

############   look at gene level:
############  
genelist <- unique(wgs_shared_data_chr$gene_symbol)
genelist <- unique(unlist(strsplit(genelist,';')))
genelist <- genelist[genelist!='NA']

out <- data.frame()
# sink(paste0("/home/bulllab/kluo/ASP/LZtest/results/LZ_WGS_chr",chr,".txt"))

for (gene in genelist){
  region <- wgs_shared_data_chr[which(!is.na(unlist(lapply(strsplit(wgs_shared_data_chr$gene_symbol,';'),match,x=gene)))),]
  ### mark rare variants by maf < 0.001
  rare_var_indi <- which(region$maf<0.001 & region$maf>0 & region$minor_allele_count<7)
  rare_var <- region[rare_var_indi,]
  if (nrow(region)==0 | (length(rare_var_indi)==0)) {
    out<- rbind(out,c(rep(NA,21),gene,nrow(region))) 
  } else {
    
    ##  family info
    f1 <- region[,colnames(region) %in% c('RD_1472_resub_Sep11_19.GT','RD_1467.GT','RD_1470.GT')]
    colnames(f1) <- c('F1_1','F1_2','F1_3')
    
    f2 <- region[,colnames(region) %in% c('RD_48033.GT','RD_13407.GT')]
    colnames(f2) <- c('F2_1','F2_2')
    
    f3 <- region[,colnames(region) %in% c('RD_1418.GT','RD_1419.GT')]
    colnames(f3) <- c('F3_1','F3_2')
    
    f5 <- region[,colnames(region) %in% c('RD97_589.GT','RD97_354.GT')]
    colnames(f5) <- c('F5_1','F5_2')
    
    f6 <- region[,colnames(region) %in% c('RD_48185.GT','RD_2256.GT')]
    colnames(f6) <- c('F6_1','F6_2')
    
    f11 <- region[,colnames(region) %in% c('RD_17861.GT','RD_21170.GT','RD_48042.GT')]
    colnames(f11) <- c('F11_1','F11_2','F11_3')
    
    f12 <- region[,colnames(region) %in% c('RD_19007.GT','RD_22873.GT')]
    colnames(f12) <- c('F12_1','F12_2')
    
    f15 <- region[,colnames(region) %in% c('RD_18010.GT','RD_19183.GT')]
    colnames(f15) <- c('F15_1','F15_2')
    
    f16 <- region[,colnames(region) %in% c('RD_48043.GT','RD_13557.GT')]
    colnames(f16) <- c('F16_1','F16_2')
    
    f17 <- region[,colnames(region) %in% c('RD_41419.GT','RD_14504.GT')]
    colnames(f17) <- c('F17_1','F17_2')
    
    f19 <- region[,colnames(region) %in% c('RD97_60.GT','RD97_56.GT')]
    colnames(f19) <- c('F19_1','F19_2')
    
    f21 <- region[,colnames(region) %in% c('RD_1199.GT','RD_1198.GT')]
    colnames(f21) <- c('F21_1','F21_2')
    
    f22 <- region[,colnames(region) %in% c('RD_31561.GT','RD97_434.GT')]
    colnames(f22) <- c('F22_1','F22_2')
    
    f23 <- region[,colnames(region) %in% c('RD_4056.GT','RD_14438.GT')]
    colnames(f23) <- c('F23_1','F23_2')
    
    f24 <- region[,colnames(region) %in% c('RD_14145.GT','RD_24058.GT','RD_24059.GT')]
    colnames(f24) <- c('F24_1','F24_2','F24_3')
    
    family_data <- cbind(f1,f2,f3,f5,f6,f11,f12,f15,f16,f17,f19,f21,f22,f23,f24)
    
    ##  IBD for each sister pair (exclude f1, f11 and f24 for now)
    region_start=region$POS[1]
    region_end=region$POS[nrow(region)]
    
    family_ibd <- data.frame()
    
    for (family in c(2,3,5,6,12,15,16,17,19,21,22,23)) {
      f_ibd <- ibd_f2_sis[ibd_f2_sis$ID1==paste0('F',family,'_F',family,'_3') & ibd_f2_sis$ID2==paste0('F',family,'_F',family,'_4'),]
      ibd_intersect <- c()
      for (i in 1:nrow(f_ibd)){
        ibd_intersect <- c(ibd_intersect,(as.numeric(f_ibd$Startbp[i]) < region_start) & (as.numeric(f_ibd$Stopbp[i])>region_end))   
      }
      ibd_f <- f_ibd$IBDType[ibd_intersect]
      if (length(ibd_f)==0) {
        ibd_f = "IBD0"
      }
      family_ibd <- rbind(family_ibd,c(family,ibd_f ))
    }
    
    ## IBD for sister trio: f1, f11 and f24
    
    for (family in c(1,11,24)) {
      f_ibd <- ibd_f2_sis[ibd_f2_sis$ID1==paste0('F',family,'_F',family,'_3') & ibd_f2_sis$ID2==paste0('F',family,'_F',family,'_4'),]
      ibd_intersect <- c()
      for (i in 1:nrow(f_ibd)){
        ibd_intersect <- c(ibd_intersect,(as.numeric(f_ibd$Startbp[i]) < region_start) & (as.numeric(f_ibd$Stopbp[i])>region_end))     
      }
      ibd_f <- f_ibd$IBDType[ibd_intersect]
      if (length(ibd_f)==0) {
        ibd_f = "IBD0"
      }
      family_ibd <- rbind(family_ibd,c(family,ibd_f ))
      
      f_ibd <- ibd_f2_sis[ibd_f2_sis$ID1==paste0('F',family,'_F',family,'_3') & ibd_f2_sis$ID2==paste0('F',family,'_F',family,'_5'),]
      ibd_intersect <- c()
      for (i in 1:nrow(f_ibd)){
        ibd_intersect <- c(ibd_intersect,(as.numeric(f_ibd$Startbp[i]) < region_start) & (as.numeric(f_ibd$Stopbp[i])>region_end))   
      }
      ibd_f <- f_ibd$IBDType[ibd_intersect]
      if (length(ibd_f)==0) {
        ibd_f = "IBD0"
      }
      family_ibd <- rbind(family_ibd,c(family,ibd_f ))
      
      f_ibd <- ibd_f2_sis[ibd_f2_sis$ID1==paste0('F',family,'_F',family,'_4') & ibd_f2_sis$ID2==paste0('F',family,'_F',family,'_5'),]
      ibd_intersect <- c()
      for (i in 1:nrow(f_ibd)){
        ibd_intersect <- c(ibd_intersect,(as.numeric(f_ibd$Startbp[i]) < region_start) & (as.numeric(f_ibd$Stopbp[i])>region_end))   
      }
      ibd_f <- f_ibd$IBDType[ibd_intersect]
      if (length(ibd_f)==0) {
        ibd_f = "IBD0"
      }
      family_ibd <- rbind(family_ibd,c(family,ibd_f ))
    }
    
    colnames(family_ibd) <- c('family',"IBD")
    # family_ibd$IBD[!startsWith(family_ibd$IBD,'IBD') ] <- 'IBD0'
    
    ## add sis id
    family_ibd$sis1 <- c("RD_48033","RD_1418","RD97_589","RD_48185","RD_19007","RD_18010",
                         "RD_48043","RD_41419","RD97_60","RD_1199","RD_31561","RD_4056",
                         "RD_1472","RD_1472","RD_1467","RD_17861","RD_17861","RD_21170",
                         "RD_14145","RD_14145","RD_24058")
    family_ibd$sis2 <- c("RD_13407","RD_1419","RD97_354","RD_2256","RD_22873","RD_19183",
                         "RD_13557","RD_14504","RD97_56","RD_1198","RD97_434","RD_14438",
                         "RD_1467","RD_1470","RD_1470","RD_21170","RD_48042","RD_48042",
                         "RD_24058","RD_24059","RD_24059")
    
    if (sum(is.na(family_ibd$IBD))>0) {
      family_ibd$IBD[which(is.na(family_ibd$IBD))] <- 'IBD0'
    }
    ###  sharing status for each pair (21)
    ###  only look at rare variants
    family_ibd$IBD[family_ibd$IBD=='IBD0'] <- 0
    family_ibd$IBD[family_ibd$IBD=='IBD1'] <- 1
    family_ibd$IBD[family_ibd$IBD=='IBD2'] <- 2
    family_ibd$IBD <- as.numeric(family_ibd$IBD)
    
    ### create parameters for EM (impute ambiguous configuration)
    n_snp <- nrow(rare_var)
    n_sample <- nrow(family_ibd)
    para <- array(NA, c(3, n_snp), list(c("u", "kn", "ks"), rare_var$POS))
    amb_sibpair <- array(FALSE, c(n_sample,n_snp))
    for(j in 1:n_snp) {
      u <- kn <- ks <- 0
      for(i in 1:n_sample) {
        index1 <- which(colnames(rare_var) %in% c('CHROM','POS','major_allele','minor_allele'))
        geno <- rare_var[j,c(index1,which(startsWith(colnames(rare_var),family_ibd$sis1[i])),which(startsWith(colnames(rare_var),family_ibd$sis2[i])))]   
        count_minor <- c(r1=sum(unlist(strsplit(as.character(geno[5]),'/'))==geno$minor_allele),r2=sum(unlist(strsplit(as.character(geno[6]),'/'))==geno$minor_allele))
        
        if (family_ibd$IBD[i]==0) {
          kn <- kn + sum(count_minor)
        }
        
        if (family_ibd$IBD[i]==2) {
          ks <- ks + count_minor[1]
        }
        
        if (family_ibd$IBD[i]==1) {
          sib1 = count_minor[1]
          sib2 = count_minor[2]
          if (sib1==1 & sib2==1) {
            u <- u + 1
            amb_sibpair[i,j] <- T
          } else {
            if(sum(count_minor)==1) kn <- kn + 1
            if(sum(count_minor)==3) {
              kn <- kn + 1
              ks <- ks + 1
            }
            if(sum(count_minor)==4) {
              kn <- kn + 2
              ks <- ks + 1
            }
          }
        } 
      }
      para[,j] <- c(u, kn, ks)
    }
    
    cn <- 4*sum(family_ibd$IBD==0) + 2*sum(family_ibd$IBD==1)  #no. of non-shared chromosomes
    cs <- sum(family_ibd$IBD==1) + 2*sum(family_ibd$IBD==2) #no. of shared chromosome
    
    prob_shared <- EM(para=para, cn=cn, cs=cs) #the probability of being a shared variant
    
    family_region <- data.frame()
    family_data <- array(NA, c(21,12), list(NULL, c('family','sis1','sis2','IBD','n_case','n_control',"s", "ns1", "ns2", "ambiguous",'xns','xs')))
    for (j in 1:21) {
      # first extract variants genotypes
      family_geno<- data.frame()
      for (fi in 1:nrow(rare_var)){
        index1 <- which(colnames(rare_var) %in% c('CHROM','POS','major_allele','minor_allele'))
        geno <- rare_var[fi,c(index1,which(startsWith(colnames(rare_var),family_ibd$sis1[j])),which(startsWith(colnames(rare_var),family_ibd$sis2[j])))]   
        count_minor <- c(r1=sum(unlist(strsplit(as.character(geno[5]),'/'))==geno$minor_allele),r2=sum(unlist(strsplit(as.character(geno[6]),'/'))==geno$minor_allele))
        family_geno <- rbind(family_geno,c(geno,count_minor))
      }
      sis1 = data.frame(t(c(family_ibd$family[j],family_ibd$sis1[j],family_geno$r1)))
      sis2 = data.frame(t(c(family_ibd$family[j],family_ibd$sis2[j],family_geno$r2)))
      family_region <- rbind(family_region,sis1,sis2)
      
      hs_i <-  (sum(family_ibd$IBD[j]==1) + 2*sum(family_ibd$IBD[j]==2))/2
      hns_i <- (2*sum(family_ibd$IBD[j]==1) + 4*sum(family_ibd$IBD[j]==0))/2
      geno1 <- family_geno$r1
      geno2 <- family_geno$r2
      
      if (family_ibd$IBD[j]==0) {
        s1 <- sum(geno1) 
        s2 <- sum(geno2)
        s <- 0
        ambig <- 0
        xns <- sum(s1>0,s2>0)
        xs <- 0
      }
      if (family_ibd$IBD[j]==2) {
        s1 <- s2 <- 0
        s <- (sum(geno1) + sum(geno2))/2
        ambig <- 0
        xns <- 0
        xs <- 1*(s>0)
      }
      if (family_ibd$IBD[j]==1) {
        s <- s1 <- s2 <- 0
        geno <- geno1 + geno2
        ambig <- which(geno1 ==1 & geno2==1)
        if (length(ambig)>0) {
          for (ai in ambig) {
            s_ns_status <- rbinom(1, 1, prob_shared[ai]) #impute if the variant is shared or not
            s <- s + 1*(s_ns_status==1)
            s1 <- s1 + 2*(s_ns_status==0)
          }
        }
        
        s <- s + sum(geno %in% c(3,4))
        s1 <- s1 + sum(geno %in% c(1,3)) + 2*sum(geno == 4)
        s2 <- 0
        ambig <- sum(geno ==2)
        xns<- 1*(s1>0)
        xs <- 0.5*(s>0)
      }
      family_data[j,] <- c(family_ibd$family[j],family_ibd$sis1[j],family_ibd$sis2[j],family_ibd$IBD[j],hs_i,hns_i,s,s1,s2,ambig,xns,xs)
    }
    # apply(family_data,2,sum)
    family_data <- as.data.frame(family_data)
    family_data[,c(4:12)] <- sapply(family_data[,c(4:12)],as.numeric)
    
    n_case <- sum(as.numeric(family_data[,'n_case']))
    n_control <- sum(as.numeric(family_data[,'n_control']))
    S1_idx <- which(family_data$IBD==1) #which sibpair is S=1 for ramdom pairing
    no_S1 <- length(which(family_data$IBD==1))
    
    # xns <- sum(as.numeric(family_data[,'xns']))
    xns <-  sum(family_data[which(family_data$IBD==0), c("ns1", "ns2")]>0) + sum(family_data[which(family_data$IBD==1), "ns1"]>0)
    # xs <- sum(as.numeric(family_data[,'xs']))
    if (no_S1>=2) {
      xs <- sum(family_data[which(family_data$IBD==2), "s"]>0) + 
        sum((family_data[S1_idx[seq(1, no_S1-(no_S1 %% 2), by=2)], "s"] + 
               family_data[S1_idx[seq(2, no_S1-(no_S1 %% 2), by=2)], "s"])>0)  
    } else {
      xs <- sum(family_data[which(family_data$IBD==2), "s"]>0) 
    }
    
    p1 <- xs/n_case
    p2 <- xns/n_control
    p <- (xs+xns)/(n_case+n_control)
    TD <- p1-p2
    var <- p*(1-p)*(1/n_case+1/n_control)
    test_stat <- TD^2/var
    p.value=pchisq(test_stat, df=1, lower=F)
    
    ##############################################
    # create regression data from family_data
    # for sib-pair with IBD=2: z=1, y=I(s>0)
    # for sib-pair with IBD=0: z=0, create 2 obs for each pair. If xns=0: 2 y's are 0
    #                                                           if xns=1: 1 y=0, 1 y=1
    #                                                           if xns=2: 2 y's are 1
    # for sib-pair with IBD=1: z=1, randomly select 2 sib-pairs with IBD=1, select 1 pair, y=I(sum(s)>0)
    #                          z=0, y=I(ns1>0)
    regression_data <- data.frame()
    
    indi2 <- which(family_data[,'IBD']==2)
    if (length(indi2)>0) {
      ibd2_case <- cbind(family_data[indi2,c('family','sis1','sis2','IBD','s')],1)
      ibd2_case$y <- 1*(ibd2_case$s>0)
      colnames(ibd2_case)  <- c('family','sis1','sis2','IBD','s','z','y')
      ibd2_case <- ibd2_case[,c('family','sis1','sis2','IBD','z','y')]
      regression_data <- rbind(regression_data,ibd2_case)
    }
    indi0 <- which(family_data[,'IBD']==0)
    if (length(indi0)>0) {
      ibd0_control <- cbind(family_data[indi0,c('family','sis1','sis2','IBD','xns')],0)
      colnames(ibd0_control)  <- c('family','sis1','sis2','IBD','xns','z')
      # xns = 0
      ibd0_xns0 <- ibd0_control[which(ibd0_control$xns==0),]
      if (nrow(ibd0_xns0)>0) {
        ibd0_xns0 <- ibd0_xns0[rep(seq_len(nrow(ibd0_xns0)), each = 2),]
        ibd0_xns0$y <- 0
      }
      # xns = 1
      ibd0_xns1 <- ibd0_control[which(ibd0_control$xns==1),]
      if (nrow(ibd0_xns1)>0) {
        ibd0_xns1 <- ibd0_xns1[rep(seq_len(nrow(ibd0_xns1)), each = 2),]
        ibd0_xns1$y <- rep(c(1,0),nrow(ibd0_xns1)/2)
      }
      # xns = 2
      ibd0_xns2 <- ibd0_control[which(ibd0_control$xns==2),]
      if (nrow(ibd0_xns2)>0) {
        ibd0_xns2 <- ibd0_xns2[rep(seq_len(nrow(ibd0_xns2)), each = 2),]
        ibd0_xns2$y <- 1
      }
      ibd0_control_all <- rbind(ibd0_xns0,ibd0_xns1,ibd0_xns2)
      ibd0_control_all <- ibd0_control_all[,c('family','sis1','sis2','IBD','z','y')]
      regression_data <- rbind(regression_data,ibd0_control_all)
    }
    indi1 <- which(family_data[,'IBD']==1)
    if (length(indi1)>0) {
      # control
      ibd1_control <- cbind(family_data[indi1,c('family','sis1','sis2','IBD','ns1','xns')],0)
      colnames(ibd1_control)  <- c('family','sis1','sis2','IBD','ns1','xns','z')
      ibd1_control$y <- 1*(ibd1_control$ns1>0)
      ibd1_control <- ibd1_control[,c('family','sis1','sis2','IBD','z','y')]
      regression_data <- rbind(regression_data,ibd1_control)
      # case
      ibd1_case <- cbind(family_data[indi1,c('family','sis1','sis2','IBD','s')],1)
      no_S1 <- length(indi1)
      if (no_S1>=2) {
        y = 1*((ibd1_case[c(1:no_S1)[seq(1, no_S1-(no_S1 %% 2), by=2)],'s'] + ibd1_case[c(1:no_S1)[seq(2, no_S1-(no_S1 %% 2), by=2)], "s"])>0) 
        ibd1_case <- ibd1_case[c(1:no_S1)[seq(1, no_S1-(no_S1 %% 2), by=2)],]
        ibd1_case$y <- y
        colnames(ibd1_case)  <- c('family','sis1','sis2','IBD','s','z','y')
        ibd1_case <- ibd1_case[,c('family','sis1','sis2','IBD','z','y')]
        regression_data <- rbind(regression_data,ibd1_case)
      }
    }
    
    # table(regression_data[,c('z','y')])
    regression_data$prs_sum <- NA
    regression_data$prs_diff <- NA
    for (i in 1:nrow(regression_data)) {
      prs1 = prs_all$prs[which(grepl(regression_data$sis1[i],rownames(prs_all)))]
      prs2 = prs_all$prs[which(grepl(regression_data$sis2[i],rownames(prs_all)))]
      regression_data$prs_sum[i] <- prs1+prs2
      regression_data$prs_diff[i] <- abs(prs1-prs2)
    }
    
    ##############################################
    # create regression data from family_data using MP's method
    # Y: # shared/non-shared allels
    # X: sharinf status: 0: non-shared; 1: shared
    # 1 obs for IBD=0/2, 2 obs (i with X=0, and 1with X=1) for IBD=1
    # N: for shared obs, N=10*IBD; for non-shared obs: N=40-20*IBD
    shared_dat <- matrix(0, nrow=n_sample, ncol=4)
    nonshared_dat <- matrix(0, nrow=n_sample,ncol=4)
    for(j in 1:n_sample){
      S_sib = family_data[j,'IBD']
      if (S_sib == 0){
        nonshared_dat[j,1] =  family_data[j,'family']
        nonshared_dat[j,2] =  family_data[j,'sis1']
        nonshared_dat[j,3] =  family_data[j,'sis2']
        nonshared_dat[j,4] = family_data[j,'ns1'] + family_data[j,'ns2']  
      }
      if (S_sib == 2) {
        shared_dat[j,1] =  family_data[j,'family']
        shared_dat[j,2] =  family_data[j,'sis1']
        shared_dat[j,3] =  family_data[j,'sis2']
        shared_dat[j,4] = family_data[j,'s']  
      }
      if (S_sib == 1) {
        nonshared_dat[j,1] =  family_data[j,'family']
        nonshared_dat[j,2] =  family_data[j,'sis1']
        nonshared_dat[j,3] =  family_data[j,'sis2']
        shared_dat[j,1] =  family_data[j,'family']
        shared_dat[j,2] =  family_data[j,'sis1']
        shared_dat[j,3] =  family_data[j,'sis2']
        nonshared_dat[j,4] = family_data[j,'ns1'] + family_data[j,'ns2']  
        shared_dat[j,4] = family_data[j,'s']  
      }
    }
    shared = as.data.frame(cbind(shared_dat, X=rep(1, n_sample),IBD=family_data$IBD,N=10*family_data$IBD))
    shared$IBD <- as.numeric(shared$IBD)
    if (length(which(shared$IBD==0))>0) {
      shared = shared[-which(shared$IBD==0),]
    }
    
    nonshared = as.data.frame(cbind(nonshared_dat, X=rep(0, n_sample),IBD=family_data$IBD, N=40-20*family_data$IBD))
    nonshared$IBD <- as.numeric(nonshared$IBD)
    if (length(which(nonshared$IBD==2))>0) {
      nonshared = nonshared[-which(nonshared$IBD==2),]
    }
    
    dat = rbind(shared, nonshared)
    colnames(dat) <- c('family','sis1','sis2','Y','X','IBD',"N")
    rownames(dat) = NULL
    dat$prs_sum <- NA
    dat$prs_diff <- NA
    for (i in 1:nrow(dat)) {
      prs1 = prs_all$prs[which(grepl(dat$sis1[i],rownames(prs_all)))]
      prs2 = prs_all$prs[which(grepl(dat$sis2[i],rownames(prs_all)))]
      dat$prs_sum[i] <- prs1+prs2
      dat$prs_diff[i] <- abs(prs1-prs2)
    }
    
    ##########################  fit regression model
    
    if (length(unique(regression_data$z))==1) {
      ELZ_Z_only_p <- ELZ_Z_only_Z <- ELZ_p <- ELZ_Z <- firth_Z_p <- firth_Z_prs_p <- lm_ELZ_Z_only_Z <- lm_ELZ_Z_only_p<- NA
    } else {
      glmfit_Z <- glm(y~z,data=regression_data,family = binomial)
      glmfit_Z_summary <- summary(glmfit_Z)
      ELZ_Z_only_p <- glmfit_Z_summary$coefficients[2,4]
      ELZ_Z_only_Z <- glmfit_Z_summary$coefficients[2,1]
      
      # firth logistic regression
      firthfit_Z <- logistf(y~z,data=regression_data)
      firthfit_Z_summary <- summary(firthfit_Z)
      firth_Z_p <- firthfit_Z_summary$prob[2]
      
      firthfit_Z_prs <- logistf(y~z+prs_sum+prs_diff,data=regression_data)
      firthfit_Z_prs_summary <- summary(firthfit_Z_prs)
      firth_Z_prs_p <- firthfit_Z_prs_summary$prob[2]
      
      glmfit_Z_prs <- glm(y~z+prs_sum+prs_diff,data=regression_data,family = binomial)
      glmfit_Z_prs_summary <- summary(glmfit_Z_prs)
      ELZ_p <- glmfit_Z_prs_summary$coefficients[2,4]
      ELZ_Z <- glmfit_Z_prs_summary$coefficients[2,1]
      
      lmfit_Z <- lm(y~z,data=regression_data)
      lmfit_Z_summary <- summary(lmfit_Z)
      lm_ELZ_Z_only_p <- lmfit_Z_summary$coefficients[2,4]
      lm_ELZ_Z_only_Z <- lmfit_Z_summary$coefficients[2,1]
    }
    if (length(unique(dat$X))==1) {
      MP_ELZ_p <- MP_ELZ_Z <- MP_ELZ_N_p <- MP_ELZ_N_Z <- NA
    } else {
      mp_reg <- glm(as.numeric(Y)~as.numeric(X)+prs_sum+prs_diff,data=dat)
      mp_reg_n <- glm(as.numeric(Y)/as.numeric(N)~as.numeric(X)+prs_sum+prs_diff,data=dat)
      MP_ELZ_Z <- mp_reg$coefficients[2]
      MP_ELZ_p <- summary(mp_reg)$coefficients[2,4]
      MP_ELZ_N_Z <- mp_reg_n$coefficients[2]
      MP_ELZ_N_p <- summary(mp_reg_n)$coefficients[2,4]
    }
    
    out<- rbind(out,c(p1,p2,test_stat,p.value,n_case,n_control,sum(regression_data$z==1),sum(regression_data$z==0),
                      ELZ_Z_only_Z,ELZ_Z_only_p,ELZ_Z,ELZ_p,lm_ELZ_Z_only_Z,lm_ELZ_Z_only_p,firth_Z_p,firth_Z_prs_p,
                      MP_ELZ_Z,MP_ELZ_p,MP_ELZ_N_p,MP_ELZ_N_Z,gene,nrow(region),length(rare_var_indi))) 
  }
  #  cat(c(LZ_p,LZ_Z_p,LZ_Z_p_p,q1,q2,p1,p2,p,c1,c2,gene,nrow(region),length(rare_var_indi)))
  # cat('',sep='"\n"')
}
# sink()
colnames(out) <- c('p1','p2','LZ_Z','LZ_p','n_case','n_control','n_case1','n_control1','ELZ_Z_only_Z','ELZ_Z_only_p',
                   'ELZ_Z','ELZ_p','lm_ELZ_Z_only_Z','lm_ELZ_Z_only_p','firth_Z_p','firth_Z_prs_p',
                   'MP_ELZ_Z','MP_ELZ_p','MP_ELZ_N_p','MP_ELZ_N_Z','gene','n_gene','n_rare_var')
saveRDS(out,paste0("/Volumes/Kexin/ASP/TRAFIC/WES_WGS_shared/results/WGS_WES_shared_chr",chr,".rds"))
