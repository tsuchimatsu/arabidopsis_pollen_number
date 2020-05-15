### Script for selection enrichment analysis 

## Input original iHS scores from data files
install.packages("rehh")
library(rehh)
## First, store the output of function "scan_hh" of library "rehh" in a object "scan_res_all"  
scan_res_all <-  read.csv("/path/to/files")
colnames(scan_res_all) <- c("CHR", "POSITION", "freq_A", "iHH_A", "iHH_D", "IES")
scan_res_all_tmp <- as.data.frame(scan_res_all)
ihs <- ihh2ihs(scan_res_all_tmp) 

## Input window-based iHS scores from data files
## Store window-based iHS scores in a object "window_ihs"
window_ihs <-  read.csv("/path/to/files")
## columns of window_ihs should consist of: 
# window_ihs[,1]: chromosome no.
# window_ihs[,2]: coordinate of start position of each 10-kb window
# window_ihs[,3]: coordinate of end position of each 10-kb window
# window_ihs[,4]: maximum value of iHS score in each 10-kb window


## Input GWAS results 
gwas_pollen <- read.csv("/path/to/GWAS files/abracadabra_pid1_P_norm_median_145_emmax_none_t76.pvals",head=T) # pollen number
gwas_ovule <- read.csv("/path/to/GWAS files/abracadabra_pid1_ovuleNumTot_emmax_none_t76.pvals",head=T) # ovule number
gwas_flower <- list(gwas_ovule,gwas_pollen)

# filename list of GWAS results of Atwell's 107 phenotypes
filenames <- read.table("/path/to/GWAS files/filenames.txt")



enrich_ratio0.05_neg <- numeric(0)
enrich_ratio0.025_neg <- numeric(0)
enrich_ratio0.1_neg <- numeric(0)
enrich_ratio0.01_neg <- numeric(0)

enrich_ratio_flower0.05_neg <- numeric(0)
enrich_ratio_flower0.025_neg <- numeric(0)
enrich_ratio_flower0.1_neg <- numeric(0)
enrich_ratio_flower0.01_neg <- numeric(0)

## enrichment of Atwell's 107 phenotypes
for(i in 1:107){
  path <- paste("/path/to/gwas_output_files", filenames[i,1],sep="")
  gwas_res <- read.csv(path,head=T)

  # if MAF > 0.1, use the followings:
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]

  # if MAF > 0.15, use the followings:
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.15),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.15),2]
  
  top <- unlist(sapply(1:length(sig_pos), function(y){
    window_ihs[(window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]),4]
  }))
  top_win <- unlist(sapply(1:length(sig_pos), function(y){
    which((window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]))
  }))

  enrich_ratio0.05_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/20)]])/length(top[!duplicated(top_win)])
  enrich_ratio0.1_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/10)]])/length(top[!duplicated(top_win)])
  enrich_ratio0.025_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/40)]])/length(top[!duplicated(top_win)])
  enrich_ratio0.01_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/100)]])/length(top[!duplicated(top_win)])
}

## enrichment of pollen and ovule phenotypes
for(i in 1:2){

  gwas_res <- gwas_flower[[i]]
  
  # if MAF > 0.1, use the followings:
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  # if MAF > 0.15, use the followings:
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.15),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.15),2]
  
  top <- unlist(sapply(1:length(sig_pos), function(y){
    window_ihs[(window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]),4]
  }))
  top_win <- unlist(sapply(1:length(sig_pos), function(y){
    which((window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]))
  }))
  enrich_ratio_flower0.05_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/20)]])/length(top[!duplicated(top_win)])
  enrich_ratio_flower0.1_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/10)]])/length(top[!duplicated(top_win)])
  enrich_ratio_flower0.025_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/40)]])/length(top[!duplicated(top_win)])
  enrich_ratio_flower0.01_neg[i] <- length(top[!duplicated(top_win)][top[!duplicated(top_win)] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/100)]])/length(top[!duplicated(top_win)])
   }

### SNP-based enrichment test

enrich_ratio_0.05_neg <- numeric(0)
enrich_ratio_0.025_neg <- numeric(0)
enrich_ratio_0.1_neg <- numeric(0)
enrich_ratio_0.01_neg <- numeric(0)

enrich_ratio_flower_0.05_neg <- numeric(0)
enrich_ratio_flower_0.025_neg <- numeric(0)
enrich_ratio_flower_0.1_neg <- numeric(0)
enrich_ratio_flower_0.01_neg <- numeric(0)

## enrichment of Atwell's 107 phenotypes (SNP-based)
for(i in 1:107){
  path <- paste("/path/to/gwas_output_files", filenames[i,1],sep="")
  gwas_res <- read.csv(path,head=T)

  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  sig_snp <- sapply(1:length(sig_chr),function(x){paste(sig_chr[x],sig_pos[x],sep="_")})
  
  top_win <- which(rownames(ihs$iHS) %in% sig_snp)
  enrich_ratio_0.05_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/20)]  ,1 ])/length(ihs$iHS[top_win,1])
  enrich_ratio_0.1_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/10)]  ,1 ])/length(ihs$iHS[top_win,1])
  enrich_ratio_0.025_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/40)]  ,1 ])/length(ihs$iHS[top_win,1])
  enrich_ratio_0.01_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/100)]  ,1 ])/length(ihs$iHS[top_win,1])
 
}

## enrichment of pollen and ovule phenotypes  (SNP-based)
for(i in 1:2){
  gwas_res <- gwas_flower[[i]]
  
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  sig_snp <- sapply(1:length(sig_chr),function(x){paste(sig_chr[x],sig_pos[x],sep="_")})
  
  top_win <- which(rownames(ihs$iHS) %in% sig_snp)
  
  enrich_ratio_flower_0.05_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/20)] ,1  ])/length(ihs$iHS[top_win,1])
  enrich_ratio_flower_0.1_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/10)] ,1  ])/length(ihs$iHS[top_win,1])
  enrich_ratio_flower_0.025_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/40)] ,1  ])/length(ihs$iHS[top_win,1])
  enrich_ratio_flower_0.01_neg[i] <- length( ihs$iHS[top_win,][ ihs$iHS[top_win,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/100)] ,1  ])/length(ihs$iHS[top_win,1])
}

### Integrate the results of enrichment test into matrices
### This is required for permutation tests and for plotting 

enrich_matrix_neg <- rbind(
  enrich_ratio0.1_neg/0.1,
  enrich_ratio0.05_neg/0.05,
  enrich_ratio0.025_neg/0.025,
  enrich_ratio0.01_neg/0.01
)

enrich_matrix_flower_neg <- rbind(
  enrich_ratio_flower0.1_neg/0.1,
  enrich_ratio_flower0.05_neg/0.05,
  enrich_ratio_flower0.025_neg/0.025,
  enrich_ratio_flower0.01_neg/0.01
)

enrich_matrix_neg[(enrich_matrix_neg==0)]  <- NA
enrich_matrix_flower_neg[(enrich_matrix_flower_neg==0)]  <- NA


## permutation test (window-based)
# enrichment of Atwell's 107 phenotypes
perm_matrix <- matrix(NA,4,107)
for(i in 1:107){
  cat(i)
  perm0.1 <- numeric(0)
  perm0.05 <- numeric(0)
  perm0.025 <- numeric(0)
  perm0.01 <- numeric(0)

  path <- paste("/path/to/gwas_output_files", filenames[i,1],sep="")
  gwas_res <- read.csv(path,head=T)
  
  # if MAF > 0.1, use the followings:
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  # if MAF > 0.15, use the followings:  
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  top <- unlist(sapply(1:length(sig_pos), function(y){
    window_ihs[(window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]),4]
  }))
  top_win <- unlist(sapply(1:length(sig_pos), function(y){
    which((window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]))
  }))
  
  relative_pos <- top_win[!duplicated(top_win)]
  for(j in 1:1000){
    rotate_ <- relative_pos + ceiling(runif(1, 0, length(window_ihs[,1])))
    rotate <- ifelse(rotate_> length(window_ihs[,1]), rotate_ - length(window_ihs[,1]), rotate_)
    
    perm0.1[j] <- length(abs(window_ihs[rotate,4])[abs(window_ihs[rotate,4]) > sort(abs(window_ihs[,4]),decreasing=T)[ceiling(length(window_ihs[,4])/10)]])/length(abs(top[!duplicated(top_win)]))/0.1
    perm0.05[j] <- length(abs(window_ihs[rotate,4])[abs(window_ihs[rotate,4]) > sort(abs(window_ihs[,4]),decreasing=T)[ceiling(length(window_ihs[,4])/20)]])/length(abs(top[!duplicated(top_win)]))/0.05
    perm0.025[j] <- length(abs(window_ihs[rotate,4])[abs(window_ihs[rotate,4]) > sort(abs(window_ihs[,4]),decreasing=T)[ceiling(length(window_ihs[,4])/40)]])/length(abs(top[!duplicated(top_win)]))/0.025
    perm0.01[j] <- length(abs(window_ihs[rotate,4])[abs(window_ihs[rotate,4]) > sort(abs(window_ihs[,4]),decreasing=T)[ceiling(length(window_ihs[,4])/100)]])/length(abs(top[!duplicated(top_win)]))/0.01
  }
  perm_matrix[1,i] <- length(perm0.1[perm0.1>enrich_matrix_neg[1,i]])/1000
  perm_matrix[2,i] <- length(perm0.05[perm0.05>enrich_matrix_neg[2,i]])/1000
  perm_matrix[3,i] <- length(perm0.025[perm0.025>enrich_matrix_neg[3,i]])/1000
  perm_matrix[4,i] <- length(perm0.01[perm0.01>enrich_matrix_neg[4,i]])/1000
}
perm_matrix_neg <- perm_matrix

## permutation test (window-based)
## enrichment of pollen and ovule phenotypes
perm_matrix <- matrix(NA,4,2)
for(i in 1:2){
  cat(i)
  perm0.1 <- numeric(0)
  perm0.05 <- numeric(0)
  perm0.025 <- numeric(0)
  perm0.01 <- numeric(0)
  perm0.005 <- numeric(0)
  
  gwas_res <- gwas_flower[[i]]
  
  # if MAF > 0.1, use the followings:
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  # if MAF > 0.15, use the followings:  
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  top <- unlist(sapply(1:length(sig_pos), function(y){
    window_ihs[(window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]),4]
  }))
  top_win <- unlist(sapply(1:length(sig_pos), function(y){
    which((window_ihs[,2] < sig_pos[y])&(window_ihs[,3] > sig_pos[y])&(window_ihs[,1]==sig_chr[y]))
  }))
  
  relative_pos <- top_win[!duplicated(top_win)]
  for(j in 1:1000){
    rotate_ <- relative_pos + ceiling(runif(1, 0, length(window_ihs[,1])))
    rotate <- ifelse(rotate_> length(window_ihs[,1]), rotate_ - length(window_ihs[,1]), rotate_)
    
    perm0.1[j] <- length(window_ihs[rotate,4][window_ihs[rotate,4] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/10)]])/length(top[!duplicated(top_win)])/0.1
    perm0.05[j] <- length(window_ihs[rotate,4][window_ihs[rotate,4] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/20)]])/length(top[!duplicated(top_win)])/0.05
    perm0.025[j] <- length(window_ihs[rotate,4][window_ihs[rotate,4] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/40)]])/length(top[!duplicated(top_win)])/0.025
    perm0.01[j] <- length(window_ihs[rotate,4][window_ihs[rotate,4] < sort(window_ihs[,4])[ceiling(length(window_ihs[,4])/100)]])/length(top[!duplicated(top_win)])/0.01
    
      }
  perm_matrix[1,i] <- length(perm0.1[perm0.1>enrich_matrix_flower_neg[1,i]])/1000
  perm_matrix[2,i] <- length(perm0.05[perm0.05>enrich_matrix_flower_neg[2,i]])/1000
  perm_matrix[3,i] <- length(perm0.025[perm0.025>enrich_matrix_flower_neg[3,i]])/1000
  perm_matrix[4,i] <- length(perm0.01[perm0.01>enrich_matrix_flower_neg[4,i]])/1000
  
}
perm_matrix_flower_neg <- perm_matrix


### permutation test (SNP-based)
# enrichment of Atwell's 107 phenotypes
perm_matrix <- matrix(NA,4,107)
for(i in 1:107){
  perm0.1 <- numeric(0)
  perm0.05 <- numeric(0)
  perm0.025 <- numeric(0)
  perm0.01 <- numeric(0)
  
  path <- paste("/path/to/gwas_output_files", filenames[i,1],sep="")
  gwas_res <- read.csv(path,head=T)
  
  # if MAF > 0.1, use the followings:
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]

  # if MAF > 0.15, use the followings:  
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  sig_snp <- sapply(1:length(sig_chr),function(x){paste(sig_chr[x],sig_pos[x],sep="_")})
  
  top_win <- which(rownames(ihs$iHS) %in% sig_snp)
  relative_pos <- top_win
  
  for(j in 1:1000){
    rotate_ <- relative_pos + ceiling(runif(1, 0, length(rownames(ihs$iHS))))
    rotate <- ifelse(rotate_> length(rownames(ihs$iHS)), rotate_ - length(rownames(ihs$iHS)), rotate_)
    
    perm0.1[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/10)]  ,1 ])/length(ihs$iHS[top_win,1])/0.1
    perm0.05[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/20)]  ,1 ])/length(ihs$iHS[top_win,1])/0.05
    perm0.025[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/40)]  ,1 ])/length(ihs$iHS[top_win,1])/0.025
    perm0.01[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/100)]  ,1 ])/length(ihs$iHS[top_win,1])/0.01
    
  }
  perm_matrix[1,i] <- length(perm0.1[perm0.1>enrich_matrix_neg[1,i]])/1000
  perm_matrix[2,i] <- length(perm0.05[perm0.05>enrich_matrix_neg[2,i]])/1000
  perm_matrix[3,i] <- length(perm0.025[perm0.025>enrich_matrix_neg[3,i]])/1000
  perm_matrix[4,i] <- length(perm0.01[perm0.01>enrich_matrix_neg[4,i]])/1000
}

## permutation test
## enrichment of pollen and ovule phenotypes (SNP-based)
perm_matrix_flower_neg <- matrix(NA,4,2)
for(i in 1:2){
  cat(i)
  perm0.1 <- numeric(0)
  perm0.05 <- numeric(0)
  perm0.025 <- numeric(0)
  perm0.01 <- numeric(0)
  
  gwas_res <- gwas_flower[[i]]
  
  # if MAF > 0.1, use the followings:
  #sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  #sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
  
  # if MAF > 0.15, use the followings:  
  sig_chr <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),1]
  sig_pos <- gwas_res[(-log10(gwas_res$scores) > 4)&(gwas_res[,4]>0.1),2]
 
  sig_snp <- sapply(1:length(sig_chr),function(x){paste(sig_chr[x],sig_pos[x],sep="_")})
  
  top_win <- which(rownames(ihs$iHS) %in% sig_snp)
  relative_pos <- top_win
  
  for(j in 1:1000){
    rotate_ <- relative_pos + ceiling(runif(1, 0, length(rownames(ihs$iHS))))
    rotate <- ifelse(rotate_> length(rownames(ihs$iHS)), rotate_ - length(rownames(ihs$iHS)), rotate_)
    
    perm0.1[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/10)]  ,1 ])/length(ihs$iHS[top_win,1])/0.1
    perm0.05[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/20)]  ,1 ])/length(ihs$iHS[top_win,1])/0.05
    perm0.025[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/40)]  ,1 ])/length(ihs$iHS[top_win,1])/0.025
    perm0.01[j] <- length( ihs$iHS[rotate,][ ihs$iHS[rotate,3] < sort(ihs$iHS[,3])[ceiling(length(ihs$iHS[,3])/100)]  ,1 ])/length(ihs$iHS[top_win,1])/0.01
    
  }

  perm_matrix_flower_neg[1,i] <- length(perm0.1[perm0.1>enrich_matrix_flower_neg[1,i]])/1000
  perm_matrix_flower_neg[2,i] <- length(perm0.05[perm0.05>enrich_matrix_flower_neg[2,i]])/1000
  perm_matrix_flower_neg[3,i] <- length(perm0.025[perm0.025>enrich_matrix_flower_neg[3,i]])/1000
  perm_matrix_flower_neg[4,i] <- length(perm0.01[perm0.01>enrich_matrix_flower_neg[4,i]])/1000
  
  }



#### generate the enrichment plot 
par(mfrow=c(1,1))
cutoff <- c(1,2,3,4)
plot(cutoff,enrich_matrix_neg[1:4,1],type="n",
     ylim=c(0,7),ylab="Fold enrichment",xlab="Top X % of selection scan",main="",cex=1,xaxt="n")
for(i in 1:107){
  points(cutoff,enrich_matrix_neg[1:4,i],type="l",col="grey")
}
points(cutoff,enrich_matrix_flower_neg[1:4,1],type="l",col="black",lwd=2)
points(cutoff,enrich_matrix_flower_neg[1:4,2],type="l",col="red",lwd=2)

abline(h=1,lty=2,lwd=2)
abline(v=1,lty=2,lwd=1)
abline(v=2,lty=2,lwd=1)
abline(v=3,lty=2,lwd=1)
abline(v=4,lty=2,lwd=1)



### generate the plot of permutation test

perm_matrix_neg[(perm_matrix_neg==1)]  <- NA
perm_matrix_flower_neg[(perm_matrix_flower_neg==1)]  <- NA

perm_matrix_neg[(perm_matrix_neg==0)]  <- 0.001
perm_matrix_flower_neg[(perm_matrix_flower_neg==0)]  <- 0.001

par(mfrow=c(1,1))
cutoff <- c(1,2,3,4)
plot(cutoff,-log10(perm_matrix_neg[1:4,1]),type="n",
     ylim=c(0,3),ylab="-log10(p value)",xlab="Top X % of selection scan",main="",cex=1,xaxt="n")
for(i in 1:107){
  points(cutoff,-log10(perm_matrix_neg[1:4,i]),type="l",col="grey")
}
points(cutoff,-log10(perm_matrix_flower_neg[1:4,1]),type="l",col="black",lwd=2)
points(cutoff,-log10(perm_matrix_flower_neg[1:4,2]),type="l",col="red",lwd=2)
abline(h=-log10(0.05),lty=2,lwd=2)
abline(v=1,lty=2,lwd=1)
abline(v=2,lty=2,lwd=1)
abline(v=3,lty=2,lwd=1)
abline(v=4,lty=2,lwd=1)

