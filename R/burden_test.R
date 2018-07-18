library(purrr)
library(dplyr)   #filter, select
library(plyr)
library(kinship2)
library(GMMAT)
library(coxme)
library(GenABEL)
library(ggplot2)
library(ggthemes)

#setwd("C:/Users/Youngjun Park/Desktop/variant_scores/")



read.table("PedFileBPDB.txt",header=T)->fam
FID=as.character(fam$FamilyID)
ID=paste(fam[,"FamilyID"],fam[,"PID"],sep=":")
Dad=ifelse(fam[,"FatherPID"]>0,paste(fam[,"FamilyID"],fam[,"FatherPID"],sep=":"),0)
Mom=ifelse(fam[,"MotherPID"]>0,paste(fam[,"FamilyID"],fam[,"MotherPID"],sep=":"),0)

ped=data.frame(FID,ID,Dad,Mom,stringsAsFactors=FALSE)
fullkins<-makekinship(ped$FID,ped$ID,ped$Dad,ped$Mom)

#read status info
ind.stat <- read.table("indiv_pheno_info.txt", header=TRUE)
ind.stat <- filter(ind.stat, ind.stat$Phenotype != -9)
colnames(ind.stat) <- c("IID","Phenotype")

#read.table("~/silo/new_imputation/oneallele/011718/snp/",header=F,stringsAsFactors = F)->filelist

filelist=c("1_1.sum", "2_1.sum", "3_1.sum", "4_1.sum", "5_1.sum", "6_1.sum", "7_1.sum", "8_1.sum", "9_1.sum", "1_1.rv.sum", "2_1.rv.sum", "3_1.rv.sum", "4_1.rv.sum", "5_1.rv.sum", "6_1.rv.sum", "7_1.rv.sum", "8_1.rv.sum", "9_1.rv.sum")

#filelist=c("10.sum", "11.sum", "12.sum", "10.rv.sum", "11.rv.sum", "12.rv.sum")

L <-matrix(ncol=8,nrow=length(filelist))
S <-matrix(ncol=8,nrow=length(filelist))

extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

for(i in 1:length(filelist)) {
  
  filename1 <- paste("allSTR",filelist[i],sep="")
  filename2 <- paste("del", filelist[i],sep="")
  tmp.sum <- read.table(filename1, stringsAsFactors = FALSE)
  allSTR.sum <- data.frame(t(tmp.sum))
  allSTR.sum <- allSTR.sum[-(1:2),]
  tmp <- data.frame(t(allSTR.sum[, -1]))
  tmp2 <- data.frame(apply(tmp, 2, function(x) as.numeric(as.character(x))))
  tmp3 <- laply(tmp2, sum)
  tmp4 <- data.frame(allSTR.sum$X1)
  tmp4[2] <- tmp3
  allSTR.sum <- tmp4
  colnames(allSTR.sum) <- c("IID","BA")
  filtered.allSTR.sum <- filter(allSTR.sum, is.element(allSTR.sum$IID, ind.stat$IID))
  fin.allSTR.sum <- inner_join(filtered.allSTR.sum, ind.stat)
  
  tmp.sum <- read.table(filename2, stringsAsFactors = FALSE)
  del.sum <- data.frame(t(tmp.sum))
  del.sum <- del.sum[-(1:2),]
  tmp <- data.frame(t(del.sum[, -1]))
  tmp2 <- data.frame(apply(tmp, 2, function(x) as.numeric(as.character(x))))
  tmp3 <- laply(tmp2, sum)
  tmp4 <- data.frame(del.sum$X1)
  tmp4[2] <- tmp3
  del.sum <- tmp4
  colnames(del.sum) <- c("IID","BD")
  filtered.del.sum <- filter(del.sum, is.element(del.sum$IID, ind.stat$IID))
  fin.sum <- inner_join(filtered.del.sum, fin.allSTR.sum)
  ##########
  ##########
  
  #infile = paste(filelist[i],".sum",sep="")
  #read.table(infile,header=T)->b
  b <- fin.sum
  
#  read.table("~/gene_lists/burden_dgv_rare/covariates.csv",header=T,sep="\t")->covariates
#  b = merge(b, covariates, all.x=T)
  
  read.table("admixing.csv",header=T,sep="\t")->admix

  b = merge(b, admix, all.x=T)
  id<-as.vector(b[,"IID"])
  tmpidx<-!is.na(match(dimnames(fullkins)[[1]],id))
  tmpkin<-fullkins[tmpidx,tmpidx]
  tmpkin2<-tmpkin*2

  # b$KBAVG[which(b$NSEG==0)] = sum(b$KB)/sum(b$NSEG)
  # do not need, almost everyone has a rare CNV at this range
  
  lmk1 = lmekin(BD ~ BA + (1|id), data=b, varlist=tmpkin2)
  b$res=lmk1$residuals
  b$R1=rntransform(lmk1$residuals)

  lmk2=lmekin(R1 ~ factor(Phenotype) + EUROPEAN + (1|id), data = b, varlist=tmpkin2)
  
  lt=extract_coxme_table(lmk2)
  p=lt$p[2]
  pADMIX=lt$p[3]
  beta=lt$beta[2]
  betaADMIX=lt$beta[3]
  
  ## pretty sure there is a better/more proper way to calculate confint
  lower=beta-1.96*lt$se[2]
  upper=beta+1.96*lt$se[2]
  x=c(filelist[i],beta,lower,upper,p,betaADMIX,pADMIX, "linear")
  L[i,] = x
  
  gmks = glmmkin(factor(Phenotype) ~ R1 + EUROPEAN, data = b, kins = as.matrix(tmpkin), family=binomial(link="logit"))
  p = (pnorm(-abs(gmks$coefficients)/sqrt(diag(gmks$cov)))*2)[2]
  pADMIX = (pnorm(-abs(gmks$coefficients)/sqrt(diag(gmks$cov)))*2)[3]
  or = exp(gmks$coeff[2])
  betaADMIX = (gmks$coeff[3])

  ## pretty sure there is a better/more proper way to calculate confint
  lower = exp(gmks$coeff[2]-1.96*sqrt(diag(gmks$cov))[2])
  upper = exp(gmks$coeff[2]+1.96*sqrt(diag(gmks$cov))[2])
  x=c(filelist[i],or,lower,upper,p,betaADMIX,pADMIX, "logistic")
  S[i,] = x
  
  graph.name <- paste(filelist[i], ".png", sep="")
  p <- ggplot(b, aes(factor(Phenotype), R1)) + geom_boxplot(color=c("steelblue", "red"))
  p <- p + coord_cartesian(ylim= c(-3, 3))
  p <- p + theme_bw() + theme(panel.grid.minor = element_blank())
  p <- p + labs(x="", y="")
  plot(p)
  ggsave(file=graph.name, height = 5, width = 2)  
}
output = as.data.frame(rbind(L,S))
names(output) = c("Geneset","Effect","Lower","Upper","P_value","Coeff_Admix","Admix_Pvalue", "Model")
write.table(output, "burden_5p1KGP_5pHIP.txt", quote=FALSE, sep="\t")
