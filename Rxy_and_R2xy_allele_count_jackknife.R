library(bootstrap)


### Rxy estimate from complete dataset
### See Do et al. (2015) Nature Genetics 47:126-131 (https://doi.org/10.1038/ng.3186) for details

## load allele count data
LOF.alleles <- read.table("OBF.BF_LOF_allele_counts.txt", sep = "\t", header = FALSE)
MOD.alleles <- read.table("OBF.BF_MOD_allele_counts.txt", sep = "\t", header = FALSE)
SYN.alleles <- read.table("OBF.BF_LOW_SYN_allele_counts.txt", sep = "\t", header = FALSE)

## create dataframes (run each variant type separately as dataframe 1)
df1<-data.frame(LOF.alleles$V2,LOF.alleles$V3)
df1<-data.frame(MOD.alleles$V2,MOD.alleles$V3)
df1<-data.frame(SYN.alleles$V2,SYN.alleles$V3)

colnames(df1)<-c("OBF_alleles","BF_alleles")
str(df1)

# Define variable
OBF_allele_count = df1$OBF_alleles
BF_allele_count = df1$BF_alleles

OBF_n = 24
BF_n = 18

#Rxy estimate from observed sample
Rxy = sum((OBF_allele_count/OBF_n)*(1-(BF_allele_count/BF_n)))/sum((BF_allele_count/BF_n)*(1-(OBF_allele_count/OBF_n)))
Rxy

##### leave-one-out sample jackknife method

# set up pseudocode for jackknife 
n = nrow(df1) # sample size
est.val = numeric(n)   #empty vector to store the jackknife estimates
for (i in 1:n) {
  est.val[i] = sum((OBF_allele_count/OBF_n)[-i] * (1-(BF_allele_count/BF_n))[-i]) / sum((BF_allele_count/BF_n)[-i] * (1-(OBF_allele_count/OBF_n))[-i])
}  ## Rxy estimate for the leave-one-out sample

## mean of the jackknife estimates
mean.jack = mean(est.val)
mean.jack

cbind(mean.jack, Rxy)  ## compare estimate from jackknife and complete dataset

# Jackknife estimate of bias
bias.jack = (n-1)*(mean.jack-Rxy)
bias.jack

# Jackknife estimate of standard error
se.jack = sqrt(((n-1)/n)*sum((est.val-mean.jack)^2))
se.jack

# Jackknife estimate of variance
#sd.jack = var(est.val)
#sd.jack
var.jack = ((n-1)^2*var(est.val/n))
var.jack

# 95% confidence interval
LL=mean.jack-1.96*sqrt(var.jack)  ## lower limit
UL=mean.jack+1.96*sqrt(var.jack)  ## upper limit
LL
UL

cat(mean.jack, bias.jack, se.jack, var.jack, LL, UL)


#################################

### R2xy estimate from complete dataset
### following methods as described in Nigenda-Morales et al.2023 Nature Communications 14:5465 (https://doi.org/10.1038/s41467-023-40052-z)
### as modified from Do et al. (2015) Nature Genetics 47:126-131 (https://doi.org/10.1038/ng.3186)

## load allele count data
LOF.alleles <- read.table("OBF.BF_LOF_allele_counts.txt", sep = "\t", header = FALSE)
MOD.alleles <- read.table("OBF.BF_MOD_allele_counts.txt", sep = "\t", header = FALSE)
SYN.alleles <- read.table("OBF.BF_LOW_SYN_allele_counts.txt", sep = "\t", header = FALSE)

## create dataframes (run each variant type separately as dataframe 1)
df1<-data.frame(LOF.alleles$V2,LOF.alleles$V3)
df1<-data.frame(MOD.alleles$V2,MOD.alleles$V3)
df1<-data.frame(SYN.alleles$V2,SYN.alleles$V3)

colnames(df1)<-c("OBF_alleles","BF_alleles")
str(df1)

# Define variable
OBF_allele_count = df1$OBF_alleles
BF_allele_count = df1$BF_alleles

OBF_n = 24
BF_n = 18

L2xnoty = sum((1 - ((2*OBF_allele_count)*(OBF_n - OBF_allele_count)) / (OBF_n*(OBF_n -1))) * ((2*BF_allele_count)*(BF_n - BF_allele_count)) / (BF_n*(BF_n -1)))
L2ynotx = sum((1 - ((2*BF_allele_count)*(BF_n - BF_allele_count)) / (BF_n*(BF_n -1))) * ((2*OBF_allele_count)*(OBF_n - OBF_allele_count)) / (OBF_n*(OBF_n -1)))

R2xy = L2xnoty/L2ynotx
R2xy

### R2xy jackknife calculations

# set up pseudocode for jackknife
n = nrow(df1) # sample size
est.val = numeric(n)
for (i in 1:n) {
  
  est.val[i] = sum((1-(((2*OBF_allele_count) * (OBF_n - OBF_allele_count))/(OBF_n *(OBF_n - 1)))[-i] ) * (((2*BF_allele_count) * (BF_n - BF_allele_count))/(BF_n *(BF_n - 1)) )[-i]) / sum((1-(((2*BF_allele_count) * (BF_n - BF_allele_count))/(BF_n *(BF_n - 1)))[-i] ) * (((2*OBF_allele_count) * (OBF_n - OBF_allele_count))/(OBF_n *(OBF_n - 1)) )[-i])
  
} ## R2xy estimate for the leave-one-out sample

mean.jack = mean(est.val)
mean.jack

cbind(mean.jack, R2xy)  ## compare estimate from jackknife and complete dataset

# Jackknife estimate of bias
bias.jack = (n-1)*(mean.jack-R2xy)
bias.jack

# Jackknife estimate of standard error
se.jack = sqrt(((n-1)/n)*sum((est.val-mean.jack)^2))
se.jack

# Jackknife estimate of variance
var.jack = ((n-1)^2*var(est.val/n))
var.jack

# 95% confidence interval
LL=mean.jack-1.96*sqrt(var.jack)  ## lower limit
UL=mean.jack+1.96*sqrt(var.jack)  ## upper limit
LL
UL

cat(mean.jack, bias.jack, se.jack, var.jack, LL, UL)