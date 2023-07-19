#Normally one would install and import the package here
library(parallel)

set.seed(42)

print(getwd())

print(parallel::detectCores())

path_to_alvm <- "C:\\Users\\lance\\Desktop\\data\\data\\ALVM_imp_maf20perc_w_Target.csv"

alvm_data <- as.matrix(read.csv(path_to_alvm))

phenotype_data <- alvm_data[,2]
feature_data <- alvm_data[,3:length(alvm_data[1,])]

phenotype_data <- matrix(as.numeric(phenotype_data),ncol=1)
dimnames(phenotype_data)[1] <- list(alvm_data[,1])
dimnames(phenotype_data)[2] <- list('Target')

feature_data <- matrix(as.numeric(feature_data),ncol=ncol(feature_data))
dimnames(feature_data)[1] <- list(alvm_data[,1])
dimnames(feature_data)[2] <- list(colnames(alvm_data)[3:length(colnames(alvm_data))])

episcan(
  geno1=feature_data,
  pheno=phenotype_data,
  phetype='case-control',
  outfile="C:\\Users\\lance\\Desktop\\data\\results\\parallel_alvm_results",
  zpthres=0.9,
  chunksize=20
)


# Omitted parameters:
# suffix defaults to .txt
# Scale defaults to true
# ncores defualts to available cores
episcan_parallelizable(
  geno1=feature_data,
  pheno=phenotype_data,
  outfile="C:\\Users\\lance\\Desktop\\data\\results\\parallel_alvm_results",
  zpthres=0.9,
  chunksize=20
)
