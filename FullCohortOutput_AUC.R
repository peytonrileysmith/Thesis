#Libraries
require(reshape2)
require(data.table)
require(glmnet)
require(ROCR)
require(ggplot2)
require(hrbrthemes)
library(kableExtra)
require(gridExtra)
library(dplyr)
library(pROC)
library(knitr)
library(webshot)

# read plink2 raw data with FID and IID, and snp names
readRaw <- function(fname) {
  dat <- fread(fname)
  snpNames <- names(dat)[-(1:6)]
  dat <- as.matrix(dat)
  print("file loaded")
  dat <- apply(dat,2,paste)
  FID <- dat[,1]
  IID <- dat[,2]
  dat <- dat[,-(1:6)]
  dat[dat=="NA"] <- "0"
  dat <- apply(dat,2,as.numeric)
  snps <- Matrix(unname(dat),sparse=T)
  return(list(FID=FID,IID=IID,snps=snps,snpNames=snpNames))
}
# prepare dataset of subjects with epi and PCs only
prepare12949 <- function(filename="pheno2_ALL12949_withJPCs.txt",plinkfile="ALL12949", hgliftfile="hglft_genome_13506_51d240.bed",pvfile="pvalues") {
  tab <- read.table(filename,header=TRUE)
  FID <- tab$FID
  IID <- tab$IID
  y <- tab$Affection.Status
  # read genotype data
  tempfile <- readRaw(paste0(plinkfile,".raw"))
  tempfile$snps <- 2-tempfile$snps
  fam <- fread(paste0(plinkfile,".fam"))
  # assemble geno matrix with correct FID & IID
  geno <- matrix(0,nrow=length(FID),ncol=ncol(tempfile$snps))
  for(i in 1:length(FID)) {
    w <- fam$V1==FID[i] & fam$V2==IID[i]
    if(sum(w)!=1) stop("matching problem")
    geno[i,] <- tempfile$snps[w,]
  }
  # initialize Xcopy with snps
  Xcopy <- geno
  if(any(is.na(Xcopy))) {
    stop("NA") # no NA here
  }
  # decode apoe status from ## to 0,1,2
  temp <- paste(tab$APOE)
  apoe <- rep(0,length(temp))
  for(i in 1:length(temp)) {
    z <- strsplit(temp[i],"")[[1]]
    if(z[1]=="4") apoe[i] <- apoe[i]+1
    if(z[2]=="4") apoe[i] <- apoe[i]+1
  }
  metaOnly <- cbind(tab$PC1,tab$PC2,tab$PC3,tab$PC4,tab$PC5,tab$PC6,tab$PC7,tab$PC8,tab$PC9,tab$PC10, apoe, tab$Sex, tab$Age)
  # population identifiers
  pop <- paste(tab$pop)# "AA", "HISP", "NHW", "NA"
  # X is Xcopy with epi covariates
  X <- cbind(Xcopy, metaOnly)
  # get NAs
  removeNA <- is.na(y) | is.na(tab$APOE) | is.na(tab$Sex) | is.na(tab$Age) | (pop=="NA")
  # remove NA
  FID <- FID[!removeNA]
  IID <- IID[!removeNA]
  y <- y[!removeNA]
  X <- X[!removeNA,]
  Xcopy <- Xcopy[!removeNA,]
  metaOnly <- metaOnly[!removeNA,]
  pop <- pop[!removeNA]
  pop2 <- as.factor(pop)
  pop2 <- relevel(pop2, ref = "NHW")
  X<-cbind(X, as.factor(pop2))
  # prepare summary stats
  bim <- fread(paste0(plinkfile,".bim"))
  pvalues <- numeric(nrow(bim))
  hglift <- fread(hgliftfile)
  pv <- fread(pvfile)
  for(i in 1:nrow(bim)) {
    pvalues[i] <- pv$V1[hglift$V2==bim$V4[i]]
  }
  return(list( FID=FID, IID=IID, y=y, X=X, Xcopy=Xcopy, metaOnly=metaOnly, pop=pop, N=sum(!removeNA), pvalues=pvalues))
}

#SNPs
raw<-readRaw("ALL12949.raw")
#raw$snpNames

#Prepare data output file
data_output<-prepare12949() #All viable subjects, 11918 total with NA values removed

set.seed(100)
#Simulate metabolites using only affection status/outcome
sim_metabolites <- function(y,num,R2=c(0.01,0.1),p=0.2,s=1) {
  result_matrix <- NULL
  n <- length(y) #number of subjects
  for(i in 1:num) {
    while(TRUE) {
      indices <- sample(1:n,size=round(n*p),replace=FALSE) #size is proportion of subjects sampled
      temp <- numeric(length(y))
      temp[indices] <- y[indices]
      temp <- temp+rnorm(n,0,s)
      corval <- cor(temp,y)**2 
      if(!is.na(corval)) {
        if(corval>=R2[1] & corval<=R2[2]) { 
          result_matrix <- cbind(result_matrix,temp)#keep only those with R2 between set bounds
          break
        }
      }
    }
  }
  return(result_matrix)
} #p is proportion of y sampled, s is sd of metabolite distribution, R2 bounded between 1% and 10%

#Create matricies from X (inputs) and y (outcome) for use in running testing and training
X<-data_output$X
X <- X[,-26]
y<- as.factor(data_output$y)

#solveLasso will return coefficient given alpha=1
solveLasso<-function(X,y, alpha=1, nfolds=10){
  cvfit <- cv.glmnet(X, y, alpha = alpha,lambda=seq(0.001, 1, by = 0.1), nfolds = nfolds, family="binomial")
  model <- glmnet(X, y, alpha = alpha,lambda = cvfit$lambda.min, standardize= T, family="binomial")
  res <- coef(model, s = cvfit$lambda.min)
  res <- as.vector(res)[-1] # Remove intercept
  return(res)
}

#solveRidge will return coefficient given alpha=0
solveRidge<-function(X,y, alpha=0, nfolds=10){
  cvfit <- cv.glmnet(X, y, alpha = alpha,lambda=seq(0.001, 1, by = 0.1), nfolds = nfolds, family="binomial")
  model <- glmnet(X, y, alpha = alpha, lambda = cvfit$lambda.min, standardize= T, family="binomial")
  res <- coef(model, s = cvfit$lambda.min)
  res <- as.vector(res)[-1] # Remove intercept
  return(res)
}

#solveEN will return coefficient given alpha=0.5
solveEN <- function(X, y, alpha_value = 0.5, lambda_values = seq(0.001, 1, by = 0.1), nfolds = 10) {
  cvfit <- cv.glmnet(X, y, alpha = alpha_value, lambda = lambda_values, nfolds = nfolds, family = "binomial")
  optimal_lambda <- cvfit$lambda.min
  model <- glmnet(X, y, alpha = alpha_value, lambda = optimal_lambda, family = "binomial")
  res <- coef(model, s = optimal_lambda)
  res<-as.vector(res)[-1]
  return(res)
}

# Function to calculate Fisher's z transformation
fisher_z <- function(r) {
  return(atanh(r))
}

# Function to compare correlations
compare_correlations <- function(r1, r2, n1, n2) {
  z1 <- fisher_z(r1)
  z2 <- fisher_z(r2)
  se_diff <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z <- (z1 - z2) / se_diff
  p_value <- 2 * (1 - pnorm(abs(z)))
  return(p_value)
}

#Initialize vectors for analysis metrics
s_vector<- numeric()
cor_vector_Lasso<-numeric()
cor_vector_EN<-numeric()
cor_vector_Ridge<-numeric()
auc_vector_Lasso<-numeric()
auc_vector_EN<-numeric()
auc_vector_Ridge<-numeric()

s_vector_sim<- numeric()
cor_vector_Lasso_sim<-numeric()
cor_vector_EN_sim<-numeric()
cor_vector_Ridge_sim<-numeric()
auc_vector_Lasso_sim<-numeric()
auc_vector_EN_sim<-numeric()
auc_vector_Ridge_sim<-numeric()

p_values_Lasso<-numeric()
p_values_EN<-numeric()
p_values_Ridge<-numeric()

p_vals_cor_Lasso<-numeric()
p_vals_cor_Ridge<-numeric()
p_vals_cor_EN<-numeric()

#Loop for with and without sim metabolites; s is sample size, using the data with epi covars, 10 PCs, sig. loci
result_matrix <- sim_metabolites(data_output$y, 10)
X_sim <- cbind(X, result_matrix)
for (s in seq(1200,10800,by=200)) { 
  #With simulated data
  sampled_indices_sim <- sample(1:nrow(X_sim), size = s, replace = FALSE)
  X_train_sim <- X_sim[sampled_indices_sim, ]
  y_train_sim <- y[sampled_indices_sim]
  Beta_train_Lasso_sim <- solveLasso(X_train_sim, y_train_sim)  # Returns beta for LASSO, alpha=1
  Beta_train_EN_sim <- solveEN(X_train_sim, y_train_sim)  # Returns beta for EN
  Beta_train_Ridge_sim <- solveRidge(X_train_sim, y_train_sim)  # Returns beta for Ridge, alpha=0
  X_test_sim <- X_sim[-sampled_indices_sim, ]
  y_test_sim <- y[-sampled_indices_sim]
  y_pred_Lasso_sim <- X_test_sim %*% Beta_train_Lasso_sim
  y_pred_EN_sim <- X_test_sim %*% Beta_train_EN_sim
  y_pred_Ridge_sim <- X_test_sim %*% Beta_train_Ridge_sim
  y_test_sim<-as.numeric(y_test_sim)
  y_pred_Lasso_sim<-as.numeric(y_pred_Lasso_sim)
  y_pred_EN_sim<-as.numeric(y_pred_EN_sim)
  y_pred_Ridge_sim<-as.numeric(y_pred_Ridge_sim)
  y_diff_Lasso_sim <- y_pred_Lasso_sim - y_test_sim
  y_diff_EN_sim <- y_pred_EN_sim - y_test_sim
  y_diff_Ridge_sim <- y_pred_Ridge_sim - y_test_sim
  corval_Lasso_sim <- cor(y_pred_Lasso_sim, y_test_sim, method="spearman")
  corval_EN_sim <- cor(y_pred_EN_sim, y_test_sim, method="spearman")
  corval_Ridge_sim <- cor(y_pred_Ridge_sim, y_test_sim, method="spearman")
  aucpred_Lasso_sim <- prediction(y_pred_Lasso_sim, y_test_sim)
  aucpred_EN_sim <- prediction(y_pred_EN_sim, y_test_sim)
  aucpred_Ridge_sim <- prediction(y_pred_Ridge_sim, y_test_sim)
  auc_Lasso_sim <- performance(aucpred_Lasso_sim, "auc")@y.values[[1]]
  auc_EN_sim <- performance(aucpred_EN_sim, "auc")@y.values[[1]]
  auc_Ridge_sim <- performance(aucpred_Ridge_sim, "auc")@y.values[[1]]
  roc_Lasso_sim <- roc(y_test_sim, y_pred_Lasso_sim)
  roc_EN_sim <- roc(y_test_sim, y_pred_EN_sim)
  roc_Ridge_sim <- roc(y_test_sim, y_pred_Ridge_sim)
  # Append results to vectors
  s_vector_sim <- c(s_vector_sim, s)
  cor_vector_Lasso_sim <- c(cor_vector_Lasso_sim, corval_Lasso_sim)
  cor_vector_EN_sim <- c(cor_vector_EN_sim, corval_EN_sim)
  cor_vector_Ridge_sim <- c(cor_vector_Ridge_sim, corval_Ridge_sim)
  auc_vector_Lasso_sim <- c(auc_vector_Lasso_sim, auc_Lasso_sim)
  auc_vector_EN_sim <- c(auc_vector_EN_sim, auc_EN_sim)
  auc_vector_Ridge_sim <- c(auc_vector_Ridge_sim, auc_Ridge_sim)
  #Without simulated data
  sampled_indicies<-sample(1:nrow(X), size=s, replace=FALSE)
  X_train<-X[+sampled_indicies,]
  y_train<-y[+sampled_indicies]
  Beta_train_Lasso <- solveLasso(X_train, y_train)  # Returns beta for LASSO, alpha is 1
  Beta_train_EN <- solveEN(X_train, y_train)  # Returns beta for EN
  Beta_train_Ridge <- solveRidge(X_train, y_train)  # Returns beta for Ridge, alpha is 0
  X_test<-X[-sampled_indicies,]
  y_test<-y[-sampled_indicies]
  y_pred_Lasso <- X_test %*% Beta_train_Lasso
  y_pred_EN <- X_test %*% Beta_train_EN
  y_pred_Ridge <- X_test %*% Beta_train_Ridge
  y_test<-as.numeric(y_test)
  y_pred_Lasso<-as.numeric(y_pred_Lasso)
  y_pred_EN<-as.numeric(y_pred_EN)
  y_pred_Ridge<-as.numeric(y_pred_Ridge)
  y_diff_Lasso <- y_pred_Lasso - y_test
  y_diff_EN <- y_pred_EN - y_test
  y_diff_Ridge <- y_pred_Ridge - y_test
  corval_Lasso <- cor(y_pred_Lasso, y_test, method="spearman")
  corval_EN <- cor(y_pred_EN, y_test, method="spearman")
  corval_Ridge <- cor(y_pred_Ridge, y_test, method="spearman")
  aucpred_Lasso <- prediction(y_pred_Lasso, y_test)
  aucpred_EN <- prediction(y_pred_EN, y_test)
  aucpred_Ridge <- prediction(y_pred_Ridge, y_test)
  auc_Lasso <- performance(aucpred_Lasso, "auc")@y.values[[1]]
  auc_EN <- performance(aucpred_EN, "auc")@y.values[[1]]
  auc_Ridge <- performance(aucpred_Ridge, "auc")@y.values[[1]]
  roc_Lasso <- roc(y_test, y_pred_Lasso)
  roc_EN <- roc(y_test, y_pred_EN)
  roc_Ridge <- roc(y_test, y_pred_Ridge)
  # Append results to vectors
  s_vector <- append(s_vector, s)
  cor_vector_Lasso<-append(cor_vector_Lasso, corval_Lasso)
  cor_vector_EN<-append(cor_vector_EN, corval_EN)
  cor_vector_Ridge<-append(cor_vector_Ridge, corval_Ridge)
  auc_vector_Lasso<-append(auc_vector_Lasso, auc_Lasso)
  auc_vector_EN<-append(auc_vector_EN, auc_EN)
  auc_vector_Ridge<-append(auc_vector_Ridge, auc_Ridge)
  #Permutation test to compare ROCs
  test_Lasso <- roc.test(roc_Lasso, roc_Lasso_sim, method = "bootstrap", boot.n = 1000)
  p_values_Lasso <- append(p_values_Lasso, test_Lasso$p.value)
  test_EN <- roc.test(roc_EN, roc_EN_sim, method = "bootstrap", boot.n = 1000)
  p_values_EN <- append(p_values_EN, test_EN$p.value)
  test_Ridge <- roc.test(roc_Ridge, roc_Ridge_sim, method = "bootstrap", boot.n = 1000)
  p_values_Ridge <- append(p_values_Ridge, test_Ridge$p.value)
  #Compare correlations
  p_val_cor_Lasso<- compare_correlations(corval_Lasso_sim, corval_Lasso, s, s)
  p_val_cor_Ridge<-compare_correlations(corval_Ridge_sim, corval_Ridge, s, s)
  p_val_cor_EN<-compare_correlations(corval_EN_sim, corval_EN, s, s)
  p_vals_cor_Lasso<-append(p_vals_cor_Lasso, p_val_cor_Lasso)
  p_vals_cor_Ridge<-append(p_vals_cor_Ridge, p_val_cor_Ridge)
  p_vals_cor_EN<-append(p_vals_cor_EN, p_val_cor_EN)
}

#Make data frames

#No simulated data included
s_vector <- data.frame(s_vector)
cor_vector_Lasso <- data.frame(cor_vector_Lasso)
auc_vector_Lasso <- data.frame(auc_vector_Lasso)
cor_vector_EN <- data.frame(cor_vector_EN)
auc_vector_EN <- data.frame(auc_vector_EN)
cor_vector_Ridge <- data.frame(cor_vector_Ridge)
auc_vector_Ridge <- data.frame(auc_vector_Ridge)
df_cor_Lasso<- data.frame(s_vector, cor_vector_Lasso)
df_auc_Lasso<- data.frame(s_vector, auc_vector_Lasso)
df_cor_EN<- data.frame(s_vector, cor_vector_EN)
df_auc_EN<- data.frame(s_vector, auc_vector_EN)
df_cor_Ridge<- data.frame(s_vector, cor_vector_Ridge)
df_auc_Ridge<- data.frame(s_vector, auc_vector_Ridge)

#Simulated data included
s_vector_sim <- data.frame(s_vector_sim)
cor_vector_Lasso_sim <- data.frame(cor_vector_Lasso_sim)
cor_vector_EN_sim <- data.frame(cor_vector_EN_sim)
cor_vector_Ridge_sim <- data.frame(cor_vector_Ridge_sim)
auc_vector_Lasso_sim <- data.frame(auc_vector_Lasso_sim)
auc_vector_EN_sim <- data.frame(auc_vector_EN_sim)
auc_vector_Ridge_sim <- data.frame(auc_vector_Ridge_sim)
df_cor_Lasso_sim<- data.frame(s_vector_sim, cor_vector_Lasso_sim)
df_auc_Lasso_sim<- data.frame(s_vector_sim, auc_vector_Lasso_sim)
df_cor_EN_sim<- data.frame(s_vector_sim, cor_vector_EN_sim)
df_auc_EN_sim<- data.frame(s_vector_sim, auc_vector_EN_sim)
df_cor_Ridge_sim<- data.frame(s_vector_sim, cor_vector_Ridge_sim)
df_auc_Ridge_sim<- data.frame(s_vector_sim, auc_vector_Ridge_sim)

#P-values
df_p_values_Lasso<-data.frame(s_vector, p_values_Lasso)
df_p_values_Ridge<-data.frame(s_vector, p_values_Ridge)
df_p_values_EN<-data.frame(s_vector, p_values_EN)
df_p_vals_cor_Lasso<-data.frame(s_vector,p_vals_cor_Lasso)
df_p_vals_cor_Ridge<-data.frame(s_vector,p_vals_cor_Ridge)
df_p_vals_cor_EN<-data.frame(s_vector,p_vals_cor_EN)

# Number of tests (p-values) for Bonferroni correction
N_Lasso_roc <- length(p_values_Lasso)
N_Ridge_roc <- length(p_values_Ridge)
N_EN_roc <- length(p_values_EN)

# Bonferroni corrected thresholds
alpha_bonf_Lasso_roc <- 0.05 / N_Lasso_roc
alpha_bonf_Ridge_roc <- 0.05 / N_Ridge_roc
alpha_bonf_EN_roc <- 0.05 / N_EN_roc

# Count significant p-values after Bonferroni correction
sig_Lasso_roc <- (sum(p_values_Lasso < alpha_bonf_Lasso_roc))/length(p_values_Lasso)*100
sig_Ridge_roc <- (sum(p_values_Ridge < alpha_bonf_Ridge_roc))/length(p_values_Ridge)*100
sig_EN_roc <- (sum(p_values_EN < alpha_bonf_EN_roc))/length(p_values_EN)*100

#Create table for full cohort and save
auc_table <- data.frame(
  Metric = c("Proportion of Significant P-Values (<0.05)"),
  LASSO = c(sprintf("%.1f%%", sig_Lasso_roc)),
  Ridge = c(sprintf("%.1f%%", sig_Ridge_roc)),
  Elastic_Net = c(sprintf("%.1f%%", sig_EN_roc))
)

# Create the HTML table for AUC
auc_table_html <- kable(auc_table, format = "html", 
                        col.names = c("Metric", "LASSO", "Ridge", "Elastic Net")) %>%
  kable_styling(full_width = F, position = "center") %>%
  add_header_above(c("Full Sample Bootstrap Results (AUC)" = 4), bold = TRUE)

# Save the AUC HTML table to a file
auc_html_file <- "fullsampletable_auc.html"
save_kable(auc_table_html, file = auc_html_file)
