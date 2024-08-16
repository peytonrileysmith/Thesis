#Libraries
require(reshape2)
require(data.table)
require(glmnet)
require(ROCR)
require(ggplot2)
require(hrbrthemes)
library(kableExtra)
require(gridExtra)
library(randomForest)
require(dplyr)
require(grid)

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

#Initialize vectors for analysis metrics
num_vector<-numeric()

cor_vector_Lasso_sim<-numeric()
cor_vector_EN_sim<-numeric()
cor_vector_Ridge_sim<-numeric()
auc_vector_Lasso_sim<-numeric()
auc_vector_EN_sim<-numeric()
auc_vector_Ridge_sim<-numeric()

#Simulated data included; Looping through 1 to 100 metabolites
for (num in 1:100) {
  result_matrix <- sim_metabolites(data_output$y, num)
  X_sim <- cbind(X, result_matrix)
  sampled_indices <- sample(1:nrow(X_sim), size = 8400)
  X_train <- X_sim[sampled_indices, ]
  y_train <- y[sampled_indices]
  Beta_train_Lasso <- solveLasso(X_train, y_train)  # Returns beta for LASSO, alpha=1
  Beta_train_EN <- solveEN(X_train, y_train)  # Returns beta for EN
  Beta_train_Ridge <- solveRidge(X_train, y_train)  # Returns beta for Ridge, alpha=0
  X_test <- X_sim[-sampled_indices, ]
  y_test <- y[-sampled_indices]
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
  corval_Lasso <- cor(y_pred_Lasso, y_test)
  corval_EN <- cor(y_pred_EN, y_test)
  corval_Ridge <- cor(y_pred_Ridge, y_test)
  aucpred_Lasso <- prediction(y_pred_Lasso, y_test)
  aucpred_EN <- prediction(y_pred_EN, y_test)
  aucpred_Ridge <- prediction(y_pred_Ridge, y_test)
  auc_Lasso <- performance(aucpred_Lasso, "auc")@y.values[[1]]
  auc_EN <- performance(aucpred_EN, "auc")@y.values[[1]]
  auc_Ridge <- performance(aucpred_Ridge, "auc")@y.values[[1]]
  # Append results to vectors
  cor_vector_Lasso_sim <- c(cor_vector_Lasso_sim, corval_Lasso)
  cor_vector_EN_sim <- c(cor_vector_EN_sim, corval_EN)
  cor_vector_Ridge_sim <- c(cor_vector_Ridge_sim, corval_Ridge)
  auc_vector_Lasso_sim <- c(auc_vector_Lasso_sim, auc_Lasso)
  auc_vector_EN_sim <- c(auc_vector_EN_sim, auc_EN)
  auc_vector_Ridge_sim <- c(auc_vector_Ridge_sim, auc_Ridge)
  num_vector<-c(num_vector, num)
}

#Make data frames
#Simulated data included
df_cor_Lasso_sim<- data.frame(num_vector, cor_vector_Lasso_sim)
df_auc_Lasso_sim<- data.frame(num_vector, auc_vector_Lasso_sim)
df_cor_EN_sim<- data.frame(num_vector, cor_vector_EN_sim)
df_auc_EN_sim<- data.frame(num_vector, auc_vector_EN_sim)
df_cor_Ridge_sim<- data.frame(num_vector, cor_vector_Ridge_sim)
df_auc_Ridge_sim<- data.frame(num_vector, auc_vector_Ridge_sim)
print(nrow(df_cor_Lasso_sim))
print(nrow(df_auc_Lasso_sim))
print(nrow(df_cor_EN_sim))
print(nrow(df_auc_EN_sim))
print(nrow(df_cor_Ridge_sim))
print(nrow(df_auc_Ridge_sim))

#Simulated
#Correlation plots
#Elastic Net
corplot_EN <- ggplot() +geom_point(data = df_cor_EN_sim, aes(num_vector, cor_vector_EN_sim), color = "green") +  labs(x="Number of Metabolites", y="Correlation",title="Correlation Trend from Elastic Net with Simulated Metabolites") + theme_ipsum() + 
  theme(plot.title = element_text(size = 14))

#LASSO
corplot_lasso <- ggplot() +geom_point(data = df_cor_Lasso_sim, aes(num_vector, cor_vector_Lasso_sim), color = "violet") +  labs(x="Number of Metabolites", y="Correlation",title="Correlation Trend from LASSO with Simulated Metabolites") + theme_ipsum() + 
  theme(plot.title = element_text(size = 14))

#Ridge
corplot_ridge <- ggplot() +geom_point(data = df_cor_Ridge_sim, aes(num_vector, cor_vector_Ridge_sim), color = "skyblue") +  labs(x="Number of Metabolites", y="Correlation",title="Correlation Trend from Ridge with Simulated Metabolites") + theme_ipsum() + 
  theme(plot.title = element_text(size = 14))

#AUC plot
#Elastic Net
aucplot_EN <- ggplot() +
  geom_point(data = df_auc_EN_sim, aes(num_vector, auc_vector_EN_sim), color = "green") +  
  labs(x = "Number of Metabolites", y = "AUC Value", 
       title = "AUC Trend from Elastic Net with Simulated Metabolites") + theme_ipsum()+
  theme(plot.title = element_text(size = 14))

#Ridge
aucplot_ridge <- ggplot() +
  geom_point(data = df_auc_Ridge_sim, aes(num_vector, auc_vector_Ridge_sim), color = "skyblue") +  
  labs(x = "Number of Metabolites", y = "AUC Value", 
       title = "AUC Trend from Ridge with Simulated Metabolites") + theme_ipsum()+
  theme(plot.title = element_text(size = 14))

#Lasso
aucplot_lasso <- ggplot() +
  geom_point(data = df_auc_Lasso_sim, aes(num_vector, auc_vector_Lasso_sim), color = "violet") +  
  labs(x = "Number of Metabolites", y = "AUC Value", 
       title = "AUC Trend from LASSO with Simulated Metabolites") + theme_ipsum()+
  theme(plot.title = element_text(size = 14)) 

ggsave(filename = "aucplotEN_testingmeta.png", plot = aucplot_EN, width = 8, height = 6, dpi = 300)
ggsave(filename = "aucplotridge_testingmeta.png", plot = aucplot_ridge, width = 8, height = 6, dpi = 300)
ggsave(filename = "aucplotlasso_testingmeta.png", plot = aucplot_lasso, width = 8, height = 6, dpi = 300)
ggsave(filename = "ENcorplot_testingmeta.png", plot = corplot_EN, width = 8, height = 6, dpi = 300)
ggsave(filename = "Ridgecorplot_testingmeta.png", plot = corplot_ridge, width = 8, height = 6, dpi = 300)
ggsave(filename = "Lassocorplot_testingmeta.png", plot = corplot_lasso, width = 8, height = 6, dpi = 300)

#Calculate slopes for diminishing returns
df_cor_Lasso_sim <- df_cor_Lasso_sim %>%
  arrange(num_vector) %>%
  mutate(slope_cor_Lasso = (cor_vector_Lasso_sim - lag(cor_vector_Lasso_sim)) / (num_vector - lag(num_vector)))

df_cor_EN_sim <- df_cor_EN_sim %>%
  arrange(num_vector) %>%
  mutate(slope_cor_EN = (cor_vector_EN_sim - lag(cor_vector_EN_sim)) / (num_vector - lag(num_vector)))

df_cor_Ridge_sim <- df_cor_Ridge_sim %>%
  arrange(num_vector) %>%
  mutate(slope_cor_Ridge = (cor_vector_Ridge_sim - lag(cor_vector_Ridge_sim)) / (num_vector - lag(num_vector)))

# Calculate slopes for AUC data
df_auc_Lasso_sim <- df_auc_Lasso_sim %>%
  arrange(num_vector) %>%
  mutate(slope_auc_Lasso = (auc_vector_Lasso_sim - lag(auc_vector_Lasso_sim)) / (num_vector - lag(num_vector)))

df_auc_EN_sim <- df_auc_EN_sim %>%
  arrange(num_vector) %>%
  mutate(slope_auc_EN = (auc_vector_EN_sim - lag(auc_vector_EN_sim)) / (num_vector - lag(num_vector)))

df_auc_Ridge_sim <- df_auc_Ridge_sim %>%
  arrange(num_vector) %>%
  mutate(slope_auc_Ridge = (auc_vector_Ridge_sim - lag(auc_vector_Ridge_sim)) / (num_vector - lag(num_vector)))

df_slope_cor <- data.frame(
  num_vector = df_cor_Lasso_sim$num_vector,
  slope_cor_Lasso = df_cor_Lasso_sim$slope_cor_Lasso,
  slope_cor_EN = df_cor_EN_sim$slope_cor_EN,
  slope_cor_Ridge = df_cor_Ridge_sim$slope_cor_Ridge
)

df_slope_auc <- data.frame(
  num_vector = df_auc_Lasso_sim$num_vector,
  slope_auc_Lasso = df_auc_Lasso_sim$slope_auc_Lasso,
  slope_auc_EN = df_auc_EN_sim$slope_auc_EN,
  slope_auc_Ridge = df_auc_Ridge_sim$slope_auc_Ridge
)

saveRDS(df_cor_Lasso_sim, file = "df_cor_Lasso_sim.rds")
saveRDS(df_auc_Lasso_sim, file = "df_auc_Lasso_sim.rds")
saveRDS(df_cor_EN_sim, file = "df_cor_EN_sim.rds")
saveRDS(df_auc_EN_sim, file = "df_auc_EN_sim.rds")
saveRDS(df_cor_Ridge_sim, file = "df_cor_Ridge_sim.rds")
saveRDS(df_auc_Ridge_sim, file = "df_auc_Ridge_sim.rds")

save_table_as_png <- function(df, file_name) {
  table <- tableGrob(df)
  png(file_name, width = 8, height = 6, units = "in", res = 300)
  grid.draw(table)
  dev.off()
}

# Save data frames as PNG
save_table_as_png(df_slope_cor, "correlation_slopes.png")
save_table_as_png(df_slope_auc, "auc_slopes.png")