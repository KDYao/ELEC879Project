# Main approach for feature selection and classification


setwd("/Users/yaokundi/Documents/Courses/Queens/ELEC879/Project/Ski")
data <- read.csv("dataset/step3_output_480.csv", header = T, sep = ',')

dataset <- data
#dataset <- read.csv("Rscripts/ski_feature_set.csv", header = T, sep = ',')



num_top_vars <- 50 # Specify the number of top variables you want to plot and print out names
library(nFactors)
library(plyr)
library(FactoMineR)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(caret)
library(parallel) # Libraries parallel and doSNOW are needed here so that we can calculate the relationship between each 
library(doSNOW)   # numerical column and the target column in parallel. Otherwise, it is slow when we have to go over 2100+
# columns
#dataset$SkillLevel <- ifelse(dataset$SkillLevel == "Pro", 1, 0) # convert string to numerical values
#dataset$SkillLevel <- ifelse(dataset$subjectId == "Intermediate1", 0, 1) # convert string to numerical values
# 1: Intermediate1, 2: Pro1, 3: Pro2
dataset$SkillLevel <- ifelse(dataset$subjectId == 1, 0, 1) # convert string to numerical values
dataset$SkillLevel<-as.factor(dataset$SkillLevel)

# We find a small portion of data has missing value, we replace these part of data with 0
################################### Zero-padding to remove missing data
dataset[is.na(dataset)] <- 0

# TODO: Difference between pro1 and pro2, maybe compare t test
# TODO: Lasso, also Random forest
# TODO: SVM see if it can be successfully classified 


# List column names which are not predictors
to_remove <- c("subjectId","experimentId","ExperimentDate","ExperimentTime","seg", "SkillLevel")

# All column names of predictors
NumericalColumns<-colnames(dataset)[!(names(dataset) %in% to_remove)]

dataset_pred <- dataset[,NumericalColumns]

# Test anova in blocks
# 1. correlation --> cor() spearman
# 2. aov for feature selection
# 3. PCA 


# Feature importance
# Wald chi-sq test
# PCA

############################# Direct spearman correlation
# Perform correlation analysis

correlations<-cor(x=dataset_pred, method = "spearman")
highCorr <- findCorrelation(correlations, cutoff = .7)

if (length(highCorr) == 0){
  low_cor_predictors_data=dataset_pred
}else{
  low_cor_names=names(dataset_pred[, -highCorr])
  low_cor_data= dataset_pred[(names(dataset_pred) %in% low_cor_names)]
  vc = varclus(as.matrix(low_cor_data), similarity = "spear")
}

############################# Anova (prev)
no_cores <- max(detectCores() - 1, 1) # Use the number of cores the computer has minus 1 to calculate the relationship in parallel.

# Leave one core out of this parallel computation to make the computer still responsive to 
# other tasks.
cluster <- makeCluster(no_cores)      # Create such mini cluster with the specified number of cores
registerDoSNOW(cluster)               # register the cluster

# Use variable clustering to remove highly correlated data
# vcobj = varclus(as.matrix(dataset[,NumericalColumns]), similarity = "spear")

aov_v <- foreach(i=1:length(NumericalColumns),
                 .packages=c("DescTools"), .combine='c') %dopar%
                 {
                   get('dataset') #distribute the dataset to each node
                   get('NumericalColumns') #distribute the NumericalColumns to each node
                   col1 <- dataset[[NumericalColumns[i]]]
                   index <- is.na(col1) # Replace NA with 0
                   col1[index] <- 0
                   if (max(col1) != min(col1)) # if it is not a constant column, conduct ANOVA to calculate the strength of relationship
                   {
                     # Use ANOVA to calculate the strength of relationships between numerical features and binary target
                     fit <- aov(col1 ~ dataset[['SkillLevel']])
                     tryCatch(EtaSq(fit)[1], error=function(e) 0)    
                   } else{ # if it is a constant column, directly assign 0 as the strength of linear relationship
                     0
                   }
                 }

stopCluster(cluster) # Release the cluster after completing the parallel computation
names(aov_v) <- NumericalColumns
aov_v <- subset(aov_v, names(aov_v)!='SkillLevel')
aov_v <- sort(aov_v, decreasing = TRUE) # sort the strength of the relationship in a descending order

barplot(head(aov_v, num_top_vars), xlab = 'Eta-squared value', beside=TRUE, # print the bar plot of the top N variables
        main = paste('Top', num_top_vars, 'Associated Numerical Variables'), 
        las=2, cex.axis = 0.5, cex.names = 0.45, space=0.5)

importantvars <- names(aov_v)[1:num_top_vars]
print(importantvars)  #print out the variable names of the top important variables

############################# PCA


# Get a dataframe of sorted pca object
getTopPCAVar <- function(pca_var_list, top_val){

  sdev = pca_var_list
  variance_propotion = sdev^2/sum(sdev^2)  
  
  return(as.vector(variance_propotion)[1:top_val])
}

# Correlation between symmetric part of body
# Symmetry data are highly correlated, so that's why we can put them into same groups
boxplot(abs(dataset[,grep('cor', tolower(colnames(dataset)))]),las=3, cex.axis=0.7, main="Correlation between symmetric part of body")

################ PCA to different groups

# Then we classify following types of predictors into different groups:
# Head LeftFoot LeftForeArm LeftHand LeftLowerLeg LeftShoulder 
# LeftUpperArm LeftUpperLeg Pelvis RightFoot RightForeArm RightHand RightLowerLeg 
# RightShoulder RightUpperArm RightUpperLeg T8

# Head:
#head_pred_cols <- grepl('head',tolower(NumericalColumns))

dataset_pca <- dataset_pred

# Scale is supported when building pca
#dataset_pca <- scale(dataset_pca)

# Shoulder
shoulder_pred_cols <- grepl('shoulder',tolower(NumericalColumns))
shoulder_data <- dataset_pca[,shoulder_pred_cols]

# Forearm
forearm_pred_cols <- grepl('forearm', tolower(NumericalColumns))
forearm_data <- dataset_pca[,forearm_pred_cols]

# Hands:
hand_pred_cols <- grepl('hand',tolower(NumericalColumns))
hand_data <- dataset_pca[,hand_pred_cols]

# Lowerleg:
lowerleg_pred_cols <- grepl('lowerleg',tolower(NumericalColumns))
lowerleg_data <- dataset_pca[,lowerleg_pred_cols]

# Foot
foot_pred_cols <- grepl('foot|feet',tolower(NumericalColumns))
foot_data <- dataset_pca[,foot_pred_cols]

# UpperLeg
upperleg_pred_cols <- grepl('upperleg',tolower(NumericalColumns))
upperleg_data <- dataset_pca[,upperleg_pred_cols]

# UpperArm
upperarm_pred_cols <- grepl('upperarm',tolower(NumericalColumns))
upperarm_data <- dataset_pca[,upperarm_pred_cols]

# Trunk (Rest part): Includes Pevis, T8 and Head
# This parts are not assymentary, and has relatively small movement scales
rest_pred_cols<- !(shoulder_pred_cols|forearm_pred_cols|hand_pred_cols|lowerleg_pred_cols|foot_pred_cols|upperleg_pred_cols|upperarm_pred_cols)
rest_data <- dataset_pred[,rest_pred_cols]

top_vals_threshold = 10

pc_shoulder = prcomp(shoulder_data, center = F, scale. = T)
pc_forearm = prcomp(forearm_data, center = F, scale. = T)
pc_hand = prcomp(hand_data, center = F, scale. = T)
pc_lowerleg = prcomp(lowerleg_data, center = F, scale. = T)
pc_foot = prcomp(foot_data, center = F, scale. = T)
pc_upperleg = prcomp(upperleg_data, center = F, scale. = T)
pc_upperarm = prcomp(upperarm_data, center = F, scale. = T)
pc_trunk = prcomp(rest_data, center = F, scale. = T)


name_list = c("Shoulder", "Forearm", "Hand", "LowerLeg", "Foot", "UpperLeg", "UpperArm", "Trunk")

topVariancePCA = as.data.frame(
  cbind(getTopPCAVar(pc_shoulder$sdev,top_vals_threshold), 
      getTopPCAVar(pc_forearm$sdev,top_vals_threshold),
      getTopPCAVar(pc_hand$sdev,top_vals_threshold),
      getTopPCAVar(pc_lowerleg$sdev,top_vals_threshold),
      getTopPCAVar(pc_foot$sdev,top_vals_threshold),
      getTopPCAVar(pc_upperleg$sdev,top_vals_threshold),
      getTopPCAVar(pc_upperarm$sdev,top_vals_threshold),
      getTopPCAVar(pc_trunk$sdev,top_vals_threshold))
)

colnames(topVariancePCA) = name_list

write.csv(topVariancePCA, file = 'temp.csv', sep=',')
library(colorspace)
matplot(x=seq(1,top_vals_threshold), y=topVariancePCA, type = 'l', lwd = 1.5,lty = 1:ncol(topVariancePCA), col = rainbow_hcl(ncol(topVariancePCA)), xlab = "Principle Components Rank", ylab="VariancePropotion")

legend("topright", legend=colnames(topVariancePCA),col=rainbow_hcl(ncol(topVariancePCA)), lty = 1:ncol(topVariancePCA), cex=0.7)


library(factoextra)
pc_shoulder$rotation[,1]

#idx = which(order(pc_shoulder$x[,1])==1)
#pc_shoulder$x[,1][order(abs(pc_shoulder$x[,1]), decreasing=TRUE)[1:10]]

sum(as.vector(pc_shoulder$sdev^2/sum(pc_shoulder$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_forearm$sdev^2/sum(pc_forearm$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_hand$sdev^2/sum(pc_hand$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_lowerleg$sdev^2/sum(pc_lowerleg$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_foot$sdev^2/sum(pc_foot$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_upperleg$sdev^2/sum(pc_upperleg$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_upperarm$sdev^2/sum(pc_upperarm$sdev^2))[1:top_vals_threshold])
sum(as.vector(pc_trunk$sdev^2/sum(pc_trunk$sdev^2))[1:top_vals_threshold])




########### We need some bootstrap and cross validation here
library(glmnet)
library(boot)
library(AUC)



###########################LASSO
bstrap_lasso <- function(df, i){
  sampledata <- df[i,]
  num_rows <- nrow(df)
  # Seperate 70% as training data, rest 30% as testing data
  training_portion <- 0.7
  training_size = round(num_rows *  training_portion)
  # Already bootstrap, not need to resample
  training_data <- sampledata[1:training_size,]
  validation_data <- sampledata[(training_size+1):num_rows,]
  X <- as.matrix(training_data[, 1:(ncol(df)-1)])
  y <- training_data[["SkillLevel"]]
  ridge_lasso_selector <- 1 #by default 1 is lasso
  #num_shrink_values <- 50
  # Build lasso regression
  lasso_fit = glmnet(X, y, alpha=ridge_lasso_selector, family="binomial")
  #plot(lasso_fit)
  # Use 10-fold cross-validation to choose the optimum shrinkage factor lambda
  cvfit = cv.glmnet(X, y, alpha=ridge_lasso_selector, family="binomial")
  #plot(cvfit)
  coef(lasso_fit,s=cvfit$lambda.min)
  min(cvfit$cvm)
  # Predict with the optimum lambda value
  yhat = predict(lasso_fit, cvfit$lambda.min, newx = as.matrix(validation_data[,1:(ncol(df)-1)]), type = "response")
  # Convert to classes according to sigmoid function
  yhat_fac = as.factor(ifelse(yhat>0.5, 1,0))
  # Accuracy of prediction result
  accuracy = round(mean(yhat_fac==validation_data[["SkillLevel"]])*100, 4)
  
  
  # Next part we calculate AUC (Area under ROC curve)
  roc <- roc(yhat, validation_data[["SkillLevel"]])
  auc = round(auc(roc),6)
  
  #print(paste("AUC=", auc, ", Accuracy=", accuracy, "%", sep=""))
  return(c(auc, accuracy))
}


set.seed(1234)
# 1. Use low correlated data
data_boot = low_cor_data
data_boot$SkillLevel = dataset$SkillLevel

# 2. Use dataset from anova
data_boot <- dataset[importantvars]
data_boot$SkillLevel = dataset$SkillLevel


# Resample dataset, cause Skill level are uniform in two groups
dataset_idx = sample(nrow(data_boot))
b <- boot(data_boot[dataset_idx,], bstrap_lasso, R=100, parallel="multicore")
auc_list = b$t[,1]
accuracy_list = b$t[,2]

mean(auc_list)
mean(accuracy_list)

########################################Logistic Regression
library(rms)


# 1. Use low correlated data
data_logit = low_cor_data
data_logit$SkillLevel = dataset$SkillLevel

# 2. Use dataset from anova
data_logit <- dataset[importantvars]
data_logit$SkillLevel = dataset$SkillLevel


data_logit = data_logit[sample(nrow(data_logit)),]
data_logit[is.na(data_logit)] <- 0

training_set_logit = data_logit[1:as.integer(nrow(data_logit)*0.7),]
testing_set_logit = data_logit[(1+as.integer(nrow(data_logit)*0.7) : nrow(data_logit)),]

training_set_logit[is.na(training_set_logit)] <- 0
testing_set_logit[is.na(testing_set_logit)] <- 0

toform = paste("SkillLevel~",paste(names(data_logit)[-ncol(data_logit)],collapse="+"))
form=as.formula(toform)

logit_fit = lrm(formula = form, data = training_set_logit,x=TRUE, y=TRUE)
logit_anova_obj = anova(logit_fit)


##########investigate previous anova result

plot(varclus(as.matrix(x),similarity = "spearman"), cex=0.95, ylab="Spearman rho-square")
abline(h=1-0,8,col="red", lwd=2)
##############################################SVM

#install.packages("e1071")
library(e1071)

data_svm <- data
data_svm$subjectId = as.factor(data_svm$subjectId)
data_svm[is.na(data_svm)] <- 0

shuffled_index = sample(nrow(data_svm))
training_size = round(nrow(data_svm) *  0.7)
training_set <- data_svm[shuffled_index[1:training_size],]
testing_set <- data_svm[shuffled_index[(training_size+1):nrow(data_svm)],]

svm_model <- svm(subjectId ~ ., data = training_set, kernal="sigmoid")
summary(svm_model)

pred <- predict(svm_model, testing_set)
table(pred, testing_set$subjectId)


# Result

# 1: Intermediate1, 2: Pro1, 3: Pro2
# pred  1  2  3
# 1 97  0  0
# 2  0 35  0
# 3  0  6 39

