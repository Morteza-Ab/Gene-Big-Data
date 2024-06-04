install.packages("RWeka")
install.packages("mice")
install.packages("Hmisc")
noYesinstall.packages("VIM")

library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(ggplot2)

dataset <- read.table("readonly/heart-ch.txt", header = TRUE, sep=",", quote="")
nrow(dataset) 
ncol(dataset) 
dim(dataset)
head(dataset, 0) 
summary(dataset)

md.pattern(dataset)  # (mice) displays all the missing values, NA for missing values

mice_plot <- aggr(dataset, col=c("green","red"),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(dataset), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern")) # (VIM) display graphically missing values
dataset[,"chol_imputed"] <- with(dataset, impute(chol, mean))  # (Hmisc) impute with mean
summary(dataset)


############################################ Data normalization


#dataset_n <- Normalize(chol_imputed ~., data = dataset) # (RWeka) normalizes all numeric except chol_imputed
#dataset_n <- Normalize(~ chol_imputed , data = dataset) # normalizes only chol_imputed
dataset_n <- Normalize(~. , data = dataset) # normalizes all variables
summary(dataset_n)


######################## Discretization has two main method, equal depth and equal width


##equal-depth.
dataset$chol_bin <- as.numeric(cut2(dataset$chol_imputed, g=3)) # (Hmisc) create 3 bins g = 3 quantile bins equal-depth
summary(dataset)
#head(dataset)

### equal-width discretization
dataset[,"chol_bin"] <- as.numeric(cut(dataset[,"chol_imputed"], 5)) # create bins of same width
summary(dataset)
#head(dataset)


########################################## Data Understanding


ggplot(dataset,aes(x=age,y=chol, color=num)) + geom_point(size = 4) # scatterplot (age, chol)

ggplot(dataset, aes(chest_pain, fill=factor(num))) + geom_bar() # stacked histogram for categorical variable

ggplot(dataset, aes(x=age, fill=factor(num))) + geom_bar() # stacked histogram for numeric variable




######################################### Predicting disease from transcriptome data

memory.limit(size=3500)
library(randomForest)
library(class)

############## Load data

setwd("/Users/mortezaabyadeh/Desktop")

mrnaNorm <- read.table("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                       header = F, fill = T, skip = 2)
mrnaIDs <- read.table("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                      header = F, fill = T, nrows = 1)
mrnaIDs <- mrnaIDs[, -1][, -1]

########### Data processing 
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
sampleType <- as.data.frame(samp)
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaNorm)
# 20531 1213 columns are patients (except the 1st for gene name) rows are expression levels for each gene
dim(mrnaIDs)
# 1 1213   the first column is the gene name, the others are one patient per row
dim(mrnaClass)
# 1 1212 one patients per row   1 = tumor, 0 = normal
table(unlist(sampClass))
#   0    1 
# 112 1100                     112 normals and 1100 tumor
sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum) 


geneNames <- mrnaNorm[1] # extract the gene names from mrnaNorm as its first column
dim(geneNames)


mrnaData = t(mrnaNorm[, -1]) # remove first column of mrnaData and transpose it to have genes as columns
rm(samp)
rm(sampClass)
rm(mrnaNorm)
gc()

##############. Classification / prediction on all features

# cell #6
#trainSet <- mrnaData
#testSet <- mrnaData
#trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
#testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
#knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
#cbr.predic = as.vector(knn.predic)
#table(cbr.predic, testClasses)
#tab <- table(cbr.predic, t(testClasses))
#error <- sum(tab) - sum(diag(tab))
#accuracy <- round(100- (error * 100 / length(testClasses)))
#print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)


bssWssFast <- function (X, givenClassArr, numClass=2)
  # between squares / within square feature selection
{
  classVec <- matrix(0, numClass, length(givenClassArr))
  for (k in 1:numClass) {
    temp <- rep(0, length(givenClassArr))
    temp[givenClassArr == (k - 1)] <- 1
    classVec[k, ] <- temp
  }
  classMeanArr <- rep(0, numClass)
  ratio <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    overallMean <- sum(X[, j]) / length(X[, j])
    for (k in 1:numClass) {
      classMeanArr[k] <- 
        sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
    }
    classMeanVec <- classMeanArr[givenClassArr + 1]
    bss <- sum((classMeanVec - overallMean)^2)
    wss <- sum((X[, j] - classMeanVec)^2)
    ratio[j] <- bss/wss
  }
  sort(ratio, decreasing = TRUE, index = TRUE)
}


# select features
dim(mrnaData)
# 1212 20531  matrix
dim(mrnaClass)
# 1 1212
dim(mrnaClassNum)
# 1 1212
dim(geneNames)
# 20531 genes
bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)


mrnaDataReduced <- mrnaData[,bss$ix[1:100]]
dim(mrnaDataReduced)
# 1212  100
trainSet <- mrnaDataReduced
testSet <- mrnaDataReduced
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
# or as.numeric(mrnaClassNum[1,])
testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)


knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1) # knn form 'class' package
knn.predic = as.vector(knn.predic)  # change knn.predic to become a vector
table(knn.predic, testClasses)      # build the confusion matrix
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))  # calculate acuracy
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)   # display acuracy after formating it as a character string


trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,])))  # concatenate gene expressions and class data
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,])))    # concatenate gene expressions and class data
colnames(trainSetClass)[101] <- "class"     # give a name to the class column
#trainSetClass$class <- as.numeric(trainSetClass$class) # for regression
trainSetClass$class <- as.factor(trainSetClass$class)  # for classification
class(trainSetClass$class)      # should be factor or categorical for classification
rf <- randomForest(class ~., trainSetClass,
                   ntree=100,
                   importance=T)      # build randomForest classifier
colnames(testSetClass)[101] <- "class"     # give a name to the class column
testSetClass$class <- as.factor(testSetClass$class)  # for classification
rf.predic <- predict(rf ,testSetClass)  # test the randomForest built model on the test set
rf.predic = as.vector(rf.predic)        # change rf.predic to become a vector
table(rf.predic, testClasses)           # build the confusion matrix
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))      # calculate acuracy
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)


nbRows <- nrow(mrnaDataReduced)
set.seed(22)       # seet random seed so that we always get same samples drawn - since they are random
trainRows <- sample(1:nbRows, .70*nbRows)
trainSet <- mrnaDataReduced[trainRows, ]
testSet <- mrnaDataReduced[-trainRows, ]
dim(trainSet)
dim(testSet)

trainClasses <- unlist(mrnaClassNum[1,trainRows], use.names=FALSE)
testClasses <- unlist(mrnaClassNum[1,-trainRows], use.names=FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
knn.predic = as.vector(knn.predic)
table(knn.predic, testClasses)
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,trainRows])))
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,-trainRows])))
colnames(trainSetClass)[101] <- "class"
trainSetClass$class <- as.factor(trainSetClass$class)  # for classification
class(trainSetClass$class)
# should be factor for classification
rf <- randomForest(class ~., trainSetClass,
                   ntree=100,
                   importance=T)
colnames(testSetClass)[101] <- "class"
testSetClass$class <- as.factor(testSetClass$class)  # for classification
rf.predic <- predict(rf ,testSetClass)
rf.predic = as.vector(rf.predic)
table(rf.predic, testClasses)
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

