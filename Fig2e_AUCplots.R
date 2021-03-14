##########################################################

rm(list=ls())
library(tidyverse)
library(caret)     
library(glmnet)
library("l0ara")
library(ggplot2)
library(caTools)
library(magrittr)
library(ROCR)
library(pROC)
library(glmnetUtils)
library(ncpen)
library(stargazer) 
library(broom)
library(ncvreg)
library(plyr)
library(pROC)

setwd("F:\\Desktop2_clone\\Donghai_Desktop2\\NCpaper_argument_discuss\\Redo_analysis\\For_testAUC")
# load(".RData")

############################# 读入数据 for data #############################

ImSig <- readRDS('ImSig.rds')

ImSig

length(ImSig)


setwd("F:\\Desktop2_clone\\Donghai_Desktop2\\NCpaper_argument_discuss\\Redo_analysis\\For_testAUC")

#### GSE91061 project
GSE91061_matrix <- readRDS('BMS038.Pre.CountTable.normalized.log.rds')
GSE91061_pData <- readRDS('BMS038_phenoData.rds')

GSE91061.pheno <- GSE91061_pData@data

range(GSE91061_matrix)


GSE91061_matrix <- t(GSE91061_matrix)

GSE91061_matrix <- data.frame(GSE91061_matrix)

GSE91061_matrix <- GSE91061_matrix[rownames(GSE91061.pheno),]


all(rownames(GSE91061_matrix)==rownames(GSE91061.pheno))

dim(GSE91061.pheno)

GSE91061.pheno$Lable <- 0
GSE91061.pheno$Lable[GSE91061.pheno$class=='Responder'] <- 1
GSE91061.pheno


GSE91061_matrix <- cbind(GSE91061.pheno,GSE91061_matrix)

GSE91061_matrix <- GSE91061_matrix[,-c(1,2)]

dim(GSE91061_matrix)

head(GSE91061_matrix[,1:20])

tail(GSE91061_matrix[,1:20])

class(GSE91061_matrix)

colnames(GSE91061_matrix)[1:30]


############################# PRJEB23709 #############################

setwd("F:\\Desktop2_clone\\Donghai_Desktop2\\NCpaper_argument_discuss\\Redo_analysis\\For_testAUC")

#### PRJEB23709 project
PRJEB23709_matrix <- read.table('DATASET-PRJEB23709_Pre_73samples.txt', sep='\t', header=T)
PRJEB23709_pData <- readRDS('PRJEB23709_Pre_73samples_phenoData.rds')

range(PRJEB23709_matrix[,-1])

PRJEB23709_matrix <- PRJEB23709_matrix[!duplicated(PRJEB23709_matrix$Symbol),]

rownames(PRJEB23709_matrix) <- PRJEB23709_matrix$Symbol

PRJEB23709_matrix <- PRJEB23709_matrix[,-1]

PRJEB23709_matrix <- t(PRJEB23709_matrix)

PRJEB23709_matrix <- data.frame(PRJEB23709_matrix)

PRJEB23709_matrix <- PRJEB23709_matrix[rownames(PRJEB23709_pData),]


all(rownames(PRJEB23709_matrix)==rownames(PRJEB23709_pData))

dim(PRJEB23709_pData)

PRJEB23709_pData$Lable <- 0
PRJEB23709_pData$Lable[PRJEB23709_pData$class=='Responder'] <- 1
PRJEB23709_pData


PRJEB23709_matrix <- cbind(PRJEB23709_pData,PRJEB23709_matrix)

PRJEB23709_matrix <- PRJEB23709_matrix[,-c(1,2)]

dim(PRJEB23709_matrix)

head(PRJEB23709_matrix[,1:10])

class(PRJEB23709_matrix)


#############################


#### GSE78220 project
GSE78220_matrix <- readRDS('GSE78220_expressionMatrix.rds')
GSE78220_pData <- readRDS('GSE78220_PhenoInfo2.rds')

rownames(GSE78220_pData) <- GSE78220_pData$sample

GSE78220_matrix <- GSE78220_matrix[,c(2,14:26,32:46)]

range(GSE78220_matrix[,-1])

GSE78220_matrix <- GSE78220_matrix[!duplicated(GSE78220_matrix$Symbol),]

rownames(GSE78220_matrix) <- GSE78220_matrix$Symbol

GSE78220_matrix <- GSE78220_matrix[,-1]

GSE78220_matrix <- t(GSE78220_matrix)

GSE78220_matrix <- data.frame(GSE78220_matrix)

GSE78220_matrix <- GSE78220_matrix[rownames(GSE78220_pData),]


all(rownames(GSE78220_matrix)==rownames(GSE78220_pData))

dim(GSE78220_pData)

GSE78220_pData$Lable <- 0
GSE78220_pData$Lable[GSE78220_pData$class=='nonPD'] <- 1
GSE78220_pData


GSE78220_matrix <- cbind(GSE78220_pData,GSE78220_matrix)

GSE78220_matrix <- GSE78220_matrix[,-c(1,2)]

dim(GSE78220_matrix)

head(GSE78220_matrix[,1:20])

tail(GSE78220_matrix[,1:20])

class(GSE78220_matrix)

colnames(GSE78220_matrix)[1:30]

######################


#### MGSP project: 103 samples
NatMed_103samples_GE_matrix <- readRDS('NatMed_103samples_GE_matrix.rds')
NatMed_103samples_pData <- readRDS('NatMed_103samples_pData.rds')
head(NatMed_103samples_GE_matrix)
range(NatMed_103samples_GE_matrix[,-1])

rownames(NatMed_103samples_GE_matrix) <- NatMed_103samples_GE_matrix$Symbol

NatMed_103samples_GE_matrix <- NatMed_103samples_GE_matrix[,-1]

NatMed_103samples_GE_matrix <- t(NatMed_103samples_GE_matrix)

NatMed_103samples_GE_matrix <- data.frame(NatMed_103samples_GE_matrix)


all(rownames(NatMed_103samples_GE_matrix)==rownames(NatMed_103samples_pData))

NatMed_103samples_pData$Lable <- 0
NatMed_103samples_pData$Lable[NatMed_103samples_pData$class=='Responder'] <- 1
NatMed_103samples_pData


NatMed_103samples_GE_matrix <- cbind(NatMed_103samples_pData,NatMed_103samples_GE_matrix)

NatMed_103samples_GE_matrix <- NatMed_103samples_GE_matrix[,-c(1,2)]

dim(NatMed_103samples_GE_matrix)

head(NatMed_103samples_GE_matrix[,1:10])

class(NatMed_103samples_GE_matrix)

MGSP_matrix <- NatMed_103samples_GE_matrix




###################################################################################################################
###################################################################################################################

############################ without batch effect correction ###########################################################


###################################################################################################################

overlap.GSE91061 <- intersect(colnames(GSE91061_matrix), ImSig)

length(overlap.GSE91061)

overlap.GSE91061


### for data3_PRJEB23709

overlap.PRJEB23709 <- intersect(colnames(PRJEB23709_matrix), ImSig)

length(overlap.PRJEB23709)

overlap.PRJEB23709


overlap.v2 <- intersect(overlap.GSE91061, overlap.PRJEB23709)

length(overlap.v2)

overlap.v2


data1.v2 <- GSE91061_matrix[,c('Lable', overlap.v2)]

dim(data1.v2)
head(data1.v2[,1:10])

dim(data1.v2)

data3 <- PRJEB23709_matrix[,c('Lable', overlap.v2)]

dim(data3)
head(data3[,1:10])
range(data3)

table(data3$Lable)


x <- model.matrix(Lable ~., data1.v2)[,-1]   

y <- data1.v2$Lable 


x.test2 <- model.matrix(Lable ~., data3)[,-1] 

y.test2 <- data3$Lable


setwd("F:/Desktop2_clone/Donghai_Desktop2/NCpaper_argument_discuss/Redo_analysis/For_testAUC/LogReg_self_fromGSE91061/")

set.seed(12345)


cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10)
lambda.min <- cv.ridge$lambda.min
lambda.1se <- cv.ridge$lambda.1se
print("***lambda.min、lambda.1se***")
print(lambda.min)
print(lambda.1se)

ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
coef_ridge <- as.matrix(coef(ridge.model))       


p2 <- predict(ridge.model, newx = x.test2, type = "response")
pred_ridge <- cbind(y.test2, p2)  


###################################################################################################################

### for data2_GSE78220

setwd("F:\\Desktop2_clone\\Donghai_Desktop2\\NCpaper_argument_discuss\\Redo_analysis\\For_testAUC")


overlap.GSE78220 <- intersect(colnames(GSE78220_matrix), ImSig)

length(overlap.GSE78220)

overlap.GSE78220


overlap.v1 <- intersect(overlap.GSE91061, overlap.GSE78220)

length(overlap.v1)

overlap.v1


data1 <- GSE91061_matrix[,c('Lable', overlap.v1)]

dim(data1)
head(data1[,1:10])

dim(data1)

data2 <- GSE78220_matrix[,c('Lable', overlap.v1)]

dim(data2)
head(data2[,1:10])
range(data2)

table(data2$Lable)


x <- model.matrix(Lable ~., data1)[,-1]   

y <- data1$Lable 


x.test <- model.matrix(Lable ~., data2)[,-1] 

y.test <- data2$Lable


setwd("F:/Desktop2_clone/Donghai_Desktop2/NCpaper_argument_discuss/Redo_analysis/For_testAUC/LogReg_self_fromGSE91061")

library(ncvreg)


set.seed(12345)


cv.MCP <- cv.ncvreg(x, y, family ="binomial", penalty="MCP") 
lambda.min <- cv.MCP$lambda.min 

print("*** lambda.min ***")
print(lambda.min)


MCP.model <- ncvreg(x, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
coef_MCP <- coef(MCP.model)


p <- predict(MCP.model, x.test, type = "response")
pred_MCP <- cbind(y.test, p)   


###################################################################################################################

### for data4_MGSP


overlap.MGSP <- intersect(colnames(MGSP_matrix), ImSig)

length(overlap.MGSP)

overlap.MGSP


overlap.v3 <- intersect(overlap.GSE91061, overlap.MGSP)

length(overlap.v3)

overlap.v3


data1.v3 <- GSE91061_matrix[,c('Lable', overlap.v3)]

dim(data1.v3)
head(data1.v3[,1:10])

dim(data1.v3)

data4 <- MGSP_matrix[,c('Lable', overlap.v3)]

dim(data4)
head(data4[,1:10])
range(data4)

table(data4$Lable)


x <- model.matrix(Lable ~., data1.v3)[,-1]   

y <- data1.v3$Lable 


x.test3 <- model.matrix(Lable ~., data4)[,-1] 

y.test3 <- data4$Lable


setwd("F:/Desktop2_clone/Donghai_Desktop2/NCpaper_argument_discuss/Redo_analysis/For_testAUC/LogReg_self_fromGSE91061/")

set.seed(12345)


cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10)
lambda.min <- cv.elastic$lambda.min
lambda.1se <- cv.elastic$lambda.1se
print("***lambda.min、lambda.1se***")
print(lambda.min)
print(lambda.1se)

elastic.model <- glmnet(x, y, alpha = 0.5, family = "binomial", lambda = cv.elastic$lambda.1se)
coef_elastic <- as.matrix(coef(elastic.model))       


p2 <- predict(elastic.model, newx = x.test3, type = "response")
pred_elastic <- cbind(y.test3, p2)  



###################################################################################################################
###################################################################################################################

############################ with batch effect correction ###########################################################

###################################################################################################################

missing.GSE91061 <- setdiff(ImSig, overlap.GSE91061)

missing.GSE91061

GSE91061.ImSig.mx <- GSE91061_matrix[,c('Lable',overlap.GSE91061)]

dim(GSE91061.ImSig.mx)

head(GSE91061.ImSig.mx[,1:10])

missing.GSE91061.mx <- matrix(NA, nrow=nrow(GSE91061.ImSig.mx),ncol=length(missing.GSE91061))

rownames(missing.GSE91061.mx) <- rownames(GSE91061.ImSig.mx)
colnames(missing.GSE91061.mx) <- missing.GSE91061

dim(missing.GSE91061.mx)
head(missing.GSE91061.mx)

GSE91061.ImSig.mx.V2 <- cbind(GSE91061.ImSig.mx, missing.GSE91061.mx)

dim(GSE91061.ImSig.mx.V2)
head(GSE91061.ImSig.mx.V2)
tail(GSE91061.ImSig.mx.V2)


###################################################################################################################

missing.PRJEB23709 <- setdiff(ImSig, overlap.PRJEB23709)

missing.PRJEB23709

PRJEB23709.ImSig.mx <- PRJEB23709_matrix[,c('Lable',overlap.PRJEB23709)]

dim(PRJEB23709.ImSig.mx)

head(PRJEB23709.ImSig.mx[,1:10])

missing.PRJEB23709.mx <- matrix(NA, nrow=nrow(PRJEB23709.ImSig.mx),ncol=length(missing.PRJEB23709))

rownames(missing.PRJEB23709.mx) <- rownames(PRJEB23709.ImSig.mx)
colnames(missing.PRJEB23709.mx) <- missing.PRJEB23709

dim(missing.PRJEB23709.mx)
head(missing.PRJEB23709.mx)

PRJEB23709.ImSig.mx.V2 <- cbind(PRJEB23709.ImSig.mx, missing.PRJEB23709.mx)
PRJEB23709.ImSig.mx.V2 <- PRJEB23709.ImSig.mx.V2[,colnames(GSE91061.ImSig.mx.V2)]
dim(PRJEB23709.ImSig.mx.V2)
head(PRJEB23709.ImSig.mx.V2)
tail(PRJEB23709.ImSig.mx.V2)


###################################################################################################################

missing.GSE78220 <- setdiff(ImSig, overlap.GSE78220)

missing.GSE78220

GSE78220.ImSig.mx <- GSE78220_matrix[,c('Lable',overlap.GSE78220)]

dim(GSE78220.ImSig.mx)

head(GSE78220.ImSig.mx[,1:10])

missing.GSE78220.mx <- matrix(NA, nrow=nrow(GSE78220.ImSig.mx),ncol=length(missing.GSE78220))

rownames(missing.GSE78220.mx) <- rownames(GSE78220.ImSig.mx)
colnames(missing.GSE78220.mx) <- missing.GSE78220

dim(missing.GSE78220.mx)
head(missing.GSE78220.mx)

GSE78220.ImSig.mx.V2 <- cbind(GSE78220.ImSig.mx, missing.GSE78220.mx)
GSE78220.ImSig.mx.V2 <- GSE78220.ImSig.mx.V2[,colnames(GSE91061.ImSig.mx.V2)]
dim(GSE78220.ImSig.mx.V2)
head(GSE78220.ImSig.mx.V2)
tail(GSE78220.ImSig.mx.V2)


###################################################################################################################

missing.MGSP <- setdiff(ImSig, overlap.MGSP)

missing.MGSP

MGSP.ImSig.mx <- MGSP_matrix[,c('Lable',overlap.MGSP)]

dim(MGSP.ImSig.mx)

head(MGSP.ImSig.mx[,1:10])

missing.MGSP.mx <- matrix(NA, nrow=nrow(MGSP.ImSig.mx),ncol=length(missing.MGSP))

rownames(missing.MGSP.mx) <- rownames(MGSP.ImSig.mx)
colnames(missing.MGSP.mx) <- missing.MGSP

dim(missing.MGSP.mx)
head(missing.MGSP.mx)

MGSP.ImSig.mx.V2 <- cbind(MGSP.ImSig.mx, missing.MGSP.mx)
MGSP.ImSig.mx.V2 <- MGSP.ImSig.mx.V2[,colnames(GSE91061.ImSig.mx.V2)]
dim(MGSP.ImSig.mx.V2)
head(MGSP.ImSig.mx.V2)
tail(MGSP.ImSig.mx.V2)



############################## PCA analysis:

library(FactoMineR)

library(factoextra)
 
 
getCommon <- function(df, common){
  df <- df[,-1]	
  df <- t(df)  
  df <- df[common,]
}


impute.fun <- function(df) {

	gene.na <- apply(apply(df,2,is.na),2,any)
	gene.na.name <- names(gene.na[gene.na])
	df.impute <- df[,c('Lable',gene.na.name)]
	
	for (i in 2:ncol(df.impute)) {
			mean1 <- mean(df.impute[df.impute[[1]]==1,i], na.rm=T)
			mean2 <- mean(df.impute[df.impute[[1]]==0,i], na.rm=T)
			df.impute[df.impute[[1]]==1 & is.na(df.impute[[i]]),i] <- mean1
			df.impute[df.impute[[1]]==0 & is.na(df.impute[[i]]),i] <- mean2
	}
	df[,c('Lable',gene.na.name)] <- df.impute
	return(df)
}

########################################################################################################
########################################################################################################
###### GSE91061 ######


common_genes <- colnames(GSE91061.ImSig.mx.V2)[-1]

GSE91061.pca <- getCommon(GSE91061.ImSig.mx.V2, common_genes)
PRJEB23709.pca <- getCommon(PRJEB23709.ImSig.mx.V2, common_genes)


dim(GSE91061.pca)
head(GSE91061.pca)

dim(PRJEB23709.pca)
head(PRJEB23709.pca)


all(rownames(GSE91061.pca)==rownames(PRJEB23709.pca))


merged <- cbind(GSE91061.pca, PRJEB23709.pca)

gene.na.status.across.samples <- apply(apply(merged,1,is.na),2,all)
(all.na.genes <- names(gene.na.status.across.samples[gene.na.status.across.samples]))

merged[all.na.genes,]

merged <- merged[!rownames(merged) %in% all.na.genes,]

dim(merged)

head(merged)


merged.t <- t(merged)

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, PRJEB23709.ImSig.mx.V2$Lable)

all(rownames(merged.t)==c(rownames(GSE91061.ImSig.mx.V2), rownames(PRJEB23709.ImSig.mx.V2)))

merged.lable.V2 <- data.frame(merged.lable)

rownames(merged.lable.V2) <- rownames(merged.t)


merged.t <- cbind(merged.lable.V2,merged.t)
colnames(merged.t)[1] <- 'Lable'
dim(merged.t)
head(merged.t)
tail(merged.t)


merged.t.impute <- impute.fun(merged.t)
dim(merged.t.impute)
head(merged.t.impute)
tail(merged.t.impute)

merged.t.impute.V2 <- t(merged.t.impute[,-1])
dim(merged.t.impute.V2)
head(merged.t.impute.V2)


####################################################################################################################
####################################################################################################################
####################################################################################################################


group_list <- c(rep('GSE91061', dim(GSE91061.pca)[2]), rep('PRJEB23709', dim(PRJEB23709.pca)[2]))
dat.pca <- PCA(t(merged.t.impute.V2), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca,repel =T,
                geom.ind = "point", # show points only (nbut not "text")
                col.ind =  group_list, # color by groups 
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, # Concentration ellipses 
                legend.title = "Dataset",
)
g

dim(dat.pca$call$X)

head(dat.pca$call$X)

any(is.na(dat.pca$call$X))


############################## remove batch effect

library(limma)

all(colnames(merged) == c(rownames(GSE91061.ImSig.mx.V2), rownames(PRJEB23709.ImSig.mx.V2)))

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, PRJEB23709.ImSig.mx.V2$Lable)

merged.lable


GSE91061_PRJEB23709_batch_corrected <- removeBatchEffect(t(dat.pca$call$X), batch=group_list)

dim(GSE91061_PRJEB23709_batch_corrected)
head(GSE91061_PRJEB23709_batch_corrected)


dat.pca.V2 <- PCA(t(GSE91061_PRJEB23709_batch_corrected), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca.V2,repel =T,
                geom.ind = "point", 
                col.ind =  group_list, 
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, 
                legend.title = "Dataset",
)
g


all(dat.pca.V2$call$X == t(GSE91061_PRJEB23709_batch_corrected))


##############################
##############################


GSE91061.with.impute <- t(GSE91061_PRJEB23709_batch_corrected)[rownames(GSE91061.ImSig.mx.V2),]

GSE91061.with.impute <- cbind(GSE91061.ImSig.mx.V2[,1],GSE91061.with.impute)

colnames(GSE91061.with.impute)[1] <- 'Lable'

GSE91061.with.impute <- data.frame(GSE91061.with.impute)

dim(GSE91061.with.impute)

head(GSE91061.with.impute)


PRJEB23709.with.impute <- t(GSE91061_PRJEB23709_batch_corrected)[rownames(PRJEB23709.ImSig.mx.V2),]

PRJEB23709.with.impute <- cbind(PRJEB23709.ImSig.mx.V2[,1],PRJEB23709.with.impute)

colnames(PRJEB23709.with.impute)[1] <- 'Lable'

PRJEB23709.with.impute <- data.frame(PRJEB23709.with.impute)

class(PRJEB23709.with.impute$Lable)

dim(PRJEB23709.with.impute)

head(PRJEB23709.with.impute)


setwd("F:/Desktop2_clone/Donghai_Desktop2/NCpaper_argument_discuss/Redo_analysis/For_testAUC/LogReg_self_fromGSE91061")

data <- GSE91061.with.impute

dim(data)
head(data[,1:10])
range(data)

table(data$Lable)


data2 <- PRJEB23709.with.impute

dim(data2)
head(data2[,1:10])
range(data2)

table(data2$Lable)


###########################################################################################################################################################################################################

x <- model.matrix(Lable ~., data)[,-1]   

y <- data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)

x.test <- model.matrix(Lable ~., data2)[,-1] 

y.test <- data2$Lable


######################### glmnet ##########################
#################################################################################

set.seed(12345)


cv.MCP <- cv.ncvreg(x, y, family ="binomial", penalty="MCP") 
lambda.min <- cv.MCP$lambda.min 

print("*** lambda.min ***")
print(lambda.min)


MCP.model <- ncvreg(x, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
coef_MCP <- coef(MCP.model)


p <- predict(MCP.model, x.test, type = "response")
pred_MCP_with_bc_forPRJ <- cbind(y.test, p)   



###################################################################################################################
###################################################################################################################
###### GSE78220 ######


common_genes <- colnames(GSE91061.ImSig.mx.V2)[-1]

GSE91061.pca <- getCommon(GSE91061.ImSig.mx.V2, common_genes)
GSE78220.pca <- getCommon(GSE78220.ImSig.mx.V2, common_genes)

dim(GSE91061.pca)
head(GSE91061.pca)


dim(GSE78220.pca)
head(GSE78220.pca)


all(rownames(GSE91061.pca)==rownames(GSE78220.pca))


merged <- cbind(GSE91061.pca, GSE78220.pca)

gene.na.status.across.samples <- apply(apply(merged,1,is.na),2,all)
(all.na.genes <- names(gene.na.status.across.samples[gene.na.status.across.samples]))

merged[all.na.genes,]

merged <- merged[!rownames(merged) %in% all.na.genes,]

dim(merged)

head(merged)


merged.t <- t(merged)

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, GSE78220.ImSig.mx.V2$Lable)

all(rownames(merged.t)==c(rownames(GSE91061.ImSig.mx.V2), rownames(GSE78220.ImSig.mx.V2)))

merged.lable.V2 <- data.frame(merged.lable)

rownames(merged.lable.V2) <- rownames(merged.t)


merged.t <- cbind(merged.lable.V2,merged.t)
colnames(merged.t)[1] <- 'Lable'
dim(merged.t)
head(merged.t)
tail(merged.t)


merged.t.impute <- impute.fun(merged.t)
dim(merged.t.impute)
head(merged.t.impute)
tail(merged.t.impute)

merged.t.impute.V2 <- t(merged.t.impute[,-1])
dim(merged.t.impute.V2)
head(merged.t.impute.V2)


group_list <- c(rep('GSE91061', dim(GSE91061.pca)[2]), rep('GSE78220', dim(GSE78220.pca)[2]))
dat.pca <- PCA(t(merged.t.impute.V2), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca,repel =T,
                geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
                col.ind =  group_list, # color by groups 颜色组
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, # Concentration ellipses 集中成椭圆
                legend.title = "Dataset",
)
g

dim(dat.pca$call$X)

head(dat.pca$call$X)

any(is.na(dat.pca$call$X))


############################## remove batch effect

library(limma)

all(colnames(merged) == c(rownames(GSE91061.ImSig.mx.V2), rownames(GSE78220.ImSig.mx.V2)))

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, GSE78220.ImSig.mx.V2$Lable)

merged.lable


GSE91061_GSE78220_batch_corrected <- removeBatchEffect(t(dat.pca$call$X), batch=group_list)

dim(GSE91061_GSE78220_batch_corrected)
head(GSE91061_GSE78220_batch_corrected)


dat.pca.V2 <- PCA(t(GSE91061_GSE78220_batch_corrected), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca.V2,repel =T,
                geom.ind = "point", 
                col.ind =  group_list, 
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, 
                legend.title = "Dataset",
)
g


all(dat.pca.V2$call$X == t(GSE91061_GSE78220_batch_corrected))


##############################

GSE91061.with.impute <- t(GSE91061_GSE78220_batch_corrected)[rownames(GSE91061.ImSig.mx.V2),]

GSE91061.with.impute <- cbind(GSE91061.ImSig.mx.V2[,1],GSE91061.with.impute)

colnames(GSE91061.with.impute)[1] <- 'Lable'

GSE91061.with.impute <- data.frame(GSE91061.with.impute)

class(GSE91061.with.impute$Lable)

dim(GSE91061.with.impute)

head(GSE91061.with.impute)


GSE78220.with.impute <- t(GSE91061_GSE78220_batch_corrected)[rownames(GSE78220.ImSig.mx.V2),]

GSE78220.with.impute <- cbind(GSE78220.ImSig.mx.V2[,1],GSE78220.with.impute)

colnames(GSE78220.with.impute)[1] <- 'Lable'

GSE78220.with.impute <- data.frame(GSE78220.with.impute)

dim(GSE78220.with.impute)

head(GSE78220.with.impute)


data <- GSE91061.with.impute

dim(data)
head(data[,1:10])
range(data)

table(data$Lable)


data2 <- GSE78220.with.impute

dim(data2)
head(data2[,1:10])
range(data2)

table(data2$Lable)


###########################################################################################################################################################################################################

x <- model.matrix(Lable ~., data)[,-1]   

y <- data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)

x.test <- model.matrix(Lable ~., data2)[,-1] 

y.test <- data2$Lable


######################### glmnet ##########################

set.seed(12345)


cv.MCP <- cv.ncvreg(x, y, family ="binomial", penalty="MCP") 
lambda.min <- cv.MCP$lambda.min 

print("*** lambda.min ***")
print(lambda.min)


MCP.model <- ncvreg(x, y, family ="binomial", lambda = cv.MCP$lambda.min, penalty="MCP")
coef_MCP <- coef(MCP.model)


p <- predict(MCP.model, x.test, type = "response")
pred_MCP_with_bc_for78 <- cbind(y.test, p)   




###################################################################################################################
###################################################################################################################
###### MGSP ######


common_genes <- colnames(GSE91061.ImSig.mx.V2)[-1]

GSE91061.pca <- getCommon(GSE91061.ImSig.mx.V2, common_genes)
MGSP.pca <- getCommon(MGSP.ImSig.mx.V2, common_genes)


dim(GSE91061.pca)
head(GSE91061.pca)

dim(MGSP.pca)
head(MGSP.pca)


all(rownames(GSE91061.pca)==rownames(MGSP.pca))


merged <- cbind(GSE91061.pca, MGSP.pca)

gene.na.status.across.samples <- apply(apply(merged,1,is.na),2,all)
(all.na.genes <- names(gene.na.status.across.samples[gene.na.status.across.samples]))

merged[all.na.genes,]

merged <- merged[!rownames(merged) %in% all.na.genes,]

dim(merged)

head(merged)


merged.t <- t(merged)

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, MGSP.ImSig.mx.V2$Lable)

all(rownames(merged.t)==c(rownames(GSE91061.ImSig.mx.V2), rownames(MGSP.ImSig.mx.V2)))

merged.lable.V2 <- data.frame(merged.lable)

rownames(merged.lable.V2) <- rownames(merged.t)


merged.t <- cbind(merged.lable.V2,merged.t)
colnames(merged.t)[1] <- 'Lable'
dim(merged.t)
head(merged.t)
tail(merged.t)


merged.t.impute <- impute.fun(merged.t)
dim(merged.t.impute)
head(merged.t.impute)
tail(merged.t.impute)

merged.t.impute.V2 <- t(merged.t.impute[,-1])
dim(merged.t.impute.V2)
head(merged.t.impute.V2)


####################################################################################################################
####################################################################################################################
####################################################################################################################


group_list <- c(rep('GSE91061', dim(GSE91061.pca)[2]), rep('MGSP', dim(MGSP.pca)[2]))
dat.pca <- PCA(t(merged.t.impute.V2), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca,repel =T,
                geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
                col.ind =  group_list, # color by groups 颜色组
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, # Concentration ellipses 集中成椭圆
                legend.title = "Dataset",
)
g

dim(dat.pca$call$X)

head(dat.pca$call$X)

any(is.na(dat.pca$call$X))


############################## remove batch effect

library(limma)

all(colnames(merged) == c(rownames(GSE91061.ImSig.mx.V2), rownames(MGSP.ImSig.mx.V2)))

merged.lable <- c(GSE91061.ImSig.mx.V2$Lable, MGSP.ImSig.mx.V2$Lable)

merged.lable


GSE91061_MGSP_batch_corrected <- removeBatchEffect(t(dat.pca$call$X), batch=group_list)

dim(GSE91061_MGSP_batch_corrected)
head(GSE91061_MGSP_batch_corrected)


dat.pca.V2 <- PCA(t(GSE91061_MGSP_batch_corrected), graph = FALSE, scale.unit = TRUE) 
g<-fviz_pca_ind(dat.pca.V2,repel =T,
                geom.ind = "point", 
                col.ind =  group_list, 
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, 
                legend.title = "Dataset",
)
g


all(dat.pca.V2$call$X == t(GSE91061_MGSP_batch_corrected))


##############################
##############################


GSE91061.with.impute <- t(GSE91061_MGSP_batch_corrected)[rownames(GSE91061.ImSig.mx.V2),]

GSE91061.with.impute <- cbind(GSE91061.ImSig.mx.V2[,1],GSE91061.with.impute)

colnames(GSE91061.with.impute)[1] <- 'Lable'

GSE91061.with.impute <- data.frame(GSE91061.with.impute)

dim(GSE91061.with.impute)

head(GSE91061.with.impute)


MGSP.with.impute <- t(GSE91061_MGSP_batch_corrected)[rownames(MGSP.ImSig.mx.V2),]

MGSP.with.impute <- cbind(MGSP.ImSig.mx.V2[,1],MGSP.with.impute)

colnames(MGSP.with.impute)[1] <- 'Lable'

MGSP.with.impute <- data.frame(MGSP.with.impute)

class(MGSP.with.impute$Lable)

dim(MGSP.with.impute)

head(MGSP.with.impute)


setwd("F:/Desktop2_clone/Donghai_Desktop2/NCpaper_argument_discuss/Redo_analysis/For_testAUC/LogReg_self_fromGSE91061/")

data <- GSE91061.with.impute

dim(data)
head(data[,1:10])
range(data)

table(data$Lable)


data2 <- MGSP.with.impute

dim(data2)
head(data2[,1:10])
range(data2)

table(data2$Lable)


###########################################################################################################################################################################################################

x <- model.matrix(Lable ~., data)[,-1]   

y <- data$Lable # y <- ifelse(train.data$Lable == "1", 1, 0)

x.test <- model.matrix(Lable ~., data2)[,-1] 

y.test <- data2$Lable


######################### glmnet ##########################
#################################################################################

set.seed(12345)


cv.ridge = cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 10)
lambda.min <- cv.ridge$lambda.min
lambda.1se <- cv.ridge$lambda.1se
print("***lambda.min、lambda.1se***")
print(lambda.min)
print(lambda.1se)


ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
coef_ridge <- as.matrix(coef(ridge.model))        


p <- predict(ridge.model, newx = x.test, type = "response")
pred_ridge_with_bc_forMGSP <- cbind(y.test, p)   


###################################################################################################################
###################################################################################################################

### multiple plots

library(RColorBrewer)

colset <- brewer.pal(n = 6, name = "Dark2")

jpeg(file = "Fig2e_from91_wo_or_with_bc_multiple_plots.jpg")

roc_wo_bc_PRJEB23709 <- plot.roc(pred_ridge[,1], pred_ridge[,2], main = "ROC curves", add =  FALSE, asp = NA, print.auc =F, col = colset[1])

auc1=round(as.numeric(gsub('Area under the curve: ','',roc_wo_bc_PRJEB23709$auc)),digits=3)

roc_wo_bc_GSE78220 <- roc(pred_MCP[,1], pred_MCP[,2], main = "Smoothing")

auc2=round(as.numeric(gsub('Area under the curve: ','',roc_wo_bc_GSE78220$auc)),digits=3)

roc_wo_bc_MGSP <- roc(pred_elastic[,1], pred_elastic[,2], main = "Smoothing")

auc3=round(as.numeric(gsub('Area under the curve: ','',roc_wo_bc_MGSP$auc)),digits=3)


roc_with_bc_PRJEB23709 <- roc(pred_MCP_with_bc_forPRJ[,1], pred_MCP_with_bc_forPRJ[,2], main = "Smoothing")

auc4=round(as.numeric(gsub('Area under the curve: ','',roc_with_bc_PRJEB23709$auc)),digits=3)

roc_with_bc_GSE78220 <- roc(pred_MCP_with_bc_for78[,1], pred_MCP_with_bc_for78[,2], main = "Smoothing")

auc5=round(as.numeric(gsub('Area under the curve: ','',roc_with_bc_GSE78220$auc)),digits=3)

roc_with_bc_MGSP <- roc(pred_ridge_with_bc_forMGSP[,1], pred_ridge_with_bc_forMGSP[,2], main = "Smoothing")

auc6=round(as.numeric(gsub('Area under the curve: ','',roc_with_bc_MGSP$auc)),digits=3)


leg1 = paste0("PRJEB23709 wo BC: AUC = ", auc1)
leg2 = paste0("GSE78220 wo BC: AUC = ", auc2)
leg3 = paste0("MGSP wo BC: AUC = ", auc3)


leg4 = paste0("PRJEB23709 with BC: AUC = ", auc4)
leg5 = paste0("GSE78220 with BC: AUC = ", auc5)
leg6 = paste0("MGSP with BC: AUC = ", auc6)



lines(roc_wo_bc_GSE78220, type = "l", lty = 3, col = colset[2])
lines(roc_wo_bc_MGSP, type = "l", lty = 4, col = colset[3])
lines(roc_with_bc_PRJEB23709, type = "l", lty = 1, col = colset[6])
lines(roc_with_bc_GSE78220, type = "l", lty = 4, col = 'blue',lwd = 4)
lines(roc_with_bc_MGSP, type = "l", lty = 1, col = 'black')

legend("bottomright", bty = "n",
       legend = c(leg1, leg2, leg3, leg4, leg5, leg6), 
       col = c(colset[c(1:3,6)],'blue','black'),lty = c(1,3,4,1,4,1), lwd = rep(4,6))

dev.off()


