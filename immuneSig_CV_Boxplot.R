rm(list=ls())

####### R version 3.6.3 (2020-02-29) -- "Holding the Windsock"

library(rocc)
library(ggplot2)
library(Biobase)
library(edgeR)
library(limma)
library(biomaRt)
library(dplyr)
library(cancerclass)
library(cowplot)
library(caret)
library("RColorBrewer")


######################## 
rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}

setwd('C:\\Donghai_Desktop2\\NCpaper_argument_discuss\\Redo_analysis')
ImSig.genes <- readRDS('ImSig.rds')

# GSE78220 28 samples
GSE78220_AltAnalyze <- readRDS('GSE78220_expressionMatrix.rds')
GSE78220_PhenoInfo2 <- readRDS('GSE78220_PhenoInfo2.rds')
GSE78220_AltAnalyze <- GSE78220_AltAnalyze[,c('Symbol',GSE78220_PhenoInfo2$sample)]
GSE78220_AltAnalyze <- GSE78220_AltAnalyze[rowSums(GSE78220_AltAnalyze[,-1]) > 1, ] # same operation with the original paper
head(GSE78220_AltAnalyze)
range(GSE78220_AltAnalyze[,-1])

# GSE91061 = BMS038 51 samples
BMS038.Pre.CountTable.normalized.log <- as.data.frame(readRDS('BMS038.Pre.CountTable.normalized.log.rds'))
BMS038_phenoData <- readRDS('BMS038_phenoData.rds')
BMS038.Pre.CountTable.normalized.log <- BMS038.Pre.CountTable.normalized.log[rowSums(BMS038.Pre.CountTable.normalized.log) > 10,] # same operation with the original paper, otherwise function fit() does not work 
BMS038.Pre.CountTable.normalized.log$Symbol <- rownames(BMS038.Pre.CountTable.normalized.log)
BMS038_PhenoInfo <- BMS038_phenoData@data
head(BMS038.Pre.CountTable.normalized.log)
range(BMS038.Pre.CountTable.normalized.log[,-ncol(BMS038.Pre.CountTable.normalized.log)])

# CC_73samples from PRJEB23709 73 samples
CC_73samples_GE <- read.table('DATASET-PRJEB23709_Pre_73samples.txt',sep="\t",header=T)
CC_73samples_pData <- readRDS('PRJEB23709_Pre_73samples_phenoData.rds')
CC_73samples_GE <- CC_73samples_GE[CC_73samples_GE$Symbol != '',]
CC_73samples_GE_matrix <- CC_73samples_GE[,c('Symbol',CC_73samples_pData$sample)]
head(CC_73samples_GE_matrix)
range(CC_73samples_GE_matrix[,-1])


#### MGSP project: 103 samples
NatMed_103samples_GE_matrix <- readRDS('NatMed_103samples_GE_matrix.rds')
NatMed_103samples_pData <- readRDS('NatMed_103samples_pData.rds')
head(NatMed_103samples_GE_matrix)
range(NatMed_103samples_GE_matrix[,-1])


################# use the same dataset for train and test with 5-fold cross-validation strategy: 
table(GSE78220_PhenoInfo2$class) # 15:13 
table(CC_73samples_pData$class) # 27:46
table(BMS038_PhenoInfo$class) # 25:26 -> GSE91061
table(NatMed_103samples_pData$class) # 56:47 -> MGSP

processData <- function(exp, phenoInfo, features, printngenes=FALSE){
  ### unify response name for cross dataset test
  phenoInfo$class <- as.character(phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'nonPD','Responder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'PD','NonResponder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'Nonresponder','NonResponder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'Progressor','NonResponder',phenoInfo$class)
  
  pData = data.frame(class=phenoInfo$class, sample=phenoInfo$sample,
                     row.names=phenoInfo$sample)
  phenoData <- new("AnnotatedDataFrame",data=pData)
  
  expdata_col_rearranged <- exp[,c("Symbol",phenoInfo$sample)]
  expdata.sig <- expdata_col_rearranged[expdata_col_rearranged$Symbol %in% features,]
  expdata.sig <- expdata.sig[,-1]
  expdata.sig <- as.matrix(expdata.sig)
  #expdata.sig <- expdata.sig[rowSums(expdata.sig) > 1,]
  rownames(expdata.sig) <- features
  
  if (printngenes){
    print(dim(expdata.sig))
  }
  ExpSet_V5 <- ExpressionSet(assayData=as.matrix(expdata.sig),phenoData=phenoData)
  return(ExpSet_V5)
  #return(list(exp=ExpSet_V5,genes=rownames(expdata.sig)))
}

trainANDtestModel <- function(traindata, testdata, features, phenoInfo_train, phenoInfo_test){
  "
  expdata: col - sample, row - gene
  features: gene list
  phenoInfo: dataframe: class, sample
  "
  common <- intersect(traindata$Symbol,testdata$Symbol)
  features <- intersect(common, features)
  
  ExpSet_train <- processData(exp = traindata, phenoInfo = phenoInfo_train, features = features)
  print(ExpSet_train@featureData@data)
  ExpSet_test <- processData(exp = testdata, phenoInfo = phenoInfo_test, features = features)
  
  predictor_V5 <- fit(ExpSet_train, method = "welch.test")
  positive.class <- unique(pData(ExpSet_test)$class)[2]
  negative.class <- unique(pData(ExpSet_test)$class)[1]
  #print(table(pData(ExpSet_test)$class))
  #print(length(features))
  print(length(features))
  prediction_V5 <- predict(predictor_V5, ExpSet_test, as.character(positive.class), ngenes=length(rownames(predictor_V5@predictor)), dist = "cor")
  
  out_V5 <- as.factor(rep(c(1,2),c(table(pData(ExpSet_test)[["class"]])[[negative.class]],table(pData(ExpSet_test)[["class"]])[[positive.class]])))
  z_V5 <- as.numeric(prediction_V5@prediction[,'z'])
  Test_V5 <- cbind(out_V5,z_V5)
  colnames(Test_V5) <- c('grp','res')
  Test_V5 <- as.data.frame(Test_V5)
  
  return(Test_V5)
}


testCV <- function(expdata, phenoInfo, features, num=5){
  set.seed(17)
  folds <- createFolds(y=phenoInfo$sample,k=num)
  test <- data.frame(label=c(), prediction=c(), nfold = c())
  for (i in 1:num){
    fold_test <- phenoInfo[folds[[i]],] #folds[[i]] for test
    fold_train <- phenoInfo[-folds[[i]],] # remaining data for train
    print(table(fold_train$class))
    print(table(fold_test$class))
    
    test0 <- trainANDtestModel(traindata = expdata, testdata = expdata, 
                               features = features, phenoInfo_train = fold_train, phenoInfo_test = fold_test)
    colnames(test0) <- c('label','prediction')
    if (i==1){
      test <- test0
    }else{
      test <- rbind(test, test0)
    }
    test0 <- na.omit(test0)
    auc_value <- as.numeric(rocdata(test0$label, test0$prediction)$stats$auc)
  }
  return(list(auc=auc_value, result=test))
}

############# compare with other signature
IFNG.Sig <- c('IFNG', 'STAT1', 'IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA')
CD8.Sig <- c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G")
PDL1.Sig <- c('PDL1','CD274','PDCD1LG2','PDCD1')
CRMA.Sig <- c('CT1.2', 'MAGEA2', 'MAGEA2A', 'MAGEA2B', 'MAGEA3', 'MAGEA6', 'MAGEA12')
IMPRES.Sig <- c("BTLA", "CD200", "CD200R1", "CD27", "CD276", "CD28", "CD40", "CD80", "CD86", "CEACAM1", "CTLA4", "IDO1",
                "IL2RB", "LAG3", "PVR", "PVRL2", "TIGIT", "TNFRSF18", "TNFRSF4", "TNFRSF9", "PDL1", "HAVCR2", "PDCD1", "PDCD1LG2", "TNFRSF14", "TNFSF4", "TNFSF9", "C10orf54")
IRG.Sig <- c('LEPR','PRLHR','NR2F2','PRL','NRP1','ANGPTL5','IGF1','TNFRSF10B','TNFRSF10A','PLAU','IFI30') # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'
LRRC15.CAF.Sig <- c('MMP11','COL11A1','C1QTNF3','CTHRC1','COL12A1','COL10A1','COL5A2','GJB2','THBS2','AEBP1','MFAP2','LRRC15','PLAU','ITGA11') # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'
T.cell.inflamed.Sig <- c('CD3D','IDO1','CIITA','CD3E','CCL5','GZMK','CD2','HLA-DRA','CXCL13','IL2RG','NKG7','HLA-E','CXCR6','LAG3','TAGAP','CXCL10','STAT1','GZMB')
IPRES.Sig <- c('ANGPT2','AXL','CCL13','CCL2','CCL7','CDH1','FAP','FLT1','1L10','LOXL2','RORS','TAGLN','TWIST2','VEGFA','VEGFC','WNT5A')
Inflammatory.Sig <- c('CCL5','CCR5','PDL1','CD3D','CD3E','CD8A','CIITA','CTLA4','CXCL10','CXCL11','CXCL13','CXCL9','GZMA','GZMB','HLA-DRA','HKA.DRB1','HLA-E',
                      'IDO1','IL2RG','ITGAL','LAG3','NKG7','PDCD1','PRF1','PTPRC','STAT1','TAGAP')
EMT.Sig <- c('CDH1','CDH3','CLDN4','EPCAM','ST14','MAL2','VIM','SNAI2','ZEB2','FN1','MMP2','AGER')
Blood.Sig <- c('ADAM17', 'CDK2', 'CDKN2A', 'DPP4', 'ERBB2', 'HLA-DRA', 'ICOS', 'ITGA4', 'LARGE', 'MYC', 'NAB2', 'NRAS', 'RHOC', 'TGFB1', 'TIMP1')

testSigs <- function(exp, pheno, dataset){
  sequence <- c('ImSig.genes','IFNG.Sig', 'CD8.Sig','CRMA.Sig','IMPRES.Sig', 'IRG.Sig', 'LRRC15.CAF.Sig',
                'T.cell.inflamed.Sig', 'IPRES.Sig', 'Inflammatory.Sig', 'EMT.Sig', 'Blood.Sig', 'PDL1.Sig')
  #sequence <- c('ImSig.genes','random.Sig')
  all_data <- list()
  roc <- data.frame(AUC=c(), Dataset=c(), sig=c())
  for (var in sequence){
    print(var)
    sig <- get(var)
    test0 <- testCV(expdata = exp, phenoInfo = pheno, features = sig)[['result']]
    Data <- data.frame(grp=test0$label, res=test0$prediction)
    Data<-na.omit(Data)
    #assign(paste0(var,'.data'), Data)
    tmp <- list(var=Data)
    names(tmp) <- var
    all_data <- append(all_data, tmp)
    auc_value <- as.numeric(rocdata(Data$grp, Data$res)$stats$auc)
    roc <- rbind(roc, data.frame(AUC=auc_value, Dataset=dataset, sig=var))
  }
  return(list(rocdata=all_data,roc=roc))
}

all_data1 <- testSigs(GSE78220_AltAnalyze, GSE78220_PhenoInfo2,'GSE78220')[['roc']]
#p1 <- rocplot.multiple.V3(all_data, title = "GSE78220 data", p.value = FALSE)

all_data2 <- testSigs(BMS038.Pre.CountTable.normalized.log, BMS038_PhenoInfo, 'GSE91061')[['roc']]
#p2 <- rocplot.multiple.V3(all_data, title = "GSE91061 data", p.value = FALSE)

all_data3 <- testSigs(CC_73samples_GE_matrix, CC_73samples_pData,'PRJEB23709')[['roc']]
#p3 <- rocplot.multiple.V3(all_data, title = "PRJEB23709 data", p.value = FALSE)

all_data4 <- testSigs(NatMed_103samples_GE_matrix, NatMed_103samples_pData,'MGSP')[['roc']]
#p4 <- rocplot.multiple.V3(all_data, title = "MGSP data", p.value = FALSE)

result <- rbind(all_data1, all_data2, all_data3, all_data4)

#plot_grid(p1, p2, p3 ,p4, labels = c("a", "b", "c", "d"),ncol = 2)

library(ggpubr)
colnames(result)[3] <- 'Signature'
result$Signature <- as.character(result$Signature)
result$Signature <- ifelse(result$Signature=='ImSig.genes','ImmuneCells.Sig',result$Signature)
ordersig <- sort(tapply(result$AUC,result[,c("Signature")],median),decreasing = TRUE)
names(ordersig)


my_comparisons <- list(c('ImmuneCells.Sig','IMPRES.Sig'),c('ImmuneCells.Sig','IFNG.Sig'),c('ImmuneCells.Sig','IPRES.Sig'),c('ImmuneCells.Sig','CRMA.Sig'),c('ImmuneCells.Sig','IRG.Sig'),c('ImmuneCells.Sig','EMT.Sig'),
c('ImmuneCells.Sig','LRRC15.CAF.Sig'),c('ImmuneCells.Sig', 'CD8.Sig'),c('ImmuneCells.Sig','Blood.Sig'),c('ImmuneCells.Sig','PDL1.Sig'),c('ImmuneCells.Sig','Inflammatory.Sig'),c('ImmuneCells.Sig','T.cell.inflamed.Sig'))


tiff("Boxplot_signatures_compare_withP.tiff", width = 11, height = 10, units = "in", res = 800)
 
ggboxplot(result, x="Signature", y="AUC", color = "Signature",
          add = "jitter", short.panel.labs = FALSE) + 
  scale_x_discrete(limits=names(ordersig)) + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5)) + 
  stat_compare_means(label = "p.format",comparisons = my_comparisons, method = 't.test',method.args = list(alternative = "greater"),paired = TRUE)
  
dev.off()

 

tiff("Boxplot_signatures_compare_noP.tiff", width = 11, height = 10, units = "in", res = 800)
 
ggboxplot(result, x="Signature", y="AUC", color = "Signature",
          add = "jitter", short.panel.labs = FALSE) + 
  scale_x_discrete(limits=names(ordersig)) + theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5)) 
  
dev.off()

 
############################## PCA analysis:
common_genes <- Reduce(intersect, list(GSE78220_AltAnalyze$Symbol, CC_73samples_GE_matrix$Symbol, BMS038.Pre.CountTable.normalized.log$Symbol, 
                                       NatMed_103samples_GE_matrix$Symbol))
common_genes <- intersect(common_genes, ImSig.genes)
getCommon <- function(df, common){
  df <- df[df$Symbol %in% common_genes,]
  rownames(df) <- df$Symbol
  df <- df[,-which(colnames(df)=='Symbol')]
}
GSE78220 <- getCommon(GSE78220_AltAnalyze, common_genes)
PRJEB23709 <- getCommon(CC_73samples_GE_matrix, common_genes)
GSE91061 <- getCommon(BMS038.Pre.CountTable.normalized.log, common_genes)
MGSP <- getCommon(NatMed_103samples_GE_matrix, common_genes)

merged <- cbind(GSE78220, PRJEB23709, GSE91061, MGSP)
group_list <- c(rep('GSE78220',dim(GSE78220)[2]), rep('PRJEB23709', dim(PRJEB23709)[2]),
                rep('GSE91061',dim(GSE91061)[2]), rep('MGSP', dim(MGSP)[2]))
dat.pca <- PCA(t(merged), graph = FALSE, scale.unit = TRUE) #'-'表示“非”
g<-fviz_pca_ind(dat.pca,repel =T,
                geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
                col.ind =  group_list, # color by groups 颜色组
                # palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, # Concentration ellipses 集中成椭圆
                legend.title = "Dataset",
)
g
