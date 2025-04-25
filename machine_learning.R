####### GSE7803 ########
#######  Subdata Analysis
rm(list = ls())
gc()
Datas=read.csv("GSE7803_Data.csv",header=T,row.names=1)
g=read.csv("GSE7803_group.csv",header=T,row.names=1)
Datav=Datas[(rownames(Datas) %in% c("201291_s_at", "202870_s_at","203755_at","204126_s_at","204444_at","204531_s_at","204709_s_at","204767_s_at","208795_s_at","210559_s_at")),]
rownames(Datav)=c("TOP2A","CDC20","BUB1B","CDC45","KIF11","BRCA1","KIF23","FEN1","MCM7","CDK1")
dim(Datav)
p=10
Tr.group=ifelse(g$Group%in%c('control'), 1,2)
Data_ML=cbind(t(Datav),Tr.group)
TrainData=as.data.frame(Data_ML)
dim(TrainData)

#######  SVM Analysis

#install.packages('caret', repos='http://cran.rstudio.com/')
library(caret)
library(e1071)
SVM_fit=svm(Tr.group~ ., data=TrainData,cost = 1,kernel = "radial",type="C-classification")
prediction=predict(SVM_fit,TrainData[,-(p+1)],type="prob")
accuracy_sv=table(prediction,TrainData$Tr.group)
accuracy_result_sv=confusionMatrix(accuracy_sv)
accuracy_result_sv

#######  Random Forest Analysis

#install.packages("randomForestSRC",repos = "http://cran.us.r-project.org")
#install.packages("randomForest")
library(randomForestSRC)
library(randomForest)
library(caret)

model.Forest<- randomForest(as.factor(Tr.group) ~ ., data=TrainData, importance=TRUE,
                            proximity=TRUE)
Forest.predTr<-predict(model.Forest,TrainData[,-(p+1)])
accuracy_rf=table(Forest.predTr,TrainData$Tr.group)
accuracy_result_rf=confusionMatrix(accuracy_rf)
accuracy_result_rf

###############  KNN Analysis
#install.packages("kknn")
library(kknn)
model.knn<- kknn(as.factor(Tr.group)~.,train=TrainData,test=TrainData)
accuracy_knn<- mean(model.knn$fit==Tr.group)
accuracy_knn

#################LOOCV##
acc<-rep(0,nrow(Data_ML))
for(i in 1:nrow(Data_ML)){
  train.data=Data_ML[-i,]
  test.data=as.vector(Data_ML[i,-(1+p)])
  SVM_fit<-svm(as.factor(Tr.group)~ ., data=train.data,cost = 1,kernel = "radial",type="C-classification")
  pred<-predict(SVM_fit,train.data)
  acc[i]<-mean(pred==Tr.group)
}
acc

#############  LOOCV K-NN  
predicted_KNN <- NULL
for(i in 1:nrow(Data_ML)){
  training<-as.data.frame(Data_ML[-i,])
  test<-data.frame(t(Data_ML[i,-(p+1)]))
  model1_KNN<-kknn(as.factor(Tr.group)~.,train=training,test=test)
  predicted_KNN[i]<- predict(model1_KNN, newdata = test)
}
AC_KNN<-mean(predicted_KNN==training[,p+1])
AC_KNN

########    LOOCV SVM    
predicted_svm <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1<-svm(as.factor(Tr.group)~ ., data=training)
  predicted_svm[i]<- predict(model1, newdata = test)
  #AC_svm[i]<-predicted_svm==training[i,p+1]
}
AC_svm<-mean(predicted_svm==training[,p+1])
AC_svm


########  LOOCV RF   
predicted_rf <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1_rf<-randomForest(as.factor(Tr.group) ~ ., data=training, importance=TRUE,proximity=TRUE)
  predicted_rf[i]<- predict(model1_rf, newdata = test)
}
AC_rf<-mean(predicted_rf==training[,p+1])
AC_rf


####### GSE9750 ########
#######  Subdata Analysis
rm(list = ls())
gc()
Datas=read.csv("GSE9750_Data.csv",header=T,row.names=1)
g=read.csv("GSE9750_group.csv",header=T,row.names=1)
Datav=Datas[(rownames(Datas) %in% c("201291_s_at", "202870_s_at","203755_at","204126_s_at","204444_at","204531_s_at","204709_s_at","204767_s_at","208795_s_at","210559_s_at")),]
rownames(Datav)=c("TOP2A","CDC20","BUB1B","CDC45","KIF11","BRCA1","KIF23","FEN1","MCM7","CDK1")
dim(Datav)
p=10
Tr.group=ifelse(g$Group%in%c('control'), 1,2)
Data_ML=cbind(t(Datav),Tr.group)
TrainData=as.data.frame(Data_ML)
dim(TrainData)

#######  SVM Analysis

library(caret)
library(e1071)
SVM_fit=svm(Tr.group~ ., data=TrainData,cost = 1,kernel = "radial",type="C-classification")
prediction=predict(SVM_fit,TrainData[,-(p+1)],type="prob")
accuracy_sv=table(prediction,TrainData$Tr.group)
accuracy_result_sv=confusionMatrix(accuracy_sv)
accuracy_result_sv

#######  Random Forest Analysis

library(randomForestSRC)
library(randomForest)
library(caret)

model.Forest<- randomForest(as.factor(Tr.group) ~ ., data=TrainData, importance=TRUE,
                            proximity=TRUE)
Forest.predTr<-predict(model.Forest,TrainData[,-(p+1)])
accuracy_rf=table(Forest.predTr,TrainData$Tr.group)
accuracy_result_rf=confusionMatrix(accuracy_rf)
accuracy_result_rf

###############  KNN Analysis
#install.packages("kknn")
library(kknn)
model.knn<- kknn(as.factor(Tr.group)~.,train=TrainData,test=TrainData)
accuracy_knn<- mean(model.knn$fit==Tr.group)
accuracy_knn

#################LOOCV##
acc<-rep(0,nrow(Data_ML))
for(i in 1:nrow(Data_ML)){
  train.data=Data_ML[-i,]
  test.data=as.vector(Data_ML[i,-(1+p)])
  SVM_fit<-svm(as.factor(Tr.group)~ ., data=train.data,cost = 1,kernel = "radial",type="C-classification")
  pred<-predict(SVM_fit,train.data)
  acc[i]<-mean(pred==Tr.group)
}
acc

#############  LOOCV K-NN  
predicted_KNN <- NULL
for(i in 1:nrow(Data_ML)){
  training<-as.data.frame(Data_ML[-i,])
  test<-data.frame(t(Data_ML[i,-(p+1)]))
  model1_KNN<-kknn(as.factor(Tr.group)~.,train=training,test=test)
  predicted_KNN[i]<- predict(model1_KNN, newdata = test)
}
AC_KNN<-mean(predicted_KNN==training[,p+1])
AC_KNN

########    LOOCV SVM    
predicted_svm <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1<-svm(as.factor(Tr.group)~ ., data=training)
  predicted_svm[i]<- predict(model1, newdata = test)
  #AC_svm[i]<-predicted_svm==training[i,p+1]
}
AC_svm<-mean(predicted_svm==training[,p+1])
AC_svm


########  LOOCV RF   
predicted_rf <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1_rf<-randomForest(as.factor(Tr.group) ~ ., data=training, importance=TRUE,proximity=TRUE)
  predicted_rf[i]<- predict(model1_rf, newdata = test)
}
AC_rf<-mean(predicted_rf==training[,p+1])
AC_rf


####### GSE168244 ########
#######  Subdata Analysis
rm(list = ls())
gc()
Datas=read.csv("GSE168244_Data.csv",header=T,row.names=1)
Datav=Datas[(rownames(Datas) %in% c("TOP2A","CDC20","BUB1B","CDC45","KIF11","BRCA1","KIF23","FEN1","MCM7","CDK1")),]
dim(Datav)
p=10
n1=18
n2=17
Tr.group <- c(rep(1,n1), rep(2,n1))
Data_ML=cbind(t(Datav),Tr.group)
TrainData=as.data.frame(Data_ML)
dim(TrainData)

#######  SVM Analysis

library(caret)
library(e1071)
SVM_fit=svm(Tr.group~ ., data=TrainData,cost = 1,kernel = "radial",type="C-classification")
prediction=predict(SVM_fit,TrainData[,-(p+1)],type="prob")
accuracy_sv=table(prediction,TrainData$Tr.group)
accuracy_result_sv=confusionMatrix(accuracy_sv)
accuracy_result_sv

#######  Random Forest Analysis

library(randomForestSRC)
library(randomForest)
library(caret)

model.Forest<- randomForest(as.factor(Tr.group) ~ ., data=TrainData, importance=TRUE,
                            proximity=TRUE)
Forest.predTr<-predict(model.Forest,TrainData[,-(p+1)])
accuracy_rf=table(Forest.predTr,TrainData$Tr.group)
accuracy_result_rf=confusionMatrix(accuracy_rf)
accuracy_result_rf

###############  KNN Analysis
#install.packages("kknn")
library(kknn)
model.knn<- kknn(as.factor(Tr.group)~.,train=TrainData,test=TrainData)
accuracy_knn<- mean(model.knn$fit==Tr.group)
accuracy_knn

#################LOOCV##
acc<-rep(0,nrow(Data_ML))
for(i in 1:nrow(Data_ML)){
  train.data=Data_ML[-i,]
  test.data=as.vector(Data_ML[i,-(1+p)])
  SVM_fit<-svm(as.factor(Tr.group)~ ., data=train.data,cost = 1,kernel = "radial",type="C-classification")
  pred<-predict(SVM_fit,train.data)
  acc[i]<-mean(pred==Tr.group)
}
acc

#############  LOOCV K-NN  
predicted_KNN <- NULL
for(i in 1:nrow(Data_ML)){
  training<-as.data.frame(Data_ML[-i,])
  test<-data.frame(t(Data_ML[i,-(p+1)]))
  model1_KNN<-kknn(as.factor(Tr.group)~.,train=training,test=test)
  predicted_KNN[i]<- predict(model1_KNN, newdata = test)
}
AC_KNN<-mean(predicted_KNN==training[,p+1])
AC_KNN

########    LOOCV SVM    
predicted_svm <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1<-svm(as.factor(Tr.group)~ ., data=training)
  predicted_svm[i]<- predict(model1, newdata = test)
  #AC_svm[i]<-predicted_svm==training[i,p+1]
}
AC_svm<-mean(predicted_svm==training[,p+1])
AC_svm


########  LOOCV RF   
predicted_rf <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1_rf<-randomForest(as.factor(Tr.group) ~ ., data=training, importance=TRUE,proximity=TRUE)
  predicted_rf[i]<- predict(model1_rf, newdata = test)
}
AC_rf<-mean(predicted_rf==training[,p+1])
AC_rf


####### GSE149450 ########
#######  Subdata Analysis
rm(list = ls())
gc()
Datas=read.csv("GSE149450_Data.csv",header=T,row.names=1)
Datav=Datas[(rownames(Datas) %in% c("ENSG00000131747.14","ENSG00000117399.13",
                                    "ENSG00000156970.12","ENSG00000093009.9",
                                    "ENSG00000138160.5","ENSG00000012048.19",
                                    "ENSG00000137807.13","ENSG00000168496.3",
                                    "ENSG00000166508.17","ENSG00000170312.15")),]
rownames(Datav)=c("KIF11","KIF23","TOP2A","FEN1","BUB1B","BRCA1","MCM7","CDK1","CDC45","CDC20")
dim(Datav)
p=10
n1=2
n2=2
Tr.group <- c(rep(1,n1), rep(2,n1))
Data_ML=cbind(t(Datav),Tr.group)
TrainData=as.data.frame(Data_ML)
dim(TrainData)

#######  SVM Analysis

library(caret)
library(e1071)
SVM_fit=svm(Tr.group~ ., data=TrainData,cost = 1,kernel = "radial",type="C-classification")
prediction=predict(SVM_fit,TrainData[,-(p+1)],type="prob")
accuracy_sv=table(prediction,TrainData$Tr.group)
accuracy_result_sv=confusionMatrix(accuracy_sv)
accuracy_result_sv

#######  Random Forest Analysis

library(randomForestSRC)
library(randomForest)
library(caret)

model.Forest<- randomForest(as.factor(Tr.group) ~ ., data=TrainData, importance=TRUE,
                            proximity=TRUE)
Forest.predTr<-predict(model.Forest,TrainData[,-(p+1)])
accuracy_rf=table(Forest.predTr,TrainData$Tr.group)
accuracy_result_rf=confusionMatrix(accuracy_rf)
accuracy_result_rf

###############  KNN Analysis
#install.packages("kknn")
library(kknn)
model.knn<- kknn(as.factor(Tr.group)~.,train=TrainData,test=TrainData)
accuracy_knn<- mean(model.knn$fit==Tr.group)
accuracy_knn

#################LOOCV##
acc<-rep(0,nrow(Data_ML))
for(i in 1:nrow(Data_ML)){
  train.data=Data_ML[-i,]
  test.data=as.vector(Data_ML[i,-(1+p)])
  SVM_fit<-svm(as.factor(Tr.group)~ ., data=train.data,cost = 1,kernel = "radial",type="C-classification")
  pred<-predict(SVM_fit,train.data)
  acc[i]<-mean(pred==Tr.group)
}
acc

#############  LOOCV K-NN  
predicted_KNN <- NULL
for(i in 1:nrow(Data_ML)){
  training<-as.data.frame(Data_ML[-i,])
  test<-data.frame(t(Data_ML[i,-(p+1)]))
  model1_KNN<-kknn(as.factor(Tr.group)~.,train=training,test=test)
  predicted_KNN[i]<- predict(model1_KNN, newdata = test)
}
AC_KNN<-mean(predicted_KNN==training[,p+1])
AC_KNN

########    LOOCV SVM    
predicted_svm <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1<-svm(as.factor(Tr.group)~ ., data=training)
  predicted_svm[i]<- predict(model1, newdata = test)
  #AC_svm[i]<-predicted_svm==training[i,p+1]
}
AC_svm<-mean(predicted_svm==training[,p+1])
AC_svm


########  LOOCV RF   
predicted_rf <- NULL
for(i in 1:nrow(Data_ML)){
  training<-Data_ML[-i,]
  test<-t(as.data.frame(Data_ML[i,-(p+1)]))
  model1_rf<-randomForest(as.factor(Tr.group) ~ ., data=training, importance=TRUE,proximity=TRUE)
  predicted_rf[i]<- predict(model1_rf, newdata = test)
}
AC_rf<-mean(predicted_rf==training[,p+1])
AC_rf
