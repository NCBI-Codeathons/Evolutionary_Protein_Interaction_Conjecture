
Matrix1 <- read.csv("mat1.csv", header = F)
colnames(Matrix1) <- c("nrdA", "nrdB", "folA", "thyA", "cydA", "cydB", "fliG", "fliM", "fliN", "trpA", "trpB", "purE", "purK", "purC", 
                       "ilvC", "ilvB","glgA", "glgC")
rownames(Matrix1) <- c("nrdA", "nrdB", "folA", "thyA", "cydA", "cydB", "fliG", "fliM", "fliN", "trpA", "trpB", "purE", "purK", "purC", 
                       "ilvC", "ilvB","glgA", "glgC")

True_Matrix <- Matrix1
True_Matrix[,]=0

True_Matrix["cydB","cydA"] = 1
True_Matrix["cydA","cydB"] = 1
True_Matrix["nrdA","nrdB"] = 1
True_Matrix["nrdB","nrdA"] = 1
True_Matrix["fliG","fliM"] = 1
True_Matrix["fliM","fliG"] = 1
True_Matrix["fliG","fliN"] = 1
True_Matrix["fliN","fliG"] = 1
True_Matrix["fliM","fliN"] = 1
True_Matrix["fliN","fliM"] = 1
True_Matrix["trpA","trpB"] = 1
True_Matrix["trpB","trpA"] = 1
True_Matrix["folA","thyA"] = 1
True_Matrix["thyA","folA"] = 1
True_Matrix["purE","purK"] = 1
True_Matrix["purK","purE"] = 1
True_Matrix["purC","purE"] = 1
True_Matrix["purE","purC"] = 1
True_Matrix["purK","purC"] = 1
True_Matrix["purC","purK"] = 1
True_Matrix["ilvC","ilvB"] = 1
True_Matrix["ilvB","ilvC"] = 1
True_Matrix["glgA","glgC"] = 1
True_Matrix["glgC","glgA"] = 1


library(heatmaply)
True_Matrix
heatmaply(True_Matrix,Rowv = FALSE, Colv = F)


Matrix3 <- as.data.frame(ifelse(Matrix1 > 0.50, 1,0))
heatmaply(Matrix3,Rowv = FALSE, Colv = F)


CM <- Confusion_Matrix(Matrix2,True_Matrix)

Confusion_Matrix <- function(pred, true){
  
  True_Positive = 0
  False_Positive= 0
  True_Negative=0
  False_Negative=0
  
  for (j in 1:length(rownames(pred))){
    print(j)
    for(i in 1:length(colnames(pred))){
      
      if(pred[i,j]== 1){
        if(true[i,j]==1){
          True_Positive = True_Positive + 1
        } else {
          False_Positive = False_Positive + 1
        }
      }
      if(pred[i,j]== 0){
        if(true[i,j]==0){
          True_Negative = True_Negative + 1
        } else {
          False_Negative = False_Negative + 1
        }
      }
    }
  }
  CM <- data.frame(True_Positive,True_Negative, False_Positive, False_Negative)
  return(CM)
}

Threshold <- c(-.5,-.4,-.3,-.2,-.1,.1,.2,.3,.4,.5)
for(i in Threshold){
  Pred <- as.data.frame(ifelse(Matrix1 > i, 1,0))
  Conf_Matrix <- Confusion_Matrix(Pred, True_Matrix)
  Conf_Matrix <- cbind(Conf_Matrix,i)
  if(!exists("ROC")){
    ROC <- Conf_Matrix
  } else {
    ROC <- rbind(ROC,Conf_Matrix)
  }
}

ggplot(ROC, aes(x=i, color))

#################################################
a <- Matrix1[upper.tri(Matrix1)]
a
b <- True_Matrix[upper.tri(True_Matrix)]
pred_vector <- Matrix1[upper.tri(Matrix1)]
true_vector <- True_Matrix[upper.tri(True_Matrix)]

library(ROCR)
 
Performance <- performance(prediction(pred_vector,true_vector), measure = "tpr", x.measure = "fpr", "auc")
Performance_AUC <- performance(prediction(pred_vector,true_vector2), measure = "auc")
Performance_AUC <-  Performance_AUC@y.values[[1]]
Performance_AUC

plot(Performance)
####################################################
Matrix2 <- read.csv("mat2.csv", header = F)
colnames(Matrix2) <- c("nrdA", "nrdB", "folA", "thyA", "cydA", "cydB", "fliG", "fliM", "fliN", "trpA", "trpB", "purE", "purK", "purC", 
                       "ilvC", "ilvB","glgA", "glgC")
rownames(Matrix2) <- c("nrdA", "nrdB", "folA", "thyA", "cydA", "cydB", "fliG", "fliM", "fliN", "trpA", "trpB", "purE", "purK", "purC", 
                       "ilvC", "ilvB","glgA", "glgC")

pred_vector2 <- Matrix2[upper.tri(Matrix2)]
true_vector2 <- True_Matrix[upper.tri(True_Matrix)]

library(ROCR)

Performance2 <- performance(prediction(pred_vector2,true_vector2), measure = "tpr", x.measure = "fpr")
Performance2_AUC <- performance(prediction(pred_vector2,true_vector2), measure = "auc")
Performance2_AUC <-  Performance2_AUC@y.values[[1]]
Performance2_AUC

plot(Performance2)

plot(Performance, col="red") 
lines(Performance2@x.values[[1]], Performance2@y.values[[1]], col = "blue", lty=2)
legend(.2,.2, legend = c("Test Protein","Housekeeping Proteins"), col = c("red", "blue"), lty = 1:2, cex = 0.8)

?performance

#################################################
library(reshape2)
library(ggplot)
MG <- melt(ROC[,-c(2,3)], id ="i")
ggplot(MG, aes(x=i, y=value, color = variable))+geom_line()

################################################
#ROC Pipeline

ROC_True_Matrix <- True_Matrix
ROC_True_Matrix[] = 0 

colnames(ROC_True_Matrix) <- c("cydA", "cydB", "fliG", "fliM", "fliN", "folA", "glgA", "glgC", "ilvB" ,"ilvC", 
                               "nrdA", "nrdB", "purC", "purE", "purK","thyA", "trpA", "trpB")


rownames(ROC_True_Matrix) <- c("cydA", "cydB", "fliG", "fliM", "fliN", "folA", "glgA", "glgC", "ilvB" ,"ilvC", 
                               "nrdA", "nrdB", "purC", "purE", "purK","thyA", "trpA", "trpB")

ROC_True_Matrix["cydB","cydA"] = 1
ROC_True_Matrix["cydA","cydB"] = 1
ROC_True_Matrix["nrdA","nrdB"] = 1
ROC_True_Matrix["nrdB","nrdA"] = 1
ROC_True_Matrix["fliG","fliM"] = 1
ROC_True_Matrix["fliM","fliG"] = 1
ROC_True_Matrix["fliG","fliN"] = 1
ROC_True_Matrix["fliN","fliG"] = 1
ROC_True_Matrix["fliM","fliN"] = 1
ROC_True_Matrix["fliN","fliM"] = 1
ROC_True_Matrix["trpA","trpB"] = 1
ROC_True_Matrix["trpB","trpA"] = 1
ROC_True_Matrix["folA","thyA"] = 1
ROC_True_Matrix["thyA","folA"] = 1
ROC_True_Matrix["purE","purK"] = 1
ROC_True_Matrix["purK","purE"] = 1
ROC_True_Matrix["purC","purE"] = 1
ROC_True_Matrix["purE","purC"] = 1
ROC_True_Matrix["purK","purC"] = 1
ROC_True_Matrix["purC","purK"] = 1
ROC_True_Matrix["ilvC","ilvB"] = 1
ROC_True_Matrix["ilvB","ilvC"] = 1
ROC_True_Matrix["glgA","glgC"] = 1
ROC_True_Matrix["glgC","glgA"] = 1

heatmaply(ROC_True_Matrix,Rowv = FALSE, Colv = F)

#########################################################


Matrix_plot <- function(x){
  vec <- unlist(x)
  vec <- unname(vec)
  # vec_ord <- vec[order(vec, decreasing = T)]
  # 
  # data=data.frame(matrix(ncol=2,nrow=length(vec_ord))) 
  # 
  # for (i in 1:length(vec_ord)){
  #   data[i,1]<- i/length(vec_ord)
  #   if(i==1){
  #     data[i,2]<- vec_ord[i]/sum(vec_ord)
  #   } else{
  #     data[i,2]<- (vec_ord[i]+vec_ord[i-1])/sum(vec_ord)
  #   }
  #   
  # }
  return(vec)
}

MT_nerd <- read.csv("MT_nrd_comp.csv")
SCA_nerd <- read.csv("SCA_nrd_comp.csv")
MT_fol <- read.csv("MT_fol_comp (1).csv")
SCA_fol <- read.csv("SCA_fol_comp.csv")
MT_CT <- read.csv("MT_fol_nrd_comp.csv")
SCA_CT <- read.csv("SCA_fol_nrd_comp.csv")

MT_nerd_data <- Matrix_plot(MT_nerd)
SCA_nrd_data <- Matrix_plot(SCA_nerd)
MT_fol_data <- Matrix_plot(MT_fol)
SCA_fol_data <- Matrix_plot(SCA_fol)
MT_CT_data <- Matrix_plot(MT_CT)
SCA_CT_data <- Matrix_plot(SCA_CT)

Dataframe <- qpcR:::cbind.na(MT_nerd_data, MT_fol_data, MT_CT_data, SCA_nrd_data, SCA_fol_data, SCA_CT_data)

library(reshape2)

melted_df <- na.omit(melt(Dataframe))


no_na <- na.omit(melted_df)

melted_df$id <- ifelse(melted_df$Var2 == "MT_nerd_data" | melted_df$Var2 == "SCA_nrd_data", "nrdA/nrdB", 
                       ifelse(melted_df$Var2 == "MT_fol_data" | melted_df$Var2 == "SCA_fol_data", "folA/thyA",
                              ifelse(melted_df$Var2 == "MT_CT_data" | melted_df$Var2 == "SCA_CT_data" , "folA/nrdB",NA)))

levels(factor(melted_df$id))

lbl = c("Mirror Tree", "Mirror Tree", "Mirror Tree", "SCA", "SCA", "SCA")

ggplot(melted_df, aes(x=Var2,y=value, fill = id))+geom_boxplot()+ ylab("Positional Interaction Score") + xlab("Model")+
  theme(axis.text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1), axis.title=element_text(size=20, face="bold"), legend.title = element_text(size=20, face="bold") ,legend.text=element_text(size=20), 
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5))+ ggtitle("Model Data Boxplot")+ scale_x_discrete(labels= lbl)

  
#############################################
library(RcppCNPy)
a <- npyLoad("SCA_Delta_0.2_Similarity_0.0.npy" , type = "interger")
dir()

