
##read the SCA and MT score data from specific directory

setwd("/project/hackathon/hackers11/shared/team3/xuexin")
options(stringsAsFactors=F)
sca<-read.table("sca_result.csv",header=F,sep=",")
sca<-as.matrix(sca)
mt<-read.table("mt_result.csv",header=F,sep=",")
mt<-as.matrix(mt)

#import package for boxplot
library(ggplot2)

group.name=c("sameProt","physInt","enzInt","negatives")

##construct the dataframe of SCA score
sca.score=c()
sca.label=c()

for(i in 1:nrow(sca)){
	dd=as.vector(sca[i,])
	sca.score=c(sca.score,na.omit(sca[i,]))
	tmp=rep(group.name[i],length(na.omit(sca[i,])))
	sca.label=c(sca.label,tmp)
}


sca.dataframe=data.frame(score=sca.score,Groups=sca.label)

library(ggplot2)

pdf("sca_boxplot.pdf")
ggplot(sca.dataframe, aes(x=Groups, y=score,fill=Groups)) + geom_boxplot()+geom_point(position = position_jitterdodge()) +theme(text = element_text(size=17), axis.text.x = element_text(angle=90, hjust=1))+ylab("SCA Score")+xlab("")
ggplot(sca.dataframe, aes(x=Groups, y=score,fill=Groups)) + geom_boxplot()+geom_point(position = position_jitterdodge()) +theme(text = element_text(size=17))+ylab("SCA Score")+xlab("")

dev.off()

##construct the dataframe of MT score
mt.score=c()
mt.label=c()

for(i in 1:nrow(mt)){
	dd=as.vector(mt[i,])
	mt.score=c(mt.score,na.omit(mt[i,]))
	tmp=rep(group.name[i],length(na.omit(mt[i,])))
	mt.label=c(mt.label,tmp)
}


mt.dataframe=data.frame(score=mt.score,Groups=mt.label)

pdf("mt_boxplot.pdf")
ggplot(mt.dataframe, aes(x=Groups, y=score,fill=Groups)) + geom_boxplot()+geom_point(position = position_jitterdodge()) +theme(text = element_text(size=17), axis.text.x = element_text(angle=90, hjust=1))+ylab("MT Score")+xlab("")
dev.off()


