#biafra ahanonu
#2013.01.09
#PCA functions

#Load PCA functions
source("../view/view.pca.R")
#Set title
titles = "USA Crime: Murder, Assault, Urban Population and Rape"
#Subtitle
subt = "Inset shows % variance captured by each component, kmeans clusters = 10, black box = 70th percentile Urban Population"
#PCA plot, see function for details
#Uncomment postscript and dev.off() to output postscript files
# postscript(imglist[count],horizontal=TRUE,onefile=FALSE,paper="a4")
# png(imglist[count],width=2000,height=1000,res=200,pointsize=10,antialias = "cleartype")
pcaplot(data=USArrests,rownames=FALSE,headers=FALSE,title=titles,subtitle=subt ,showlabels=TRUE,pointsize=1.5)
# dev.off()