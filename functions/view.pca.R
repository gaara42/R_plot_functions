#biafra ahanonu
#2013.01.09
#PCA functions

#Required packages
library(calibrate)
#For score(), to get PCA points
library(vegan)
#For subplots
library(Hmisc)

pcaplot <- function(data,filesep=",",headers=TRUE,rownames=1,analysisRows="ALL",correlation=TRUE,categorize=FALSE,categories=c(),categoryColumn="",labelsize=1.2,axisthickness=2,title="",biplotGraph=FALSE,showlabels=FALSE,clusters=10,subtitle="",sublist="none",pointsize=2) {
	#Data must have headers and rownames!

	# data = read.table(data,sep=filesep,header=FALSE)
	# df = as.data.frame(t(data),row.names=t(data[1,]))
	# names(df) = row.names(t(data))

	#If given a link to a file, get the data, assume format is already correct
	if(class(data)=="character"){
		if(rownames==FALSE){data = read.table(data,sep=filesep,header=headers)
		}else{data = read.table(data,sep=filesep,header=headers,row.names=rownames)
		}
	}

	# If showlabels=TRUE, use this rownames of the data.frame
	# labels = paste("+",row.names(data[,analysisRows]))
	if(analysisRows=="ALL"){
		labels = row.names(data)
	}else{
		labels = row.names(data[,analysisRows])
	}

	#Remove NaNs
	data = cleanPCAInput(data)
	if(analysisRows=="ALL"){
		pca = princomp(data,cor=correlation,na.rm=TRUE)
	}else{
		pca = princomp(data[,analysisRows],cor=correlation,na.rm=TRUE)
	}
	#Sites gives you the x,y coordinates of the pca analysis
	pca.species = scores(pca,display="species")
	pca.sites = scores(pca,display="sites")
	#Get names of rows...
	dataRowNames = row.names(pca.sites)

	#Use largest magnitude values, better looking plot
	pct = 0.2
	xmax = max(abs(pca.sites[,1]))+max(abs(pca.sites[,1]))*pct
	ymax = max(abs(pca.sites[,2]))+max(abs(pca.sites[,2]))*pct
	xlims = c(-xmax,xmax)
	ylims = c(-ymax,ymax)

	if(categorize==TRUE){
		plot(pca.sites,pch=row.names(pca.sites),bty="l",main=title,ylab="Component 2",xlab="Component 1",cex=labelsize,cex.lab=labelsize,cex.axis=labelsize)
		pca.points.categorize(pca.sites,categories,categoryColumn)
	}else{
		if(biplotGraph==FALSE){
			plotdata = pca.sites
			plotdata = ""
			plot(plotdata,pch=22,col=rgb(1,1,1),bg=rgb(1,0,0),bty="n",main=title,ylab="Component 2",xlab="Component 1",cex=labelsize,cex.lab=labelsize,cex.axis=labelsize,xlim=xlims,ylim=ylims,xaxt="n",yaxt="n")
			#Thicken lines, plot subtitle
			mtext(subtitle)
			box(lwd=axisthickness,bty="l")
			axis(1,lwd=axisthickness)
			axis(2,lwd=axisthickness)
			# axis(2,lwd=2,bty="l")

			#Subplot of amount of variance captured in each component
			variance = pca$sdev^2/sum(pca$sdev^2)
			subplot(pca.plot.pvar(pca),x=xmax,y=ymax,size=c(1,1),hadj=1,vadj=1)

			#Plot points into clusters
			pca.points.cluster(pca.sites,clusters,pointsize)

			#Show labels
			if(showlabels==TRUE){textxy(pca.sites[,1],pca.sites[,2],labs=labels,cx=1)}	

			#Highlight points of interest
			icol = 3
			high = quantile(data[,icol],prob=.7)
			pointsO = pca.sites[data[,icol]>high,]
			points(pointsO,pch=22,col="black",lwd=2,bg=rgb(.5,.5,.5,.3),cex=pointsize)

			if(sublist=="none"){
			}else{
				filename = sublist
				subList = read.table(filename,header=FALSE)
				subList = as.vector(subList[,1])
				# filename = "../../data/processed/2013_01_09/yeastALL.all.TE.r6.exit45.v1"
				# data1 = read.table(filename,header=FALSE,sep=",")
				#Find all points in the data where subRow is in subList
				result = apply(data, 1, function(x) all(x[1] %in% subList)) 
				subData = pca.sites[result,]
				points(subData,pch=22,col="black",lwd=2,bg=rgb(.5,.5,.5,.3),cex=1.2)
			}
		}else{
			biplot(pca,xlabs=labels,main=title,col="blue",pch=22)
		}
	}
}
pca.plot.pvar <- function(pca){
	#Plot the variance captured by each PCA component
	pvar = pca$sdev^2/sum(pca$sdev^2)
	title = "PCA Analysis: %Var of Components"
	title = ""
	# "Variance (%)"
	# "Component"
	print(pvar)
	bardata = barplot(pvar,col=rainbow(length(pvar)),border=rgb(1,1,1),xaxt="n",yaxt="n",main=title,ylab="",xlab="",bg="white")
	axisthickness=1.2
	box(lwd=axisthickness,bty="l")
	axis(1,lwd=axisthickness,at=bardata,labels=c(1:length(pvar)))
	axis(2,lwd=axisthickness)
}
pca.points.categorize <- function(data,categories,categoryColumn){
	#Vector with colors
	color = rainbow(clusters)
	catNum = 1
	for (category in categories) {
		cdata = data[data[[categoryColumn]]==category,]
		cdata = cdata[,1:2]
		points(cdata,pch=22,col=rgb(1,1,1),bg=color[catNum])
		catNum = catNum + 1
	}
}
pca.points.cluster <- function(data,clusters,pointsize){
	#Only cluster first two components
	data = data[,1:2]
	#Get kmeans data
	algorithms = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
	pca.clusters = kmeans(x=data,centers=clusters,algorithm=algorithms[4])	
	#Vector with colors
	color = rainbow(clusters)
	#For each cluster, plot points in the cluster with different color
	for (cluster in c(1:clusters)){
		clusterData = data[pca.clusters$cluster==cluster,]
		points(clusterData,pch=22,col=rgb(1,1,1),bg=color[cluster],cex=pointsize)
	}
}
cleanPCAInput <- function(data){
	#Replace NaNs with column mean, i.e. non-informative point in PCA
	columns = c(1:ncol(data))
	for (column in columns) {
		cdata = data[,column]
		#Skip columns that aren't numeric/integer
		if(class(cdata[1])=="factor"){next}
		cmean=mean(cdata,na.rm=TRUE)
		cdata[is.na(cdata)]=cmean
		#Set data to mean=0, var=1
		cdata = (cdata - cmean)/sd(cdata,na.rm=TRUE)
		data[,column]=cdata
	}
	return(data)
}