#biafra ahanonu
#2013.01.09
#PCA functions

#Required packages
library(calibrate)
#For score(), to get PCA points
library(vegan)
#For subplots
library(Hmisc)

#Notes on inputs
# data
	#Input data, can either be a data.frame organized with variables in columns and rows containing datapoints or a string that references a properly organized file. No default, must supply something.
# filesep=","
	#Separating character used to read in data if data is a list.
# headers=TRUE
	#Decide if imported data has headers
# rownames=1
	#Make data.frame rows have rownames from the corresponding column in the data. Set to FALSE if you don't want rownames.
# analysisRows="ALL"
	#Either run PCA on all rows or give a character vector corresponding to columns that should be analyzed. e.g. c("Homes","Cars","FamilyMembers")
# correlation=TRUE
	#Should princomp use correlations
# categorize=FALSE
	#Color code data points based on which category they belong to of a particular categoryColumn
# categories=c()
	#Character vector, Give the categories to separate on
# categoryColumn=""
	#Character string, data column to categorize
# labelsize=1.2
	#Label size for points on the plot
# axisthickness=2
	#Thickness of the axis, default is slightly larger to make it look nice.
# title=""
	#Character string, title of the graph
# subtitle=""
	#Character string, subtitle of the graph
# footnote=format(Sys.time(), "%d %b %Y")
	#Character string, gray footnote at the bottom of the graph, always will have the date even if you change the default.
# biplotGraph=FALSE
	#Use biplot to show data
# showlabels=FALSE
	#Show row label identifiers next to each data point, set to TRUE if want labels
# sublist=FALSE
	#Character string, give the file path to a sublist. If a string is provided, compares data in file to sublistColumn data. File should have list organized as a single column.
# sublistColumn=1
	#Column to use for comparison against the sublist
# pointsize=2
	#Size of points on the graph
# quantileColumn=FALSE
	#Character string, if column given, color codes data points by quantile
# kmeans.show=TRUE
	#Show kmean clustering of data points
# clusters=10
	#Number of kmeans clusters to have

pcaplot <- function(data,filesep=",",headers=TRUE,rownames=1,analysisRows="ALL",correlation=TRUE,categorize=FALSE,categories=c(),categoryColumn="",labelsize=1.2,axisthickness=2,title="",biplotGraph=FALSE,showlabels=FALSE,clusters=10,subtitle="",sublist=FALSE,pointsize=2,footnote=format(Sys.time(), "%d %b %Y"),quantileColumn=FALSE,kmeans.show=TRUE,sublistColumn=1) {

	#If given a link to a file, get the data, assume format is already correct
	if(class(data)=="character"){
		if(rownames==FALSE){data = read.table(data,sep=filesep,header=headers)
		}else{data = read.table(data,sep=filesep,header=headers,row.names=rownames)
		}
	}

	# If showlabels=TRUE, use this rownames of the data.frame
	if(analysisRows=="ALL"){
		labels = row.names(data)
	}else{
		labels = row.names(data[,analysisRows])
	}

	#Make a copy of the original data for later use (e.g. with boxplots)
	datacp = data
	#NaNs are set to column mean, i.e. non-informative point in PCA
	data = cleanPCAInput(data)

	#Decide whether to analyze all rows
	if(analysisRows=="ALL"){
		pca = princomp(data,cor=correlation,na.rm=TRUE)
	}else{
		#Is a integer or character vector given for analysisRows
		if(class(analysisRows)=="character"){
			data_pca = data[names(data) %in% analysisRows]
			# data_pca = cleanPCAInput(data_pca)
			print(summary(data_pca))
			pca = princomp(data_pca,cor=correlation,na.rm=TRUE)
		}else{
			pca = princomp(data[,analysisRows],cor=correlation,na.rm=TRUE)
		}
	}

	#Sites gives you the x,y coordinates of PCA
	pca.species = scores(pca,display="species")
	pca.sites = scores(pca,display="sites")
	#Get PCA row names
	dataRowNames = row.names(pca.sites)

	#Use largest magnitude values and add a pct, better looking plot
	pct = 0.2
	xmax = max(abs(pca.sites[,1]))+max(abs(pca.sites[,1]))*pct
	ymax = max(abs(pca.sites[,2]))+max(abs(pca.sites[,2]))*pct
	xlims = c(-xmax,xmax)
	ylims = c(-ymax,ymax)

	if(categorize==TRUE){
		#Plot function as normal
		plot(pca.sites,pch=row.names(pca.sites),bty="l",main=title,ylab="Component 2",xlab="Component 1",cex=labelsize,cex.lab=labelsize,cex.axis=labelsize)
		pca.points.categorize(pca.sites,categories,categoryColumn)
	}else{
		if(biplotGraph==FALSE){
			#Only make plot if not showing clusters to avoid redundant plotting
			if(kmeans.show==TRUE){
				plotdata = ""
			}else{
				plotdata = pca.sites
			}

			plot(plotdata,pch=22,col=rgb(1,1,1),bg=rgb(1,0,0),bty="n",main=title,ylab="Component 2",xlab="Component 1",cex=labelsize,cex.lab=labelsize,cex.axis=labelsize,xlim=xlims,ylim=ylims,xaxt="n",yaxt="n")

			#Thicken lines, plot subtitle
			makeFootnote(footnote)
			mtext(subtitle)
			box(lwd=axisthickness,bty="l")
			axis(1,lwd=axisthickness)
			axis(2,lwd=axisthickness)
			# axis(2,lwd=2,bty="l")

			#Subplot of amount of variance captured in each component
			variance = pca$sdev^2/sum(pca$sdev^2)
			subplot(pca.plot.pvar(pca),x=xmax,y=ymax,size=c(1,1),hadj=1,vadj=1)

			#Plot points into clusters
			if(kmeans.show==TRUE){
				pca.points.cluster(pca.sites,clusters,pointsize)
			}

			#Show labels
			if(showlabels==TRUE){
				textxy(pca.sites[,1],pca.sites[,2],labs=labels,cx=1)
			}

			#Highlight points of interest based on quantile grouping
			if(quantileColumn!=FALSE&&class(quantileColumn)=="character"){
				highlight.quantile(quantileColumn,data,datacp,pca.sites,pointsize,pct,xmax,ymax)
			}

			#Highlight set of proteins from a list
			if(sublist!=FALSE&&class(sublist)=="character"){
				highlight.sublist(sublist,data,pca.sites,pointsize,pct,sublistColumn)
			}
		}else{
			biplot(pca,xlabs=labels,main=title,col="blue",pch=22)
		}
	}
}

highlight.sublist <- function(sublist,data,pca.sites,pointsize,pct,sublistColumn){
	#Highlight a subset of the points based on values from a list
	filename = sublist
	subList = read.table(filename,header=FALSE)
	subList = as.vector(subList[,1])
	#Find all points in the data where subRow is in subList
	result = apply(data, 1, function(x) all(x[sublistColumn] %in% subList)) 
	subData = pca.sites[result,]
	points(subData,pch=22,col="black",lwd=2,bg=rgb(.5,.5,.5,.3),cex=pointsize)
}

highlight.quantile <- function(separatingCol,data,datacp,pca.sites,pointsize,pct,xmax,ymax){
	#Highlight points of interest
	icol = separatingCol
	#Get quantile data points
	dataQuantiles = quantile(data[[icol]],na.rm=TRUE)
	colors = rainbow(length(dataQuantiles))
	boxdata = list()
	for(q in c(2:length(dataQuantiles))){
		high = dataQuantiles[q]
		low = dataQuantiles[q-1]
		#Get idx for values inside the current quantile
		pointsO = pca.sites[data[[icol]]<high&data[[icol]]>low,]
		points(pointsO,pch=22,col=colors[q],lwd=2,cex=pointsize)
		tempData = as.vector(data[[icol]])
		result = tempData<high&tempData>low
		dataQuant = datacp[result,]	
		#Store in boxdata for use in subplot
		boxdata[[q]] = dataQuant[[icol]]
	}

	#Plot barplot showing distribution of highlighted points
	subplot(subplot.box(boxdata[2:length(boxdata)],icol),x=xmax,y=-ymax+ymax*pct,size=c(1,1),hadj=1,vadj=0)
}

makeFootnote <- function(footnoteText=format(Sys.time(), "%d %b %Y"),size= .7, color= grey(.5)){
   require(grid)
   pushViewport(viewport())
   footnoteText = paste(footnoteText,format(Sys.time(), "%d %b %Y"),sep=" | ")
   grid.text(label= footnoteText ,
             x = unit(1,"npc") - unit(2, "mm"),
             y= unit(2, "mm"),
             just=c("right", "bottom"),
             gp=gpar(cex= size, col=color))
   popViewport()
}
subplot.box <- function(data,icol){
	#Plot the variance captured by each PCA component
	# pvar = pca$sdev^2/sum(pca$sdev^2)
	title = "icol"
	#Check if last column is order of magnitude larger than previous, if so, log y-axis
	comp = mean(data[[length(data)]])/mean(data[[length(data)-1]])
	if(abs(comp) > 1 ){logval="y"
	}else{logval=""}
	colors=rainbow(length(data)+1)
	#Plot boxplot
	bardata = boxplot(data,col=colors[2:length(colors)],border="grey",xaxt="n",yaxt="n",main=title,ylab="",xlab="",bg="white",pch=".",lwd=.5,frame=F,log=logval)
	#Alter axis
	axisthickness=1.2
	box(lwd=axisthickness,bty="l")
	axis(1,lwd=axisthickness)
	axis(2,lwd=axisthickness)
}
pca.plot.pvar <- function(pca){
	#Plot the variance captured by each PCA component
	pvar = pca$sdev^2/sum(pca$sdev^2)
	title = "PCA Analysis: %Var of Components"
	#Plot barplot of variance
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