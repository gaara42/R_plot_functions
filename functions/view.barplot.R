# biafra ahanonu
# updated: 2013.10.03 [13:36:36]
# implements barplot via ggplot2

#load ggplot
require(ggplot2)
require(plyr)
ggplotErrorBarplot <-function(data,x,y,fill,addPoints=FALSE,...){
	result = tryCatch({
		# makes a barplot with errorbars using ggplot, variable inputs are strings
		# biafra ahanonu, updated: 2013.12.28

		data$fill = data[,names(data) %in% c(fill)]
		data  <- ddply(data,c(x,fill),
		            transform,
		            types = paste(as.character(fill)," - (n=",length(fill),")",sep = ""))
		types = "types"

		# create functions to get the lower and upper bounds of the error bars
		stderr <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
		lowsd <- function(x){
			val=median(x)-stderr(x);
			if(val<0){
				val=median(x)-median(x)*0.01;
			}
			return(val)
		}
		highsd <- function(x){return(median(x)+stderr(x))}

		# create a ggplot
		thisPlot = ggplot(data,aes_string(x=x,y=y,fill=types), color="white")+
		# first layer is barplot with means
		stat_summary(fun.y=median, geom="bar", position=position_dodge(), colour='white')+
		# second layer overlays the error bars using the functions defined above
		stat_summary(fun.y=median, fun.ymin=lowsd, fun.ymax=highsd, geom="errorbar", position=position_dodge(.9),color = 'black', size=1, width=0.2)+
		# stat_bin(geom="text", aes_string(x=x,y = 0.5, label="..count..", vjust=2),position=position_dodge())+
		theme(line = element_blank(),panel.background = element_rect(fill = "white", colour = NA), text = element_text(size=15))
		# stat_bin(data=data, aes_string(x), geom="text", color="white", inherit.aes=FALSE)
		if(addPoints==TRUE){
			thisPlot = thisPlot+geom_point()
		}
		print(thisPlot)
		return(thisPlot)
	}, error = function(err) {
		print(err)
		print(traceback())
		return(FALSE)
	}, finally = {
		return(thisPlot)
		# print(Sys.time()-startTime); flush.console();
		# stop the cluster
		# stopCluster(cl)
		# return(data.frame())
	})
}