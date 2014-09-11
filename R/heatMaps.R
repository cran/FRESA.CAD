heatMaps <-
function(varList,varRank=NULL,outcome,idnames=NULL,dataframe,title="Heat Map",colCluster="NA",prediction=NULL,outcomeGain=0) 
{

if (!require(gplots)) {
  heatmap.2 <- function(...) {}
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

	# creates a own color palette from red to blue
	my_palette <- colorRampPalette(c("blue","cyan","black","yellow","red"))(n = 224)

	# defines the color breaks manually for a "skewed" color transition
	col_breaks = c(
	seq(-2.0,-1.10,length=50),  		# for blue
	seq(-1.1,-0.20,length=50),          # for cyan
	seq(-0.20,0.20,length=25),          # for black
	seq(0.20,1.1,length=50),            # for yelos
	seq(1.1,2.0,length=50))             # for red

	vnames <- as.vector(varList[,1]);
	if (is.null(varRank)) 
	{
		topvarID <- seq(1, length(vnames), 1)
		hits <- seq(1, length(vnames), 1)
	}
	else
	{
		topvarID <- as.numeric(rownames(varRank));
		hits <- as.vector(varRank)
	}
	if (outcomeGain==0)
	{
		orderData <- cbind(dataframe[,outcome]);
	}
	else
	{
		orderData <- cbind(outcomeGain*dataframe[,outcome]-outcomeGain/2);
	}
	if (!is.null(prediction))
	{
		orderData <- cbind(orderData,as.vector(prediction));
	}

	cnames <- vector();
#	attach(dataframe,name="FRESA:HeatMapDataFrame")
#	on.exit(detach("FRESA:HeatMapDataFrame"))
	for ( i in 1:length(topvarID))
	{
		if (hits[i]>0)
		{
			s <- toString(varList[topvarID[i],1]);
			s <- str_replace_all(s,"I\\(","\\(");
			orderData <- cbind(orderData,with(dataframe,eval(parse(text=s))))
			s <- str_replace_all(s," ","");
			s <- str_replace_all(s,"\\*","");
			if (!is.na(str_locate(s,")")[1]))
			{			
				s <- str_replace_all(s,")","");
				s <- str_replace_all(s,"\\(","");
				s <- str_replace_all(s,"<","L");
				s <- str_replace_all(s,">","G");
				s <- str_replace_all(s,"\\-","");
			}
			cnames <- append(cnames,s);
		}
	}
#	detach("FRESA:HeatMapDataFrame")

	orderData <- as.data.frame(orderData);
	if (!is.null(idnames))
	{
		rownames(orderData) <- dataframe[,idnames];
	}
	if (is.null(prediction))
	{
		colnames(orderData) <- append(outcome,cnames);
	}
	else
	{	
		nm <- append(outcome,"Prediction");
		colnames(orderData) <- append(nm,cnames);
	}
	dataMat <- as.matrix(orderData);

	orderData <- eval(parse(text=paste("with(orderData,orderData[order(",outcome,",orderData[,2]),])")))
	
	orderData <- as.matrix(orderData);
	
	if (is.null(prediction))
	{
		if (colCluster=="NA")
		{
			heatMap <- heatmap.2(orderData, 
			main = title, 			# heat map title
			notecol="black",      	# change font color of cell labels to black
			density.info="none",  	# turns off density plot inside color legend
			trace="none",         	# turns off trace lines inside the heat map
			margins =c(12,9),     	# widens margins around plot
			col=my_palette,       	# use on color palette defined earlier 
			breaks=col_breaks,    	# enable color transition at specified limits
			dendrogram="row", 		# only draw a row dendrogram
			Colv="NA") 				# turn off column clustering
		}
		else
		{
			heatMap <-heatmap.2(orderData, 
			main = title, 			# heat map title
			notecol="black",      	# change font color of cell labels to black
			density.info="none",  	# turns off density plot inside color legend
			trace="none",         	# turns off trace lines inside the heat map
			margins =c(12,9),     	# widens margins around plot
			col=my_palette,       	# use on color palette defined earlier 
			breaks=col_breaks,    	# enable color transition at specified limits
			) 				
		}
	}
	else
	{
			heatMap <- heatmap.2(orderData, 
			main = title, 			# heat map title
			notecol="black",      	# change font color of cell labels to black
			density.info="none",  	# turns off density plot inside color legend
			trace="none",         	# turns off trace lines inside the heat map
			margins =c(12,9),     	# widens margins around plot
			col=my_palette,       	# use on color palette defined earlier 
			breaks=col_breaks,    	# enable color transition at specified limits
			dendrogram="none", 		# sorted by prediction
			Colv="NA",				# turn off column clustering
			Rowv="NA") 				# turn off row clustering
	}

	result <- list(dataMatrix=dataMat,
	orderMatrix=orderData,
	heatMap=heatMap)
	return (result);
}
