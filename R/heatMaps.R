heatMaps <-
function(variableList,varRank=NULL,Outcome,idnames=NULL,data,title="Heat Map",colCluster=FALSE,prediction=NULL,outcomeGain=0) 
{


if (!requireNamespace("gplots", quietly = TRUE)) {
   install.packages("gplots", dependencies = TRUE)
} 
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
   install.packages("RColorBrewer", dependencies = TRUE)
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


	vnames <- as.vector(variableList[,1]);
	frm <- paste(Outcome," ~ ",1);
	for (i in 1:nrow(variableList))
	{
		frm <- paste(frm," + ",vnames[i])
	}
	modelFrame <- model.frame(formula(frm),data);

	if (is.null(varRank)) 
	{
		topvarID <- seq(1, length(vnames), 1)
		hits <- seq(1, length(vnames), 1)
	}
	else
	{
		topvarID <- as.numeric(rownames(varRank));
		hits <- as.vector(as.matrix(varRank))
	}
	if (outcomeGain==0)
	{
		orderData <- cbind(data[,Outcome]);
	}
	else
	{
		orderData <- cbind(outcomeGain*data[,Outcome]-outcomeGain/2);
	}
	orderData <- as.data.frame(orderData);
	colnames(orderData) <- c(Outcome);
	cn = 1; 
	if (!is.null(prediction))
	{
		cn = cn + 1;
		orderData <- cbind(orderData,as.vector(prediction));
		colnames(orderData)[cn] <- "prediction";
	}
	for ( i in 1:length(topvarID))
	{
		if (hits[i]>0)
		{
			orderData <- cbind(orderData,modelFrame[,topvarID[i]+1])
			cn = cn + 1;
			colnames(orderData)[cn] <- colnames(modelFrame)[topvarID[i]+1];
		}
	}

	dataMat <- as.matrix(orderData);
	

	orderData <- eval(parse(text=paste("with(orderData,orderData[order(",Outcome,",orderData[,2]),])")))
	
	orderData <- as.matrix(orderData);
	
	if (is.null(prediction))
	{
		if (colCluster==FALSE)
		{
			heatMap <- gplots::heatmap.2(orderData, 
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
			heatMap <- gplots::heatmap.2(orderData, 
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
			heatMap <- gplots::heatmap.2(orderData, 
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
