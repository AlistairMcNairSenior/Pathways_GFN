
inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) { 
	
	library(geometry)
	
	# https://tolstoy.newcastle.edu.au/R/e8/help/09/12/8784.html
	calpts <- as.matrix(calpts) 
	testpts <- as.matrix(testpts) 
	p <- dim(calpts)[2] 
	cx <- dim(testpts)[1] # rows in testpts
	nt <- dim(hull)[1] # number of simplexes in hull 
	nrmls <- matrix(NA, nt, p)
	
	degenflag <- matrix(TRUE, nt, 1) 

	for (i in 1:nt) { 
	
		nullsp<-t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))

		if (dim(nullsp)[1] == 1){
			nrmls[i,]<-nullsp
		  degenflag[i]<-FALSE
		}
	}

	if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
		nrmls <- nrmls[!degenflag,] 
		nt <- dim(nrmls)[1] 
		
		center = apply(calpts, 2, mean) 
		a <- calpts[hull[!degenflag,1],] 
		
		nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
		
		dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
		nrmls <- nrmls*matrix(dp, nt, p)
		
		aN <- diag(a %*% t(nrmls)) 
		val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 
		
		val[abs(val) < tol] <- 0 
		as.integer(sign(val)) 
}



# A function created to find the outer perimeter over which the surface should be fitted
findConvex<-function(x, y, rgnames, res=101, x.limits=NA, y.limits=NA){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	
	# Either use the specifiec limits for the grid, or use pretty
	if(length(which(is.na(x.limits) == T)) > 1){
		px<-pretty(x)	
	}else{
		px<-x.limits	
	}
	if((length(which(is.na(y.limits) == T)) > 1)){
		py<-pretty(y)	
	}else{
		py<-y.limits	
	}	
	
	# Get the matrix
	x.new<-seq(min(px, na.rm=T),max(px, na.rm=T),len=res)
	y.new<-seq(min(py, na.rm=T),max(py, na.rm=T),len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}



# Function to run LRTs on GAMs with nutrient smoothing by treatment for outcomes
# returns a list of dim 2
# element 1 is a list of dim treatments of LRT stats for each outcome
# element 2 is a list of dim treatments of models for each outcome
# Takes a dataset (data), a vector of column names corresponding to outcomes (outcomes), the name of the column containing the treatments (treatment.col), a vector of labels for the Treatments (Treatments), the label for the control (Control), a vector of GAM families of length outcomes (families), k value for the GAM (k), labels for intake columns (intake.P, intake.C and intake.F)

GAM.LRT<-function(data, outcomes, treatment.col, Treatments, Control, families, k, intake.P="intake.P", intake.C="intake.C", intake.F="intake.F", save=F){
	
	# Formatting the progress bar
	pb <- txtProgressBar(min = 0, max = length(outcomes) * length(Treatments), style = 3)
	progress<-0
	
	# Ensure we have the intake and Drug columns
	data$intake.P<-data[,intake.P]
	data$intake.C<-data[,intake.C]
	data$intake.F<-data[,intake.F]
	data$Drug<-data[,treatment.col]
	
	# List to hold the models
	models<-list()
	LRT.tests<-list()
	
	
	for(d in 1:length(Treatments)){
		
		# Create the dth dataset
		data.d<-rbind(data[which(data$Drug == Treatments[d]), ], data[which(data$Drug == Control), ])
		data.d<-droplevels(data.d)
		
		# List to hold the models for drug d
		models.d<-list()
		
		# Table to hold results of dth set of LRTs
		LRT.results.d<-data.frame(outcome=outcomes, Res.DF.add=NA, Res.Dev.add=NA, DF.add=NA, Dev.add=NA, p.add=NA, Res.DF.int=NA, Res.Dev.int=NA, DF.int=NA, Dev.int=NA, p.int=NA)
		
		# Loop to do outcomes for dth drug treatment
		for(o in 1:length(outcomes)){
			
			# Get an oth copy
			data.d.o<-data.d
			
			# Specify the outcome
			data.d.o$outcome<-data.d.o[,outcomes[o]]
			
			# Drop any missing data
			drop<-which(is.na(data.d.o$outcome) == T)
			if(length(drop) > 0){
				data.d.o<-data.d.o[-drop,]
			}
			
			# Fit the GAM with smoothing by group and without
			gamma<-log(dim(data.d.o)[1])/2
			GAM.int<-gam(outcome ~ s(intake.P, intake.C, intake.F, by = Drug, k = k) + Drug, data=data.d.o, method="REML", family=families[[o]], gamma=gamma)
			GAM.add<-gam(outcome ~ s(intake.P, intake.C, intake.F, k = k) + Drug, data=data.d.o, method="REML", family=families[[o]], gamma=gamma)
			GAM<-gam(outcome ~ s(intake.P, intake.C, intake.F, k = k), data=data.d.o, method="REML", family=families[[o]], gamma=gamma)
			test<-anova(GAM, GAM.add, GAM.int, test="Chisq")	
			
			# Save the reuslts of test and the GAM with int
			LRT.results.d[o,c(2:6)]<-test[2,]
			LRT.results.d[o,c(7:11)]<-test[3,]				
			models.d[[o]]<-GAM.int
			names(models.d)[[o]]<-outcomes[o]
			
			# remove this set of results
			rm(GAM.int)
			rm(GAM)
			rm(test)
			rm(data.d.o)
			
			# Update the progress
			progress<-progress + 1
			setTxtProgressBar(pb, progress)
		}
		
		# Save externally if we are doing that
		if(save == T){
			# Table to hold results of LRT
			file.name<-paste0(Treatments[d], "_LRT_GAM_Results.csv")
			write.table(LRT.results.d, file=file.name, sep=",", row.names=F, col.names=names(LRT.results.d))
			file.name<-paste0(Treatments[d], "GAMs.Rdata")
			save(models.d, file=file.name)
		}
		
		# Save the models and the LRTs
		models[[d]]<-models.d
		names(models)[[d]]<-Treatments[d]
		LRT.tests[[d]]<-LRT.results.d
		names(LRT.tests)[[d]]<-Treatments[d]
		
		# Remove dth set of results
		rm(models.d)
		rm(LRT.results.d)
		
	}
	
	output<-list()
	output[[1]]<-LRT.tests
	output[[2]]<-models
	return(output)	
}

# Function to make nice surfaces using ggplot

ggSurface<-function(GAM, data, XYZ, labels, predict_val=NA, surf_min=NA, surf_max=NA, x.limits=NA, y.limits=NA, z.val=NA, exclude, subtitle="", lab_size=3, nlevels=5, contour_at=NA, skip=0, axis_size=15, axis_lab_size=15){
	
	require(ggplot2)
	require(sp)
	require(geometry)
	require(mgcv)
	require(metR)
	
	# This specifies the color scheme for surface
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(256)
	
	# List to hold the plots
	plots_list<-list()

	# List for the order of plots
	nutrient.order<-XYZ[c(1,2,3)]
	
	# List for the labels
	labels.order<-labels[c(1,2,3)]
				
	# Values to predict over, if they are unspecified
	if(is.na(x.limits)[1] == T){
		x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
	}
	if(is.na(y.limits)[1] == T){
		y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
	}		
			
	# If we do not specify values to slice at, use the 25, 50, and 75 %ile
	if(is.na(z.val) == T){
		z.val<-round(median(data[,nutrient.order[3]]))
	}
			
	# Fitted list to hold some results for later
	x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=501)
	y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=501)
	z.new<-z.val
	predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
	names(predictors)<-nutrient.order
	in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
			
	# Add the predictors for the additional 'confounders'
	predictors<-cbind(predictors, predict_val)
	predictors<-predictors[-which(in.poly == 0),]
			
	# Do the predictions
	predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, newdata.guaranteed=T)
		
	# Get the proedictions for the kth trait					
	predictions_k<-predictions
							
	# Find the min and max values across all predictions
	mn<-surf_min
	mx<-surf_max
	if(is.na(mn)==T){
		mn<-min(predictions_k, na.rm=T)
	}
	if(is.na(mx)==T){
		mx<-max(predictions_k, na.rm=T)
	}
	locs<-(range(predictions_k, na.rm=TRUE) - mn) / (mx-mn) * 256	
	
	plot_data<-predictors
	plot_data$fit<-predictions_k
	plot_data$x<-plot_data[,nutrient.order[1]]
	plot_data$y<-plot_data[,nutrient.order[2]]
	
	# Set the contour
	if(is.na(contour_at)[1] == T){
		contour_use<-signif((max(predictions_k, na.rm=T)-min(predictions_k, na.rm=T))/nlevels, 1)
	}else{
		contour_use<-contour_at	
	}
	
	# Make the plot
	plot<-ggplot(plot_data, aes(x=x, y=y)) +
			geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
			scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
			geom_contour(data=plot_data, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +	
			geom_label_contour(data=plot_data, aes(x=x, y=y, z=fit), size=lab_size, binwidth=contour_use, skip=skip) +
			theme_bw() +
			labs(x = labels.order[1], y = labels.order[2], subtitle=subtitle) +
			theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_lab_size)) +
			theme(title=element_text(size=15)) + 
			xlim(x.limits) +
			ylim(y.limits)
				
	return(plot)

}


# Write GAM

write.GAM<-function(GAM, csv.file){
	
	# create the csv file to write to, if you want one
	write.table(Sys.time(), file=csv.file, sep=",", row.names=F, col.names=F)
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Add the formula
	write.table(as.character(GAM$formula), file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Get the summary
	summary_GAM<-summary(GAM)
	
	# Write the linear terms
	p.table<-round(summary_GAM$p.table, 4)
	p.table<-as.data.frame(cbind(row.names(p.table), p.table))
	names(p.table)[1]<-"Coef."
	suppressWarnings(write.table(p.table, file=csv.file, sep=",", row.names=F, col.names=names(p.table), append=T))
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Write the smooth terms
	if(length(summary_GAM$s.table) > 0){
		s.table<-round(summary_GAM$s.table, 4)
		s.table<-as.data.frame(cbind(row.names(s.table), s.table))
		names(s.table)[1]<-"Coef."
		suppressWarnings(write.table(s.table, file=csv.file, sep=",", row.names=F, col.names=names(s.table), append=T))	
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	}
	
	# Write the n and deviance explained
	dev.expl<-paste0("n = ", summary_GAM$n, ": % Dev. Explained = ", round(summary_GAM$dev.expl * 100, 2), ": AIC = ", round(AIC(GAM)))
	write.table(dev.expl, file=csv.file, sep=",", row.names=F, col.names=F, append=T)		
	write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)


}
