

rm(list=ls())
detach(data)

directory<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Pathways Study/Cell Metab Resubmission"
#directory<-"/Users/asenior/Dropbox (Sydney Uni)/Pathways Study/Cell Metab Submission"
setwd(directory)

# Load  the libraries
library(mgcv)
library(sp)
library(lattice)
library(ellipse)
library(gplots)
library(geometry)
library(MASS)
source("Headers.R")
library(arm)
library(plotly)
library(gridExtra)

# Changing the settings of one of the packages (not sure which one), but this is an essential step
options(na.action = NULL)

#################################################################
######### READ IN THE DATA AND WORK OUT THE TREATMENTS ##########
#################################################################

# Read in the data
data<-read.csv("phenotypealldata0312020.csv")
dim(data)
data$Drug<-as.factor(data$Drug)
data$percent.fat.10m<-data$percent.fat.10m/100

# Make sure the intakes of PCF and are right
data$protein.intake<-data$intake.E * (data$P/100)
data$carbohydrate.intake<-data$intake.E * (data$C/100)
data$fat.intake<-data$intake.E * (data$F/100)

# Drop any missing intake data
intakes<-c("protein.intake", "carbohydrate.intake", "fat.intake")
for(i in 1:length(intakes)){
	tag<-which(is.na(data[, intakes[i]]) == T)
	if(length(tag) > 0){
		data<-data[-tag,]
	}
}

# There was one animal with a fglu.ins product of -0.04 - the value cannot be negative, and is likely a measurment reading from a very low animal OR the machine read through with missing data. Amanda B votes missing so I added as NA
data$glu.ins.10m[which(data$glu.ins.10m < 0)]<-NA

# There is also one animal with IGF1 conc of 0 - I will set this as missing data too
data$IGF1.conc.10m[which(data$IGF1.conc.10m == 0)]<-NA

# treatments and controls
Treatments<-c("metformin", "resveratrol", "rapamycin")
Control<-"control"

# Set the WD for Fig 1 
setwd(paste0(directory, "/Fig_1"))

#################################################################
## DEFINE THE OUTCOMES AND TRIM OUT THOSE WITH MISSING DATA #####
#################################################################

# Which outcomes to analyse
outcomes<-c("g.body.wt.10m", "percent.fat.10m", "FGF21.conc.10m", "IGF1.conc.10m", "glu.ins.10m")
titles<-c("Body Weight (g)", "% Body Fat", "FGF21 Conc.", "IGF1 Conc.", "HOMA")

# Screen for missing data - drop missing rows for each outcome and create a histogram
pdf("histograms.pdf")
data_list<-list()
for(i in 1:length(outcomes)){
	tag<-which(is.na(data[,outcomes[i]]) == T)
	if(length(tag) > 0){
		data_list[[i]]<-data[-tag,]
	}else{
		data_list[[i]]<-data
	}
	hist(data_list[[i]][,outcomes[i]], main=outcomes[i])
}
dev.off()

# Model all as gaussian for now and specify k
families<-list(gaussian(link="identity"), betar(link="logit"), Gamma(link="log"), Gamma(link="log"), Gamma(link="identity"))

################################################################
################### Surfaces for each group ####################
################################################################

# Groups
groups<-c(Control, Treatments)

# Lets do PCF, PFC, CFP
XYZ_list<-list(c("protein.intake", "carbohydrate.intake", "fat.intake"), c("protein.intake", "fat.intake", "carbohydrate.intake"), c("carbohydrate.intake", "fat.intake", "protein.intake"))
labels_list<-list(c("Protein kJ/day", "Carbohydrate kJ/day", "Fat kJ/day"), c("Protein kJ/day", "Fat kJ/day", "Carbohydrate kJ/day"), c("Carbohydrate kJ/day", "Fat kJ/day", "Protein kJ/day"))

# List the limits for xs, ys - determined by trial and error
x_lims<-list(c(0, 30), c(0, 35), c(5, 45))
y_lims<-list(c(5, 37.5), c(5, 45), c(5, 17.5))

# Based on the previous submission I believe a fat to protein ratio of 1.4 maximises BW so will draw a rail on this on the central panels
FP_ratio<-1.4

# Run through for each outcome, and for each fit the gam, comparisons between the control and treatment surfaces are taken from the LRT above

for(o in 1:length(outcomes)){

	# Get the control surfaces
	bw_surfaces<-list()
	
	# Counter for the surfaces
	counter<-1
	
	# Find the right data and fit the model
	data_g<-data_list[[o]]
	data_g$outcome<-data_g[,outcomes[o]]
	gamma<-log(mean(table(data_g$Drug)))/2
	k<-80
	GAM<-gam(outcome ~ s(protein.intake, carbohydrate.intake, fat.intake, k = k, by=Drug) + Drug, data=data_g, family=families[[o]], gamma=gamma)
	write.GAM(GAM, csv.file=paste0(outcomes[o], "output.csv"))
	
	# pdf for the gam check
	pdf(paste0(outcomes[o], "_GAM_Check.pdf"))
	
	# Run the gam check on the GAM
	print(outcomes[o])
	print(gam.check(GAM))

	# Close the pdf for the gam checks
	dev.off()

	# Do a few tests to see if smoothing  by drug group improves model fit
	GAM_null<-gam(outcome ~ s(protein.intake, carbohydrate.intake, fat.intake, k = k) + Drug, data=data_g, family=families[[o]], gamma=gamma)
	
	print(anova(GAM_null, GAM, test="Chisq"))
	print(AIC(GAM_null, GAM))
	
	# Loop for each groups
	for(j in 1:length(groups)){
		
		# Just use the data for the jth group
		data_gj<-data[which(data$Drug == groups[j]),]
				
		# Make the three surfaces and find the min and max
		surf_min<-NA
		surf_max<-NA
		ranges_list<-list()
		for(i in 1:3){
			surf<-ggSurface(GAM=GAM, data=data_gj, XYZ=XYZ_list[[i]], labels=labels_list[[i]], exclude=NA, lab_size=5, axis_lab_size=22, predict_val=data.frame(Drug=groups[j]))
			surf_min<-c(surf_min, min(surf$data$fit))
			surf_max<-c(surf_max, max(surf$data$fit))
			ranges_list[[i]]<-c(range(pretty(range(surf$data[,1]))), range(pretty(range(surf$data[,2]))))
		}
		
		# Adjust the surface colour scale mins and maxs by 2% to ensure range is good
		surf_min<-min(surf_min, na.rm=T)
		if(surf_min > 0){
			surf_min<-surf_min * 0.98
		}
		if(surf_min < 0){
			surf_min<-surf_min * 1.02
		}
		surf_max<-max(surf_max, na.rm=T)
		if(surf_max > 0){
			surf_max<-surf_max * 1.02
		}
		if(surf_max < 0){
			surf_max<-surf_max * 0.98
		}
			
		# Now make again and actually use the min and max
		for(i in 1:3){
			p<-ggSurface(GAM=GAM, data=data_gj, XYZ=XYZ_list[[i]], labels=labels_list[[i]], exclude=NA, surf_min=surf_min, surf_max=surf_max, lab_size=5, axis_lab_size=22, x.limits=ranges_list[[i]][c(1,2)], y.limits=ranges_list[[i]][c(3,4)], predict_val=data.frame(Drug=groups[j]))
			
			# If we are on the central panel add the rail
			if(i == 2){
				p<-p+geom_abline(intercept=0, slope=FP_ratio, col="darkgrey", size=2)
			}
			
			bw_surfaces[[counter]]<-p
			counter<-counter+1
		}
	
	}

	# Now lets arrange all those plots
	pdf(paste0(outcomes[o], ".pdf"), height=20, width=15)
	
	# Lay them out nicely
	grid.arrange(bw_surfaces[[1]], 
					bw_surfaces[[2]],
					bw_surfaces[[3]], 
					bw_surfaces[[4]], 
					bw_surfaces[[5]], 
					bw_surfaces[[6]], 
					bw_surfaces[[7]], 
					bw_surfaces[[8]], 
					bw_surfaces[[9]], 
					bw_surfaces[[10]], 
					bw_surfaces[[11]], 
					bw_surfaces[[12]],
					layout_matrix=rbind(c(1,2,3),
										c(4,5,6),
				  					  	c(7,8,9),
				  					  	c(10,11,12)))
	
	dev.off()
	
	# Now I will make a plot along the grey vector on each drug
	
	# Get the P and F along the vector
	prediction_data<-data.frame(protein.intake=seq(5, 12.5, 0.001))
	prediction_data$fat.intake<-prediction_data$protein.intake * FP_ratio
	
	# Add in the C values and drug groups
	prediction_data$carbohydrate.intake<-median(data$carbohydrate.intake[which(data$Drug == groups[1])])
	prediction_data$Drug<-groups[1]
	for(j in 2:length(groups)){
		prediction_j<-prediction_data
		prediction_j$carbohydrate.intake<-median(data$carbohydrate.intake[which(data$Drug == groups[j])])
		prediction_j$Drug<-groups[j]
		prediction_data<-rbind(prediction_data, prediction_j)
	}
	
	# Now do the prediction from the GAM
	out<-predict(GAM, newdata=prediction_data, type="response", se.fit=T)
	prediction_data$y<-out$fit
	prediction_data$se<-out$se.fit
	prediction_data$intake<-prediction_data$protein.intake + prediction_data$fat.intake
	prediction_data$lci<-prediction_data$y - prediction_data$se
	prediction_data$uci<-prediction_data$y + prediction_data$se
	
	# Now do the plot
	p<-ggplot(data=prediction_data, aes(x=intake, y=y, col=Drug, fill=Drug)) +
		geom_ribbon(aes(x=intake, ymin=lci, ymax=uci), alpha=0.2, linetype=0) +
		geom_line(size=1.5) +
		theme_bw() + 
		xlab("Intake kJ/day") + ylab(titles[o]) +
		theme(legend.position=c(0.2,0.8), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15), legend.title=element_text(size=15))
	
	# Add to a pdf
	pdf(paste0(outcomes[o], "_vector_plot.pdf"), height=6, width=6)
		grid.arrange(p)
	dev.off()
}


################################################################
################## GAM with smoothing by group LRT #############
################################################################

output<-GAM.LRT(data=data, outcomes=outcomes, treatment.col="Drug", Treatments=Treatments, Control=Control, families=families, k=40, save=F, intake.P="protein.intake", intake.C="carbohydrate.intake", intake.F="fat.intake")

# Pull out the results of the LRT and the models themselves
LRT_results<-output[[1]]

LRT_results[[1]]$p.int < 0.05
LRT_results[[2]]$p.int < 0.05
LRT_results[[3]]$p.int < 0.05