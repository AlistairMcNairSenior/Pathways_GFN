

rm(list=ls())
detach(data)

directory<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Pathways Study/Cell Metab Submission"
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
data<-read.csv("LipidPeroxidation master.csv")
dim(data)
data$Drug<-as.factor(data$Drug)

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

# In this analysis we are only interested in control animals
data<-data[which(data$Drug == "control"),]

# Set the WD for Fig 4 
setwd(paste0(directory, "/Fig_4"))

#################################################################
## DEFINE THE OUTCOMES AND TRIM OUT THOSE WITH MISSING DATA #####
#################################################################

# Which outcomes to analyse
outcomes<-c("MDAnmol.mg")

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
families<-list(gaussian(link="log"))

################################################################
################### Surfaces for each group ####################
################################################################

# Lets do PCF, PFC, CFP
XYZ_list<-list(c("protein.intake", "carbohydrate.intake", "fat.intake"), c("protein.intake", "fat.intake", "carbohydrate.intake"), c("carbohydrate.intake", "fat.intake", "protein.intake"))
labels_list<-list(c("Protein kJ/day", "Carbohydrate kJ/day", "Fat kJ/day"), c("Protein kJ/day", "Fat kJ/day", "Carbohydrate kJ/day"), c("Carbohydrate kJ/day", "Fat kJ/day", "Protein kJ/day"))

# List the limits for xs, ys - determined by trial and error
x_lims<-list(c(0, 30), c(0, 35), c(5, 45))
y_lims<-list(c(5, 37.5), c(5, 45), c(5, 17.5))

# Run through for each outcome, and for each fit the gam, comparisons between the control and treatment surfaces are taken from the LRT above

for(o in 1:length(outcomes)){

	# Get the control surfaces
	bw_surfaces<-list()
	
	# Counter for the surfaces
	counter<-1
	
	# Find the right data and fit the model
	data_g<-data_list[[o]]
	data_g$outcome<-data_g[,outcomes[o]]
	gamma<-max(1.5, log(mean(table(data_g$Drug)))/2)
	k<-21
	GAM<-gam(outcome ~ s(protein.intake, carbohydrate.intake, fat.intake, k = k), data=data_g, family=families[[o]], gamma=gamma)
	write.GAM(GAM, csv.file=paste0(outcomes[o], "output.csv"))
	
	# pdf for the gam check
	pdf(paste0(outcomes[o], "_GAM_Check.pdf"))
	
	# Run the gam check on the GAM
	print(outcomes[o])
	print(gam.check(GAM))

	# Close the pdf for the gam checks
	dev.off()
					
	# Make the three surfaces and find the min and max
	surf_min<-NA
	surf_max<-NA
	ranges_list<-list()
	for(i in 1:3){
		surf<-ggSurface(GAM=GAM, data=data_g, XYZ=XYZ_list[[i]], labels=labels_list[[i]], exclude=NA, lab_size=5, axis_lab_size=22)
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
		bw_surfaces[[counter]]<-ggSurface(GAM=GAM, data=data_g, XYZ=XYZ_list[[i]], labels=labels_list[[i]], exclude=NA, surf_min=surf_min, surf_max=surf_max, lab_size=5, axis_lab_size=22, x.limits=ranges_list[[i]][c(1,2)], y.limits=ranges_list[[i]][c(3,4)], contour_at=0.07, skip=1)
		counter<-counter+1
	}

	# Now lets arrange all those plots
	pdf(paste0(outcomes[o], ".pdf"), height=5, width=15)
	
	# Lay them out nicely
	grid.arrange(bw_surfaces[[1]], 
					bw_surfaces[[2]],
					bw_surfaces[[3]], 
					layout_matrix=rbind(c(1,2,3)))
	
	dev.off()

}

