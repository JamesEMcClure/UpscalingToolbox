# Functions to read average logs for two-fluid simulations

library(ggplot2)
#library(extrafont)
library(reshape2)
library(pspline)

ReadTwoPhaseAverages<-function(PATH){
	# Read two phase averages from 'timelog.tcat' file and non-dimensionalize
	# simluation input and output files must be present at PATH
	#    Color.in
	#    Domain.in
	#    timelog.tcat

        # Read the input parameters
        FILENAME=paste0(PATH,"/Color.in")
        PRM<-scan(file=FILENAME)        
        IFT<-PRM[2]*5.796
        VISC<-(PRM[1]-0.5)/3
        Fx<-PRM[6]
        Fy<-PRM[7]
        Fz<-PRM[8]
	PI=3.14159265358979

        FILENAME=paste0(PATH,"/Domain.in")
        DM<-scan(file=FILENAME)        
        Px<-DM[1]
        Py<-DM[2]
        Pz<-DM[3]
        Nx<-DM[4]
        Ny<-DM[5]
       	Nz<-DM[6]

	NX=Px*Nx
	NY=Py*Ny
	NZ=Pz*Nz

        # Read the time history of averages 
        FILENAME=paste0(PATH,"/timelog.tcat")
        DATA<-read.csv(file=FILENAME,head=TRUE,sep=" ")

	# Compute the Sauter mean for non-dimensionalization
	#D32=

        # interfacial tension
        DATA$IFT<-IFT
        # dynamic viscosity
        DATA$viscosity<-VISC
	
	DATA$Source<-PATH
		
	return(DATA)
}

