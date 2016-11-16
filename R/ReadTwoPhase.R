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
	# external force
	DATA$Fx<-Fx
	DATA$Fy<-Fy
	DATA$Fz<-Fz
	
	DATA$Source<-PATH
		
	return(DATA)
}

ReadComponent<-function(FILENAME){
    DATA<-read.csv(file=FILENAME,head=TRUE,sep=" ")
    
    # Remove non-sequential rows from WP, DATA (based on time)
    # (redundant information from restarts)
    DATA$dt<-append(diff(DATA$time,lag=1),0,after=0)
    DATA$TAG<-"KEEP"
    MAXLABEL=0
    for (i in 1:nrow(DATA)){
        #if (DATA$label[i] > MAXLABEL){
        #   MAXLABEL=DATA$label[i]
        #}
        if (DATA$dt[i]<0){
            TIME=DATA$time[i-1]
            SKIP=DATA$dt[i]/1000
            while(DATA$time[i]<TIME+1000){
                #skip these rows
                DATA$TAG[i]="DISCARD"
                i=i+1
            }
            # update all ensuing labels so that none are repeated
            #for (j in i:nrow(DATA)){
            #   DATA$label[j] = DATA$label[j] + MAXLABEL
            #}
        }
    }
    DATA<-subset(DATA,DATA$TAG=="KEEP")
    DATA$TAG<-NULL
    return(DATA)
}

LoadGangliaTrajectories<-function(PATH){
    
    # Read the non-wetting phase component averages
    FILENAME=paste0(PATH,"/components.NWP.tcat")
    NWP<-ReadComponent(FILENAME)
    
    # Read the input parameters
    FILENAME=paste0(PATH,"/Color.in")
    PRM<-scan(file=FILENAME)
    IFT<-PRM[2]*5.796
    VISC<-(PRM[1]-0.5)/3
    Fx<-PRM[6]
    Fy<-PRM[7]
    Fz<-PRM[8]
    
    FILENAME=paste0(PATH,"/Domain.in")
    DM<-scan(file=FILENAME)
    Px<-DM[1]
    Py<-DM[2]
    Pz<-DM[3]
    Nx<-DM[4]
    Ny<-DM[5]
    Nz<-DM[6]
    
    # compute the domain length and volume
    Lx=Px*Nx
    Ly=Py*Ny
    Lz=Pz*Nz
    VOLUME=Lx*Ly*Lz
    
    NWP$IFT<-IFT	
    NWP$VOL<-VOLUME
    NWP$mu<-VISC
    return(NWP)
    
}

