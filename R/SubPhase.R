library(ggplot2)
library(pspline)
library(reshape2)
library(data.table)

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
	return(DATA)	
}


ExtractSubphase<-function(PATH){

	# Read the non-wetting phase component averages	
	FILENAME=paste0(PATH,"/components.NWP.tcat")
	NWP<-ReadComponent(FILENAME)

	# Read the wetting phase component averages
	FILENAME=paste0(PATH,"/components.WP.tcat")
	WP<-ReadComponent(FILENAME)
	
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

	# Get the biggest blob for all time values
	BMAX<-NWP[as.logical(with(NWP,ave(vol,time,FUN=function(x) x==max(x)))),]
	
	# Map all the non-wetting components to integrals (instead of averages)
	NWP$pwn<-NWP$pn*NWP$awn
 	NWP$pn<-NWP$pn*NWP$vol
	NWP$vx<-NWP$vx*NWP$vol
	NWP$vy<-NWP$vy*NWP$vol
	NWP$vz<-NWP$vz*NWP$vol
	NWP$vsq<-NWP$vsq*NWP$vol
	NWP$trJwn<-NWP$trJwn*NWP$awn
	NWP$count<-1
	
	# Generate blob sum totals for selected quantities
	NWP.table<- data.table(NWP,key='time')

	BSUMS<-NWP.table[,list(vol=sum(vol),pn=sum(pn),awn<-sum(awn),ans<-sum(ans),vnx<-sum(vx),vny<-sum(vy),
			  vnz<-sum(vz),vsq<-sum(vsq),
			  Euler=sum(Eulern),pwn<-sum(pwn),count<-sum(count),trJwn<-sum(trJwn)),by='time']

	setnames(BSUMS,c("time","Vn","pn","awn","ans","vnx","vny","vnz","vsq","Euler","pwn","count","trJwn"))

	BSUMS$pn<-BSUMS$pn/BSUMS$Vn
	BSUMS$vnx<-BSUMS$vnx/BSUMS$Vn
	BSUMS$vny<-BSUMS$vny/BSUMS$Vn
	BSUMS$vnz<-BSUMS$vnz/BSUMS$Vn
	BSUMS$vsq<-BSUMS$vsq/BSUMS$Vn
	BSUMS$pwn<-BSUMS$pwn/BSUMS$awn
	BSUMS$trJwn<-BSUMS$trJwn/BSUMS$awn

	setnames(BSUMS,c("time","Vn.global","pn.global","awn.global","ans.global",
		"vnx.global","vny.global","vnz.global","vsq.global","Euler.global","pwn.global","B0","trJwn.global"))

	BSUMS$IFT<-IFT
	BSUMS$Viscosity<-VISC
	
	TMP<-merge(BSUMS,WP[,c("time","vol","pw","vx","vy","vz")],by="time")

	BMAX<-BMAX[,c("time","vol","pn","awn","ans","vx","vy","vz","vsq","Eulern","trJwn")]
	setnames(BMAX,c("time","Vn.c","pn.c","awn.c","ans.c","vnx.c","vny.c","vnz.c","vsq.c","Euler.c","trJwn.c"))

	DATA<-merge(TMP,BMAX,by="time",allow.cartesian=TRUE)

	# estimate the total volume of the porespace from the data
	PoreVolume=mean(DATA$vol+DATA$Vn.global)
	porosity=1.0 - PoreVolume/VOLUME

	DATA$Vn.d<-DATA$Vn.global-DATA$Vn.c
	DATA$awn.d<-DATA$awn.global-DATA$awn.c
	DATA$trJwn.d<-(DATA$trJwn.global*DATA$awn.global-DATA$trJwn.c*DATA$awn.c)/DATA$awn.d
	DATA$pn.d<-(DATA$pn.global*DATA$Vn.global-DATA$pn.c*DATA$Vn.c)/DATA$Vn.d
	DATA$vnx.d<-(DATA$vnx.global*DATA$Vn.global-DATA$vnx.c*DATA$Vn.c)/DATA$Vn.d
	DATA$vny.d<-(DATA$vny.global*DATA$Vn.global-DATA$vny.c*DATA$Vn.c)/DATA$Vn.d
	DATA$vnz.d<-(DATA$vnz.global*DATA$Vn.global-DATA$vnz.c*DATA$Vn.c)/DATA$Vn.d
	DATA$vsq.d<-(DATA$vsq.global*DATA$Vn.global-DATA$vsq.c*DATA$Vn.c)/DATA$Vn.d


  	DATA$sw<-1-DATA$Vn.global/(VOLUME*porosity)
	DATA$sn.d<-DATA$Vn.d/(VOLUME*porosity)
	DATA$pwn.d<-DATA$pw-DATA$pn.d
	DATA$pwn.c<-DATA$pw-DATA$pn.c
	DATA$pwn.global<-DATA$pw-DATA$pn.global

	FILENAME=paste0(PATH,"/blobsum.tcat")	
	write.table(DATA,file=FILENAME,row.names=FALSE,sep=" ",quote=FALSE)

	return(DATA)
}

PlotSubPhase<-function(DATA,BASENAME){

#interfacial tension in kPa-cm
#IFTPHYS=0.0024
# same thing in Pa-mm
IFTPHYS=24.0
# viscosity in Pa s (N s / m^2) x 10^{-3}
MUPHYS=1.0
#bead diameter in mm (DPHYS) and voxels (DSIM)
DPHYS=1.0
DSIM=60

PWN=expression(p^n~-~p^w)
PC=expression(gamma^{wn}~J[w]^{wn})
SW=expression(s^{bar(bar(w))})
SNC=expression(s^{bar(bar(n[c]))})
AWN=expression(epsilon^{bar(bar(wn))})
eni=expression(epsilon^{n[i]})

p<-ggplot()+
	geom_line(data=DATA,aes(sw,sn.d,colour="DATA"))+
	geom_line(data=DATA,aes(sw,sn.d,colour="DATA"))+
	xlab(SW)+ylab(SNC)+
	scale_colour_manual(name="Source",values=c("red","blue"))+
        theme_bw()


FILENAME=paste0(BASENAME,"-snr.pdf")
ggsave(FILENAME,p,width=5.5,height=4.0)

p<-ggplot()+
	geom_line(data=DATA,aes(sw,IFT*trJwn.global,colour="global"))+
	geom_line(data=subset(DATA,DATA$sn.d>0.01),aes(sw,IFT*trJwn.d,colour="trapped"))+
	geom_line(data=DATA,aes(sw,IFT*trJwn.c,colour="connected"))+
	xlab(SW)+ylab(PC)+
	scale_x_continuous(limits=c(0,1))+
	scale_colour_manual(name="Region",values=c("red","black","dodgerblue1"))+
        theme_bw()


FILENAME=paste0(BASENAME,"-Jwn.pdf")
ggsave(FILENAME,p,width=5.5,height=4.0)

p<-ggplot()+
	geom_line(data=DATA,aes(sw,pwn.global,colour="global"))+
	geom_line(data=subset(DATA,DATA$sn.d>0.01),aes(sw,pwn.d,colour="trapped"))+
	geom_line(data=DATA,aes(sw,pwn.c,colour="connected"))+
	xlab(SW)+ylab(PWN)+
	scale_colour_manual(name="Region",values=c("red","black","dodgerblue1"))+
	scale_x_continuous(limits=c(0,1))+
        theme_bw()

FILENAME=paste0(BASENAME,"-pwn.pdf")
ggsave(FILENAME,p,width=5.5,height=4.0)

}