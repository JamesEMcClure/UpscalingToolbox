# Functions to read average logs for two-fluid simulations

library(ggplot2)
library(extrafont)
library(reshape2)
library(pspline)

ReadPistonTimelog<-function(PATH, RADIUS, LENGTH, VOLUME){

        # Read the input parameters
        FILENAME=paste0(PATH,"/Color.in")
        PRM<-scan(file=FILENAME)        
        IFT<-PRM[2]*5.796
        VISC<-(PRM[1]-0.5)/3
        Fx<-PRM[6]
        Fy<-PRM[7]
        Fz<-PRM[8]
	PI=3.14159265358979
	
        # Read the time history of averages 
        FILENAME=paste0(PATH,"/timelog.tcat")
        DATA<-read.csv(file=FILENAME,head=TRUE,sep=" ")

        # interfacial tension
        DATA$IFT<-IFT
        # dynamic viscosity
        DATA$viscosity<-VISC
	
	DATA$Case<-PATH

	DATA$awn<-DATA$awn*VOLUME
	DATA$aws<-DATA$aws*VOLUME
	DATA$ans<-DATA$ans*VOLUME

	DATA$U <- predict(sm.spline(DATA$time,DATA$sw),DATA$time,1)*LENGTH
	DATA$dswdt <- predict(sm.spline(DATA$time,DATA$sw),DATA$time,1)
	DATA$Ca<-(-1)*DATA$vawnsz*VISC/IFT
	DATA$dawndt <-predict(sm.spline(DATA$time,DATA$awn),DATA$time,1)
#	DATA$dansdt <-predict(sm.spline(DATA$time,DATA$ans),DATA$time,1)*VOLUME/(2*PI*RADIUS)
	DATA$dansdt <-predict(sm.spline(DATA$time,DATA$ans),DATA$time,1)/(2*PI*RADIUS)
	DATA$pc<-(DATA$pn-DATA$pw)*RADIUS/IFT	
	
	return(DATA)
}

mainfont="Helvetica"

# additional definitions
dswdt=expression(L~frac(ds^{bar(bar(w))},dt))
dawndt=expression(frac(d~epsilon^{bar(bar(wn))},dt))
dansdt=expression(frac(V,2~pi~R)~frac(d~epsilon^{bar(bar(ws))},dt))
LHSA=expression(frac((p^{n}-p^{w})~R,gamma^{wn})-J[w]^{wn}~R)
RHSAa=expression(frac(eta~R,gamma^{wn})~frac(d~s^{bar(bar(w))},dt))
RESA=expression(frac((p^{n}-p^{w})~R,gamma^{wn})-J[w]^{wn}~R-c^{wn}~frac(eta~R,gamma^{wn})~frac(d~s^{bar(bar(w))},dt))
RHSAb=expression(frac(gamma^{wn}~(epsilon^{bar(bar(wn))}~R-epsilon[eq]^{bar(bar(wn))}~R),p^{wn}~R))
RHSB=expression(frac(1,2)~G[zz]^{wn}~(~v[z]^{bar(w)}~+~v[z]^{bar(n)}))

# define the TCAT variables, length scale and IFT
source("~/Programs/LBPM-WIA/example/R/DefsTCAT.R")
gamma = 0.058
D=18
L=400
VOLUME=40*40*400
PI=3.14159265
Jeq=1.711

C1<-ReadPistonTimelog("SET1/Case1",D,L,VOLUME)
C1$Set<-"A"
C1$cwnseq<-0
C2<-ReadPistonTimelog("SET1/Case2",D,L,VOLUME)
C2$Set<-"A"
C2$cwnseq<-0
C3<-ReadPistonTimelog("SET1/Case3",D,L,VOLUME)
C3$Set<-"A"
C3$cwnseq<-0
C4<-ReadPistonTimelog("SET1/Case4",D,L,VOLUME)
C4$Set<-"A"
C4$cwnseq<-0
C5<-ReadPistonTimelog("SET1/Case5",D,L,VOLUME)
C5$Set<-"A"
C5$cwnseq<-0
C6<-ReadPistonTimelog("SET1/Case6",D,L,VOLUME)
C6$Set<-"A"
C6$cwnseq<-0

C7<-ReadPistonTimelog("SET2/Case1",D,L,VOLUME)
C7$Set<-"B"
C7$cwnseq<-(-0.4328)
C8<-ReadPistonTimelog("SET2/Case2",D,L,VOLUME)
C8$Set<-"B"
C8$cwnseq<-(-0.4328)
C9<-ReadPistonTimelog("SET2/Case3",D,L,VOLUME)
C9$Set<-"B"
C9$cwnseq<-(-0.4328)
C10<-ReadPistonTimelog("SET2/Case4",D,L,VOLUME)
C10$Set<-"B"
C10$cwnseq<-(-0.4328)
C11<-ReadPistonTimelog("SET2/Case5",D,L,VOLUME)
C11$Set<-"B"
C11$cwnseq<-(-0.4328)
C12<-ReadPistonTimelog("SET2/Case6",D,L,VOLUME)
C12$Set<-"B"
C12$cwnseq<-(-0.4328)

C13<-ReadPistonTimelog("SET3/Case1",D,L,VOLUME)
C13$Set<-"C"
C13$cwnseq<-(-0.83656)
C14<-ReadPistonTimelog("SET3/Case2",D,L,VOLUME)
C14$Set<-"C"
C14$cwnseq<-(-0.83656)
C15<-ReadPistonTimelog("SET3/Case3",D,L,VOLUME)
C15$Set<-"C"
C15$cwnseq<-(-0.83656)
C16<-ReadPistonTimelog("SET3/Case4",D,L,VOLUME)
C16$Set<-"C"
C16$cwnseq<-(-0.83656)
C17<-ReadPistonTimelog("SET3/Case5",D,L,VOLUME)
C17$Set<-"C"
C17$cwnseq<-(-0.83656)
C18<-ReadPistonTimelog("SET3/Case6",D,L,VOLUME)
C18$Set<-"C"
C18$cwnseq<-(-0.83656)

# Read the equilibrium cases
#   Eq1 matches the contact angle for SET1
#   Eq2 matches the contact angle for SET2
#   Eq3 matches the contact angle for SET3
# Note the trapped bubble has two interfaces, piston case has one! 
#   (2x interfacial area / common curve length)
Eq1<-ReadPistonTimelog("EQ/Case1",D,160,160*40*40)
Eq2<-ReadPistonTimelog("EQ/Case2",D,160,160*40*40)
Eq3<-ReadPistonTimelog("EQ/Case3",D,160,160*40*40)
Eq4<-ReadPistonTimelog("EQ/Case4",D,160,160*40*40)

#Bind all of the data into frames by case
SET1<-rbind(C1,C2,C3,C4,C5,C6)
SET2<-rbind(C7,C8,C9,C10,C11,C12)
SET3<-rbind(C13,C14,C15,C16,C17,C18)

# Assign equilibrium variables for each set
SET1$Set<-"A"
SET1$awneq<-Eq1$awn[NROW(Eq1$awn)]*0.5
SET1$Jwneq<-Eq1$trJwn[NROW(Eq1$trJwn)]
SET1$cwnseq<-Eq1$cwns[NROW(Eq1$cwns)]

SET2$Set<-"B"
SET2$awneq<-Eq2$awn[NROW(Eq2$awn)]*0.5
SET2$Jwneq<-Eq2$trJwn[NROW(Eq2$trJwn)]
SET2$cwnseq<-Eq2$cwns[NROW(Eq2$cwns)]

SET3$Set<-"C"
SET3$awneq<-Eq3$awn[NROW(Eq3$awn)]*0.5
SET3$Jwneq<-Eq3$trJwn[NROW(Eq3$trJwn)]
SET3$cwnseq<-Eq3$cwns[NROW(Eq3$cwns)]

# Crop out measurements that are:
#   (1) early in time due to non-equilibrium starting condition
#   (2) close to the exit boundary

model<-lm(formula = cwns - cwnseq ~ Ca, data = SET2)
slope<-summary(model)$coefficients[2,1]

SET1<-subset(SET1,SET1$time>20000)
SET1<-subset(SET1,SET1$sw>0.1)

SET2<-subset(SET2,SET2$time>20000)
SET2<-subset(SET2,SET2$sw>0.1)

SET3<-subset(SET3,SET3$time>20000)
SET3<-subset(SET3,SET3$sw>0.1)

# Bind the full data set
Full<-rbind(SET1,SET2,SET3)

Full$Ca<-Full$Ca*(-1) 
Full$vawnz<-Full$vawnz*(-1)
Full$vawnsz<-Full$vawnsz*(-1)
Full$trJwn<-D*Full$trJwn
Full$U<-Full$U*(-1)
#Full$awn<-Full$awn/D
#Full$aws<-Full$aws/D
#Full$ans<-Full$ans/D
#Full$lwns<-Full$lwns/D/D
Full$RHSA<-Full$pc*400/D-Full$trJwn*400/D

Full$LHSeq<-Full$pc-Full$trJwn
Full$LHS<-Full$pc-Full$trJwn
Full$NonDimA<-Full$Ca*D/L
Full$NonDimB<-(Full$awn-Full$awneq)*D/Full$pc
Full$ResA<-Full$LHS-2808.1693*Full$NonDimA+0.1038

myfitA<-lm(NonDimA~LHS+NonDimB,data=Full)
summary.lm(myfitA)
Err <- summary(myfitA)$coefficients[1,1]
cwn <- summary(myfitA)$coefficients[2,1]
kwn <- summary(myfitA)$coefficients[3,1]

Full$RHSB<-0.5*Full$Gwnzz*(Full$vawz+Full$vanz)

# Dynamic contact angle - linear function of the capillary number
MOD1<-lm(formula = cwns ~ Ca, data = SET1)
INT1<-summary(MOD1)$coefficients[1,1]
SLP1<-summary(MOD1)$coefficients[2,1]
MOD2<-lm(formula = cwns ~ Ca, data = SET2)
INT2<-summary(MOD2)$coefficients[1,1]
SLP2<-summary(MOD2)$coefficients[2,1]
MOD3<-lm(formula = cwns  ~ Ca, data = SET3)
INT3<-summary(MOD3)$coefficients[1,1]
SLP3<-summary(MOD3)$coefficients[2,1]

Eq1<-subset(Eq1,Eq1$time>40000)
Eq2<-subset(Eq2,Eq2$time>40000)
Eq3<-subset(Eq3,Eq3$time>40000)

p<-ggplot()+geom_point(data=SET1,aes(Ca,cwns,colour="SET1",shape="Dynamic"))+
	geom_point(data=SET2,aes(Ca,cwns,colour="SET2",shape="Dynamic"))+
	geom_point(data=SET3,aes(Ca,cwns,colour="SET3",shape="Dynamic"))+
	geom_abline(intercept=INT1,slope=SLP1,aes(colour="SET1"))+
	geom_abline(intercept=INT2,slope=SLP2,aes(colour="SET2"))+
	geom_abline(intercept=INT3,slope=SLP3,aes(colour="SET3"))+
	ylab(expression(cos~varphi^{bar(bar(wns))}))+
	guides(colour=FALSE)+
	geom_point(data=Eq1,aes(Ca,cwns,colour="SET1",shape="Equilibrium"),size=3)+
	geom_point(data=Eq2,aes(Ca,cwns,colour="SET2",shape="Equilibrium"),size=3)+
	geom_point(data=Eq3,aes(Ca,cwns,colour="SET3",shape="Equilibrium"),size=3)+
	scale_shape_manual(name="Simulation",values=c(20,4))

ggsave("cwns-dynamic.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(LHS,NonDimA,colour=Set))+
	geom_point() + ylab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + 
	xlab(C^{wn}~Delta~P^{c}) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=cwn,intercept=Err,colour="gray")
						      
ggsave("piston-dswdt-nondim.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,NonDimA,colour=Case))+
	geom_line() + xlab("time") + ylab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + 
	theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
 
ggsave("piston-dswdt-nondim-A.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,Err+cwn*LHS,colour=Case))+
	geom_line() + xlab("time") + 
	ylab(expression(c^{wn~paste("*")}~(p^{wn~paste("*")}-p^{c~paste("*")})+epsilon^{~paste("*")})) + 
	theme_bw() + theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
 
ggsave("piston-dswdt-nondim-B.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,kwn*NonDimB,colour=Case))+
	geom_line() + xlab("time") + 
	ylab(expression(Delta~epsilon^{bar(bar(wn))~paste("*")}/p^{c~paste("*")})) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
  
ggsave("piston-dswdt-nondim-C.pdf",p,width=4.5,height=3.5)


Steady<-subset(Full,time>37000 & time<39000,select=c(RHSA,Ca,vawnz,RHSB))
fitSteady<-lm(Ca~RHSA,data=Steady)

# Plot the saturation against time
p<-ggplot(Full,aes(time,sw,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(sw) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-sw.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,awn,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewn) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-ewn.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,ans,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ens) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-ens.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,aws,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ews) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-ews.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,lwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-ewns.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,vawz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))
					      
ggsave("piston-vaw.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,vanz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vanz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
  
ggsave("piston-van.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,vawnz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-vawn.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,vawnsz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnsz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-vawns.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,cwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(cwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-cwns.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,KNwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KNwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-KNwns.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,KGwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KGwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-KGwns.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,U,colour=Case)) +
	geom_line() + xlab("time") + ylab(dswdt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))

ggsave("piston-dswdt.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,Gwnxx,colour=Case))+
	geom_line() + xlab("time") + ylab(Gwnxx) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-Gwnxx.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,Gwnzz,colour=Case))+
	geom_line() + xlab("time") + ylab(Gwnzz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-Gwnzz.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,Ca,colour=Case))+
	geom_line() + xlab("time") + ylab("Ca") + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-Ca.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,trJwn,colour=Case))+
	geom_line() + xlab("time") + ylab(JwnD) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))

						      
ggsave("piston-Jwn.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,dansdt,colour=Case))+
	geom_line() + xlab("time") + ylab(dansdt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-dansdt.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(time,dawndt,colour=Case))+
	geom_line() + xlab("time") + ylab(dawndt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("piston-dawndt.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(NonDimA,LHS,colour=Case))+
	geom_point() + xlab(RHSAa) + ylab(LHSA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3))) 
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("piston-dswdt-nondim-A.pdf",p,width=4.5,height=3.5)

#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("piston-dswdt-nondim-A-eq.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(NonDimA,LHSeq,colour=Case))+
	geom_point() + xlab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + ylab(expression(p^{paste(wn,"*")}-p^{paste(c,"*")})) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=2917.1011,intercept=-0.1527,colour="gray")
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("piston-dswdt-nondim.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(NonDimB,LHS,colour=Case))+
	geom_point() + xlab(RHSAb) + ylab(LHSA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3))) 
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("piston-dswdt-nondim-B.pdf",p,width=4.5,height=3.5)

p<-ggplot(Full,aes(NonDimB,ResA,colour=Case))+
	geom_point() + xlab(RHSAb) + ylab(RESA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))

#ggsave("piston-dswdt-residual-A.pdf",p,width=4.5,height=3.5)

#p<-ggplot(subset(Full,time>37000 & time<39000),aes(RHSB,vawnz,colour=Case))+
p<-ggplot(Full,aes(RHSB,vawnz,colour=Set))+
	geom_point() + xlab(RHSB) + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=1,intercept=0,colour="gray")
#	stat_smooth(aes(group=1),method=lm,fullrange=TRUE,se=FALSE,colour="gray")

ggsave("piston-Gwn-vawn.pdf",p,width=4.5,height=3.5)

#p<-ggplot(SSA,aes(Ca,RHSA))+
#	geom_point() + xlab(LHSA) + ylab(RHSA) + theme_bw() +
#	theme(text=element_text(family=mainfont,size=rel(3))) +
#	theme(legend.text=element_text(size=rel(4))) +
#	coord_cartesian(xlim=c(-0,52)) +
#	theme(legend.title=element_text(size=rel(3)))
						      
ggplot(data=Full,aes(x=NonDimA,y=cwn*LHS+kwn*NonDimB+Err,colour=Set))+
	geom_point()+
	xlab("Measured Value")+
	ylab("Predicted Value")