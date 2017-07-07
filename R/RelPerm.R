library(ggplot2)
library(pspline)
library(reshape2)
require(RColorBrewer)
require(mgcv)

ComputeRelPerm<-function(DATA, PERM, POROSITY, D32){
    
    #DATA<-data.frame(as.list(colMeans(DATA)))
    DATA<-tail(DATA,1)
    
    # retain the assigned sauter mean and porosity
#    D32<-POROSITY*D32/(1-POROSITY)
 #   PERM<-PERM*(1-POROSITY/POROSITY)^2
    DATA$D32<-D32
    DATA$Porosity<-POROSITY
    # LBM force and dimensionless force
    DATA$Fo<-DATA$Fz*D32*D32*D32/(DATA$viscosity*DATA$viscosity)
    
    # Reynolds number computed for each phase
    DATA$Rew<-DATA$vawz*D32/DATA$viscosity
    DATA$Ren<-DATA$vanz*D32/DATA$viscosity
    
    #Relative permeabilities
    DATA$knr<-(DATA$Ren*(1-DATA$sw)/DATA$Fo)/PERM
    DATA$kwr<-(DATA$Rew*DATA$sw/DATA$Fo)/PERM
    
    # Volumetric flow rate
    #DATA$Qwz<-DATA$vawz*DATA$sw*CrossSectionArea*Porosity
    #DATA$Qnz<-DATA$vanz*(1-DATA$sw)*CrossSectionArea*Porosity
    #DATA$Qz<-(DATA$vawz*DATA$sw+DATA$vanz*(1-DATA$sw))*CrossSectionArea*Porosity
    
    # Capillary number computed from volumetric flow rate
    DATA$Ca<-DATA$vawz*DATA$viscosity/DATA$IFT
    
    # fractional flow rate
    #DATA$fw<-DATA$vawz*DATA$sw/(DATA$vawz*DATA$sw+DATA$vanz*(1-DATA$sw))
    
    # capillary pressure and phase pressure difference
    DATA$pc<- D32*DATA$trJwn
    DATA$pcs<-D32/((1/DATA$trJwn)-1.0)
    DATA$pwn<-(DATA$pn-DATA$pw)*D32/DATA$IFT
    DATA$Rawn<-(2.0/DATA$trJwn)
    DATA$awnr<-DATA$awn*(DATA$Rawn-2.0)^2/(DATA$Rawn^2)

    DATA$An<-DATA$An*D32
    
    #	RESULT<-DATA[,c("time","sw","pw","pn","awn","ans","aws","Jwn","lwns","vawz",
    #			"vanz","vawnz","vawnsz","Euler","IFT","Image","viscosity",
    #			"Fz","Fo","Rew","Ren","knr","kwr","Qwz","Qnz","Qz",
    #			"Ca","fw","pc","pwn")]
    
    return(DATA)

}

RelPerm<-function(PATTERN){

    # Function to return rel perm from each path
    D32<-scan(file="D32.txt")
    # Media porosity
    Porosity<-scan(file="Porosity.txt")
    #dimionsionless permeability
    Permeability<-scan(file="K.txt")

    #############################################################
    # READ IN THE CASES
    #############################################################

    list<-list.files(pattern=PATTERN)
    Data=ReadTwoPhaseAverages(list[1])
    Data<-ComputeRelPerm(Data,Permeability,Porosity,D32)
    for (k in 2:length(list)){
        print(list[k])
        tmp=ReadTwoPhaseAverages(list[k])
        tmp<-ComputeRelPerm(tmp,Permeability,Porosity,D32)
        Data<-rbind(Data,tmp)
    }

    return(Data)
}

PlotData<-function(Data){
    p<-ggplot()+
    geom_point(data=Data,aes(x=sw,y=awn*D32,shape=Process),size=3)+
    geom_line(data=Data,aes(x=sw,y=awn*D32,linetype=Process))+
    scale_x_continuous(limits=c(0,1))+
    xlab(expression(s^{bar(bar(w))})) +
    ylab(expression(epsilon^{bar(bar(wn))}~D))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


    ggsave("awn-sw-Ca.pdf",p,width=6.5,height=4.5)

    p<-ggplot()+
    xlab(expression(s^{bar(bar(w))})) + ylab("Relative Permeability") +
    geom_point(data=Data,aes(x=sw,y=knr,colour="Oil",shape=Process),size=3) +
    geom_point(data=Data,aes(x=sw,y=kwr,colour="Water",shape=Process),size=3) +
    geom_line(data=Data,aes(x=sw,y=knr,linetype=Process,colour="Oil")) +
    geom_line(data=Data,aes(x=sw,y=kwr,linetype=Process,colour="Water")) +
    scale_x_continuous(limits=c(0,1))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_colour_manual(name="Fluid",values=c("#FF3300","#3399FF"))

    ggsave("RelPerm.pdf",p,width=6.5,height=4.5)

    p<-ggplot()+
    geom_point(data=Data,aes(x=sw,y=Euler*(D32)^3,shape=Process),size=3)+
    geom_line(data=Data,aes(x=sw,y=Euler*D32^3,linetype=Process))+
    xlab(expression(s^{bar(bar(w))})) +
    ylab(expression(chi[n]~D^3))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_x_continuous(limits=c(0,1))

    ggsave("Euler-sw.pdf",p,width=6.5,height=4.5)

    p<-ggplot()+
    geom_point(data=Data,aes(x=sw,y=pwn,shape=Process),size=3)+
    geom_line(data=Data,aes(x=sw,y=pwn,linetype=Process))+
    scale_x_continuous(limits=c(0,1))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    xlab(expression(s^{bar(bar(w))})) +
    ylab(expression(p[n]~-~p[w]~paste(")")~D~paste("/")~gamma[wn]))

    ggsave("pwn-sw.pdf",p,width=6.5,height=4.5)

    p<-ggplot()+
    geom_point(data=Data,aes(x=sw,y=pc,shape=Process),size=3)+
    geom_line(data=Data,aes(x=sw,y=pc,linetype=Process))+
    scale_x_continuous(limits=c(0,1))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    xlab(expression(s^{bar(bar(w))})) +
    ylab(expression(J[w]^{wn}~paste(")")~D))

    ggsave("pc-sw.pdf",p,width=6.5,height=4.5)

}