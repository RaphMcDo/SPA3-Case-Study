library(TMB)
library(dplyr)
logit<-function(x) x/(1-x)
NY<-22

#Running the TLM model first, as it is much faster

dyn.unload(dynlib("SPA3_TLM"))
compile("SPA3_TLM.cpp")
dyn.load(dynlib("SPA3_TLM"))

load("spa3_unid_TLM_data.RData")

#If don't want to use prior:
#data$prior_q<-0

par<-list()
par$log_sigma_tau<--1
par$log_sigma_phi<--1
par$log_sigma_m<--1
par$log_sigma_epsilon<--1
par$log_sigma_upsilon<--1
par$log_q_R<--1
par$log_q_I<--1
par$log_S<--1
par$logit_p_I<-logit(0.5)
par$logit_p_IR<-logit(0.5)

par$log_B<-rep(log(max(data$I,na.rm=T)*10),NY+1)
par$log_R<-rep(log(max(data$IR,na.rm=T)*10),NY+1)
par$log_m<-rep(log(0.3),NY+1)

random<-c("log_R","log_B","log_m")
maps = list()

obj_temp <- MakeADFun( data=data, parameters=par,random=random, map = maps ,DLL="SPA3_TLM")
Opt_temp <- try(nlminb(start=obj_temp$par,obj=obj_temp$fn,gr=obj_temp$gr,
                  control=list(eval.max=10000,iter.max=10000)),T)
Report_temp = obj_temp$report()
rep_temp <- sdreport(obj_temp,bias.correct=F)
Opt_temp$message

#Running the SEBDAM model

nknot<-25
mesh_n<-282

dyn.unload(dynlib("SPA_3_SEBDAM"))
compile("SPA_3_SEBDAM.cpp")
dyn.load(dynlib("SPA_3_SEBDAM"))

load("spa3_unid_data.RData")

data$prior_q<-1
#If don't want to use prior:
#data$prior_q<-0

par<-list()
par$log_sigma_epsilon<--1
par$log_sigma_upsilon<--1
par$log_S<--1
par$log_kappa_B<--1
par$log_tau_B<-1
par$log_kappa_R<--1
par$log_tau_R<-1
par$log_kappa_m<-0
par$log_tau_m<-0
# par$log_qI<--1
par$log_qI<-rep(-1,data$n_s)
par$log_qR<--1
par$log_H_input_B<-c(0,0)
par$log_H_input_R<-c(0,0)
par$log_H_input_m<-c(0,0)

par$log_R0<-log(150)
par$log_B0<-log(3000)
par$log_m0<-log(0.1)
par$logit_p_I<-logit(0.5)
par$logit_p_IR<-logit(0.5)

#Random effects
par$omega_B<-matrix(rep(0,mesh_n*((data$n_t+1))),ncol=(data$n_t+1))
par$omega_R<-matrix(rep(0,mesh_n*(data$n_t)),ncol=(data$n_t))
par$omega_m<-matrix(rep(0,mesh_n*((data$n_t+1))),ncol=(data$n_t+1))

random<-c("omega_R","omega_B","omega_m")

maps = list()

obj <- MakeADFun( data=data, parameters=par,random=random, map = maps, DLL = "SPA_3_SEBDAM")
Opt <- try(nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                  control=list(eval.max=10000,iter.max=10000)),T)
Report = obj$report()
rep <- sdreport(obj,bias.correct=F)
rep
Opt$message

matYear<-c(rep(1997:2019,each=nknot))
matYear1<-c(rep(1997:2018,each=nknot))
knots<-rep(1:nknot,NY+1)
knots1<-rep(1:nknot,NY)

library(sf)
library(shape)
library(ggplot2)

B<-data.frame(B=as.vector(Report$B),Year=matYear,
              knotID=knots,
              se=rep$sd[which(names(rep$value)=="B")])
B_plot<-left_join(gridded_bound,B,by=c("knotID"))

R<-data.frame(R=as.vector(Report$R),Year=matYear1,
              knotID=knots1,
              se=rep$sd[which(names(rep$value)=="R")])
R_plot<-left_join(gridded_bound,R,by=c("knotID"))

m<-data.frame(m=as.vector(Report$m),Year=matYear,
              knotID=knots,
              se=rep$sd[which(names(rep$value)=="m")])
m_plot<-left_join(gridded_bound,m,by=c("knotID"))

st_crs(ns_map)<-NA
ns_map$geometry<-ns_map$geometry/1000

#Spatial predictions
#B
spatial_B.plot<-ggplot()+
  geom_sf(data=B_plot,aes(fill=B),col=NA)+
  facet_wrap(~Year)+
  scale_fill_continuous(name="Predicted Biomass \nDensity (kg/km^2)",high="red",low="white")+
  geom_sf(data=ns_map)+theme_bw()+
  coord_sf(xlim=c(670,765),ylim=c(4835,4940))+
  xlab("Easting")+ylab("Northing")

#R
spatial_R.plot<-ggplot()+
  geom_sf(data=R_plot,aes(fill=R),col=NA)+
  facet_wrap(~Year)+
  scale_fill_continuous(name="Predicted Recruit \nDensity (kg/km^2)",high="red",low="white")+
  geom_sf(data=ns_map)+theme_bw()+
  coord_sf(xlim=c(670,765),ylim=c(4835,4940))+xlab("Easting")+ylab("Northing")

#m
spatial_m.plot<-ggplot()+
  geom_sf(data=m_plot,aes(fill=m),col=NA)+
  facet_wrap(~Year)+
  scale_fill_continuous(name="Predicted Natural \nMortality",high="red",low="white")+
  geom_sf(data=ns_map)+theme_bw()+
  coord_sf(xlim=c(670,765),ylim=c(4835,4940))+xlab("Easting")+ylab("Northing")

rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
Eigen_B<-eigen(Report$H_B)
Eigen_m<-eigen(Report$H_m)
Major_B<-Eigen_B$vectors[,1]*Eigen_B$values[1] * Report$Range_B
Minor_B<-Eigen_B$vectors[,2]*Eigen_B$values[2] * Report$Range_B
Major_R<-Eigen_B$vectors[,1]*Eigen_B$values[1] * Report$Range_R
Minor_R<-Eigen_B$vectors[,2]*Eigen_B$values[2] * Report$Range_R
Major_m<-Eigen_m$vectors[,1]*Eigen_m$values[1] * Report$Range_m
Minor_m<-Eigen_m$vectors[,2]*Eigen_m$values[2] * Report$Range_m
Range_B = 1.1 * c(-1,1) * max(abs( cbind(Major_B,Minor_B) ))
Range_R = 1.1 * c(-1,1) * max(abs( cbind(Major_R,Minor_R) ))
Range_m = 1.1 * c(-1,1) * max(abs( cbind(Major_m,Minor_m) ))
plot( 1, type="n", xlim=c(-450,450), ylim=c(-450,450), xlab="", ylab="")
shape::plotellipse( rx=rss(Major_B), ry=rss(Minor_B),
                    angle=-1*(atan(Major_B[1]/Major_B[2])/(2*pi)*360-90),
                    lcol=c("black")[1], lty=c("solid")[1])
shape::plotellipse( rx=rss(Major_R), ry=rss(Minor_R),
                    angle=-1*(atan(Major_R[1]/Major_R[2])/(2*pi)*360-90), 
                    lcol=c("black")[1], lty=c("dotted")[1])
shape::plotellipse( rx=rss(Major_m), ry=rss(Minor_m),
                    angle=-1*(atan(Major_m[1]/Major_m[2])/(2*pi)*360-90), 
                    lcol=c("black")[1], lty=c("twodash")[1])
mtext(side=1, outer=FALSE, line=2, text="Eastings (km)")
mtext(side=2, outer=FALSE, line=2, text="Northings (km)")

#Set up to compare total predictions with temporal model
CIB<-data.frame(B=Report_temp$B,
                SE=rep_temp$sd[which(names(rep_temp$value)=="B")]) %>% 
  mutate(low=B-(1.96*SE),high=B+(1.96*SE))
CIR<-data.frame(R=Report_temp$R,
                SE=rep_temp$sd[which(names(rep_temp$value)=="R")]) %>% 
  mutate(low=R-(1.96*SE),high=R+(1.96*SE))
CIR$low[CIR$low<0]<-0
CIm<-data.frame(m=Report_temp$m,
                SE=rep_temp$sd[which(names(rep_temp$value)=="m")]) %>% 
  mutate(low=m-(1.96*SE),high=m+(1.96*SE))
CIm$low[CIm$low<0]<-0

totals<-data.frame(B=Report$totB, m=Report$mean_m, 
                   se_B=rep$sd[which(names(rep$value)=="totB")], 
                   se_m=rep$sd[which(names(rep$value)=="mean_m")])
totR<-data.frame(R=Report$totR, 
                 se_R=rep$sd[which(names(rep$value)=="totR")])
totR[23,]<-c(NA,NA)
totals<-totals %>% cbind(totR) %>% mutate(lowB=B-(1.96*se_B),
                                          highB=B+(1.96*se_B),
                                          lowm=m-(1.96*se_m),
                                          highm=m+(1.96*se_m),
                                          lowR=R-(1.96*se_R),
                                          highR=R+(1.96*se_R))

comp_B<-ggplot()+geom_line(data=CIB,aes(x=1997:2019,y=B,col="TLM"))+
  geom_point(data=CIB,aes(x=1997:2019,y=B,col="TLM"))+
  geom_ribbon(data=CIB,aes(x=1997:2019,ymin=low,ymax=high),fill="red",alpha=0.2)+
  geom_line(data=totals,aes(x=1997:2019,y=B,col="SEBDAM"))+
  geom_point(data=totals,aes(x=1997:2019,y=B,col="SEBDAM"))+
  geom_ribbon(data=totals,aes(x=1997:2019,ymin=lowB,ymax=highB),fill="blue",alpha=0.2)+
  scale_color_manual(values=c("blue","red"),name="Model")+
  xlab("Year")+ylab("Total Predicted Biomass (metric tonnes)")

comp_R<-ggplot()+geom_line(data=CIR[-23,],aes(x=1997:2018,y=R,col="TLM"))+
  geom_point(data=CIR[-23,],aes(x=1997:2018,y=R,col="TLM"))+
  geom_ribbon(data=CIR[-23,],aes(x=1997:2018,ymin=low,ymax=high),fill="red",alpha=0.2)+
  geom_line(data=totals[-23,],aes(x=1997:2018,y=R,col="SEBDAM"))+
  geom_point(data=totals[-23,],aes(x=1997:2018,y=R,col="SEBDAM"))+
  geom_ribbon(data=totals[-23,],aes(x=1997:2018,ymin=lowR,ymax=highR),fill="blue",alpha=0.2)+
  scale_color_manual(values=c("blue","red"),name="Model")+
  xlim(c(1997,2018))+
  xlab("Year")+ylab("Predicted Total Recruitment (metric tonnes)")

comp_m<-ggplot()+geom_line(data=CIm,aes(x=1997:2019,y=m,col="TLM"))+
  geom_point(data=CIm,aes(x=1997:2019,y=m,col="TLM"))+
  geom_ribbon(data=CIm,aes(x=1997:2019,ymin=low,ymax=high),fill="red",alpha=0.2)+
  geom_line(data=totals,aes(x=1997:2019,y=m,col="SEBDAM"))+
  geom_point(data=totals,aes(x=1997:2019,y=m,col="SEBDAM"))+
  geom_ribbon(data=totals,aes(x=1997:2019,ymin=lowm,ymax=highm),fill="blue",alpha=0.2)+
  scale_color_manual(values=c("blue","red"),name="Model")+
  xlab("Year")+ylab("Mean Predicted Natural Mortality")
