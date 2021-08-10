#Simulations for spatial model
rm(list=ls())
library(TMB)
library(sp)
library(INLA)
library(raster)
library(ggplot2)
library(dplyr)
library(grid)
logit<-function(x) log(x/(1-x))
invlogit<-function(x) exp(x)/(1+exp(x))

nknot<-25

matYear<-c(rep(1999,nknot),rep(2000,nknot),
           rep(2001,nknot),rep(2002,nknot),rep(2003,nknot),
           rep(2004,nknot),rep(2005,nknot),rep(2006,nknot),
           rep(2007,nknot),rep(2008,nknot),rep(2009,nknot),
           rep(2010,nknot),rep(2011,nknot),rep(2012,nknot),
           rep(2013,nknot),rep(2014,nknot),rep(2015,nknot),
           rep(2016,nknot),rep(2017,nknot),rep(2018,nknot),rep(2019,nknot))

matYear1<-c(rep(1999,nknot),rep(2000,nknot),
            rep(2001,nknot),rep(2002,nknot),rep(2003,nknot),
            rep(2004,nknot),rep(2005,nknot),rep(2006,nknot),
            rep(2007,nknot),rep(2008,nknot),rep(2009,nknot),
            rep(2010,nknot),rep(2011,nknot),rep(2012,nknot),
            rep(2013,nknot),rep(2014,nknot),rep(2015,nknot),
            rep(2016,nknot),rep(2017,nknot),rep(2018,nknot))

#Create a 100X100km square
x_coord<-c(175,225,225,175)
y_coord<-c(175,175,225,225)
xy<-cbind(x_coord,y_coord)
p = Polygon(xy)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

#Randomly sample number of locations for tow locations
NY<-20
yearly_tows<-125
set.seed(30)
random_loc<-spsample(sps,NY*yearly_tows,"random")
#Setting up knots based on these locations
# nknot<-25
set.seed(26)
# set.seed(110) #for 40 years
knots<-kmeans(cbind(random_loc@coords[,1],random_loc@coords[,2]),nknot,nstart=25)
knots.loc<-as.data.frame(knots[[2]])
knotID<-knots$cluster

#Create the mesh
utm.prj4s <- CRS("+init=epsg:32619")
sp_knots<-SpatialPoints(knots.loc,proj4string = utm.prj4s)
mesh<- inla.mesh.2d(sp_knots, 
                    max.edge=c(5,30),cutoff=2, #original is 30 and 1.5
                    boundary=inla.sp2segment(sps))
#Create anisotropic spde object for model
spde<-inla.spde2.matern(mesh)
# Triangle info
Dset = 1:2
TV = mesh$graph$tv           # Triangle to vertex indexing
V0 = mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = mesh$loc[TV[,2],Dset]
V2 = mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0
# Calculate Areas
TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
spde_aniso <- list(
  "n_s"      = spde$n.spde,
  "n_tri"    = nrow(TV),
  "Tri_Area" = Tri_Area,
  "E0"       = E0,
  "E1"       = E1,
  "E2"       = E2,
  "TV"       = TV - 1,
  "G0"       = spde$param.inla$M0,
  "G0_inv"   = as(diag(1/diag(spde$param.inla$M0)), "dgTMatrix"))

#To create area, must create grid and attribute each cell to closest knot
grid<-raster(extent(sps),resolution=c(1,1))
grid<-raster::extend(grid,c(1,1))
gridpoly<-rasterToPolygons(grid)
clipgrid<-raster::intersect(gridpoly,sps)

centlon<-c()
centlat<-c()
centsize<-c()
for (i in 1:length(clipgrid@polygons)) {
  centlon<-c(centlon,clipgrid@polygons[[i]]@labpt[1])
  centlat<-c(centlat,clipgrid@polygons[[i]]@labpt[2])
  centsize<-c(centsize,clipgrid@polygons[[i]]@area)
}
centroids<-cbind(centlon,centlat)
distmat<-pointDistance(centroids,sp_knots@coords,lonlat=F)
polyknotID<-c()
for (i in 1:length(clipgrid@polygons)) {
  polyknotID<-c(polyknotID,which(distmat[i,]==min(distmat[i,])))
}
centroids<-cbind(centroids,polyknotID,centsize,c(as.character(1:length(clipgrid@polygons))))
centroids<-as.data.frame(centroids)
centroids$knotID<-as.vector(centroids$knotID)
colnames(centroids)<-c("lon","lat","knotID","area","id")
centroidsID<-centroids[,c(3,5)]
centroids$area<-as.numeric(as.character(centroids$area))
stratarea<-aggregate(area~knotID,sum,data=centroids)

#Setup for plotting purposes
clipfort<-fortify(clipgrid)
plotloc<-as.data.frame(sp_knots@coords)
plotloc$group<-rep(1,length(plotloc[,1]))
plotloc$knotID<-as.factor(1:nknot)
plotloc<-left_join(plotloc,stratarea,by="knotID")
stratplot<-left_join(clipfort,centroidsID,by="id")

ggplot() +geom_polygon(data=stratplot,aes(x=long,y=lat,group=group,fill=knotID),
                       color=NA)+
  scale_fill_grey()+
  xlab("Easting")+ylab("Northing")+
  theme(legend.position="none")

#Now, can go into model since the objects are created
dyn.unload(dynlib("SPA_3_SEBDAM"))
compile("SPA_3_SEBDAM.cpp")
dyn.load(dynlib("SPA_3_SEBDAM"))

data<-list()
data$I<-rep(1,NY*yearly_tows)
data$IR<-rep(1,NY*yearly_tows)
data$C<-matrix(rep(1,(NY+1)*nknot),ncol=(NY+1),nrow=nknot)
data$n_bin<-rep(1,NY*yearly_tows)
data$L<-rep(1,NY*yearly_tows)
data$area<-stratarea$area

#Indices
data$n_i<-length(data$I)
data$n_t<-NY
data$n_s<-nknot
data$n_m<-mesh$n

#Factors
#If even sampling
data$s_i<-rep(rep((1:nknot)-1,yearly_tows/nknot),NY)
#If random, so doesn't have to be the same effort everywhere
# data$s_i<-sample(1:nknot,length(data$I),replace=T)
data$t_i<-rep(0:(NY-1),each=yearly_tows)
data$v_i =mesh$idx$loc-1

data$n_tows<-rep(1,NY)
data$pos_tows_I<-rep(1,NY)
data$pos_tows_IR<-rep(1,NY)

#Growth rates
data$gI<-rep(1,NY)
data$gR<-rep(1,NY)

#Mesh and spatial structure
data$spde_aniso<-spde_aniso

#Include prior
data$prior_q<-1

truepar<-list()
truepar$log_sigma_epsilon<-log(0.5)
truepar$log_sigma_upsilon<-log(0.5)
truepar$log_S<-log(0.6)
set.seed(35)
truepar$log_qI<-log(rbeta(data$n_s,10,12))
# truepar$log_qI<-log(0.4)
truepar$log_qR<-log(0.2)
truepar$logit_p_I<-rep(logit(0.95),1)
truepar$logit_p_IR<-rep(logit(0.5),1)
truepar$log_B0<-log(600)
truepar$log_R0<-log(100)
truepar$log_m0<-log(0.3)

#Spatial true parameters
Range_B<-40
Range_R<-60
Range_m<-20
truepar$log_kappa_B<-log(sqrt(8)/Range_B)
truepar$log_kappa_R<-log(sqrt(8)/Range_R)
truepar$log_kappa_m<-log(sqrt(8)/Range_m)

Sigma_B<-0.2
Sigma_R<-0.2
Sigma_m<-0.1

truepar$log_tau_B<-log(1/((Sigma_B)*sqrt(4*pi)*exp(truepar$log_kappa_B)))
truepar$log_tau_R<-log(1/((Sigma_R)*sqrt(4*pi)*exp(truepar$log_kappa_R)))
truepar$log_tau_m<-log(1/((Sigma_m)*sqrt(4*pi)*exp(truepar$log_kappa_m)))

#For anisotropy, simply use the ones obtained from model fit for simplicity
truepar$log_H_input_B<-c(-1.2,-0.5)
truepar$log_H_input_m<-c(-0.5,0.4)

#Random effects
truepar$omega_B<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))
truepar$omega_R<-matrix(rep(0,mesh$n*(data$n_t)),ncol=data$n_t)
truepar$omega_m<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))

truepar$pred_B<-rep(log(3000),data$n_s)
truepar$pred_m<-rep(log(0.1),data$n_s)

random<-c("omega_R","omega_m","omega_B")
maps = list()

obj <- MakeADFun( data=data, parameters=truepar,random=random, map = maps, DLL="SPA_3_SEBDAM" )
simdata <- obj $ simulate(complete=T)
max(simdata$B)
# check<-checkConsistency(obj,n=150)
# check
# checksum<-summary(check)
# checksum

# dyn.unload(dynlib("SPA_3_SEBDAM"))
# compile("SPA_3_SEBDAM.cpp")
# dyn.load(dynlib("SPA_3_SEBDAM"))

par<-list()
par$log_sigma_epsilon<--1
par$log_sigma_upsilon<--1
par$log_S<--1
# par$log_qI<--1
par$log_qI<-rep(-1,data$n_s)
par$log_qR<--1
par$logit_p_I<-rep(logit(0.5),1)
par$logit_p_IR<-rep(logit(0.5),1)
par$log_B0<-log(3000)
par$log_R0<-log(350)
par$log_m0<-log(0.1)
par$log_kappa_B<-0
par$log_kappa_R<-0
par$log_kappa_m<-0
par$log_tau_B<-0
par$log_tau_R<-0
par$log_tau_m<-0
par$log_H_input_B<-c(0,0)
par$log_H_input_m<-c(0,0)
par$omega_B<-matrix(rep(0,mesh$n*(data$n_t+1)),ncol=data$n_t+1)
par$omega_R<-matrix(rep(0,mesh$n*(data$n_t)),ncol=data$n_t)
par$omega_m<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))
par$pred_B<-rep(log(3000),data$n_s)
par$pred_m<-rep(log(0.1),data$n_s)

nrep<-200
Report_list<-list()
sim_save<-list()
diff_list<-list()
diff_list2<-list()
rep_list<-list()
msg_list<-list()
for (i in 1:nrep) {
  tryCatch({
    simdata <- obj $ simulate(complete=T)
    sim_save[[i]]<-simdata
    par$log_B0<-log(max(simdata$I,na.rm=T)*5)
    par$pred_B<-rep(log(max(simdata$I,na.rm=T)*5),data$n_s)
    obj2 <- MakeADFun ( simdata , par , random=random, map=maps )
    Opt2 <- try(nlminb(start=obj2$par,obj=obj2$fn,gr=obj2$gr,
                       control=list(eval.max=1000,iter.max=1000),silent=T),T)
    rep_list[[i]] <- sdreport(obj2,bias.correct=F)
    if (!is.null(rep_list[[i]])){
      Report_list[[i]] = obj2$report()
      diff_list[[i]]<-data.frame(diff_B=Report_list[[i]]$areaB-simdata$areaB,
                                 diff_R=Report_list[[i]]$areaR-simdata$areaR,
                                 diff_m=Report_list[[i]]$m-simdata$m)
      msg_list[[i]] <- Opt2$message
    }
  }, error=function(e) {})
}

converge<-0
false<-0
x_conv<-0
singular<-0
null_vec<-nrep
chosen_ones<-c()
for (i in 1:nrep){
  if (!is.null(msg_list[[i]])){null_vec<-null_vec-1
  if (msg_list[[i]]=="relative convergence (4)") {converge<-converge+1
  chosen_ones<-c(chosen_ones,i)}
  if (msg_list[[i]]=="false convergence (8)") {false<-false+1}
  if(msg_list[[i]]=="singular convergence (7)"){singular<-singular+1}
  if(msg_list[[i]]=="both X-convergence and relative convergence (5)"){x_conv<-x_conv+1
  chosen_ones<-c(chosen_ones,i)}
  }
}
conv_frame<-data.frame(converge=converge,false=false,singular=singular,x=x_conv)

#parameters
sigma_epsilon<-c(rep(NA,nrep))
sigma_upsilon<-c(rep(NA,nrep))
S<-c(rep(NA,nrep))
qI<-matrix(rep(NA,nrep*nknot),ncol=nknot)
qR<-c(rep(NA,nrep))
p_I<-c(rep(NA,nrep))
p_IR<-c(rep(NA,nrep))
B0<-c(rep(NA,nrep))
R0<-c(rep(NA,nrep))
m0<-c(rep(NA,nrep))
kappa_B<-c(rep(NA,nrep))
kappa_R<-c(rep(NA,nrep))
kappa_m<-c(rep(NA,nrep))
tau_B<-c(rep(NA,nrep))
tau_m<-c(rep(NA,nrep))
tau_R<-c(rep(NA,nrep))
log_H_input_B1<-c(rep(NA,nrep))
log_H_input_B2<-c(rep(NA,nrep))
log_H_input_m1<-c(rep(NA,nrep))
log_H_input_m2<-c(rep(NA,nrep))
SigmaO_B<-rep(NA,nrep)
SigmaO_R<-rep(NA,nrep)
SigmaO_m<-rep(NA,nrep)
# so its summary(rep_list[[i]])[(25+3*(data$n_m*data$n_t)+2*(data$n_s*(data$n_t+1))),]
#to start the parameters
for (i in chosen_ones){
  sigma_epsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_epsilon")]
  sigma_upsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_upsilon")]
  S[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="S")]
  qR[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qR")]
  p_I[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="p_I")]
  p_IR[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="p_IR")]
  B0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="B0")]
  R0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="R0")]
  m0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="m0")]
  kappa_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_B")]
  kappa_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_R")]
  kappa_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_m")]
  tau_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_B")]
  tau_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_R")]
  tau_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_m")]
  log_H_input_B1[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="log_H_input_B")][1]
  log_H_input_B2[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="log_H_input_B")][2]
  log_H_input_m1[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="log_H_input_m")][1]
  log_H_input_m2[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="log_H_input_m")][2]
  SigmaO_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_B")]
  SigmaO_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_R")]
  SigmaO_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_m")]
  # } 
}

for (i in chosen_ones){
  for (j in 1:nknot){
    qI[i,j]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qI")][j]
  }
}

qI_diffs<-matrix(rep(NA,nrep*nknot),ncol=nknot)
for (i in chosen_ones){
  qI_diffs[i,]<-qI[i,]-exp(truepar$log_qI)
}

all_params<-data.frame(sigma_epsilon,sigma_upsilon,
                       S,
                       qR,p_I,
                       p_IR,B0,R0,m0,kappa_B,
                       kappa_R,kappa_m,tau_B,tau_m,tau_R,
                       log_H_input_B1,log_H_input_B2,
                       log_H_input_m1,log_H_input_m2,
                       SigmaO_B,SigmaO_R,SigmaO_m)

#For estimq, necessary to do this to actually see things
plot_kappa_B<-kappa_B
plot_kappa_B[which(kappa_B>0.5)]<-NA
length(which(is.na(plot_kappa_B)))
plot_kappa_R<-kappa_R
plot_kappa_R[which(kappa_R>0.5)]<-NA
length(which(is.na(plot_kappa_R)))
plot_kappa_m<-kappa_m
plot_kappa_m[which(kappa_m>0.5)]<-NA
length(which(is.na(plot_kappa_m)))
plot_log_H_input_B2<-log_H_input_B2
plot_log_H_input_B2[which(abs(log_H_input_B2)>5)]<-NA
length(which(is.na(plot_log_H_input_B2)))
plot_log_H_input_m2<-log_H_input_m2
plot_log_H_input_m2[which(abs(log_H_input_m2)>5)]<-NA
length(which(is.na(plot_log_H_input_m2)))
plot_log_H_input_B1<-log_H_input_B1
plot_log_H_input_B1[which(abs(log_H_input_B1)>5)]<-NA
length(which(is.na(plot_log_H_input_B1)))
plot_log_H_input_m1<-log_H_input_m1
plot_log_H_input_m1[which(abs(log_H_input_m1)>5)]<-NA
length(which(is.na(plot_log_H_input_m1)))
plot_SigmaO_B<-SigmaO_B
plot_SigmaO_B[which(SigmaO_B>1)]<-NA
length(which(is.na(plot_SigmaO_B)))
plot_SigmaO_R<-SigmaO_R
plot_SigmaO_R[which(SigmaO_R>1)]<-NA
length(which(is.na(plot_SigmaO_R)))
plot_SigmaO_m<-SigmaO_m
plot_SigmaO_m[which(SigmaO_m>1)]<-NA
length(which(is.na(plot_SigmaO_m)))
hist_plot<-data.frame(value=c(sigma_epsilon,sigma_upsilon,
                              S,
                              qR,p_I,
                              p_IR,B0,R0,m0,
                              plot_kappa_B,plot_kappa_R,plot_kappa_m,
                              # kappa_B,kappa_R,kappa_m,
                              tau_B,tau_m,tau_R,
                              plot_log_H_input_B1,plot_log_H_input_B2,
                              plot_log_H_input_m1,plot_log_H_input_m2,
                              plot_SigmaO_B,plot_SigmaO_R,plot_SigmaO_m),
                      # log_H_input_B1,log_H_input_B2,
                      # log_H_input_m1,log_H_input_m2,
                      # SigmaO_B,SigmaO_R,SigmaO_m),
                      parameter=rep(c("sigma_epsilon","sigma_upsilon",
                                      "S",
                                      "q_R",
                                      "pI",
                                      "pIR","B0",
                                      "R0","m0",
                                      "kappa_B",
                                      "kappa_R",
                                      "kappa_m","tau_B",
                                      "tau_m","tauR",
                                      "H_input_B1","H_input_B2",
                                      "H_input_m1","H_input_m2",
                                      "SigmaO_B","SigmaO_R","SigmaO_m")
                                    ,each=nrep),
                      true=rep(c(exp(truepar$log_sigma_epsilon),
                                 exp(truepar$log_sigma_upsilon),exp(truepar$log_S),
                                 exp(truepar$log_qR),
                                 invlogit(truepar$logit_p_I),invlogit(truepar$logit_p_IR),
                                 exp(truepar$log_B0),
                                 exp(truepar$log_R0),exp(truepar$log_m0),
                                 exp(truepar$log_kappa_B),exp(truepar$log_kappa_R),
                                 exp(truepar$log_kappa_m),
                                 exp(truepar$log_tau_B),exp(truepar$log_tau_m),exp(truepar$log_tau_R),
                                 truepar$log_H_input_B[1],truepar$log_H_input_B[2],
                                 truepar$log_H_input_m[1],truepar$log_H_input_m[2],
                                 1/((exp(truepar$log_tau_B))*sqrt(4*pi)*exp(truepar$log_kappa_B)),
                                 1/((exp(truepar$log_tau_R))*sqrt(4*pi)*exp(truepar$log_kappa_R)),
                                 1/((exp(truepar$log_tau_m))*sqrt(4*pi)*exp(truepar$log_kappa_m))),each=nrep))

hist_plot$parameter<-factor(hist_plot$parameter,
                            labels=c(expression(B [0]),expression("H" ["input"]^"B1"),
                                     expression("H" ["input"]^"B2"),
                                     expression("H" ["input"]^"m1"),expression("H" ["input"]^"m2"),
                                     expression(kappa ["B"]),expression(kappa ["m"]),
                                     expression(kappa ["R"]),expression("m" [0]),
                                     expression("p" ["I"]),
                                     expression("p" ["I"]^"R"),
                                     expression("q" ["R"]),
                                     expression("R" [0]),expression("S"),
                                     expression(sigma [epsilon]),
                                     expression(sigma[upsilon]),
                                     expression(sigma[B]^Omega),
                                     expression(sigma[m]^Omega),
                                     expression(sigma[R]^Omega),
                                     expression(tau ["B"]),
                                     expression(tau ["m"]),expression(tau["R"])))


all_plot_params<-data.frame(sigma_epsilon,sigma_upsilon,
                            S,
                            qR,p_I,
                            p_IR,B0,R0,m0,
                            plot_kappa_B,plot_kappa_R,plot_kappa_m,
                            tau_B,tau_m,tau_R,
                            plot_log_H_input_B1,plot_log_H_input_B2,
                            plot_log_H_input_m1,plot_log_H_input_m2,
                            plot_SigmaO_B,plot_SigmaO_R,plot_SigmaO_m)


all_params<-data.frame(sigma_epsilon,sigma_upsilon,
                       S,
                       qR,p_I,
                       p_IR,B0,R0,m0,
                       kappa_B,kappa_R,kappa_m,
                       tau_B,tau_m,tau_R,
                       log_H_input_B1,log_H_input_B2,
                       log_H_input_m1,log_H_input_m2,
                       SigmaO_B,SigmaO_R,SigmaO_m)

theme_replace(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.5, linetype = "solid",
                                       colour = "black"),
              panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"))

plot_qI_diff<-data.frame(value=as.vector(qI_diffs),
                         knot=rep(1:25,each=nrep),
                         true=rep(exp(truepar$log_qI),each=nrep))

#Look by knot for pattern
ggplot(plot_qI_diff)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0))+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~knot)+
  xlab("Estimated Value")+ylab("Frequency")

#If no pattern, which is the case when informed, then
ggplot(plot_qI_diff)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept=0))+
  xlab("Difference between estimated and simulation value")+ylab("Frequency")

ggplot(hist_plot)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = true,col="True Value"))+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  xlab("Estimated Value")+ylab("Frequency")

diff_B<-matrix(ncol=nrep,nrow=(NY+1))
diff_R<-matrix(ncol=nrep,nrow=NY)
diff_m<-matrix(ncol=nrep,nrow=(NY+1))
perc_diff_m<-matrix(ncol=nrep,nrow=(NY+1))
for (i in chosen_ones) {
  diff_B[,i]<-(Report_list[[i]]$totB-sim_save[[i]]$totB)/sim_save[[i]]$totB
  diff_R[,i]<-(Report_list[[i]]$totR-sim_save[[i]]$totR)/sim_save[[i]]$totR
  diff_m[,i]<-(Report_list[[i]]$mean_m-sim_save[[i]]$mean_m)/sim_save[[i]]$mean_m
}
Bquant<-matrix(ncol=3,nrow=(NY+1))
Rquant<-matrix(ncol=3,nrow=NY)
mquant<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY)){
  Bquant[i,]<-quantile(diff_B[i,],probs=c(0.25,0.5,0.75),na.rm=T)  
  if (i <= NY){Rquant[i,]<-quantile(diff_R[i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant[i,]<-quantile(diff_m[i,],probs=c(0.25,0.5,0.75),na.rm=T) 
}
Bquant<-as.data.frame(Bquant)
colnames(Bquant)<-c("t","med","n")
Rquant<-as.data.frame(Rquant)
colnames(Rquant)<-c("t","med","n")
mquant<-as.data.frame(mquant)
colnames(mquant)<-c("t","med","n")

ggplot(Bquant)+geom_ribbon(aes(x=1:21,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:21,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  xlab("Year")+ylab("(Predicted - True Total Biomass)/True Total Biomass")

ggplot(Rquant)+geom_ribbon(aes(x=1:20,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:20,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - True Recruitment")
  xlab("Year")+ylab("(Predicted - True Total Recruitment)/True Total Recruitment")

ggplot(mquant)+geom_ribbon(aes(x=1:21,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:21,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  xlab("Year")+ylab("Predicted - True Mean Mortality")

knot_diff_B<-array(dim=c(nrep,NY+1,nknot))
knot_diff_R<-array(dim=c(nrep,NY,nknot))
knot_diff_m<-array(dim=c(nrep,NY+1,nknot))
for (i in which(!is.na(plot_SigmaO_B))) {
  for (j in 1:(NY+1)){
    for (k in 1:nknot){
      knot_diff_B[i,j,k]<-(Report_list[[i]]$B[k,j]-sim_save[[i]]$B[k,j])/sim_save[[i]]$B[k,j]
      knot_diff_m[i,j,k]<-(Report_list[[i]]$m[k,j]-sim_save[[i]]$m[k,j])/sim_save[[i]]$m[k,j]
      if (j <= NY) knot_diff_R[i,j,k]<-(Report_list[[i]]$R[k,j]-sim_save[[i]]$R[k,j])/sim_save[[i]]$R[k,j]
    }
  }
}
Bquant_big<-matrix(ncol=3,nrow=(NY+1))
Rquant_big<-matrix(ncol=3,nrow=NY)
mquant_big<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY+1)){
  Bquant_big[i,]<-quantile(knot_diff_B[,i,],probs=c(0.25,0.5,0.75),na.rm=T)  
  if (i <= NY){Rquant_big[i,]<-quantile(knot_diff_R[,i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant_big[i,]<-quantile(knot_diff_m[,i,],probs=c(0.25,0.5,0.75),na.rm=T) 
}
Bquant_big<-as.data.frame(Bquant_big)
colnames(Bquant_big)<-c("t","med","n")
Rquant_big<-as.data.frame(Rquant_big)
colnames(Rquant_big)<-c("t","med","n")
mquant_big<-as.data.frame(mquant_big)
colnames(mquant_big)<-c("t","med","n")

#Temporally over all knots
ggplot(Bquant_big)+geom_ribbon(aes(x=1:21,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:21,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  xlab("Year")+ylab("(Predicted - True Biomass Density)/True Biomass Density")

ggplot(Rquant_big)+geom_ribbon(aes(x=1:20,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:20,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - True Recruitment")
  xlab("Year")+ylab("(Predicted - True Recruitment)/True Recruitment")

ggplot(mquant_big)+geom_ribbon(aes(x=1:21,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:21,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  xlab("Year")+ylab("Predicted - True Mean Mortality")

load("boxplot_spat.RData")
#For total predictions
ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_B),x="1A"))+
  geom_boxplot(aes(y=as.vector(diff_B),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.3,0.4))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Biomass Density)/Simulated Biomass Density")

ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_R),x="1A"))+
  geom_boxplot(aes(y=as.vector(diff_R),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,10))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Total Recruitment)/Simulated Total Recruitment")

ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_m),x="1A"))+
  geom_boxplot(aes(y=as.vector(diff_m),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,1.5))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("Predicted - Simulated Mean Natural Mortality")

#Knot-based ones
ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_B),x="1A"))+
  geom_boxplot(aes(y=as.vector(knot_diff_B),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,8))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Biomass Density)/Simulated Biomass Density")

ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_R),x="1A"))+
  geom_boxplot(aes(y=as.vector(knot_diff_R),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,10))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Recruit Density)/Simulated Recruit Density")

ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_m_estimq),x="1A"))+
  geom_boxplot(aes(y=as.vector(knot_diff_m_estimq),x="1A"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.15,0.15))+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("Predicted - Simulated Natural Mortality")

#Plot Spatial difference median

B_mat<-array(dim=c(nrep,NY+1,nknot))
R_mat<-array(dim=c(nrep,NY,nknot))
m_mat<-array(dim=c(nrep,NY+1,nknot))
for (i in which(!is.na(qR))) {
  for (j in 1:(NY+1)){
    for (k in 1:nknot){
      B_mat[i,j,k]<-(Report_list[[i]]$B[k,j]-sim_save[[i]]$B[k,j])/sim_save[[i]]$B[k,j]
      m_mat[i,j,k]<-Report_list[[i]]$m[k,j]-sim_save[[i]]$m[k,j]
      if (j <= NY) R_mat[i,j,k]<-(Report_list[[i]]$R[k,j]-sim_save[[i]]$R[k,j])/sim_save[[i]]$R[k,j]
    }
  }
}
B_mat<-B_mat[!is.na(qR),,]
m_mat<-m_mat[!is.na(qR),,]
R_mat<-R_mat[!is.na(qR),,]
quant_spec<-function(x){quantile(x,prob=0.5,na.rm=T)}
median_Bs<-apply(B_mat, c(2,3), quant_spec)
median_Rs<-apply(R_mat,c(2,3),quant_spec)
median_ms<-apply(m_mat,c(2,3),quant_spec)

plot_medB<-data.frame(B=as.vector(median_Bs),Year=matYear,knotID=as.factor(rep(1:nknot,NY+1)))
plot_medB<-left_join(stratplot,plot_medB,by="knotID")

ggplot() +geom_polygon(data=subset(plot_medB,Year==2010),
                       aes(x=long,y=lat,group=group,fill=B),
                       color=NA)+
  scale_fill_continuous(name="Median Percent \nDifference in Predicted \nBiomass Density",high="red",low="white")+
  facet_wrap(~Year)+
  xlab("Easting")+ylab("Northing")

plot_medR<-data.frame(R=as.vector(median_Rs),Year=matYear1,knotID=as.factor(rep(1:nknot,NY)))
plot_medR<-left_join(stratplot,plot_medR,by="knotID")

ggplot() +geom_polygon(data=plot_medR,aes(x=long,y=lat,group=group,fill=R),
                       color=NA)+
  scale_fill_continuous(name="Median Percent \nDifference in Predicted \nRecruit Biomass \nDensity",high="red",low="white")+
  facet_wrap(~Year)+
  xlab("Easting")+ylab("Northing")

plot_medm<-data.frame(m=as.vector(median_ms),Year=matYear,knotID=as.factor(rep(1:nknot,NY+1)))
plot_medm<-left_join(stratplot,plot_medm,by="knotID")

ggplot() +geom_polygon(data=plot_medm,aes(x=long,y=lat,group=group,fill=m),
                       color=NA)+
  scale_fill_continuous(name="Median Difference in \nNatural Mortality",high="red",low="white")+
  facet_wrap(~Year)+
  xlab("Easting")+ylab("Northing")

