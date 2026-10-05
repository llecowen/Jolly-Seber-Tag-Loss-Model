#*************************************************************************************
######################################################################################
#Paramater-Specific Design Matrix Specification:
#User must identify "model"
######################################################################################
#*************************************************************************************
model=3;

nstd_survival=(nsample-1)*ngroup;
nstd_capture=nsample*ngroup;
nstd_tagloss=sum((nsample-1):1)*ngroup;
nstd_entry=nsample*ngroup;
nstd_bstar=nstd_entry;
#nstd_release=nsample*ngroup

nstd=sum(nstd_survival,nstd_capture,nstd_tagloss,nstd_entry)#nstd_release is also summed into nstd if it is accounted for by the model


if(model==1 || model==4)# parameters vary across groups (except bstars, they vary across sample times only)
{ survival.design<- survival.design_1 
  survival.design<- as.matrix(survival.design);
  
  capture.design<- capture.design_1
  capture.design<- as.matrix(capture.design);
  
  tagloss.design<- tag.design_1
  tagloss.design<- as.matrix(tagloss.design);
  
  bstar.design<- bstar_1
  bstar.design<- as.matrix(bstar.design);
  
  nbeta_survival=ncol(survival.design);
  nbeta_capture=ncol(capture.design);
  nbeta_tagloss=ncol(tagloss.design);
  nbeta_bstar=ncol(bstar.design);
  nbeta_entry=nbeta_bstar;
}

if(model==2 || model==5)#parameters vary across groups and sample times (bstars vary only across sample times)
{ survival.design<- survival.design_2 
  survival.design<- as.matrix(survival.design);
  
  capture.design<- capture.design_2
  capture.design<- as.matrix(capture.design);
  
  tagloss.design<- tag.design_2
  tagloss.design<- as.matrix(tagloss.design);
  
  bstar.design<- bstar_1
  bstar.design<- as.matrix(bstar.design);
 
  nbeta_survival=ncol(survival.design);
  nbeta_capture=ncol(capture.design);
  nbeta_tagloss=ncol(tagloss.design);
  nbeta_bstar=ncol(bstar.design);
  nbeta_entry=nbeta_bstar;
}

if(model==3 || model==6)#do NOT vary across groups or sample times, (except bstars vary across sample times)
{ survival.design<- survival.design_3 
  survival.design<- as.matrix(survival.design);
  
  capture.design<- capture.design_3
  capture.design<- as.matrix(capture.design);
  
  tagloss.design<- tag.design_3
  tagloss.design<- as.matrix(tagloss.design);
  
  bstar.design<- bstar_1
  bstar.design<- as.matrix(bstar.design);
  
  nbeta_survival=ncol(survival.design);
  nbeta_capture=ncol(capture.design);
  nbeta_tagloss=ncol(tagloss.design);
  nbeta_bstar=ncol(bstar.design);
  nbeta_entry=nbeta_bstar;
}

if(model==7)#do NOT vary across groups, DO vary across sample times (for testing data#2 and data#5 w/o considering group heterogeneity)
{ survival.design<- survival.design_7; 
  survival.design<- as.matrix(survival.design);
  
  capture.design<- capture.design_7
  capture.design<- as.matrix(capture.design);
  
  tagloss.design<- tag.design_7
  tagloss.design<- as.matrix(tagloss.design);
  
  bstar.design<- bstar_1
  bstar.design<- as.matrix(bstar.design);
  
  nbeta_survival=ncol(survival.design);
  nbeta_capture=ncol(capture.design);
  nbeta_tagloss=ncol(tagloss.design);
  nbeta_bstar=ncol(bstar.design);
  nbeta_entry=nbeta_bstar;
}


nbeta= sum(nbeta_survival,nbeta_capture,nbeta_tagloss,nbeta_bstar)

######################################################################################
#Creation of the DESIGN MATRIX (from the parameter-specific design matrices)
######################################################################################
design_matrix=make_design_matrix(survival.design,capture.design,tagloss.design, bstar.design,nstd,nbeta,nstd_survival,nbeta_survival,nstd_capture,nbeta_capture,nstd_tagloss,nbeta_tagloss)


#*************************************************************************************
######################################################################################
#Initialization of standard parameters:
#std_survival= matrix (ngroup x (nsample-1))
#std_capture= matrix (ngroup x nsample)
#std_tagloss= 3D array (ngroup x (nsample-1) x (nsample-1)) -- the 2nd dimension corresponds to the first sample time at which the animal was captured and the 3rd dimension corresponds to the last sample time at which the animal was captured
######################################################################################
#*************************************************************************************
std_survival= matrix(0.85,nrow=ngroup,ncol=(nsample-1));
std_capture= matrix(0.85,nrow=ngroup,ncol=nsample);
#for(j in 1:nsample) # Don't need for JSTL
#{ std_capture[2,j]=0.8;
#}
std_tagloss= array(0.25,dim=c(ngroup,(nsample-1),(nsample-1)));
for(g in 1:ngroup)
{ for(f in 1:(nsample-1))
  { for(j in 1:(nsample-1))
  	{ if(j<f){std_tagloss[g,f,j]=NA}
  	}
  }
}
std_entry= matrix(nrow=ngroup,ncol=nsample);
for(g in 1:ngroup)
{ std_entry[g,1]=0.2;
  std_entry[g,2]=0.2;
  std_entry[g,3]=0.6;
}
std_bstar=matrix(nrow=ngroup,ncol=nsample);
for(g in 1:ngroup)
{ for(j in 1:nsample)
	{ std_bstar[g,j]=bstar(j,g,std_entry,nsample)}
}

######################################## FOR DATA#5 #######################################
if(model==2 || model == 5)
{
std_survival[1,1] = 0.8
std_survival[1,2] = 0.7
std_survival[2,1] = 0.6
std_survival[2,2] = 0.5
std_capture[1,1] = 0.75
std_capture[2,1] = 0.9
}

###########################################################################################
stdparms_initial= pack_std(std_survival,std_capture,std_tagloss,std_bstar, nstd_survival,nstd_capture,nstd_tagloss,nstd_bstar,ngroup,nsample)


######################################################################################
#Transformation of initial standard parameters to initial beta parameters (packed in one long vector)
######################################################################################
beta_initial= transform_std(design_matrix,stdparms_initial);

#*************************************************************************************
#Maximizing the likelihood:
#*************************************************************************************
flag=0
if(model==2 || model==7)
 {flag=1
 }
out=nlm(f=negloglik,p=beta_initial,first=first,last=last,capt.hist=capt.hist,tag.hist=tag.hist,numtag=numtag,ngroup=ngroup,nobs_total=nobs_total,nobs=nobs,nsample=nsample,ntag=ntag,nstd_survival=nstd_survival,nstd_capture=nstd_capture,nstd_tagloss=nstd_tagloss,flag=flag,Nhat=Nhat,iterlim=3000,hessian=TRUE);


######################################################################################
#lmax is the maximized log-likelihood for the model
#Hess is the hessian matrix required to calculate the covariance
######################################################################################
lmax=-out$minimum;
Hess=out$hessian;
######################################################################################

#The following are the estimated standard parameters
std_parms=transform_beta(design_matrix,out$estimate)

std_survival=matrix(std_parms[1:nstd_survival],nrow=ngroup,ncol=(nsample-1),byrow=TRUE);
std_survival

std_capture=matrix(std_parms[(nstd_survival+1):(nstd_survival+nstd_capture)],nrow=ngroup,ncol=nsample,byrow=TRUE);
std_capture

std_tagloss=array(NA,dim=c(ngroup,(nsample-1),(nsample-1)));
  index=(nstd_survival+nstd_capture+1);
  for(g in 1:ngroup)
  { for(f in 1:(nsample-1))
  	{ for(j in 1:(nsample-1))
  	  { if(j>=f)
  	  	{ std_tagloss[g,f,j]=std_parms[index]
  	      index=(index+1);
  	    }
  	  }
  	}
  }
std_tagloss

std_bstar=matrix(std_parms[(nstd_survival+nstd_capture+nstd_tagloss+1):length(std_parms)],nrow=ngroup,ncol=nsample,byrow=TRUE);
std_bstar

std_entry=matrix(nrow=ngroup,ncol=nsample);
  for(g in 1:ngroup)
  { for(j in 1:nsample)
  	{ std_entry[g,j]= bstar_to_b(g,j,std_bstar,nsample,ngroup)}
  }
std_entry
######################################################################################
Cov_std= COV(design_matrix,Hess,std_parms) #Covariance
Cov_std1=matrix(data=0,nrow=nstd+1,ncol=nstd+1)
Cov_std1[1:nstd,1:nstd]=Cov_std[,1:nstd]   ##Add Nhat into the parameters
Var_std= VAR(design_matrix,Hess,std_parms)
######################################################################################
std.err=sqrt(Var_std)
std.err
######################################################################################
AIC=calc_AIC(design_matrix,lmax,nbeta)
AIC
lmax
nbeta
######################################################################################
N_est = matrix(data=NA,nrow=ngroup,ncol=nsample);
for(g in 1:ngroup){
	for(j in 1:nsample){
	  N_est[g,j] = N_est_func(std_survival,std_capture,std_entry,j,g,ngroup,nsample,nobs)
      
	}
}
N_est

P_0 = P_0_func(std_survival,std_capture, std_tagloss,std_entry,ngroup,nsample)
Nhat = vector (mode="double",length=ngroup)
for (g in 1:ngroup)
 {Nhat[g]=nobs[g]/(1-P_0[g])
 }
Nhat


