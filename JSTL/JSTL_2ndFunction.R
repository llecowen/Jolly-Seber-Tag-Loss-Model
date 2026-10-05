######################################################################################
#The make_design_matrix function uses the (above) parameter-specific design matrices to create one large design matrix
######################################################################################
make_design_matrix<- function(survival.design,capture.design,tagloss.design, bstar.design,nstd,nbeta,nstd_survival,nbeta_survival,nstd_capture,nbeta_capture, nstd_tagloss,nbeta_tagloss)
{ design_matrix<- matrix(0,nrow=nstd,ncol=nbeta);
  design_matrix[1:nstd_survival,1:nbeta_survival]=survival.design;
  design_matrix[(nstd_survival+1):(nstd_survival+nstd_capture),(nbeta_survival+1):(nbeta_survival+nbeta_capture)]=capture.design;
  design_matrix[(nstd_survival+nstd_capture+1):(nstd_survival+nstd_capture+nstd_tagloss),(nbeta_survival+nbeta_capture+1):(nbeta_survival+nbeta_capture+nbeta_tagloss)]=tagloss.design;
  design_matrix[(nstd_survival+nstd_capture+nstd_tagloss+1):nstd,(nbeta_survival+nbeta_capture+nbeta_tagloss+1):nbeta]=bstar.design;
  design_matrix
}

######################################################################################
#The initial standard parameters are contained in matrices to begin: Each of the survival, capture, tagloss, entry (and perhaps "release"/"loss on capture") parameters are read in as matrices where the rows correspond to the population groups and the columns correspond to the sample times. The "pack_std" function uses these matrices to create a single vector of initial std parameters which can be used with the "mylogit" function. This vector contains first all the initial "survival" parameters for the 1st group, then all the survival parameters for the 2nd group, then for the 3rd group, etc. The vector then contains all the initial "capture" parms for the 1st group, then those for the 2nd group, and so on.  Following the capture parms are the "tagloss" parms (similarly listed by sample times and groups) and finally the "entry" parms, which follow the same listing order.
#std_survival= matrix (ngroup x (nsample-1)) that specifies initial survival parms
#std_capture= matrix (ngroup x nsample) that specifies initial capture parms
#std_tagloss= array (ngroup x (nsample-1) x (nsample-1)) that specifies initial tagloss parms
#std_entry= matrix (ngroup x nsample) that specifies initial entry parms
#std_bstar= matrix (ngroup x nsample) that specifies initial bstar parms
######################################################################################
pack_std<- function(std_survival,std_capture,std_tagloss,std_bstar, nstd_survival,nstd_capture,nstd_tagloss,nstd_bstar,ngroup,nsample)
{ stdparms.vect<- vector(mode="double");
  index=1;
  for(g in 1:ngroup)
  {	for(j in 1:(nstd_survival/ngroup))
  	{ stdparms.vect[index]=std_survival[g,j];
  	  index=index+1;
  	}
  }
  for(g in 1:ngroup)
  { for(j in 1:(nstd_capture/ngroup))
  	{ stdparms.vect[index]=std_capture[g,j];
  	  index=index+1;
  	}
  }
  for(g in 1:ngroup)
  { for(f in 1:(nsample-1))
    { for(j in 1:(nsample-1))
  	  { if(j>=f)
  	  	{ stdparms.vect[index]=std_tagloss[g,f,j];
  	  	  index= index+1;
  	    }
  	  }
  	}
  }
  for(g in 1:ngroup)
  { for(j in 1:(nstd_bstar/ngroup))
  	{ stdparms.vect[index]=std_bstar[g,j];
  	  index=index+1;
  	}
  }
  stdparms.vect
}


######################################################################################
#The mylogit function is the logit function, logit(x) = ln(x/(1-x))
#This function is called upon by the "transform_std" function
######################################################################################
mylogit<- function(std_parms)
{ ans<- vector(mode="double",length=length(std_parms))
  for(h in 1:length(std_parms))
  { if(std_parms[h]<=0){ans[h]=-100}
     else{ if(std_parms[h]>=1){ans[h]=100}
   		    else{ ans[h]=log(std_parms[h])-log(1-std_parms[h])}
   	     }
  }
  ans
}


######################################################################################
#The transform_std function is used once at the beginning of the program.  It uses the mylogit function to transform the initial standard parameters to beta parameters so the latter may be estimated using maximum likelihood techniques.
######################################################################################
transform_std<- function(design_matrix, stdparms_initial)
{ y= mylogit(stdparms_initial); #y = logit(stdparms_initial) = X*beta
  ans= lsfit(design_matrix, y, intercept=0);
  beta_initial= as.vector(ans$coefficients, mode="double")
}


######################################################################################
#The transform_beta function: Once the beta parameters have been estimated, the transform_beta function uses an anti-logit function to obtain the original standard parameters from the beta estimates.  Due to the invariance property of MLEs, the estimated beta parameters can, through the transform_beta function, yield MLEs of the original (standard) parameters.
######################################################################################
transform_beta<- function(design_matrix,beta)
{ std_parms = as.vector(exp(design_matrix %*% beta)/(1+exp(design_matrix %*% beta)))}


######################################################################################
#The chi0 function calculates the necessary ?? values for those animals that are never caught.  Here ?? = P(not observing an animal after time j given that the animal is not tagged). All the necessary ?? values are calculated here through a recursive process. These will be necessary in the identification of the likelihood function.
######################################################################################
chi0<-function(samp_time,group,std_survival,std_capture, std_tagloss,nsample)
{ if(samp_time==nsample){ans=1}
   else{ ans= (1-std_survival[group,samp_time]+std_survival[group,samp_time]*(1-std_capture[group,(samp_time+1)])*chi0((samp_time+1),group,std_survival,std_capture, std_tagloss,nsample))
  	   }
}


######################################################################################
#The chi1 function calculates the necessary ?? values if an animal has 1 tags at the final sample time it is caught.  Here, ?? = P(not observing an animal after time j given animal is tagged). This function calculates all the necessary ?? values through a recursive process. These will be necessary in the identification of the likelihood function.
######################################################################################
chi1<-function(first_time,samp_time,group,std_survival,std_capture, std_tagloss,nsample)
{ if(samp_time==nsample){ans=1}
   else{ ans= 1-std_survival[group,samp_time]+std_survival[group,samp_time]*(1-std_capture[group,samp_time+1])*std_tagloss[group,first_time,samp_time]*chi1(first_time,(samp_time+1),group,std_survival,std_capture,std_tagloss,nsample)+std_survival[group,samp_time]*(1-std_tagloss[group,first_time,samp_time])
       }
  ans 
}


######################################################################################
#The chi2 function calculates the necessary ?? values if an animal has 2 tags at the final sample time it is caught.  Here, ?? = P(not observing an animal after time j given animal is tagged). This function calculates all the necessary ?? values through a recursive process. These will be necessary in the identification of the likelihood function.
######################################################################################
chi2<-function(first_time,samp_time,group,std_survival,std_capture,std_tagloss, nsample)
{ if(samp_time==nsample){ans=1}
   else{ ans= 1-std_survival[group,samp_time]+std_survival[group,samp_time]*(1-std_capture[group,samp_time+1])*(std_tagloss[group,first_time,samp_time])^2*chi2(first_time,(samp_time+1),group,std_survival,std_capture,std_tagloss,nsample)+2*std_survival[group,samp_time]*(1-std_capture[group,samp_time+1])*std_tagloss[group,first_time,samp_time]*(1-std_tagloss[group,first_time,samp_time])*chi1(first_time,(samp_time+1), group,std_survival,std_capture,std_tagloss,nsample)+std_survival[group,samp_time]*(1-std_tagloss[group,first_time,samp_time])^2
	   }
  ans
}


######################################################################################
#The bstar function calculates the bstar paramaters from the b (std_entry) parameters.  It is called upon by the Psi function
######################################################################################
bstar<-function(samp_time,group,std_entry,nsample)
{ if(samp_time==1){ ans=std_entry[group,samp_time]}
   else{ ans=0;
   		 for(j in samp_time:nsample)
   	     { ans= ans+std_entry[group,j];}
   	ans= std_entry[group,samp_time]/ans;
	   }
}


######################################################################################
#The bstar_to_b function is the opposite of the bstar function.  It takes a bstar and turns it back to it's corresponding b (std_entry) value.  It is called upon by the Psi function and later the P_0 function.
######################################################################################
bstar_to_b<- function(group,samp_time,std_bstar,nsample,ngroup)
{ if(samp_time==1){b=std_bstar[group,samp_time]}
	else{ b=1;
          for(u in 1:(samp_time-1))
          { b=b*(1-std_bstar[group,u]);}
          b=b*std_bstar[group,samp_time];
        }
}


######################################################################################
#The function Psi:  ??????¡ì?? = P(animal enters population, is still alive, and is not seen before time i). The Psi function calculates the necessary ??????¡ì?? values through a recursive process. These will be necessary in the identification of the likelihood function
######################################################################################
Psi<- function(samp_time,group,std_bstar,std_survival,std_capture,nsample,ngroup)
{ if(samp_time==1){ans= std_bstar[group,1];}
	else{ ans= Psi((samp_time-1),group,std_bstar,std_survival,std_capture,nsample, ngroup)*(1-std_capture[group,samp_time-1])*std_survival[group,samp_time-1] +bstar_to_b(group,samp_time,std_bstar,nsample,ngroup);
		}
  ans
}	


######################################################################################
#The main function creates the likelihoods for each animal from when it is first seen to when it is last seen
######################################################################################
main<-function(group,obs_num,first,last,std_survival,std_capture,std_tagloss, capt.hist,tag.hist,ngroup,nobs,nsample,ntag)
{ lik_main=1;
  g=group;
  i=obs_num;
  if(first[g,i]<last[g,i])
  { for(j in first[g,i]:(last[g,i]-1))#loop for the survival (phi) parameters
    { lik_main= lik_main*std_survival[g,j];
  	}
  	  
    for(j in first[g,i]:last[g,i])#loop for the recapture (p) parameters
    { if(capt.hist[g,i,j]==1){ lik_main= lik_main*std_capture[g,j]}
       else{ lik_main= lik_main*(1-std_capture[g,j])}
    }
   	  
    for(d in 1:ntag)#loop for the tag-retention (Lambda) parameters
    { j=first[g,i];
      while(j<last[g,i])
      { nextt= j+1;
        while(capt.hist[g,i,nextt]==0)
   	    { nextt=nextt+1;}
   	    if(tag.hist[g,i,j,d]==1 & tag.hist[g,i,nextt,d]==1)
   	    { lik_main= lik_main*prod(std_tagloss[g,first[g,i],j:(nextt-1)]);}
   	    if(tag.hist[g,i,j,d]==1 & tag.hist[g,i,nextt,d]==0)
   	    { lik_main= lik_main*(1-prod(std_tagloss[g,first[g,i],j:(nextt-1)]));}
   	    j=nextt;
      }
    }
  }
  if(first[g,i]==last[g,i])
  { lik_main=lik_main*std_capture[g,first[g,i]];}
  lik_main
}


######################################################################################
#The P_0 value is used in the calculation of Nhat.initial and within the log-likelihood specification; here, P_0= probability of never being seen
######################################################################################
P_0_func<- function(std_survival,std_capture,std_tagloss,std_entry,ngroup,nsample)
{ P_0=vector(mode="double",length=ngroup);
  for(g in 1:ngroup)
  { P0=0;
    for(j in 1:nsample)
    {#print(std_entry)
    P0= P0+std_entry[g,j]*(1-std_capture[g,j])*chi0(j,g,std_survival,std_capture, std_tagloss,nsample);
    }
    P_0[g]=P0;
  }
  P_0
}


######################################################################################
#Nhat.initial is the initial estimate of N, the super-population size; it is necessary for calculationg L_1^A
######################################################################################
Nhat_func<- function(nobs,P_0,ngroup)
{ Nhat=vector(mode="double",length=ngroup);
  for(g in 1:ngroup)
  { Nhat[g]= nobs[g]/(1-P_0[g])}
  Nhat
}


######################################################################################
#The negloglik function creates the negative-log-likelihood function for the model
######################################################################################
negloglik<- function(first,last,beta,capt.hist,tag.hist,numtag,ngroup, nobs_total,nobs,nsample,ntag,nstd_survival, nstd_capture,nstd_tagloss,flag,Nhat)
{ std_parms= transform_beta(design_matrix,beta);
#	print(std_parms)#TESTING!!.........
  std_survival=matrix(std_parms[1:nstd_survival],nrow=ngroup,ncol=(nsample-1),byrow=TRUE);
  std_capture=matrix(std_parms[(nstd_survival+1):(nstd_survival+nstd_capture)],nrow=ngroup,ncol=nsample,byrow=TRUE);
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
  std_bstar=matrix(std_parms[(nstd_survival+nstd_capture+nstd_tagloss+1):length(std_parms)],nrow=ngroup,ncol=nsample,byrow=TRUE);
#  print("bstar")
#  print(std_bstar)
  std_entry=matrix(nrow=ngroup,ncol=nsample);
  for(g in 1:ngroup)
  { for(j in 1:nsample)
  	{ std_entry[g,j]= bstar_to_b(g,j,std_bstar,nsample,ngroup)}
  }
#    print("std_entry")
#    print(std_entry)
#  print(std_bstar);##TESTING
#  print(std_entry);##TESTING 
#  print(std_survival);

if (flag==1)
  {std_capture[,1]=1
   std_capture[,nsample]=1
  }

  lik<- matrix(1, nrow=ngroup,ncol=max(nobs[1:ngroup]));
  for(g in 1:ngroup)
  { for(i in 1:nobs[g])
  	{ lik[g,i]= lik[g,i]*Psi(first[g,i],g,std_bstar,std_survival,std_capture, nsample,ngroup);
  	  lik[g,i]= lik[g,i]*main(g,i,first,last,std_survival,std_capture,std_tagloss, capt.hist,tag.hist,ngroup,nobs,nsample,ntag);
  	  if(numtag[g,i]==1)
  	  { lik[g,i]= lik[g,i]*chi1(first[g,i],last[g,i],g,std_survival,std_capture, std_tagloss,nsample)
   	  }
   	  if(numtag[g,i]==2)
   	  { lik[g,i]= lik[g,i]*chi2(first[g,i],last[g,i],g,std_survival,std_capture, std_tagloss,nsample)
   	  }
   }
  }
  
  P_0=P_0_func(std_survival,std_capture,std_tagloss,std_entry,ngroup,nsample);
  for(g in 1:ngroup){
  	#print(P_0)
  	P_0[g]=max(P_0[g],0.0001)   #ensures no log(0) will occur
  	#CHECK is "0.0001" appropriate?
  }
  
  Nhat=Nhat_func(nobs,P_0,ngroup);

  negloglik.vect=vector(mode="double",length=ngroup);
  ##print(sum(log(lik[g,])));#TESTING
##  print(P_0)##TESTING
  #print(Nhat)#TESTING
  for(g in 1:ngroup)
  { negloglik.vect[g]=(sum(log(lik[g,]))+max(0,Nhat[g]-nobs[g])*log(P_0[g])+lgamma(Nhat[g]+1)-lgamma(Nhat[g]-nobs[g]+1));
  }	
##    print(negloglik.vect)#TESTING
  negloglik=-sum(negloglik.vect)
}


######################################################################################
#The COV function calculates the covariance 
######################################################################################
COV<- function(design_matrix,Hess,stdparms)
{ library(MASS)
  inv_Hess=ginv(Hess,tol=sqrt(.Machine$double.eps));
  Cov_XB=design_matrix%*%inv_Hess%*%t(design_matrix);
  std.vect=stdparms*(1-stdparms);
  COV=diag(std.vect)%*%Cov_XB%*%diag(std.vect);
  for(r in 1:nrow(design_matrix)){ print(COV[r,r])}
  COV
}

######################################################################################
VAR<- function(design_matrix,Hess,stdparms)
{ library(MASS)
  inv_Hess=ginv(Hess,tol=sqrt(.Machine$double.eps));
  Cov_XB=design_matrix%*%inv_Hess%*%t(design_matrix);
  std.vect=stdparms*(1-stdparms);
  COV=diag(std.vect)%*%Cov_XB%*%diag(std.vect);
  VAR.vect=vector(mode="double");
  i=1;
  for(r in 1:nrow(design_matrix)){ VAR.vect[i]=COV[r,r]; i=i+1;}
  VAR.vect
}



######################################################################################
#The AIC function calculates the Akaike Information Criterion (AIC), which is an index that is used in selecting between competing models.  Higher values of AIC indicate a less favorable model.
######################################################################################
calc_AIC<- function(design_matrix,lmax,nstd)#lmax is the maximized log-likelihood
{  ##nstd= nrow(design_matrix);
   AIC= 2*nstd-2*lmax# returns AIC
}
#####################################################################################
N_est_func<-function(std_survival,std_capture,std_entry,samp_time,group,ngroup, nsample,nobs){
	P_0 = P_0_func(std_survival,std_capture, std_tagloss,std_entry,ngroup,nsample);
	Nhat = Nhat_func(nobs,P_0,ngroup);
	if(samp_time==1){ N_est=Nhat[group]*std_entry[group,samp_time]}
	if(samp_time>1 && samp_time<=nsample){ 
		N_est = N_est_func(std_survival,std_capture,std_entry,samp_time-1,group,ngroup,nsample,nobs)*std_survival[group,samp_time-1] + Nhat[group]*std_entry[group,samp_time-1]
	}
	N_est
}

