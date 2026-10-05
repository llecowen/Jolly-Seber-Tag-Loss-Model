###########################################################################
#The following GOF code is specific to 2 groups and 2 tags
###########################################################################

#Saturated Model: P_hist is a vector 
P_hist = nobs/nobs_total;
P_hist

#Number of unique tag histories in each group -- currently works only for two groups (ngroup=2)
#groupNum= as.vector(as.numeric(strsplit(data.inp$V4,split=";")))+1
#numUnique= vector(length=ngroup)
#for(g in 1:ngroup)
#{	numUnique[g]=0;
#	for(h in 1: length(groupNum))
#	{	if(groupNum[h]==g)
#		{	numUnique[g]=numUnique[g]+1;
#		}
#	}	
#}




######################################################################
# Estimate the probability of double-tagging at the 1st capture time #
######################################################################

#for each group, calculate the number tagged of each type (double or single) at each sample occassion
n_tagged<-array(data=0,dim=c(ngroup,nsample,ntag))
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { if(sum(tag.hist[g,i,first[g,i],])==1)#single-tagged
  	{ n_tagged[g,first[g,i],1]=n_tagged[g,first[g,i],1]+1
  	}
  	else if(sum(tag.hist[g,i,first[g,i],])==2)#double-tagged
  	{ n_tagged[g,first[g,i],2]=n_tagged[g,first[g,i],2]+1
  	}
  }
}
n_tagged


# For each group, calculate the probability of being double-tagged at each sample occassion

#probability of double tagging at each sample time (& for each group)
p_double<-array(data=0,dim=c(ngroup,nsample))
for(g in 1:ngroup)
{ for(j in 1:nsample)
  { p_double[g,j]=n_tagged[g,j,2]/sum(n_tagged[g,j,])
  }
}
p_double



################################################################################
#The loglik function produces the log-likelihood for a specific observation i in group g
################################################################################

loglik_func<- function(group,observation_num,first,last,beta,capt.hist,tag.hist, numtag,ngroup,nobs_total,nobs,nsample,ntag,nstd_survival,nstd_capture, nstd_tagloss)
{ 
  std_parms= transform_beta(design_matrix,beta)
  #std_parms=c(0.6,0.24,0.8, 0.627, 1, 0.7, 1, 1, 0.55, 1, 0.7, 0.5, 0.5, 0.9, 0.7, 0.7, 0.16,0.375,1, 0.13, 0.375, 1) 
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
  std_entry=matrix(nrow=ngroup,ncol=nsample);
  for(g in 1:ngroup)
  { for(j in 1:nsample)
  	{ std_entry[g,j]= bstar_to_b(g,j,std_bstar,nsample,ngroup)}
  }
  
  g=group;
  i=observation_num;
  
  lik=1 
  lik= lik*Psi(first[g,i],g,std_bstar,std_survival,std_capture, nsample,ngroup);
  lik= lik*main(g,i,first,last,std_survival,std_capture,std_tagloss, capt.hist,tag.hist,ngroup,nobs,nsample,ntag);
  if(numtag[g,i]==1)
  { lik= lik*chi1(first[g,i],last[g,i],g,std_survival,std_capture,std_tagloss, nsample)
  }
  if(numtag[g,i]==2)
  { lik= lik*chi2(first[g,i],last[g,i],g,std_survival,std_capture,std_tagloss, nsample)
  }
  P_0=P_0_func(std_survival,std_capture,std_tagloss,std_entry,ngroup,nsample);
  ## hat=Nhat_func(nobs,P_0,ngroup);
  loglik=log(lik)
  loglik
}





################################################################################
#Develop the conditional likelihood
################################################################################

loglik_cond=vector(length=ngroup,mode="numeric")
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { if(sum(tag.hist[g,i,first[g,i],])==1) #single tag at first history
  	{P_0=P_0_func(std_survival,std_capture,std_tagloss,std_entry,ngroup,nsample)
       a=loglik_func(g,i,first,last,out$estimate,capt.hist,tag.hist,numtag,ngroup,nobs_total,nobs,nsample,ntag,nstd_survival,nstd_capture,nstd_tagloss)
       loglik_cond[g] = loglik_cond[g]+ a - log(max(1-P_0[g], 0.001)) + log(max(1 - p_double[g,first[g,i]],0.001))
  	}
  	
  	else #double tag at first history
  	{ P_0=P_0_func(std_survival,std_capture,std_tagloss,std_entry,ngroup,nsample)
        a=loglik_func(g,i,first,last,out$estimate,capt.hist,tag.hist,numtag,ngroup,nobs_total,nobs,nsample,ntag,nstd_survival,nstd_capture,nstd_tagloss)
        loglik_cond[g]=loglik_cond[g]+ a - log(max(1-P_0[g], 0.001)) +log(max(p_double[g,first[g,i]],0.001))
# print(loglik_cond[g])
	}
#a=loglik_cond
#  print(c(g, a[g]))
  }
}




################################################################################
## Saturated log-likelihood (ignoring factorial terms)
#################################################################################
freq= as.vector(as.numeric(strsplit(data.inp$V2,split=";")))
#groupNum= as.vector(as.numeric(strsplit(data.inp$V4,split=";")))+1

#numUnique= vector(length=ngroup)
#for(g in 1:ngroup)
#{	numUnique[g]=0;
#	for(h in 1: length(groupNum))
#	{	if(groupNum[h]==g)
#		{	numUnique[g]=numUnique[g]+1;
#		}
#	}	
#}


N_unique = vector (length=ngroup)
#for(g in 1:ngroup)
#{ N_unique[g]=numUnique[g+1]
# }
N_unique[ngroup]=length(data.inp$V3)


index=0
loglik_sat<- vector(length=ngroup,mode="numeric")
for(g in 1:ngroup)
{ for(i in 1:N_unique[g])
  {   
  index=index+1
  loglik_sat[g] = loglik_sat[g] + freq[index]*log(freq[index]/nobs[g]) 
  }
}

Tsat=sum(loglik_sat)
Tcum=sum(loglik_cond)

T = sum(loglik_sat) - sum(loglik_cond)
T

df<-sum(N_unique)-1-nbeta-ngroup*nsample
df

pv=1-pchisq(2*T,df)
pv
