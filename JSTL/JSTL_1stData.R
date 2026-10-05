##########################################################################
##1. Create the numeric vector of counts (frequency of capture history)
##2. Indentify ngroup (the number of groups) and nobs_total (the number of observed animals)
##3. Create a matrix of the individual tag histories
##4. Extract the semicolon and expand the last numeric vector of groups
##5. Expand the list of histories to get one single history for each animal
##6. The count vector is used to expand the history vector
############################################################################
ntag=2;


data.inp=read.table(file=paste("Data1000b.2",42, "txt", sep="."), colClasses="character")


ngroup= dim(data.inp)[2]-2;#number of groups in model
ngroup
count<- as.numeric(unlist(data.inp$V2));
count
nobs_total= sum(count);#number of animals observed in experiment
nobs_total
nsample= nchar(data.inp$V1[1])/ntag; #number of sample times used in experiment
nsample

history <- matrix(as.numeric(unlist(strsplit(rep(data.inp$V1,times=count), ""))),ncol=nchar(data.inp$V1[1]), byrow=TRUE)  ##the tag history for each individual


# Jun 2012 added code to make sure the data contains only unique tag histories where 11 10 10 is not different from 11 10 10 
# (previously we assumed these were different)
# check the tag histories by summing tag history at each sample time.  As long as these are equal then the tag histor is deemed unique.  

flag=numeric(length(history[,1])) # flag set to 1 if the history is unique
nfreq=numeric(length(history[,1])) # frequency of each unique tag history
flag2=0          				# temporary flag set to 1 as long as the history is unique at time j
flag[1]=0 # first tag history is unique
for(i in 1:(length(history[,1])-1)){
	if(flag[i]==0){
		nfreq[i]=1}
   for(k in (i+1):length(history[,1])){
      for(j in seq(1,ntag*nsample,by=2)){
         if(sum(history[i,j:(j+1)])==sum(history[k,j:(j+1)] & flag[k]==0)){            
         	flag2=1 # tag histories are the same so far
         }
         else{
         	  flag2=0 # tag histories are not equal
            break # As soon as tag histories are deemed not equal, break the loop
                
         }
      }
      if(flag2==1){ # When a replicate tag history is found, we add one to the frequency of the unique history.
      nfreq[i]=nfreq[i]+1
      flag[k]=1
   }
 }
}
if(flag[i+1]==0){nfreq[i+1]=1}

all.data=cbind(flag, nfreq, history)
new.data=all.data[which(all.data[,1]==0),]



count<- new.data[,2]
nobs_total= sum(count);#number of animals observed in experiment
nobs_total

unique.history <- new.data[,3:length(new.data[1,])]


group.vect<- vector(mode="double",length=(ngroup*nobs_total));#this intermediary-step group.vect is later used to create the intermediary-step group matrix below

vect<- as.numeric(unlist(strsplit(rep(data.inp[,(ngroup+2)],times=count),";")));
#vect here is the last column of data.inp, which specifies those individuals in the last group, minus the ";"

for(i in 1:nobs_total)#this loop fills the group.vect with the last column of 1's and 0's in vect
{ g=(ngroup-1)*nobs_total+i;
  group.vect[g]=vect[i]; 
}

#Only for GJSTL when there are more than 1 group
#g=1;
#for(v in 3:(3+ngroup-1))#this loop fills the group.vect with the #other (not last) columns of 1's and 0's in vect, which specify #which group each observation belongs to
#{ vect<- as.numeric(unlist(rep(data.inp[,v], times=count)))
#  i=1;
#  while(i<=nobs_total)
#  { group.vect[g]=vect[i];
#  	g=g+1;
#  	i=i+1; 
#  }
#}

group.mat<- matrix(group.vect,ncol=ngroup, nrow=nobs_total);#intermediary tool: matrix which identifies which observations correspond to which group

######################################################################################
#The nobs vector specifies the number of observated animals belong to each group
######################################################################################
nobs<- vector(mode="double",length=ngroup)
for(g in 1:ngroup){nobs[g]=sum(group.mat[,g]);}
nobs

######################################################################################
#The group.vect is a vector corresponding to the original "history" matrix; it identifies which group each observed animal in the "history" belongs to.
######################################################################################
group.vect<-vector(mode="double",length=nobs_total)
for(i in 1:nobs_total)
{ for(g in 1:ngroup)
  { if(group.mat[i,g]==1){group.vect[i]=g}
  }
}

######################################################################################
#The hist.list is a list containing the history matrix but separated into groups, that is, all the observations belonging to 
##group 1 are in a matrix contained as the first element in the list; the observations belonging to group 2 are in a matrix contained 
##as the second element in the list, and so on.
######################################################################################
hist.list<- list()
for(g in 1:ngroup)
{ hist.mat<- matrix(nrow=nobs[g],ncol=(nsample*ntag));
  index=1;
  for(i in 1:nobs_total)
  { if(group.vect[i]==g)
  	{ hist.mat[index,]=history[i,];
  	  index=index+1;
  	}
  }
  hist.list[[g]]=hist.mat
}

######################################################################################
#The tag.hist array is 4-dimensional and uses the hist.list to create group-specific tag histories
######################################################################################
tag.hist<- array(data=NA,dim=c(ngroup,max(nobs),nsample,ntag))
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { index=1;
    for(j in 1:nsample) 
    { for(d in 1:ntag)
	  { tag.hist[g,i,j,d]= hist.list[[g]][i,index];
  	    index=index+1;
  	  }
  	}
  }
}


######################################################################################
#The capt.hist array uses the tag history (created above) to identify the capture history of each group-specific observation
######################################################################################
capt.hist<- array(NA, dim=c(ngroup, max(nobs),nsample));
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { for(j in 1:nsample)
 	{ if(sum(tag.hist[g,i,j,1:ntag])>=1){capt.hist[g,i,j]=1}
   	   else{capt.hist[g,i,j]=0}
   	}
  }
}


######################################################################################
#The "first" matrix specifies the first sample time at which each animal was observed 
######################################################################################
first<- array(NA,dim=c(ngroup,max(nobs)));
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { for(j in 1:nsample)
    { if(capt.hist[g,i,j]==1){first[g,i]= j; break}
    }
  }
}


######################################################################################
#The "last" matrix specifies the last sample time at which each animal was observed
######################################################################################
last<- array(NA,dim=c(ngroup,max(nobs)));
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { for(j in nsample:1)
    { if(capt.hist[g,i,j]==1){last[g,i]= j; break}
    }
  }
}


######################################################################################
#The numtag matrix specifies the number of tags at the last sample time for each animal observed 
######################################################################################
numtag<- array(NA,dim=c(ngroup,max(nobs)));
for(g in 1:ngroup)
{ for(i in 1:nobs[g])
  { numtag[g,i]= sum(tag.hist[g,i,last[g,i],1:ntag]);
  }
}

