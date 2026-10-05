# Tag history generation for the JS model
# January 22, 2013

# statistics to be fixed before generation
n=10 #number of individuals in the population
nsample=5 #number of sample times

#creation of vectors
entrypoint=numeric(length=n)
b=numeric(length=nsample) #probability of birth/immigration (sum to 1)
bstar=numeric(length=nsample) # fraction entering the population
phi=numeric(length=(nsample-1)) #probability of survival
p=numeric(length=nsample)  #probability of capture
first=numeric(length=n) #first entry time of individual i
history=matrix(rep(0,n*nsample), nrow=n, ncol=nsample)
freq=numeric(length=n) #frequency of each capture history.  Normally this is equal to 1 but if lost on capture this is -1

#initialization of parameters
b[1:nsample]=1/nsample
for(i in 1:nsample){
	bstar[i]=b[i]/sum(b[i:nsample])
	}
p[1:nsample]=0.9
phi[1:nsample-1]=0.8
loss_parm=0.1
freq[1:n]=0

for(i in 1:n){
# Determine when the individual enters the population (before entrypoint[i])
	j=1
	alive=0
	while(alive==0){
		if(runif(1,0,1)<bstar[j]){
			alive=1
			entrypoint[i]=j
			break
		}#endif
		j=j+1
	}#endwhile
# Determine when the individual was first captured
	l=entrypoint[i]
	while((l<=nsample) & (alive==1)){
		if(runif(1,0,1)<p[l]){
			first[i]=l
			freq[i]=1
			break
			}#endif
		if((l<nsample) & (first[i]==0)){
			if(runif(1,0,1)>phi[l]){
				alive=0
			}#endif
		}#endif
		l=l+1
	}#endwhile
	
# Determine tagging history
	if(first[i]>0 & (first[i]<nsample)){
		position=first[i]
		history[i,position]=1 
		# Determine if lost on capture
		loss=0
		if(runif(1,0,1)<loss_parm){ 
			loss=1
			freq[i]=-1
			} #endif
		alive=alive*(1-loss)
		for(k  in (first[i]+1):nsample){
		# Does the individual survive?
			if(runif(1,0,1)>phi[k-1]){
				alive=0
			}#endif
	
		# Is the animal captured and if so are they lost on capture?	
			if(runif(1,0,1) <p[k] & (alive==1)){ #animal is captured, alive and has at least one tag
				loss=0
				if(runif(1,0,1)<loss_parm){ #are they lost on capture?
					loss=1
					freq[i]=-1
				}#endif
				alive=alive*(1-loss)
				history[i,k]=1
			}#endif
		}#endfor
	}#endif
	if(first[i]==nsample){
		history[i,first[i]]=1
		if(runif(1,0,1)<loss_parm){ #are they lost on capture?
			freq[i]=-1
		} #endif
	}#endif	
}#endfor


# Delete zero histories
obs_history=matrix(0,nrow=n,ncol=nsample)
index=0
for(i in 1:length(first)){
	if(first[i]!=0){
		index=index+1
		obs_history[index,]=history[i,]
	}
}
obs_history=obs_history[1:index,]

data<- cbind(obs_history," 1"," 1;")
write(t(data),file="Data1000.txt",ncolumns=nsample+2, sep="")

