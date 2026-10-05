# Tag history generation for the JSTL model
# January 21, 2013

# statistics to be fixed before generation
n=10 #number of individuals in the population
nsample=5 #number of sample times
frac_double=0.50 # fraction of individuals double tagged

#creation of vectors
entrypoint=numeric(length=n)
b=numeric(length=nsample) #probability of birth/immigration (sum to 1)
bstar=numeric(length=nsample) # fraction entering the population
phi=numeric(length=(nsample-1)) #probability of survival
p=numeric(length=nsample)  #probability of capture
lambda=matrix(rep(0, (nsample-1)*(nsample-1)),nrow=nsample-1, ncol=nsample-1) #probability of tag retention from release group i between periods j and j+1 (lambda[i,j])
first=numeric(length=n) #first entry time of individual i
a=matrix(rep(0,n*nsample), nrow=n, ncol=nsample) #alive status of individual i at time j
history=matrix(rep(0,n*nsample*2), nrow=n, ncol=(nsample*2))

#initialization of parameters
b[1:nsample]=1/nsample
for(i in 1:nsample){
	bstar[i]=b[i]/sum(b[i:nsample])
	}
p[1:nsample]=1
phi[1:nsample-1]=1
loss_parm=0
lambda[1:nsample-1, 1:nsample-1]=1

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
	a[i,entrypoint[i]]=1
	while((l<=nsample) & (alive==1)){
		if(runif(1,0,1)<p[l]){
			first[i]=l
			break
			}#endif
		if((l<nsample) & (first[i]==0)){
			if(runif(1,0,1)>phi[l]){
				alive=0
	#			a[i,l+1]=alive
			}#endif
		}#endif
		l=l+1
	}#endwhile
	
# Determine tagging history
	if(first[i]>0 & (first[i]<nsample)){
		position=first[i]*2-1
		if(runif(1,0,1)<frac_double){ #is the individual double tagged?
			history[i,position]=1 #double tagged gets history '11'
			history[i,position+1]=1
			tag1=1 # tag 1 status
			tag2=1 #tag 2 status
			}
		else{
			 history[i,position]=1 #single tag gets history'10'
			 tag1=1
			 tag2=0
			}#endif
			print(history)
		# Determine if lost on capture
		loss=0
		if(runif(1,0,1)<loss_parm){ 
			loss=1
			} #endif
		alive=alive*(1-loss)
		for(k  in (first[i]+1):nsample){
			print(k)
		# Does the individual survive?
			if(runif(1,0,1)>phi[k-1]){
				alive=0
			}#endif
	
		# Does the individual retain its tags?	
			if(runif(1,0,1)>lambda[first[i],k-1]){
				tag1=0
			}#endif
			if(runif(1,0,1)>lambda[first[i],k-1]){
				tag2=0
			}#endif
		# Is the animal captured and if so are they lost on capture?	
			if(runif(1,0,1) <p[k] & (alive==1) & ((tag1+tag2)>0)){ #animal is captured, alive and has at least one tag
				loss=0
				if(runif(1,0,1)<loss_parm){ #are they lost on capture?
					loss=1
				}#endif
				alive=alive*(1-loss)
				position=2*k-1
				history[i,position]=tag1
				history[i,position+1]=tag2	
			}#endif
		}#endfor
	}#endif
	if(first[i]==nsample){
	position=first[i]*2-1
		if(runif(1,0,1)<frac_double){ #is the individual double tagged?
			history[i,position]=1 #double tagged gets history '11'
			history[i,position+1]=1
			tag1=1 # tag 1 status
			tag2=1 #tag 2 status
			}
		else{
			 history[i,position]=1 #single tag gets history'10'
			 tag1=1
			}#endif	
	}
}#endfor


# Delete zero histories
obs_history=matrix(0,nrow=n,ncol=2*nsample)
index=0
for(i in 1:length(first)){
	if(first[i]!=0){
		index=index+1
		obs_history[index,]=history[i,]
	}
}
obs_history=obs_history[1:index,]

%data<- cbind(obs_history," 1"," 1;")
%write(t(data),file="Data1000.txt",ncolumns=2*nsample+2, sep="")


obs_history
