#!/home/lcowen/software/R-2.14.0/bin/Rscript

library(RMark);
#setwd('/Users/lcowen/Documents/Batch Marks/Takis Sims')

# JS data generated in matlab analysed in Rmark
# March 28, 2013


# Functions that are called on
run.models=function(){
	model.JS=mark(model.process,model.ddl, model="POPAN")
}
all.results=list()


# Read in the history matrix

sim_num=1

for(sim_num in 1:1){
   history=read.table(file=paste("history", sim_num, "out", "txt", sep="."), colClasses="character")


# Determine when the individual was first captured
   n=length(history[,1]) # number of individuals in the population including those not seen
   nsample=length(history[1,]) # number of sample times

   first=numeric(length=n) #first capture time of individual i
   for(i in 1:n){
	   for(j in 1:nsample){
          if(history[i,j]==1){
       	     first[i]<-j
       	     break
          }		
	   }
   }
# remove 000 histories
	temp=cbind(history,first)
	temp2=subset(temp,first>0)
	history2=temp2[,1:nsample]
	



# Save the data to a file ready for Rmark- remove 000 histories first
	temp=cbind(history,first)

	temp2=subset(temp,first>0)
	history2=as.matrix(temp2[,1:nsample])

	cat(history2[1, 1:nsample], " 1;", sep="", file=paste("JShist",sim_num, "txt", sep="."), append=FALSE, fill=TRUE)
	    for (i in 2:length(history2[,1])) 
	       {cat(history2[i, 1:nsample], " 1;", sep="", file=paste("JShist",sim_num, "txt", sep="."), append=TRUE, fill=TRUE)
	   	   }
 

# Convert data to Rmark dataframe	
	model1=convert.inp(paste("JShist",sim_num, "txt", sep="."))
	model.process=process.data(model1, model="POPAN")
	model.ddl=make.design.data(model.process)

	all.results[[sim_num]]=run.models()
	Nderiv=popan.derived(model.process, all.results[[sim_num]])

# Print resulting estimates and standard errors to results files
	if(sim_num==1){
        cat(Nderiv$Nbyocc$N,	"\n", file="popsize.txt",append=FALSE)
	cat(Nderiv$Nbyocc$se, "\n", file="popsize.se.txt", append=FALSE)
	}
	else{
	cat(Nderiv$Nbyocc$N,  "\n", file="popsize.txt",append=TRUE)
        cat(Nderiv$Nbyocc$se, "\n", file="popsize.se.txt", append=TRUE)
	}
}#end of simulation	
#cleanup(ask=FALSE)

