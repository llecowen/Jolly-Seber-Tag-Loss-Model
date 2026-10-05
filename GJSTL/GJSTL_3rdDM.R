#############################################################
##No group or time: model 3##
#############################################################
survival.design_3<-matrix(data=1,nrow=(nsample-1)*ngroup,ncol=1)

capture.design_3<-matrix(data=1,nrow=(nsample)*ngroup,ncol=1)

bstar_0<-matrix(data=0,nrow=nsample,ncol=nsample)
index=0
for (i in 1:nsample)
 {index=index+1
   {bstar_0[index,i]=1}
}
bstar_1<-matrix(data=0,nrow=ngroup*nsample,ncol=nsample)
for (g in 1:ngroup)
{
 bstar_1[((g-1)*nsample+1):(g*nsample),]=bstar_0
}

R = nsample * (nsample-1) / 2 
tag.design_3<-matrix(data=1,nrow=ngroup*R,ncol=1)

#############################################################
##Only group: model 1########
#############################################################
survival.design_1<-matrix(data=0,nrow=ngroup*(nsample-1),ncol=ngroup)
for(g in 1:ngroup)
 {
     survival.design_1[((g-1)*(nsample-1)+1):((nsample-1)*g),g]=1
 }

capture.design_1<-matrix(data=0,nrow=ngroup*nsample,ncol=ngroup)
for(g in 1:ngroup)
 {
     capture.design_1[((g-1)*nsample+1):(nsample*g),g]=1
 }

tag.design_1<-matrix(data=0,ncol=ngroup,nrow=ngroup*R)
for(g in 1:ngroup)
 {tag.design_1[((g-1)*R+1):(g*R),g]=1
 }

################################################################
##Only time: model 7###########
################################################################
survival.design<-matrix(data=0,nrow=(nsample-1),ncol=(nsample-1))
index=0
for (i in 1:(nsample-1))
 {index=index+1
   {survival.design[index,i]=1}
}
survival.design_7 <-matrix(data=0,nrow=ngroup*(nsample-1),ncol=(nsample-1))
for (g in 1:ngroup)
{
 survival.design_7[((g-1)*(nsample-1)+1):(g*(nsample-1)),]=survival.design
}

capture.design<-matrix(data=0,nrow=nsample,ncol=nsample)
index=0
for(i in 1:nsample)
 {index=index+1
  {capture.design[index,i]=1}
 }
capture.design_7<-matrix(data=0,nrow=ngroup*nsample,ncol=nsample)
for (g in 1:ngroup)
{
 capture.design_7[((g-1)*nsample+1):(g*nsample),]=capture.design
}

tag.design_7<-matrix(data=0,ncol=(nsample-1),nrow=(ngroup*R))
shift=0
index=0
for (g in 1:ngroup)
 {
  for (i in 1:(nsample-1))
 {
    for (j in 1:(nsample-i))
 {
      index = index + 1
      
      tag.design_7[index, ((j+shift-1) %% (nsample-1))+1] = 1
    }   
   shift = shift + 1  }
  
} 
 

###########################################################################
###Both group and time: model 2#############
#####################################################################
survival.design_2<-matrix(data=0,ncol=ngroup*(nsample-1),nrow=ngroup*(nsample-1))
for (g in 1:ngroup)
{
 survival.design_2[((g-1)*(nsample-1)+1):(g*(nsample-1)),((g-1)*(nsample-1)+1):(g*(nsample-1))]=survival.design
}

capture.design_2<-matrix(data=0,ncol=ngroup*nsample,nrow=ngroup*nsample)
for (g in 1:ngroup)
{
 capture.design_2[((g-1)*nsample+1):(g*nsample),((g-1)*nsample+1):(g*nsample)]=capture.design
}

tag.design_2<-matrix(data=0,ncol=ngroup*(nsample-1),nrow=ngroup*R)
for (g in 1:ngroup)
 { tag.design_2[((g-1)*R+1):(g*R),((g-1)*(nsample-1)+1):(g*(nsample-1))]=tag.design_7[1:R,1:(nsample-1)]
 }



