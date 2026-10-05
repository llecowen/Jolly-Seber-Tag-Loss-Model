birth=matrix(data=NA,nrow=ngroup,ncol=nsample)
for (g in 1:ngroup)  ###compute the net births at each period
 {for (j in 1:nsample)
   {
    birth[g,j]=Nhat[g]*std_entry[g,j]
   }
 }
B=vector(length=(ngroup*nsample))
for (g in 1:ngroup)
 {for (j in 1:nsample)
   {B[((g-1)*nsample+j)]=Nhat[g]*std_entry[g,j]}
 }
###B=vector(length=(ngroup*(nsample+1)))
###B[1:nsample]=birth[1,]
###B[(nsample+2):(ngroup*(nsample+1)-1)]=birth[2,]
B


b=vector(length=(ngroup*nsample))
for (g in 1:ngroup)
 {for (j in 1:nsample)
   {b[((g-1)*nsample+j)]=std_bstar[g,j]}
 }
###b[1:nsample]=std_bstar[1,]
###b[(nsample+2):(ngroup*(nsample+1)-1)]=std_bstar[2,]
b


partial=matrix(data=0,nrow=(ngroup*nsample),ncol=(nsample+1))
  
for (g in 1:ngroup)
 {for (i in ((g-1)*nsample+1):(g*nsample))
    {
     partial[i,(nsample+1)]=B[i]/Nhat[g]
    }
 }
for (i in 1:(nsample-1))
   {partial[i,i]=B[i]/b[i]}
#for (i in (nsample+1):(ngroup*nsample))
#   {for (j in 1:(nsample-1)) 
#      {if ((i-j)==nsample)
#         {partial[i,j]=B[i]/b[i]}
#      }
#  }
for (i in 1:nsample)
      {for (j in 1:nsample)
        {if (j>i)
         {partial[j,i]=-B[j]/(1-b[i])}
        }
      }
#for (i in (nsample+1):(ngroup*nsample))
#     {for (j in (nsample+1):(ngroup*nsample))
#       {if (j>i)
#         {partial[j,j-i]=-B[j]/(1-b[j-i])}
#       }
#     }
H=matrix(data=0,nrow=(ngroup*nsample),ncol=nstd+1)
H[,(nstd+1-(nsample+1)+1):(nstd+1)]=partial[,1:(nsample+1)]

cov_net_birth=matrix(data=0,nrow=(ngroup*nsample),ncol=(ngroup*nsample))

for(g in 1:ngroup)
 {
  for (i in ((g-1)*nsample+1):(nsample*g))
    {for (j in ((g-1)*nsample+1):(nsample*g))
      {for (k in 1:(nstd+1))
        {for (l in 1:(nstd+1))
         {cov_net_birth[i,j]=cov_net_birth[i,j]+H[i,k]*H[j,l]*Cov_std1[k,l]}
        }
      }
    }
  }

se_net_birth=vector(mode="numeric",length=(nsample*ngroup))
for(i in 1:(ngroup*nsample))
  {se_net_birth[i]=sqrt(max(0,cov_net_birth[i,i]))}
se_net_birth

cov_pop_size=matrix(data=0,nrow=(ngroup*nsample),ncol=(ngroup*nsample))
se_pop_size=vector(mode="numeric",length=(nsample*ngroup))

H1=matrix(data=0,nrow=(ngroup*nsample),ncol=nstd+1)
phi=vector(length=(ngroup*(nsample-1)))
for (g in 1:ngroup)
 {for (j in 1:(nsample-1))
   {phi[((g-1)*(nsample-1)+j)]=std_survival[g,j]}
 }
phi
pop_size=vector(length=(ngroup*nsample))
for (g in 1:ngroup)
 {for (j in 1:nsample)
   {pop_size[((g-1)*nsample+j)]=N_est[g,j]}
 }
pop_size


H1[1,]=H[1,]
#H1[nsample+1,]=H[nsample+1,]
for (i in 2:nsample)
  {H1[i,]=H[i,]
   for (j in 1:(nstd+1))
     {H1[i,j]=H1[i,j]+H1[i-1,j]*phi[i-1]
     }
   H1[i,i-1]=H1[i,i-1]+pop_size[i-1]
  }
#for (i in (nsample+2):(ngroup*nsample))
#  {H1[i,]=H[i,]
#   for (j in 1:(nstd+1))
#     {H1[i,j]=H1[i,j]+H1[i-1,j]*phi[i-2]
#     }
#   H1[i,i-2]=H1[i,i-2]+pop_size[i-1]
#  }


for(g in 1:ngroup)
 {
  for (i in ((g-1)*nsample+1):(nsample*g))
    {for (j in ((g-1)*nsample+1):(nsample*g))
      {for (k in 1:(nstd+1))
        {for (l in 1:(nstd+1))
         {cov_pop_size[i,j]=cov_pop_size[i,j]+H1[i,k]*H1[j,l]*Cov_std1[k,l]}
        }
      }
    }
  }

for(i in 1:(ngroup*nsample))
  {se_pop_size[i]=sqrt(max(0,cov_pop_size[i,i]))}
se_pop_size


