optstrat<-function(X,n=2,nhmin=2,L=3,type_aloc="N",TIMECPU=Inf)
{
Round_sample<-function(Nh,nh,n,nhmin,L,rgn,nt,Vx,tx)
{
rg<-rep(list(-rgn:rgn),L)
cnh<-as.matrix(expand.grid(rg))
nr<-round(nh)
nhr=t(apply(cnh,1,function(x) nr+x))
nhr<-nhr[which(apply(nhr,1,min)>=nhmin),]
nhr<-nhr[which(apply(nhr,1,function(xn) min(Nh-xn))>=0),]
nhr<-nhr[which(abs(apply(nhr,1,sum)-n)<=nt),]
dfn<-apply(nhr,1,function(nx) sum((nh-nx)^2))
ix<-which.min(dfn)
nhr<-as.numeric(nhr[ix,])
cvtn<-100*sqrt(sum(Nh^2/nhr*Vx*(1-nhr/Nh)))/tx
return(list(nhr=nhr,cvtn=cvtn))

}   
   
   

###Calculates Variance considering sampling allocation#####################   
Calc_Sigma_Sigma2<-function(Nh,X,L,type_aloc,n,nhmin,ARRED=FALSE)
{
   N<-sum(Nh)   
   s2<-matrix(0,ncol=L,nrow=2)
   tx<-sum(X)
   for(h in 1:L)
   {
      s2[1,h]<-var(X[1:Nh[h]])
      s2[2,h]<-sd(X[1:Nh[h]])
      X<-X[-c(1:Nh[h])]
   }
   if (type_aloc=="N") {nh<-n*Nh*s2[2,]/sum(Nh*s2[2,])}
   if (type_aloc=="P") {nh<-n*Nh/N}
   if (type_aloc=="U") {nh<-rep(n/L,L)}
   df.real<-1-nh/Nh
   df.real<-ifelse(df.real<0,0,df.real)
   cvtx.real.fator<-100*sqrt(sum(Nh^2/nh*s2[1,]*df.real))/tx
   cvtx.real<-100*sqrt(sum(Nh^2/nh*s2[1,]))/tx
   ni=MultAlloc::BSSM_FD(Nh,s2[1,],tx,Cust=n,nmin=nhmin)$nh
   #ni=MultAlloc::BSSM_FC(Nh,s2[1,],tx,cvt=cvtx.real.fator/100,nmin=nhmin)$nh
   df.integer<-1-ni/Nh
   df.integer<-ifelse(df.integer<0,0,df.integer)
   cvtx.integer.fator<-100*sqrt(sum(Nh^2/ni*s2[1,]*df.integer))/tx
   cvtx.integer<-100*sqrt(sum(Nh^2/ni*s2[1,]))/tx
  
   if (ARRED==TRUE)
      {nia<-Round_sample(Nh,nh,n,nhmin,L,4,0,s2[1,],tx)
       ni.arred<-nia$nhr
       if (length(ni.arred)==0) {ni.arred<-round(nh)}
       cvtx.arred.fator<-nia$cvtn
       cvtx.arred<-100*sqrt(sum(Nh^2/ni.arred*s2[1,]))/tx
      } 
   else {nhr=NULL;cvtr=NULL}
   
   
   return(list(s2=s2,nh=nh,nho=ni,nha=ni.arred,
                 cvtx.real.fator=cvtx.real.fator,cvtx.real=cvtx.real,
                 cvtx.integer.fator=cvtx.integer.fator,
                 cvtx.integer=cvtx.integer,
                 cvtx.arred.fator=cvtx.arred.fator,
                 cvtx.arred=cvtx.arred))
                 
}   
####################################################################################

##############################################################################
Minimize_Variance_NPU<-function(type_aloc,nhmin)
{   
if (type_aloc=="N")
  {VijNij<-parApply(cl=clust,combine,2,function(x) sd(X[which(X>=vx[x[1]] & X<=vx[x[2]])]))
   W<-Nij/n
  }
if (type_aloc=="P")
  {VijNij<-parApply(cl=clust,combine,2,function(x) var(X[which(X>=vx[x[1]] & X<=vx[x[2]])]))
   W<-Nij/n  
  }
if (type_aloc=="U")
 {VijNij<-parApply(cl=clust,combine,2,function(x) var(X[which(X>=vx[x[1]] & X<=vx[x[2]])]))
  W<-(L*Nij^2)/n  
 }
   
t3<-proc.time()   
Sigma<-matrix(0,ncol=Nsr,nrow=Nsr)
Sigma[t(combine)]<-VijNij
Sigma<-Sigma+t(Sigma)
##################################################################################
###Generate model - matrix A, fobj, b and c

#################################################################################

tc<-choose(Nsr,2)

################################### Fobj ########################################
fobj<-rep(parApply(cl=clust,combine,2,function(x) W[x[1],x[2]]*Sigma[x[1],x[2]]),L)
#################### Constraint 20 - F1##########################################

R1<-matrix(0,ncol=L*tc,nrow=L)
li=1;ls=tc
cat("Generating Constraint 1   \n")
for(h in 1:L)
{R1[h,li:ls]<-rep(1,tc); li<-ls+1; ls<-ls+tc}
b<-rep(1,L+2)
desig<-rep("=",L+2)
#################################################################################

#################### Constraint 21 - F1##########################################
cat("Generating Constraint 2   \n")
R2<-apply(combine,2,function(x) ifelse(x[1]==1 & x[2]<Nsr-2*(L-2),1,0))
R2<-c(R2,rep(0,tc*(L-1)))

#################### Constraint 22 - F1##########################################
cat("Generating Constraint 3   \n")
R3<-apply(combine,2,function(x) ifelse(x[1]>=(2*(L-1)+1) & x[2]==Nsr,1,0))
R3<-c(rep(0,tc*(L-1)),R3)

#################### Constraint 23 - F1##########################################
cat("Generating Constraint 4   \n")
xh=-combine[2,]
xhm=combine[1,]
Ru=NULL
for(h in 1:(L-1))
 { b=c(b,1)
   desig=c(desig,"=")
   R<-c(rep(0,tc*(h-1)),xh,xhm,rep(0,tc*(L-1-h)))
   Ru=rbind(Ru,R)
}

cat("  \n")
###############################################################################
A<-rbind(R1,R2,R3,Ru)
t3<-(proc.time()-t3)[3]

model <- list()
model$obj<-fobj
model$modelsense<-"min"
model$rhs<-b
model$sense<-desig
model$vtype<-"B"
model$A<-A
params<-list()
params$NumericFocus <-1
params$TimeLimit<-TIMECPU
params$Heuristics=HE

cat("Applying Gurobi - Solver \n\n")
result <- gurobi(model, params)

xg=round(result$x)

t2<-result$runtime
y=matrix(xg,ncol=choose(Nsr,2),byrow=TRUE)
index<-apply(y,1,function(w) which(w==1))
cuts<-combine[,index]
Nh<-Nij[t(cuts)]
bh<-as.numeric(names(Xr[cuts]))
bh<-bh[seq(2,2*(L-1),2)]
S2<-Calc_Sigma_Sigma2(Nh,X,L,type_aloc,n,nhmin,ARRED=TRUE)
tmodel<-t1+t3

##############################################################################
return(list(nvars=length(xg),nconstraints=nrow(A),
            tgurobi=t2,bh=bh,Nh=Nh,vhx=S2$s2[1,],
            nh=S2$nh,nho=S2$nho,nha=S2$nha,
            cvtx.r.f=S2$cvtx.real.fator,cvtx.r=S2$cvtx.real,
            cvtx.i.f=S2$cvtx.integer.fator,cvtx.i=S2$cvtx.integer,
            cvtx.a.f=S2$cvtx.arred.fator,cvtx.a=S2$cvtx.arred
        ))

}





##################################################################################
################Global ###########################################################
Minimize_Variance_global<-function(n,L,nhmin)
{
   
t1<-proc.time()

tc<-choose(Nsr,2)
b<-n
desig<-"="
###################Restrição 29#######################################
cat("Generating Constraint 1   \n")
R1<-NULL
for(h in 1:L)  {for(j in nhmin:n) {R1<-c(R1,rep(j,tc))}}
#####################################################################

###################Restrição 30######################################
cat("Generating Constraint 2   \n")
R2<-matrix(0,ncol=(n-nhmin+1)*tc*L,nrow=L)
for(h in 1:L) {li=1+(h-1)*(tc*(n-nhmin+1));ls=tc*(n-nhmin+1)*h;R2[h,li:ls]<-1}
b<-c(b,rep(1,L))
desig=c(desig,rep("=",L))
#####################################################################

###################Restrição 31######################################
cat("Generating Constraint 3  \n")
R3<-rep(0,(n-nhmin+1)*tc*L)
li<-1
ls<-Nsr-2*(L-1)-1
for(j in 1:(n-nhmin+1))
 { R3[li:ls]<-1
   li<-j*tc+1
   ls<-j*tc+Nsr-1
}      
b<-c(b,1)
desig<-c(desig,"=")
####################################################################
###################Restrição 32#####################################
cat("Generating Constraint 4   \n")
Xc<-rep(0,tc)
Xc[which(combine[2,]==Nsr)]<-1
R4<-NULL
for(j in nhmin:n) {R4<-c(R4,Xc)}
b<-c(b,1)
desig<-c(desig,"=")


###################Restrição 33#####################################
####################################################################
cat("Generating Constraint 5  \n")
R5<-matrix(0,ncol=(n-nhmin+1)*L*tc,nrow=L-1)
ls<-(n-nhmin+1)*tc
li<-1
for(h in 1:(L-1)) 
 { 
   R5[h,li:(h*ls)]=-rep(combine[2,],(n-nhmin+1))
   R5[h,(ls*h+1):((h+1)*ls)]<-rep(combine[1,],(n-nhmin+1))
   li<-h*ls+1
}   
b<-c(b,rep(1,L-1))
desig<-c(desig,rep("=",L-1))


###################Restrição 34#####################################
####################################################################
cat("Generating Constraint 6   \n")
RX<-NULL
R6<-matrix(0,ncol=tc*(n-nhmin+1)*L,nrow=L)
li<-1;ls<-(n-nhmin+1)*tc;
for(h in 1:L) 
 {
 for(ns in nhmin:n) {Rt=c(Nij[t(combine)]-rep(ns,times=tc));RX=c(RX,Rt)}
 R6[h,li:ls]<-RX
 RX<-NULL
 li<-ls+1
 ls<-ls+((n-nhmin+1)*tc)
 }
b<-c(b,rep(0,L))
desig<-c(desig,rep(">",L))


A<-rbind(R1,R2,R3,R4,R5,R6)


fracn<-1/rep(rep(nhmin:n,each=tc),L)
nij<-rep(rep(Nij[t(combine)],(n-nhmin+1)),L)

Sij<-rep(rep(sij,(n-nhmin+1)),L)

fobj<-(nij^2*fracn)*Sij

   
model <- list()
model$obj<-fobj
model$modelsense<-"min"
model$rhs<-b
model$sense<-desig
model$vtype<- "B"
model$A <- A
params  <- list()
params$NumericFocus <-1
params$ObjScale=OS
params$TimeLimit<-TIMECPU
params$Heuristics=HE

t1<-(proc.time()-t1)[3]

cat("Applying Gurobi - Solver \n\n")
result <- gurobi(model, params)
t2<-result$runtime
nh<-rep(0,L)


xm<-round(result$x)
li=1
ls=tc*(n-nhmin+1)
Ni<-matrix(0,ncol=2,nrow=L)
for(h in 1:L)
 {u<-which(xm[li:ls]==1)
  p<-u%%tc
  if (p==0) {p=tc}
  Ni[h,]<-combine[,p]
  xm<-xm[-c(li:ls)]
 }   

Nh<-Nij[Ni]
bh<-rep(0,L-1)
for(h in 1:(L-1)) {bh[h]<-vx[Ni[h,2]]}
varx<-rep(0,L)
bz<-c(-Inf,bh,Inf)
for(h in 1:L) {varx[h]<-var(X[which(X>bz[h] & X<=bz[h+1])])}


xm<-round(result$x*R1)
nh<-xm[which(xm>0)]

cvf<-100*sqrt(sum(varx*Nh^2/nh*(1-nh/Nh)))/Tx
cvx<-100*sqrt(sum(varx*Nh^2/nh))/Tx


return(list(nvars=length(result$x),nconstraints=nrow(A),
            nh=nh,Nh=Nh,bh=bh,
            vhx=varx,cvtx.i.f=cvf,cvtx.i=cvx,
            tgurobi=t2)) 

}




OS=-1
HE=0.05
NF=1
############################Main Program###################################
library(MultAlloc)
library(gurobi)
library(parallel)
nucleos<-detectCores(logical=F)
clust<-makeCluster(nucleos)
###########################################################################
t1<-proc.time()
Tx<-sum(X)    
N<-length(X)  
X<-sort(X)   
Xr<-table(X) 
vx<-sort(unique(X))  
fx<-as.numeric(Xr) 
Nsr<-length(vx)
cat("Population Size = ",N,"\n")
cat("Number of different values for X = ",Nsr,"\n")
combine<-combn(Nsr,2)
Fij<-apply(combine,2,function(x) sum(fx[x[1]:x[2]])) 
Nij<-matrix(0,ncol=Nsr,nrow=Nsr)
Nij[t(combine)]<-Fij
Nij<-Nij+t(Nij)
sij<-parApply(cl=clust,combine,2,function(x) var(X[which(X>=vx[x[1]] & X<=vx[x[2]])]))
t1<-(proc.time()-t1)[3]
###############################################################################
if (type_aloc=="G")  
   {sg<-Minimize_Variance_global(n,L,nhmin)
    stopCluster(clust)
    return(sg)
   }
if (type_aloc %in% c("N","P","U"))
   {sg<-Minimize_Variance_NPU(type_aloc,nhmin)
    stopCluster(clust)
    return(sg)
   }
}






