# Likelihood Function for MS(k) Regression

MS_Var_Likelihood_Output<-function(params){
k<-2
paramstmp<-params[1:(length(params)-k*k)]
paramstmp1<-params[1:(length(params)-k*k)]
ifelse(abs(paramstmp<.001),paramstmp<-params*10000,paramstmp<-paramstmp1)

#   Calculation of some preliminary variables

nrdm<-dim(dep);
nr<-nrdm[1]
n_dep<-nrdm[2]

Param_Coeffs<-params[((k*n_dep*k)+1):((k*n_dep*k)+k*n_dep)]
Indep_S_param<-array(1,dim=c(nr,length(Param_Coeffs)))
for (i in 1:length(Param_Coeffs)) {
Indep_S_param[,i]<-Indep_S_param[,i]*Param_Coeffs[i]
}
Param_VarCoVar<-params[1:(k*n_dep*k)]
Param_VarCoVar<-matrix(Param_VarCoVar,n_dep,n_dep*k)
Param_TransProbs<-matrix(params[(length(params)-(k*k)+1):length(params)],k,k)
Param_TransProbs[1,2]<-1-Param_TransProbs[2,2]
Param_TransProbs[2,1]<-1-Param_TransProbs[1,1]
Param_TransProbstmp<-1-Param_TransProbs
Coeff_Max<-array(max(Param_TransProbs),dim=c(k,k))
Coeff_Min<-array(min(Param_TransProbs),dim=c(k,k))
Coeff_Init<-rbind(c(.5,.5),c(.5,.5))
ifelse(Coeff_Max>.99,Param_TransProbs<-Coeff_Init,Param_TransProbs<-Param_TransProbstmp)
ifelse(Coeff_Min<0,Param_TransProbs<-Coeff_Init,Param_TransProbs<-Param_TransProbstmp)

# Constrain the variance Covariance Matrix.
for (i in 1:k){
Cov_Mat<-Param_VarCoVar[,((i-1)*n_dep+1):((i-1)*n_dep+n_dep)]
for (r in 1:n_dep){
for (s in 1:n_dep){
Cov_Mat[r,s]<-Cov_Mat[s,r]
}
}
Param_VarCoVar[,((i-1)*n_dep+1):((i-1)*n_dep+n_dep)]<-Cov_Mat
}


Mean=array(0,dim=c(nr,k*n_dep));
e=array(0,dim=c(nr,k*n_dep));
n=array(0,dim=c(nr,k));

# Vectorized main engine

for (i in 0:(k-1)){
    Mean[,(((i*n_dep)+1):((i*n_dep)+n_dep))]=Indep_S_param[,(((i*n_dep)+1):((i*n_dep)+n_dep))]; # Conditional Mean for each state
    e[,(((i*n_dep)+1):((i*n_dep)+n_dep))]=dep-Mean[,(((i*n_dep)+1):((i*n_dep)+n_dep))]; # F_Probsrror for each state
} 

CovDep<-cov(dep)
DMult<-1
TMult<-array(0,dim=c(nr,k))
Det_M<-array(0,dim=c(1,k))

for (p in 1:nr) {
Res<-e[p,]

for (i in 1:k) {
inv_mat<-(Param_VarCoVar[,((i-1)*i+1):((i-1)*i+k)])
inv_mattmp<-(Param_VarCoVar[,((i-1)*i+1):((i-1)*i+k)])
d_inv_mat<-det(Param_VarCoVar[,((i-1)*i+1):((i-1)*i+k)])
d_inv_mat<-array(d_inv_mat,dim=c(k,k))
ifelse(d_inv_mat==0,inv_mat<-CovDep,inv_mat<-inv_mattmp)
ifelse(is.nan(inv_mat),inv_mat<-CovDep,inv_mat<-inv_mattmp)
inv_mat<-solve(inv_mat)
Det_M[i]<-det(inv_mat)
Res_S<-Res[((i-1)*i+1):((i-1)*i+k)]
MultA<-inv_mat*Res_S
MultB<-1
for (u in 1:k) {
s<-sum(MultA[,u])
MultB[u]<-s
}
Mult<-sum(MultB*Res_S)
DMult[i]<-Mult;
}
eta<-(1/sqrt(2*pi))^(n_dep)*sqrt(abs(Det_M)) * exp(-0.5*DMult);
n[p,]<-eta
}




    F_Probs=array(0,dim=c(nr,k));
    f=array(0,dim=c(nr,1));
    
    # Setting up first probs of F_Probs
    
    F_Probs[1,]=array(1/k,dim=c(k,1));

   
 for (i in 2:nr){
    F_Probsadj<-matrix(F_Probs[i-1,],k,1)
    AdjDen<-Param_TransProbs;
    densum<-array(1,dim=c(1,k))
    nadj<-matrix(n[i],k,1)
    for (r in 1:k){
    AdjDen[,r]<-Param_TransProbs[,r]*F_Probsadj[r]
}
   for (r in 1:k){
    densum[r]<-sum(AdjDen[r,])
}


        f[i]=sum(densum*n[i,]) # MS Filter equation
        F_Probs[i,]=(densum*n[i,])/f[i];   # MS Filter equation for probabilities

}

# Negative sum of log likelihood for fmincon (fmincon minimizes the function)

LL<-(-(sum(log(f[2:nr]))));
LL[is.nan(LL)]<-Inf
LL1<-LL
ifelse(max(Param_TransProbs)>.99,LL<-Inf,LL1)
ifelse(min(Param_TransProbs)<0,LL<-Inf,LL1)
LL
Coeff.Sparam<-matrix(Param_Coeffs,k,k)

# Format Output Structure. 
# Filtered Probabilities.
state1<-F_Probs[,1]
state2<-F_Probs[,2]
F_Probs<-cbind(state1,state2)

# Transition Matrix.
state1<-Param_TransProbs[,1]
state2<-Param_TransProbs[,2]
Param_TransProbs<-cbind(state1,state2)

# Coefficients.
state1<-Coeff.Sparam[,1]
state2<-Coeff.Sparam[,2]
Coeff.Sparam<-cbind(state1,state2)

# State Densities.
state1<-n[,1]
state2<-n[,2]
n<-cbind(state1,state2)

ModelOutput<-list("filtered probs"=F_Probs,residuals=e,"transition matrix"=Param_TransProbs,Coefficients=Coeff.Sparam,"log-likelihood"=LL,densities=n)
return(ModelOutput)
}