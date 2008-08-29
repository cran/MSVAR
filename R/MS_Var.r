MS_Var<-function(dep){
# Set-up Dependent Variables.

dim_data<-dim(dep)
dep_n<-matrix(1,dim_data[1],2)
dep_n[,1]<-dep[,1]
dep_n[,2]<-dep[,2]

dep<<-dep_n
indep<-dep/dep
indep[56,2]<-1
k<-2
S<-c(2)

nrdm<-dim(dep);
nr<-nrdm[1];
n_dep<-nrdm[2]

n_indep<-length(indep);
n_S<-sum(S);
n_nS<-n_indep-n_S;
S_S<-dep/dep;

count<-0;
countS<-0;

# Calculate Initial Parameters.

S_S<-1
indep_S<-indep

indep_ols<-1
for (i in 1:n_dep) {
param_ols_S<-lm(dep[,i]~indep[,i]); # simple Ols for param0 of switcing variables
param_ols_S<-coefficients(param_ols_S)
param_ols_S<-param_ols_S[1]
indep_ols[i]<-param_ols_S;
}

param0_indep_STmp<-c(indep_ols)
param0_indep_S<-c(indep_ols)
for (i in 2:k) {
c<-
p<-i-1
param0_indep_S<-c(param0_indep_S,-param0_indep_STmp)
}

# Calculate Initial Sigma and Transition Matrix Values.
sigma<-cov(dep)
sigma<-rep(sigma,k)
tranmatrix<-array(.5,dim=c(1,k*k))

        param0<-c(sigma,param0_indep_S,tranmatrix)
# Lower Bounds for the Variance Covariance Matrix.
sigmalB<-array(-Inf,dim=c(1,n_dep*n_dep))
for (i in 1:n_dep) {
SigmaTmp<-sigmalB[((i-1)*i+1):((i-1)*i+n_dep)]
SigmaTmp[i]<-0
sigmalB[((i-1)*i+1):((i-1)*i+n_dep)]<-SigmaTmp
}
  
# lower and upper bounds
uB<-c(array(Inf,dim=c(1,length(sigma))), array(Inf,dim=c(1,length(param0_indep_S))), array(.99,dim=c(1,length(tranmatrix))))       
lB<-c(rep(sigmalB,2), array(-Inf,dim=c(1,length(param0_indep_S))), array(0,dim=c(.01,length(tranmatrix))))       

# Run nlminb.

print('Running Optimization with nlminb')
res<-nlminb(param0,MS_Var_Likelihood,gradient=NULL,hessian=NULL,scale=1,control=list(iter.max=5000),lower=lB, upper = uB)
print('Optimization Completed')

# Calculate Output for Optimal Parameters.
ModelOutput<-MS_Var_Likelihood_Output(res$par);

# calculating smoothed probabilities

Prob_t_1<-array(0,dim=c(nr,k));
Prob_t_1[1,]=array(1/k,dim=c(k,1));

for (i in 2:nr){
   Eadj<-matrix(ModelOutput$"filtered probs"[i-1,],k,1)
    AdjDen<-ModelOutput$"transition matrix";
    densum<-array(1,dim=c(1,k))
    nadj<-matrix(ModelOutput$densities[i],k,1)
    for (r in 1:k){
    AdjDen[,r]<-ModelOutput$"transition matrix"[,r]*Eadj[r]
}
   for (r in 1:k){
    densum[r]<-sum(AdjDen[r,])
}
Prob_t_1[i,]<-densum
}

smooth_value<-array(0,dim=c(1,k))
filtProb<-ModelOutput$"filtered probs";

P<-(ModelOutput$"transition matrix");

smoothProb<-array(0,dim=c(nr,k));
smoothProb<-filtProb;  # last observation for starting filter

for (i in (nr-1):1)  {   # work backwards in time for smoothed probs
    
    for (j1 in 1:k){
        for (j2 in 1:k){
            
            smooth_value[1,j2]<-smoothProb[i+1,j2]*filtProb[i,j1]*P[j2,j1]/Prob_t_1[i+1,j2];
            
}
        smoothProb[i,j1]<-sum(smooth_value);
    }
    
}

# Plot Filtered Probabilities and Return Output.
state1<-matrix(ModelOutput$"filtered probs"[,2],nr,1)
plot(state1,type='l',main="filtered probabilities for state 2")

# Plot Smoothed Probabilities and Return Output.
state1<-(smoothProb[,2])
plot(state1,type='l',main="smoothed probabilities for state 2")

# Smoothed Probabilities.
state1<-smoothProb[,1]
state2<-smoothProb[,2]
smoothProb<-cbind(state1,state2)

ModelOutput$"smoothed Probs"<-smoothProb
return(ModelOutput)
}