# PERFORM JOINT SIMULATION STUDY. 

# Generate K = 3 classes
K=3
n.vals = c(150,200,300)
p=100
N=100 # Number of simulations
fracs.disagreements = c(0,0.05,0.1,0.2,0.5,1)


# CASE 1: 100% edge agreement ------------------------------------------
set.seed(12)
seeds.1=sample(1:1000,N)
res.1 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=0,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.1)

# CASE 2: 95% edge agreement ------------------------------------------
set.seed(123)
seeds.2=sample(1:1000,N)
res.2 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=0.05,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.2)

# CASE 3: 90% edge agreement ------------------------------------------
set.seed(12344)
seeds.3=sample(1:1000,N)
res.3 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=0.1,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.3)

# CASE 4: 80% edge agreement ------------------------------------------
#set.seed(12344)
set.seed(1234567)
seeds.4=sample(1:1000,N)
res.4 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=0.2,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.4)

# CASE 5: 50% edge agreement ------------------------------------------
set.seed(123456)
seeds.5=sample(1:1000,N)
res.5 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=0.5,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.5)

# CASE 6: 0% edge agreement ------------------------------------------
set.seed(1234567)
seeds.6=sample(1:1000,N)
res.6 = JoStARS_simulation(K,n.vals,p,N=N,frac.disagreement=1,ebic.gamma=0,scale=T,stars.thresh = 0.1,lambda2.max=0.1,lambda2.init=0.01,seeds=seeds.6)

# PRINT RESULTS -------------------------------------------------------
res.list=list(res.1,res.2,res.3,res.4,res.5,res.6)
print_results_JGL(res.list,fracs.mutated=fracs.disagreements,show.specificity = T)





