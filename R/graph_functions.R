# SHOULD BE RUN BEFORE OTHER SCRIPTS
library(igraph)
library(huge)
library(glasso)
library(MASS)
library(mvtnorm)
library(JGL)
library(Matrix)
library(foreach)
library(doParallel)

sparsity = function(g){
  # Function for calculating the sparsity of a graph. Checks if diagonal elements are zero.
  g=g+0 # Added this later
  p=dim(g)[1]
  if(g[1,1]==0) return(sum(g!=0)/(p^2-p))
  else return((sum(g!=0)-p)/(p^2-p) )
}

# Function for assigning colors to each vertex according to their degree. Returns vector of color codes. 
color.by.degree = function(adj.mat){
  # adj.mat is the igraph adjacency matrix object.
  # Luminance is in [0,100], where low degrees have high luminance. 
  luminance= 100-(degree(adj.mat)-min(degree(adj.mat)))/(max(degree(adj.mat))-min(degree(adj.mat)))*100 
  return( hcl(c=80,l=luminance,h=0) ) # Pink colors, darker means higher degree. 
}

gaussianloglik = function(sample.cov,theta,n){
  # Find the loglikelihood. Sample.cov is empirical covariance matrix, theta is the estimated precision matrix. n is number of observations
  p=nrow(theta)
  return(-p*n*log(2*pi)/2 + n*log(det(theta))/2 -n*sum(diag(sample.cov%*%theta))/2)
}


eBIC = function(sample.cov,theta,n,gamma){
  # Find eBIC score.
  # sample.cov is the empirical covariance matrix, theta is the estimated precision matrix, n no of observations
  # gamma is to be tuned. 
  p=nrow(theta)
  theta2=theta
  diag(theta2) = rep(0,p)
  d=sum(theta2!=0)/2 # Number of edges
  loglik = -p*n*log(2*pi)/2 + n*log(det(theta))/2 -n*sum(diag(sample.cov%*%theta))/2
  ebic = -2*loglik + d*log(n) + 4*d*gamma*log(p)
  return(ebic)
}


# Find Matthews correlation coefficient for estimated graph
MCC = function(g,g.hat){
  p = nrow(g[,])
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat) = rep(0,p) 
  tp = sum(g.hat ==1 & g ==1)/10 # True positives. Divide by 10 to avoid integer overflow. 
  fp = sum(g.hat ==1 & g ==0)/10 # False positives
  tn = (sum(g.hat == 0 & g == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat == 0 & g == 1)/10 # False negatives
  return((tp*tn - fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
}

confusion.matrix = function(g,g.hat){
  # Newly changed: all must be halved! Since we look at the whole adjacency matrix. 
  p = nrow(g[,])
  g = g[,]
  g.hat = g.hat[,]
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat[,]) = rep(0,p) 
  tp = sum(g.hat[,] ==1 & g[,] ==1)/10 # True positives. Divide by 10 to avoid integer overflow. 
  fp = sum(g.hat[,] ==1 & g[,] ==0)/10 # False positives
  tn = (sum(g.hat[,] == 0 & g[,] == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat[,] == 0 & g[,] == 1)/10 # False negatives
  return(matrix(10*c(tp,fp,fn,tn),nrow=2,byrow=T)/2)
}

recall = function(g,g.hat){
  # Sensitivity or True Positive Rate. What fraction of true edges are discovered?
  conf.mat = confusion.matrix(g,g.hat)
  if(conf.mat[1,1]==0 & conf.mat[2,1]==0) return(1) # Avoid zero division
  else return(conf.mat[1,1]/(conf.mat[1,1]+conf.mat[2,1]))
}

precision = function(g,g.hat){
  # Positive Predictive Value. What fraction of the estimated edges are actually true?
  # Equal to 1-FDR
  conf.mat = confusion.matrix(g,g.hat)
  if(conf.mat[1,1]==0 & conf.mat[1,2]==0) return(1) # Avoid zero division
  else return(conf.mat[1,1]/(conf.mat[1,1]+conf.mat[1,2]))
}

specificity = function(g,g.hat){
  # Specificity. How many of the edges that should have been negative are? TN/(TN+FP)
  conf.mat = confusion.matrix(g,g.hat)
  if(conf.mat[1,2]==0 & conf.mat[2,2]==0) return(1) # Avoid zero division
  else return(conf.mat[2,2]/(conf.mat[1,2]+conf.mat[2,2]))
}

FDR = function(g,g.hat){
  # False Discovery Rate
  return(1-precision(g,g.hat))
}




mutate.graph= function(graph,fraction){
  # Mutate a given fraction of the edges of a graph. 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  prec.mat = graph$omega
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = graph$sigma
  adj.mat = graph$theta
  data=graph$data
  p = ncol(graph$omega)
  
  adj.mat.upper = adj.mat
  adj.mat.upper[lower.tri(adj.mat.upper)]=0
  diag(adj.mat.upper) =0
  edges = which(adj.mat.upper==1,arr.ind=T) # Edge pairs.
  n.mutations = floor(nrow(edges)*fraction)
  
  if(n.mutations==0 | is.na(n.mutations)){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    ans$data = data
    return(ans)
  }
  
  edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
  edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] # The nodes to 'leave out'
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    id.stay = edges.to.change[i,1]
    id.remove = edges.to.change[i,2]
    id.add=nodes.add[i]
    # Swap prec mat elements in rows. Then cols, the order does not matter!
    prec.mat[id.stay,id.add] = tmp.prec[id.stay,id.remove]
    prec.mat[id.stay,id.remove] = tmp.prec[id.stay,id.add]
    prec.mat[id.add,id.stay] = tmp.prec[id.remove,id.stay]
    prec.mat[id.remove,id.stay] = tmp.prec[id.add,id.stay]
    # swap adj mat rows
    adj.mat[id.stay,id.add] = tmp.adj[id.stay,id.remove]
    adj.mat[id.stay,id.remove] = tmp.adj[id.stay,id.add]
    adj.mat[id.add,id.stay] = tmp.adj[id.remove,id.stay]
    adj.mat[id.remove,id.stay] = tmp.adj[id.add,id.stay]
  }
  ans = list()
  ans$cov.mat=solve(prec.mat)
  ans$cov.mat[which(abs(ans$cov.mat)<1e-4,arr.ind = T)] = 0
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  # Generate new data
  ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
  return(ans)
}




make.df.of.changed.edges = function(mat1,mat2,colnam){
  changes.ind = which(mat1!=0 & mat2 ==0,arr.ind=T)
  changes= data.frame(as.character(colnam[changes.ind[,1]]),as.character(colnam[changes.ind[,2]]))
  df=data.frame(t(apply(changes,1,sort)))
  df = unique(df)
  rownames(df) = 1:nrow(df)
  colnames(df) = c('Gene1','Gene2')
  return(df)
}

matrix.distance = function(mat1,mat2){
  # Function for evaluating distance between two matrices. 
  # Simplified: disregards the possibility of overlap. 
  p = nrow(mat1)
  mat1[which(abs(mat1)<10^(-4),arr.ind=T)] = 0 # Disregard almost-zero elements
  mat2[which(abs(mat2)<10^(-4),arr.ind=T)] = 0 # Disregard almost-zero elements
  m1 = cov2cor(as.matrix(forceSymmetric(mat1)))
  m2 = cov2cor(as.matrix(forceSymmetric(mat2)))
  observed.dist = sum(abs(abs(m1)-abs(m2)))
  expected.dist = (sum(abs(m1))+sum(abs(m2))-2*p)
  if(expected.dist==0) return(0)
  else return(observed.dist/expected.dist)
}
