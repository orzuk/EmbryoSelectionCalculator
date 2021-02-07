// Utilities for chromosomal selection 
#include <iostream>
#include <math.h>
// #include <map>
#include <random>

// #include <Rcpp.h> // for including R with Rcpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// #include "utilities_ToR.h"


using namespace std;
using namespace Rcpp;  // for including R with Rcpp 


// cpp code: 

// [[Rcpp::export]]
long is_pareto_optimal_rcpp(NumericVector x, NumericMatrix X_mat)
{
    double epsilon = 0.0000000000000000000000001;

    if(X_mat.nrow() == 0)
        return(TRUE);
    long n_row = X_mat.nrow();
  
    for(long i=0; i<n_row; i++)
    {
        if(max(x - X_mat(i,_)) < -epsilon)
            return(FALSE);
    }
    return(TRUE);

  // check if row maxs are good: 
//  return( min(rowMaxs(t(replicate(n.row, x))+epsilon - X.mat, value=TRUE)) >= 0 )
}

// Get only the pareto optimal vectors out of X_mat 
List get_pareto_optimal_vecs_rcpp(NumericMatrix X_mat)
{
    long n = X_mat.nrow();
    long i, ctr=0;
    NumericVector is_pareto(n);
    for(i=0; i < n; i++)
        is_pareto[i] = is_pareto_optimal_rcpp(X_mat(i,_), X_mat); //   max(colMins(t(replicate(n, X[i,]))+epsilon - X, value=TRUE)) >= 0
    long n_pareto = sum(is_pareto);
    NumericVector pareto_inds(n_pareto);
    for(i=0; i < n; i++)
        if(is_pareto[i])
            pareto_inds[ctr++] = i;    

    NumericMatrix X_pareto(n_pareto, X_mat.ncol());
    for(i=0; i<n_pareto; i++)
        X_pareto.row(i) = X_mat.row(pareto_inds[i]);
    List pareto;
	pareto["pareto_inds"] = pareto_inds;
	pareto["X_pareto"] = X_pareto;
 	return(pareto); 
}



/**
// R code: 
# Check if vector x is not dominated by any vector in X.mat
is_pareto_optimal <- function(x, X.mat)
{
  epsilon = 0.0000000000000000000000001
  if(isempty(X.mat))
    return(TRUE)
  if(is.null(dim(X.mat)))
    n.row <- 1
  else
    n.row <- dim(X.mat)[1]
  return( min(rowMaxs(t(replicate(n.row, x))+epsilon - X.mat, value=TRUE)) >= 0 )
}

# Extract only pareto-optimal vectors in a matrix
get_pareto_optimal_vecs <- function(X.mat)
{
  n = dim(X.mat)[1]
  if(is.null(n)) # here X.mat is a vector - only one vector 
    return(list(pareto.X=X.mat, pareto.inds=1))
  is.pareto = rep(0,n)
  for(i in 1:n)
    is.pareto[i] = is_pareto_optimal(X.mat[i,], X.mat) #   max(colMins(t(replicate(n, X[i,]))+epsilon - X, value=TRUE)) >= 0
  pareto.inds <- which(as.logical(is.pareto))
  return(list(pareto.X=X.mat[pareto.inds,], pareto.inds=pareto.inds))
}
**/
