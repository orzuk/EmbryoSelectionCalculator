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
// [[Rcpp::export]]
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
	pareto["pareto.inds"] = pareto_inds;
	pareto["pareto.X"] = X_pareto;
 	return(pareto); 
}


// Add a vector to a matrix of already pareto optimal vectors 
// [[Rcpp::export]]
NumericMatrix update_pareto_optimal_vecs_rcpp(NumericVector x, NumericMatrix X_mat)
{
    if(!is_pareto_optimal_rcpp(x, X_mat))  // new vector is dominated 
        return X_mat; // unchanged 

    // new vector is pareto. Check which vectors are dominated by it:
    double epsilon = 0.0000000000000000000000001;    
    long n = X_mat.nrow();
    long i, ctr=0;
    NumericVector is_pareto(n);
    for(i=0; i<n; i++)
        is_pareto[i] = (max(X_mat(i,_)-x) < -epsilon); 
    long n_pareto = sum(is_pareto)+1; // include new vector 
    NumericVector pareto_inds(n_pareto-1);
    for(i=0; i < n; i++)
        if(is_pareto[i])
            pareto_inds[ctr++] = i;    

    NumericMatrix X_pareto(n_pareto, X_mat.ncol());
    for(i=0; i<n_pareto-1; i++)
        X_pareto.row(i) = X_mat.row(pareto_inds[i]);
    X_pareto.row(n_pareto-1) = x;
    return X_pareto;
}


// Unite two lists of pareto-optinal vectors. Keep only pareto optimals in the joint list - 
// [[Rcpp::export]]
List union_pareto_optimal_vecs_rcpp(NumericMatrix X_mat1, NumericMatrix X_mat2)
{
    long n1 = X_mat1.nrow();
    long n2 = X_mat2.nrow();
    NumericVector is_pareto1(n1);
    NumericVector is_pareto2(n2);
    long i; 
    long ctr = 0;

    // First determine who is pareto-optimal from the first list 
    for(i=0; i<n1; i++)
        is_pareto1[i] = is_pareto_optimal_rcpp(X_mat1(i,_), X_mat2); 
    for(i=0; i<n2; i++)
        is_pareto2[i] = is_pareto_optimal_rcpp(X_mat2(i,_), X_mat1); 
    long n_pareto1 = sum(is_pareto1); // include new vector 
    long n_pareto2 = sum(is_pareto2); // include new vector 

    NumericVector pareto_inds1(n_pareto1);
    NumericVector pareto_inds2(n_pareto2);
    for(i=0; i < n1; i++)
        if(is_pareto1[i])
            pareto_inds1[ctr++] = i;    
    ctr = 0; 
    for(i=0; i < n2; i++)
        if(is_pareto2[i])
            pareto_inds2[ctr++] = i;    
    NumericMatrix X_pareto(n_pareto1+n_pareto2, X_mat1.ncol());  // get pareto-optimal vectors from both sets 
    for(i=0; i<n_pareto1; i++)
        X_pareto.row(i) = X_mat1.row(pareto_inds1[i]);
    for(i=0; i<n_pareto2; i++)
        X_pareto.row(i+n_pareto1) = X_mat2.row(pareto_inds2[i]);
  
    List pareto;
	pareto["pareto.inds1"] = pareto_inds1;
	pareto["pareto.inds2"] = pareto_inds2;
	pareto["pareto.X"] = X_pareto;
 	return pareto; 
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
