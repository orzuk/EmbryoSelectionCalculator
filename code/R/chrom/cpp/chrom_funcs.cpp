// Utilities for chromosomal selection 
// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <iostream>
#include <math.h>
// #include <map>
#include <random>
#include <boost/math/distributions/normal.hpp>

// C:\Program Files\boost\boost_1_75_0\boost\math\distributions
// #include "C:/Program Files/boost/boost_1_75_0/boost/math/distributions/normal.hpp"

// #include <Rcpp.h> // for including R with Rcpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// #include "utilities_ToR.h"

/**
namespace boost{ namespace math{ 
   
template <class RealType = double, 
          class Policy   = policies::policy<> >
class normal_distribution;

typedef normal_distribution<> normal;

template <class RealType, class Policy>
class normal_distribution
{
public:
   typedef RealType value_type;
   typedef Policy   policy_type;
   // Construct:
   normal_distribution(RealType mean = 0, RealType sd = 1);
   // Accessors:
   RealType mean()const; // location.
   RealType standard_deviation()const; // scale.
   // Synonyms, provided to allow generic use of find_location and find_scale.
   RealType location()const;
   RealType scale()const;
};

}} // namespaces


template <class RealType, class Policy>
RealType quantile(const Distribution-Type<RealType, Policy>& dist, const RealType& p);
**/

using namespace std;
using namespace Rcpp;  // for including R with Rcpp 
using namespace arma; // for armadillo
using namespace boost; // for boost mathematical packages 

// cpp code: 

// Inner product of two vectors  
// [[Rcpp::export]]
double my_inner_mult(NumericVector U, NumericVector V)
{
    double r = 0; 
	for (long i = 0; i < U.length(); i++)
		r += U[i]*V[i]; 
	return r;
}


// Multiply matrix by a vector 
// [[Rcpp::export]]
NumericVector my_mv_mult(NumericMatrix A, NumericVector V)
{
	NumericVector r(A.nrow());
	long i, j;

	for (i = 0; i < A.nrow(); i++)
	{
		r[i] = 0.0;
		for (j = 0; j < A.ncol(); j++)
			r[i] += A(i, j) * V[j];
	}
	return r;
}


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
            pareto_inds[ctr++] = i + 1; // 1-based indexing (R)    

    NumericMatrix X_pareto(n_pareto, X_mat.ncol());
    for(i=0; i<n_pareto; i++)
        X_pareto.row(i) = X_mat.row(pareto_inds[i]-1);
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
            pareto_inds1[ctr++] = i;    // NO one-based indices (for R)
    ctr = 0; 
    for(i=0; i < n2; i++)
        if(is_pareto2[i])
            pareto_inds2[ctr++] = i;      // NO one-based indices (for R) 
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



// The loss for vectors 
// Input: 
// X_c_mat - a matrix with vectors of length C
// loss_type - what loss to compute 
// loss_params - parameters of loss 
//
// Output: 
// loss.vec - a vector of losses for each row of X
// [[Rcpp::export]]
NumericVector loss_PS_mat_rcpp(NumericMatrix X_c_mat, string loss_type, List loss_params)
{
  long C = X_c_mat.nrow();
  long T = X_c_mat.ncol();
  long i, j;

  Rcout << "Inside Compute loss PS mat" << endl; 
  NumericVector loss_vec(C);
  if(loss_type == "quant")
    loss_vec = my_mv_mult(X_c_mat, as<NumericVector>(loss_params["theta"])); // scalar product 
  
  if((loss_type == "stabilizing") || (loss_type == "balancing"))
  {
    NumericMatrix X_c_mat2(C,T);
    for(i=0; i<C; i++)
        for(j=0; j<T; j++)
            X_c_mat2(i,j) = X_c_mat(i,j)*X_c_mat(i,j);
    loss_vec = my_mv_mult(X_c_mat2, as<NumericVector>(loss_params["theta"])); // scalar product 
  }
  
  if(loss_type == "disease")  // weighted disease probability 
  {
      Rcout << "Inside Compute loss PS mat Disease" << endl; 

    boost::math::normal s; 
    NumericVector z_K = clone(as<NumericVector>(loss_params["K"])); //  = as<double>(loss_params["K"]);
    NumericVector one_over_sqrt_h_ps = clone(as<NumericVector>(loss_params["h.ps"]));
    Rcout << "Prevalences are: " << z_K << endl; 
    for(i = 0; i < T; i++)  
    {
//      Rcout << "Prev i = " << z_K[i] << endl; 
      z_K[i] = quantile(s, z_K[i]); // as<double>(loss_params["K"][i]));
 //     Rcout << "Quant i = " << z_K[i] << endl; 
      one_over_sqrt_h_ps[i] = 1.0 / sqrt(1-one_over_sqrt_h_ps[i]);
    }
  //  Rcout << "z_K are now: " << z_K << endl; 
  //  Rcout << "one_over_sqrt_h_ps: " << one_over_sqrt_h_ps << endl; 
//    if(min(z_K) < 0)
//   {
//      Rcout << "Error! Z_k changed!!";
//      return(loss_vec);
//    }


//    Rcout << "K = " << K << endl; 
//    double z_K = quantile(s, as<double>(loss_params["K"]));
    Rcout << "Loop Compute loss PS mat Disease z_k = " << z_K << endl; 
    for(i = 0; i < C; i++) // loop on chromosomes 
    {
      NumericVector r = (z_K - X_c_mat(i,_)); //  / sqrt(1-as<double>(loss_params["h.ps"]));
//      Rcout << "i = " << i << ", r = " << r << endl;   
//      Rcout << "Vec: " <<  as<NumericVector>(loss_params["theta"]) << endl;
      NumericVector v = as<NumericVector>(loss_params["theta"]);
//      Rcout << "v = " << v << endl;
      NumericVector w = r * one_over_sqrt_h_ps;
      Rcout << "w = " << w << endl;
      for(j = 0; j < T; j++) // apply normal cdf 
        w[j] = cdf(s, w[j]);
      Rcout << "pnorm(w) = " << w << endl; 
//      double u = my_inner_mult(    (z_K - X_c_mat(i,_)) * one_over_sqrt_h_ps, as<NumericVector>(loss_params["theta"]) ); //  / sqrt(1-as<double>(loss_params["h.ps"])), v ); 
//      Rcout << "u = " << u << endl; 
      loss_vec[i] = my_inner_mult(w,  as<NumericVector>(loss_params["theta"])); // my_inner_mult(  (z_K - X_c_mat(i,_)) / sqrt(1-as<double>(loss_params["h.ps"])), as<NumericVector>(loss_params["theta"]) )); 
    }
  }
  return loss_vec;
}  

/** R code below: 
loss_PS_mat <- function(X.c.mat, loss.type, loss.params)
{
  n <- dim(X.c.mat)[1]
  
  #  X.c <- compute_X_C_mat(X, C)
  if(loss.type == "quant")
    loss.vec <- X.c.mat %*% loss.params$theta # scalar product 
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    loss.vec <- (X.c.mat*X.c.mat) %*% loss.params$theta    
  }
  
  if(loss.type == 'disease')  # weighted disease probability 
  {
    z.K <- qnorm(loss.params$K)
    loss.vec <- t(pnorm( (z.K-t(X.c.mat))/sqrt(1-loss.params$h.ps)  )) %*% loss.params$theta 
  }
  return(loss.vec)
}  
**/




/*** Set al B&B In comment !!! ***/

// New: implement all B&B inside cpp to save time 
// [[Rcpp::export]]
List optimize_C_branch_and_bound_rcpp(arma::cube X, string loss_type, List loss_params)
{
    long i, j, c; 
    long L1, L2;

//    Dimension dim = X.attr("dim");
    long M = X.n_rows; // dim[0];
    long C = X.n_cols; // dim[1];
    long T = X.n_slices; // dim[2];

    Rcout << "START B&B CPP: " << endl;
//    NumericMatrix X0(wrap(X(span(0), span(), span())));
    NumericMatrix X0(C, T); // X0(wrap(X.row(0)));

    // copy using a loop: 
    for(i = 0; i < C; i++)
        for(j = 0; j < T; j++)
            X0(i,j) = X(0,i,j);

    X.row(0).print(); 
    Rcout << "printed 0 row" << endl; 
//    mat M(wrap(X.row(0))); 
 //   Rcout << "Converted to matrix " << endl; 
  //  Rcout << M << endl; 
//    NumericMatrix X0(wrap(X.row(0)));


    Rcout << "START get pareto CPP: " << endl;
    Rcout << "X0 = " << endl << X0 << endl; 
    Rcout << "X0 Dims: " << X0.nrow() << ", " << X0.ncol() << endl; 

//    long T = X.n_slices; // dim[2];
    List par_X = get_pareto_optimal_vecs_rcpp(X0); // wrap(X.slice(0)));  // Save only Pareto-optimal vectors . Needs fixing 
    Rcout << "Copy list outputs c:" << endl;
    NumericVector c_vec = as<NumericVector>(par_X["pareto.inds"]);
    Rcout << "Copied vector!" << endl; 
    NumericMatrix cur_c(c_vec.length(), 1);
    cur_c(_,0) = c_vec; // par_X["pareto.inds"]);
    Rcout << "Copy list outputs X:" << endl;
    NumericMatrix cur_X = as<NumericMatrix>(par_X["pareto.X"]);
    long L = cur_X.nrow();

    Rcout << "Compute pareto, cur.c dims: " << cur_c.nrow() << endl; 
    Rcout << "Compute pareto, cur.X dims: " << cur_X.nrow() << ", " << cur_X.ncol() << endl; 

    NumericVector L_vec(M); // vector of zeros 
    L_vec[0] = L;
    NumericMatrix new_X; 
    NumericMatrix new_c;
    for(i=1; i<M; i++) //  in 2:M) # loop on block
    {
    //    Rcout << "Loop on blocks: " << i << endl;
        L = cur_X.nrow();
        new_X = cur_X;
        for(j = 0; j<L; j++)
            new_X(j,_) = new_X(j,_) + as<NumericVector>(wrap(X.subcube( span(i), span(0), span() )));  
    //    Rcout << "Copied X tensor blocks: " << i << endl;
        new_c = NumericMatrix(L, i+1);
        for(j = 0; j<i; j++)    
            new_c(_, j) = cur_c(_, j);
        new_c(_,i) = NumericVector(L, 1); // add 1 to last column 

        NumericMatrix temp_X(L, i+1);
        List temp_X_list;
        List union_X;
        for(c = 1; c < C; c++) // next add other vectors 
        {
    //        Rcout << "Copied X tensor blocks: c=" << c << endl;
            temp_X = cur_X; 
            for(j = 0; j < L; j++)
                temp_X(j,_) = temp_X(j,_)  + as<NumericVector>(wrap(X.subcube( span(i), span(c), span() )));  
            temp_X_list = get_pareto_optimal_vecs_rcpp(temp_X);
            union_X = union_pareto_optimal_vecs_rcpp(new_X, temp_X_list["pareto.X"]);
            L1 = as<NumericVector>(union_X["pareto.inds1"]).length();
            L2 = as<NumericVector>(union_X["pareto.inds2"]).length();
            new_X = as<NumericMatrix>(union_X["pareto.X"]);
            Rcout << "add c=" << c << endl;
            NumericMatrix add_c(L, i+1);
            for(j = 0; j < i; j++)
                add_c(_, j) = cur_c(_, j);
            add_c(_, i) = NumericVector(L, 1); // add 1 to last column 
    //        Rcout << "New c=" << c << endl;
    //        Rcout << "L1=" << L1 << " L2=" << L2 << endl;
            cur_c = clone(new_c); // save previous
            new_c = NumericMatrix(L1+L2, i+1);
    //        Rcout << "new_c_dim=" << new_c.nrow() << " , " << new_c.ncol() << endl;

    //        Rcout << "cur_c_dim_for_L1=" << cur_c.nrow() << " , " << cur_c.ncol() << endl;
    //        Rcout << "Pareto inds1: " << as<NumericVector>(union_X["pareto.inds1"]) << endl;
    //        Rcout << "Pareto inds2: " << as<NumericVector>(union_X["pareto.inds2"]) << endl; // 1 based indices for R!! 

            for(j = 0; j < L1; j++)
                new_c(j,_) = cur_c(as<NumericVector>(union_X["pareto.inds1"])[j],_); // self copying
            Rcout << "Finished first loop: L1=" << L1 << " L2=" << L2 << endl;
            for(j = 0; j < L2; j++)
            { 
    //          Rcout << " Loop j=" << j << endl;
              NumericVector vv = as<NumericVector>(union_X["pareto.inds2"]);
    //          Rcout << " Nooma j=" << j << endl;
              NumericVector nnn = as<NumericVector>(union_X["pareto.inds2"])[j];
    //         Rcout << " Nooma2 j=" << j << " L1+j=" << j << endl;
    //          Rcout << "new_c_dim_now=" << new_c.nrow() << " , " << new_c.ncol() << endl;
              new_c(L1+j,_) = nnn; 
    //          Rcout << "new_c_dim_again=" << new_c.nrow() << " , " << new_c.ncol() << endl;
              new_c(L1+j,_) = add_c(as<NumericVector>(union_X["pareto.inds2"])[j],_); // assignment from the same variable! wtf?
    //          Rcout << " Hey j=" << j << endl;
            }
    //        Rcout << "Finished second loop: L1=" << L1 << " L2=" << L2 << endl;
    
        }
        cur_X = new_X;
        cur_c = new_c;
        L_vec[i] = new_X.nrow();          
    } // end loop on blocks 

    Rcout << "Finished big loop" << endl; 
    // Finally find the cost-minimizer out of the Pareto-optimal vectors
    L = cur_X.nrow();
    NumericVector loss_vec = loss_PS_mat_rcpp(cur_X, loss_type, loss_params);
    Rcout << "Calculated loss vec" << endl; 
    long i_min = which_min(loss_vec);  // find vector minimizing loss 



//    (list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), 
//             loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X, pareto.opt.c = cur.c))
    List sol;

    sol["opt.X"] = cur_X(i_min,_);    
    sol["opt.c"] = cur_c(i_min,_);
    double opt_loss = min(loss_vec);
    sol["opt.loss"] = opt_loss; // min(loss_vec);
/***/    
    sol["loss.vec"] = loss_vec;
    sol["L.vec"] = L_vec; 
    sol["pareto.opt.X"] = cur_X; 
    sol["pareto.opt.c"] = cur_c; 
/***/   
    return sol;

}

/*** END OF B&B ***/

/** # R code:  A branch and bound algorithm for finding the X combination with minimal loss 
optimize_C_branch_and_bound <- function(X, loss.C, loss.params)
{
    

  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  
  par.X <- get_pareto_optimal_vecs_rcpp(X[1,,]) # Save only Pareto-optimal vectors . Needs fixing 
  cur.c <- t(t(par.X$pareto.inds))
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  L.vec <- rep(0, M)
  L.vec[1] = L
  for(i in 2:M) # loop on block
  {
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    if(i == M)
      print(paste0("B&B i=", i, " L=", L))
    new.X <- c()
    new.c <- c()
    
    # We know that the first vectors are pareto optimal
    if(L>1)
    {
      new.X <- sweep(cur.X, 2, X[i,1,], "+")
    } else
      new.X <- matrix(cur.X + X[i,1,], nrow=1) # check that it doesn't flip x
    new.c <- cbind(cur.c, rep(1, L) )

    # new version: create sums and take union
    for(c in 2:C)
    {
      temp.X <- sweep(cur.X, 2, X[i,c,], "+")
      temp.X <- get_pareto_optimal_vecs_rcpp(temp.X)
      union.X <- union_pareto_optimal_vecs_rcpp(new.X, temp.X$pareto.X)
      new.X <- union.X$pareto.X
      add.c <- cbind(cur.c, rep(c, L) )
      new.c <- rbind(new.c[union.X$pareto.inds1,], add.c[union.X$pareto.inds2,]) # need to modify here indices 
    }
    cur.X <- new.X
    cur.c <- new.c
    if(is.null(dim(new.X)[1]))
      L.vec[i] = 1
    else
      L.vec[i] = dim(new.X)[1]  
        
  } # end loop on blocks 
  
  # Finally find the cost-minimizer out of the Pareto-optimal vectors
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
  {
    L = 1
    loss.vec = loss_PS(cur.X, loss.C, loss.params)
  }  else
  {
    loss.vec <- loss_PS_mat(cur.X, loss.C, loss.params)
  }
  i.min <- which.min(loss.vec) # find vector minimizing loss 
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), 
              loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X, pareto.opt.c = cur.c))
}
**/ 
