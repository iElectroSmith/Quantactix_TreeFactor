


//#include <RcppArmadillo.h>

#include <csignal>

#include "common.h"
#include "state.h"
#include "APTree.h"
#include "model.h"
#include "json_io.h"

 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
//Rcpp::List 
   void    
           TreeFactor_APTree_cpp( 
  /*                                  arma::vec R , 
                                  arma::vec Y , 
                                  arma::mat X , 
                                  arma::mat Z , 
                                  arma::mat H , 
                                  arma::vec portfolio_weight , 
                                  arma::vec loss_weight , 
                                  arma::vec stocks , 
                                  arma::vec months , 
                                  arma::vec unique_months , 
                                  arma::vec first_split_var , 
                                  arma::vec second_split_var , 
                                  size_t num_stocks , 
                                  size_t num_months , 
                                  size_t min_leaf_size = 100 , 
                                  size_t max_depth = 5 , 
                                  size_t num_iter = 30 , 
                                  size_t num_cutpoints = 4 , 
                                  double eta = 1.0 , 
                                  bool equal_weight = false , 
                                  bool no_H = false , 
                                  bool abs_normalize = false , 
                                  bool weighted_loss = false , 
                                  bool stop_no_gain = false , 
                                  double lambda_mean = 0 , 
                                  double lambda_cov = 0 */
           ) ;
