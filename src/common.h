#ifndef GUARD_common_h
#define GUARD_common_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <omp.h>


#include <iostream>
#include <armadillo>
 

//#include "RcppArmadillo.h"
//#include "Rcpp.h"
//#include "omp.h"

#include <csignal>

using namespace std ;
using namespace arma ;
//using namespace Rcpp ;

#define LTPI 1.83787706640934536

std::ostream& operator<<( std::ostream& out , const std::vector<double>& v ) ;
std::ostream& operator<<( std::ostream& out , const std::vector<size_t>& v ) ;
std::ostream& operator<<( std::ostream& out , const std::vector<bool>& v ) ;
std::ostream& operator<<( std::ostream& out , const std::vector<std::vector<double>>& v ) ;
std::ostream& operator<<( std::ostream& out , const std::vector<std::vector<size_t>>& v ) ;

double fastLm( const arma::vec& y , const arma::mat& X ) ;

double fastLm_weighted( const arma::vec& y , const arma::mat& X , const arma::vec& weight ) ;

bool sum( std::vector<bool>& v ) ;

class leaf_data
{

public:
    std::vector<double> R ;
    std::vector<size_t> vec_months ;
    std::vector<size_t> vec_stocks ;
    std::vector<double> vec_weight ;

    leaf_data( size_t N ) : R( N , 0.0 ) , vec_months( N , 0 ) , vec_stocks( N , 0 ) , vec_weight( N , 0 ) { }

} ;

struct node_info
{

    std::size_t id ; //node id
    std::size_t var ;  //variable
    double cutPoint ;       //cut point // different from BART

    std::vector<double> vec_theta ;

} ;

double log_normal_density( arma::vec& R , arma::mat& cov ) ;

// functions below are for Lasso regression
double soft_c( double a , double lambda ) ;

double lasso_loss( const arma::mat& X , 
                   const arma::mat& Y , 
                   const arma::vec& beta , 
                   double lambda ) ;

arma::vec lasso_fit_standardized( const arma::mat& X , 
                                  const arma::mat& Y , double lambda ,
                                  const arma::vec& beta_ini , 
                                  double eps ) ;

//// indepenent sampler of univariate regression model with conjugate prior
//Rcpp::List runireg_rcpp_loop( arma::vec const& y , arma::mat const& X , arma::vec const& betabar ,
//    arma::mat const& A , double nu , double ssq , size_t R , size_t keep ) ;

void int_to_bin( size_t num , std::vector<size_t>& s ) ;


// 用于打印调试信息的宏定义
#define DEBUG_PRINT(msg) std::cout << "DEBUG: " << __FUNCTION__ << ": " << msg << std::endl;
// 用于打印调试信息的宏定义
#define DEBUG_PRINT_SPACE  std::cout << "  "  << std::endl ;


#include <sstream>
#include <iomanip> // For std::fixed and std::setprecision


void printMat(const arma::mat& m , size_t start_row = 0, size_t start_col = 0, 
                                   size_t num_rows = 0, size_t num_cols = 0 ) ;
 
void printMat(const arma::umat& m , size_t start_row  , size_t start_col , 
                                   size_t num_rows  , size_t num_cols  = 0   ) ;		 
 
 
void  printVec(const arma::vec& v , size_t start = 0, size_t length = 0 )  ;
 


#endif


