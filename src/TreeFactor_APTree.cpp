

//#ifdef _WIN32
//// Windows doesn't have sys/time.h, so we provide a dummy implementation or exclude it
//#include <windows.h>
//#else
//#include <sys/time.h>
//#endif

//#include <RcppArmadillo.h>
//#include <RInside.h>

#include <csignal>

#include "common.h"
#include "state.h"
#include "APTree.h"
#include "model.h"
#include "json_io.h"

#include "TreeFactor_APTree.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

//// [[Rcpp::export]]
//Rcpp::List read_params(const std::string& file) {
//    return Rcpp::readRDS(file);
//}



arma::vec read_vec(const std::string& filename) {
    arma::vec vec;
    vec.load(filename, arma::raw_ascii);
    return vec;
}

arma::mat read_mat(const std::string& filename) {
    arma::mat mat;
    mat.load(filename, arma::raw_ascii);
    return mat;
}

template <typename T>
T read_scalar(const std::string& filename) {
    std::ifstream file(filename);
    T value;
    file >> value;
    return value;
}



//template<typename T>
//T read_scalar(const std::string& filename) {
//    T value;
//    std::ifstream file(filename);
//    if (file.is_open()) {
//        file >> value;
//        file.close();
//    }
//    else {
//        throw std::runtime_error("Unable to open file: " + filename);
//    }
//    return value;
//}

// Specialization for reading boolean values
template<>
bool read_scalar(const std::string& filename) {
    std::string value;
    std::ifstream file(filename);
    if (file.is_open()) 
    {
        file >> value;
        file.close();
        if (value == "TRUE" || value == "true")
        {
            return true;
        }
        else 
        { 
            return  false;  
        }
    }
}

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
           )
{


    DEBUG_PRINT_SPACE;

    DEBUG_PRINT( "" );

    //arma::vec R ;   
    //arma::vec Y ;   
    //arma::mat X ;   
    //arma::mat Z ;   
    //arma::mat H ;   
    //arma::vec portfolio_weight ;   
    //arma::vec loss_weight ;   
    //arma::vec stocks ;   
    //arma::vec months ;   
    //arma::vec unique_months ;   
    //arma::vec first_split_var ;   
    //arma::vec second_split_var ;   
    //size_t num_stocks ;   
    //size_t num_months ;   
    //size_t min_leaf_size = 100 ;   
    //size_t max_depth = 5 ;   
    //size_t num_iter = 30 ;   
    //size_t num_cutpoints = 4 ;   
    //double eta = 1.0 ;   
    //bool equal_weight = false ;   
    //bool no_H = false ;   
    //bool abs_normalize = false ;   
    //bool weighted_loss = false ;   
    //bool stop_no_gain = false ;   
    //double lambda_mean = 0 ;   
    //double lambda_cov = 0  ;


    // 读取参数
    arma::vec vec_R_train = read_vec("params/R_train.txt");
    arma::vec vec_Y_train = read_vec("params/Y_train1.txt");
    arma::mat matx_X_train = read_mat("params/X_train.txt");
    arma::mat matx_Z_train = read_mat("params/Z_train.txt");
    arma::mat matx_H_train = read_mat("params/H_train1.txt");

    arma::vec portfolio_weight = read_vec("params/portfolio_weight_train.txt");
    arma::vec loss_weight = read_vec("params/loss_weight_train.txt");
    arma::vec stocks = read_vec("params/stocks_train.txt");
    arma::vec months = read_vec("params/months_train.txt");
    arma::vec unique_months = read_vec("params/unique_months_train.txt");
    arma::vec first_split_var = read_vec("params/first_split_var.txt");
    arma::vec second_split_var = read_vec("params/second_split_var.txt");

    size_t num_stocks = read_scalar<size_t>("params/num_stocks.txt");
    size_t num_months = read_scalar<size_t>("params/num_months.txt");
    size_t min_leaf_size = read_scalar<size_t>("params/min_leaf_size.txt");
    size_t max_depth = read_scalar<size_t>("params/max_depth.txt");
    size_t num_iter = read_scalar<size_t>("params/num_iter.txt");
    size_t num_cutpoints = read_scalar<size_t>("params/num_cutpoints.txt");
    double eta = read_scalar<double>("params/eta.txt");
    bool equal_weight = read_scalar<bool>("params/equal_weight.txt");
    bool no_H = read_scalar<bool>("params/no_H.txt");
    bool abs_normalize = read_scalar<bool>("params/abs_normalize.txt");
    bool weighted_loss = read_scalar<bool>("params/weighted_loss.txt");
    bool stop_no_gain = read_scalar<bool>("params/stop_no_gain.txt");


    double lambda_mean = read_scalar<double>("params/lambda_mean.txt");
    double lambda_cov = read_scalar<double>("params/lambda_cov.txt");
    lambda_mean = 0;
    lambda_cov = 0;


    no_H = false ;

    DEBUG_PRINT_SPACE;

    // 打印调试信息
    std::cout << "DEBUG: " << __FUNCTION__ << ":\n";

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT( "matx_Z_train : "   );
    printMat( matx_Z_train , 0 , 0 , 20 , matx_Z_train.n_cols ) ;
    DEBUG_PRINT_SPACE;
    DEBUG_PRINT( "matx_H_train : "  );
    printMat( matx_H_train , 0 , 0 , 20 , matx_H_train.n_cols ) ;


    std::cout << "num_stocks: " << num_stocks << "\n";
    std::cout << "num_months: " << num_months << "\n";
    std::cout << "min_leaf_size: " << min_leaf_size << "\n";
    std::cout << "max_depth: " << max_depth << "\n";
    std::cout << "num_iter: " << num_iter << "\n";
    std::cout << "num_cutpoints: " << num_cutpoints << "\n";
    std::cout << "eta: " << eta << "\n";
    std::cout << "equal_weight: " << equal_weight << "\n";
    std::cout << "no_H: " << no_H << "\n";
    std::cout << "abs_normalize: " << abs_normalize << "\n";
    std::cout << "weighted_loss: " << weighted_loss << "\n";
    std::cout << "stop_no_gain: " << stop_no_gain << "\n";
    
    DEBUG_PRINT_SPACE;

    std::cout << "lambda_mean: " << lambda_mean << "\n";
    std::cout << "lambda_cov: " << lambda_cov << "\n";

    DEBUG_PRINT_SPACE;

    //std::raise(SIGTRAP)
    //Rcpp::stop("Breakpoint reached.");
    //Rcpp::browser( );
    //Rcpp::Rcout << "______Debug point reached." << std::endl; 
    //Rcpp::Rcerr << "______Debug point reached______." << std::endl;
    //std::cout<< "______Debug point reached." << std::endl; 
    //Rprintf("______Debug point reached.") ;

    // we assume the number of months is continuous
    std::map<size_t , size_t> months_list ;
    assert( num_months == unique_months.n_elem ) ;

    // a mapping from month to index from zero to num_months - 1
    // it is not necessary to normalize months, adjust from zero in the input
    for( size_t i = 0 ; i < num_months ; i++ )
    {
        // count from zero
        months_list[ unique_months( i ) ] = i ;
    }


    // initialize state class to save data objects
    CState state(matx_X_train , 
                 vec_Y_train , 
                 vec_R_train , 
                 matx_Z_train , 
                 matx_H_train , 
                 portfolio_weight , 
                 loss_weight , 
                 stocks , 
                 months , 
                 first_split_var , 
                 second_split_var , 
                 num_months , 
                 months_list , 
                 num_stocks , 
                 min_leaf_size , 
                 max_depth , 
                 num_cutpoints , 
                 equal_weight , 
                 no_H , 
                 abs_normalize , 
                 weighted_loss , 
                 stop_no_gain , 
                 eta , 
                 lambda_mean , 
                 lambda_cov ) ;


    CAPTreeModel model( lambda_cov ) ;

    // calculate Xorder matrix, each index is row index of the data in the X matrix, but sorted from low to high
    // matx_X_train 排序后的 matx_Xsort ，在原matx_X_train中的序号
    // matx_Xsort(0 ,0 )= -1.3926670684 , matx_Xorder(0 ,0 )= 68165 为 对应 matx_X_train( 68165 , 0 ) = -1.392667068458
    arma::umat matx_Xorder( matx_X_train.n_rows , matx_X_train.n_cols , arma::fill::zeros ) ;

    arma::mat matx_Xsort( matx_X_train.n_rows , matx_X_train.n_cols , arma::fill::zeros ) ;
    matx_Xsort  = arma::sort( matx_X_train  ) ;
    matx_Xsort.save( "./tmpParams/matx_Xsort.txt" , arma::raw_ascii )  ;

     for( size_t i = 0 ; i < matx_X_train.n_cols ; i++ )
    {
        matx_Xorder.col( i ) = arma::sort_index( matx_X_train.col( i ) ) ;

       

    }

 
   printMat(  matx_Xorder ,  0 , 0  , 20 , matx_Xorder.n_cols   );
   
   printMat(  matx_X_train ,  0 , 0  , 20 , matx_X_train.n_cols   );



    // initialize tree class
    CAPTree treeRoot( state.num_months , 
                        1 , 
                        state.num_obs_all , 
                        1 , 
                        0 , 
                        &matx_Xorder ) ;

    treeRoot.set_numData_inNode( matx_X_train.n_rows ) ;

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT( "set_numData_inNode : " << matx_X_train.n_rows );
    DEBUG_PRINT_SPACE;

    // initialize the portfolio at the root node
    model.initialize_portfolio( state , &treeRoot ) ;

    // initialize the proper regressor matrix for the criterion
    // Rt ~ Zt * Ft + Ht
    // create a matrix of Zt * Ft + Ht
    model.initialize_regressor_matrix( state ) ;

    bool break_flag = false ;

    std::vector<double> vec_criterion_values ;

    //// out put of split criterion evaluations
    //// a list (number of iters), each one has length of all possible candidates
    //Rcpp::List all_criterion = Rcpp::List::create( ) ;

    arma::vec temp_vec ;

    int stop = 0 ;
    for( size_t iter = 0 ; iter < num_iter ; iter++ )
    {

        DEBUG_PRINT_SPACE ;
        DEBUG_PRINT_SPACE ;
        DEBUG_PRINT( "===================================================================="  << " iter = "  <<  iter    ) ;
        DEBUG_PRINT( " iter = "  <<  iter    ) ;
        DEBUG_PRINT_SPACE ;  

        // main function that grows the tree
        treeRoot.grow( break_flag , model , state , iter , vec_criterion_values ) ;

        temp_vec.set_size( vec_criterion_values.size( ) ) ;

        for( size_t i = 0 ; i < vec_criterion_values.size( ) ; i++ )
        {
            temp_vec( i ) = vec_criterion_values[ i ] ;
        }

        //// output vector of all split criterion for debugging
        //all_criterion.push_back( temp_vec , to_string( iter ) ) ;

        if( break_flag )
        {
            break ;
        }

        stop++ ;
        if (  stop  == 3 )
        {
            int aaa = 0 ;
        }

    }

    arma::vec vec_leaf_node_index ;
    arma::mat matx_all_leaf_portfolio , matx_leaf_weight , matx_factor ;

    model.calculate_factor( treeRoot , 
                            vec_leaf_node_index , 
                            matx_all_leaf_portfolio , 
                            matx_leaf_weight , 
                            matx_factor , 
                            state ) ;

    cout << "fitted tree " << endl ;
    DEBUG_PRINT_SPACE ;
    DEBUG_PRINT_SPACE ;

    cout.precision( 3 ) ;

    cout << treeRoot << endl ;


    //std::stringstream trees ;
    //Rcpp::StringVector output_tree( 1 ) ;
    //trees.precision( 10 ) ;
    //trees.str( std::string( ) ) ;
    //trees << root ;
    //output_tree( 0 ) = trees.str( ) ;

    //// return pointer to the tree structure, cannot be restored if saving the environment in R
    //// APTree *root_pnt = &root ;
    //// Rcpp::XPtr<APTree> tree_pnt(root_pnt, true) ;
    //Rcpp::StringVector json_output( 1 ) ;
    //json j = tree_to_json( root ) ;
    //json_output[ 0 ] = j.dump( 4 ) ;

    //// calculating the pricing error of the factor, run regression
    //double loss = model.calculate_R2( state , ft ) ;

    //return Rcpp::List::create(  Rcpp::Named( "R" ) = R ,
    //                            Rcpp::Named( "X" ) = X ,
    //                            Rcpp::Named( "Xorder" ) = Xorder ,
    //                            Rcpp::Named( "tree" ) = output_tree ,
    //                            Rcpp::Named( "leaf_weight" ) = leaf_weight ,
    //                            Rcpp::Named( "leaf_id" ) = leaf_node_index ,
    //                            Rcpp::Named( "ft" ) = ft ,
    //                            Rcpp::Named( "portfolio" ) = all_leaf_portfolio ,
    //                            Rcpp::Named( "json" ) = json_output ,
    //                            Rcpp::Named( "R2" ) = loss ,
    //                            Rcpp::Named( "all_criterion" ) = all_criterion ) ;

    return ; 

}