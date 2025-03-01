

#ifndef GUARD_state_h
#define GUARD_state_h

#include "common.h"

class CState
{

public:
    arma::mat* m_matx_X_train ; // pointer to the charateristics matrix
    arma::vec* m_matx_Y_train ;
    arma::vec* m_vec_R_train ; // pointer to the return vector
    arma::mat* R_mat ; // pointer to the return matrix, for TSTree only
    arma::mat* m_matx_Z_train ; // placeholder
    arma::mat* F ; // for Bayes tree
    arma::mat* m_matx_regressor ; // for Bayes tree
    arma::mat* m_matx_H_train ; // placeholder
    arma::vec* m_vec_portfolio_weight ;
    arma::vec* m_vec_loss_weight ;
    arma::vec* m_vec_stocks ; // pointer to the index of stocks, same number of rows as X
    arma::vec* m_vec_months ; // months indicator
    arma::vec* m_vec_first_split_var ;
    arma::vec* m_vec_second_split_var ;
    // the two vectors below are for the APTree model 2, first cut at Macro variable
    arma::vec* m_vec_third_split_var ;
    arma::vec* m_vec_deep_split_var ;
    arma::mat* m_matx_first_split_mat ; // for APTree model 2 only
    arma::mat* m_matx_split_candidate_mat ;
    std::map<size_t , size_t>* m_map_months_list ; // list of UNIQUE months

    size_t num_obs_all ;
    size_t num_stocks ;
    size_t num_months ;
    size_t min_leaf_size ;
    size_t max_depth ;
    size_t num_cutPoints ;
    size_t num_regressors ; // for Bayes tree
    size_t numOfCharitisc ;              // p  //number of charateristics
    std::vector<double> m_vec_split_candidates ;

    bool m_b_equal_weight ;
    bool m_b_no_H ;
    bool m_b_abs_normalize ;
    bool m_b_weighted_loss ;
    bool m_b_stop_no_gain ;
    double m_d_overall_loss ;
    double m_d_sigma ;
    double m_d_tau ;
    double m_d_lambda ;
    double m_d_lambda_mean ;
    double m_d_lambda_cov ;
    double m_d_eta ;
    bool m_b_flag_first_cut ;

    // prior parameters for the Bayes tree
    double a ;
    double b ;
    double xi_normal ;
    double xi_spike ;
    double xi_slab ;
    size_t p_normal_prior ;
    size_t p_spike_slab ;

    // state for APTree model
    CState( arma::mat& X , 
           arma::vec& Y , 
           arma::vec& R , 
           arma::mat& Z , 
           arma::mat& H , 
           arma::vec& portfolio_weight , 
           arma::vec& loss_weight , 
           arma::vec& stocks , 
           arma::vec& months , 
           arma::vec& first_split_var , 
           arma::vec& second_split_var , 
           size_t& num_months , 
           std::map<size_t , size_t>& months_list , 
           size_t& num_stocks , 
           size_t& min_leaf_size , 
           size_t& max_depth , 
           size_t& num_cutpoints , 
           bool& equal_weight , 
           bool& no_H , 
           bool& abs_normalize , 
           bool& weighted_loss , 
           bool& stop_no_gain , 
           double& eta , 
           double& lambda_mean , 
           double& lambda_cov )
    {

        DEBUG_PRINT_SPACE;
        DEBUG_PRINT("");

        this->m_matx_X_train = &X ;
        this->m_matx_Y_train = &Y ;
        this->m_vec_R_train = &R ;
        this->m_matx_Z_train = &Z ;
        this->m_matx_H_train = &H ;
        this->m_vec_portfolio_weight = &portfolio_weight ;
        this->m_vec_loss_weight = &loss_weight ;
        this->m_vec_stocks = &stocks ;
        this->m_vec_months = &months ;
        this->m_map_months_list = &months_list ;
        this->m_vec_first_split_var = &first_split_var ;
        this->m_vec_second_split_var = &second_split_var ;
        this->m_vec_third_split_var = 0 ;
        this->m_vec_deep_split_var = 0 ;
        this->num_months = num_months ;
        this->num_stocks = num_stocks ;
        this->min_leaf_size = min_leaf_size ;
        this->max_depth = max_depth ;
        this->num_cutPoints = num_cutpoints ;
        this->m_vec_split_candidates.resize( num_cutpoints ) ;
        this->numOfCharitisc = X.n_cols ;
        this->num_obs_all = X.n_rows ;
        this->m_b_equal_weight = equal_weight ;
        this->m_b_no_H = no_H ;
        this->m_b_abs_normalize = abs_normalize ;
        this->m_b_weighted_loss = weighted_loss ;
        this->m_b_stop_no_gain = stop_no_gain ;

        this->m_d_overall_loss = std::numeric_limits<double>::max( ) ;
        this->m_d_sigma = 0.0 ;
        this->m_d_tau = 0.0 ;
        this->m_d_lambda = 0.0 ;
        this->m_d_eta = eta ;
        this->m_matx_first_split_mat = 0 ;
        this->num_regressors = 0 ;
        this->m_d_lambda_mean = lambda_mean ;
        this->m_d_lambda_cov = lambda_cov ;

        for( size_t i = 0 ; i < num_cutpoints ; i++ )
        {
            m_vec_split_candidates[ i ] = 2.0 / ( num_cutpoints + 1 ) * ( i + 1 ) - 1 ;
        }

        DEBUG_PRINT_SPACE;
        cout << "The split value candidates are " << m_vec_split_candidates << endl ;
        DEBUG_PRINT_SPACE;

    }

    // state for APTree model2
    CState( arma::mat& X , 
            arma::vec& Y , 
            arma::vec& R , 
            arma::mat& Z , 
            arma::mat& H , 
            arma::vec& portfolio_weight , 
            arma::vec& loss_weight , 
            arma::vec& stocks , 
            arma::vec& months , 
            arma::vec& first_split_var , 
            arma::vec& second_split_var , 
            arma::vec& third_split_var , 
            arma::vec& deep_split_var , 
            size_t& num_months , 
            std::map<size_t , size_t>& months_list , 
            size_t& num_stocks , 
            size_t& min_leaf_size , 
            size_t& max_depth , 
            size_t& num_cutpoints , 
            bool& equal_weight , 
            bool& no_H , 
            bool& abs_normalize , 
            bool& weighted_loss , 
            bool& stop_no_gain , 
            double& lambda , 
            size_t& num_obs_all , 
            arma::mat& first_split_mat )
    {

        this->m_matx_X_train = &X ;
        this->m_matx_Y_train = &Y ;
        this->m_vec_R_train = &R ;
        this->m_matx_Z_train = &Z ;
        this->m_matx_H_train = &H ;
        this->m_vec_portfolio_weight = &portfolio_weight ;
        this->m_vec_loss_weight = &loss_weight ;
        this->m_vec_stocks = &stocks ;
        this->m_vec_months = &months ;
        this->m_map_months_list = &months_list ;
        this->m_vec_first_split_var = &first_split_var ;
        this->m_vec_second_split_var = &second_split_var ;
        this->m_vec_third_split_var = &third_split_var ;
        this->m_vec_deep_split_var = &deep_split_var ;

        this->num_months = num_months ;
        this->num_stocks = num_stocks ;
        this->min_leaf_size = min_leaf_size ;
        this->max_depth = max_depth ;
        this->num_cutPoints = num_cutpoints ;
        this->m_vec_split_candidates.resize( num_cutpoints ) ;
        this->numOfCharitisc = X.n_cols ;
        this->m_b_equal_weight = equal_weight ;
        this->m_b_no_H = no_H ;
        this->m_b_abs_normalize = abs_normalize ;
        this->m_b_weighted_loss = weighted_loss ;
        this->m_b_stop_no_gain = stop_no_gain ;
        this->m_d_overall_loss = std::numeric_limits<double>::max( ) ;
        this->m_d_sigma = 0.0 ;
        this->m_d_tau = 0.0 ;
        this->m_d_lambda = lambda ;
        this->m_d_eta = 0.0 ;
        this->num_obs_all = num_obs_all ;
        this->m_matx_first_split_mat = &first_split_mat ;
        this->num_regressors = 0 ;

        for( size_t i = 0 ; i < num_cutpoints ; i++ )
        {
            m_vec_split_candidates[ i ] = 2.0 / ( num_cutpoints + 1 ) * ( i + 1 ) - 1 ;
        }

        cout << "The split value candidates are " << m_vec_split_candidates << endl ;

    }

} ;

#endif