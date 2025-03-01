#ifndef GUARD_model_h
#define GUARD_model_h
#include "state.h"

class CTree ;
class CAPTree ;

class CModel
{

public:
    double lambda ;

    CModel( double lambda ) { this->lambda = lambda ; }

    virtual double criterion( CState& state , 
                              leaf_data& data ) { return 0.0 ; } ;

    virtual void update_leaf_theta( CState& state , 
                                    arma::umat& matx_Xorder ,
                                    CTree* leaf_node ) { return ; } ;

    virtual void calculate_criterion( CState& state , 
                                      arma::umat& matx_Xorder ,
                                      size_t& split_var , 
                                      size_t& split_point , 
                                      size_t& num_obs_left , 
                                      size_t& num_obs_right , 
                                      CTree* tree_pointer , 
                                      bool& splitable ) { return ; } ;

    virtual double calculate_criterion_one_candidate( CState& state , 
                                                      arma::umat& matx_Xorder ,
                                                      size_t var , 
                                                      size_t ind , 
                                                      size_t num_obs ) { return 0.0 ; } ;

} ;


class CAPTreeModel : public CModel
{

public:
    arma::mat m_matx_regressor ;

    CAPTreeModel( double lambda ) : CModel( 1.0 ) { this->lambda = lambda ; }

    void check_node_splitability( CState& state , 
                                  std::vector<CAPTree*>& vec_bottom_nodes ,
                                  std::vector<bool>& vec_node_splitability ) ;

    void calculate_criterion( CState& state , 
                              std::vector<CAPTree*>& vec_bottom_nodes ,
                              std::vector<bool>& vec_node_splitability ,
                              size_t& split_node , 
                              size_t& split_var , 
                              size_t& split_point , 
                              bool& splitable , 
                              std::vector<double>& vec_criterion_values ) ;

    void calculate_criterion_APTree_TS( CState& state , 
                                        std::vector<CAPTree*>& vec_bottom_nodes ,
                                        std::vector<bool>& vec_node_splitability ,
                                        size_t& split_node , 
                                        size_t& split_var , 
                                        size_t& split_point , 
                                        bool& splitable ) ;

    void split_node( CState& state , 
                     CAPTree* node , 
                     size_t split_var , 
                     size_t split_point ) ;

    void split_node_APTree_TS( CState& state , 
                               CAPTree* node , 
                               size_t split_var , 
                               size_t split_point ) ;

    void initialize_portfolio( CState& state , CAPTree* node ) ;

    void initialize_regressor_matrix( CState& state ) ;

    void predict_AP( arma::mat& matx_X ,
                     CAPTree& root , 
                     arma::vec& vec_months ,
                     arma::vec& vec_leaf_index ) ;

    void calculate_criterion_one_variable( CState& state , 
                                           size_t var , 
                                           std::vector<CAPTree*>& vec_bottom_nodes ,
                                           size_t node_ind , 
                                           std::vector<double>& vec_output ,
                                           arma::vec& vec_weighted_return_all ,
                                           arma::vec& vec_cumu_weight_all ,
                                           arma::vec& vec_num_stocks_all ) ;

    void calculate_criterion_one_variable_APTree_TS( CState& state , 
                                                     size_t var , 
                                                     std::vector<CAPTree*>& vec_bottom_nodes,
                                                     size_t node_ind , 
                                                     std::vector<double>& vec_output ,
                                                     arma::vec& vec_weighted_return_all ,
                                                     arma::vec& vec_cumu_weight_all ,
                                                     arma::vec& vec_num_stocks_all ,
                                                     size_t var_ind ) ;

    void node_sufficient_stat( CState& state , 
                               arma::umat& matx_Xorder ,
                               arma::vec& vec_weighted_return_all , 
                               arma::vec& vec_cumu_weight_all , 
                               arma::vec& vec_num_stocks_all ) ;

    void calculate_factor( CAPTree& root , 
                           arma::vec& vec_leaf_node_index ,
                           arma::mat& matx_all_leaf_portfolio ,
                           arma::mat& matx_leaf_weight ,
                           arma::mat& matx_fator ,
                           CState& state ) ;

    double calculate_R2( CState& state , 
                         arma::mat& matx_fator) ;

} ;

#endif

