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

    virtual double criterion( State& state , 
                              leaf_data& data ) { return 0.0 ; } ;

    virtual void update_leaf_theta( State& state , 
                                    arma::umat& Xorder , 
                                    CTree* leaf_node ) { return ; } ;

    virtual void calculate_criterion( State& state , 
                                      arma::umat& Xorder , 
                                      size_t& split_var , 
                                      size_t& split_point , 
                                      size_t& num_obs_left , 
                                      size_t& num_obs_right , 
                                      CTree* tree_pointer , 
                                      bool& splitable ) { return ; } ;

    virtual double calculate_criterion_one_candidate( State& state , 
                                                      arma::umat& Xorder , 
                                                      size_t var , 
                                                      size_t ind , 
                                                      size_t num_obs ) { return 0.0 ; } ;

} ;


class CAPTreeModel : public CModel
{

public:
    arma::mat regressor ;

    CAPTreeModel( double lambda ) : CModel( 1.0 ) { this->lambda = lambda ; }

    void check_node_splitability( State& state , 
                                  std::vector<CAPTree*>& bottom_nodes_vec , 
                                  std::vector<bool>& node_splitability ) ;

    void calculate_criterion( State& state , 
                              std::vector<CAPTree*>& bottom_nodes_vec , 
                              std::vector<bool>& node_splitability , 
                              size_t& split_node , 
                              size_t& split_var , 
                              size_t& split_point , 
                              bool& splitable , 
                              std::vector<double>& criterion_values ) ;

    void calculate_criterion_APTree_TS( State& state , 
                                        std::vector<CAPTree*>& bottom_nodes_vec , 
                                        std::vector<bool>& node_splitability , 
                                        size_t& split_node , 
                                        size_t& split_var , 
                                        size_t& split_point , 
                                        bool& splitable ) ;

    void split_node( State& state , 
                     CAPTree* node , 
                     size_t split_var , 
                     size_t split_point ) ;

    void split_node_APTree_TS( State& state , 
                               CAPTree* node , 
                               size_t split_var , 
                               size_t split_point ) ;

    void initialize_portfolio( State& state , CAPTree* node ) ;

    void initialize_regressor_matrix( State& state ) ;

    void predict_AP( arma::mat& X , 
                     CAPTree& root , 
                     arma::vec& months , 
                     arma::vec& leaf_index ) ;

    void calculate_criterion_one_variable( State& state , 
                                           size_t var , 
                                           std::vector<CAPTree*>& bottom_nodes_vec , 
                                           size_t node_ind , 
                                           std::vector<double>& output , 
                                           arma::vec& weighted_return_all , 
                                           arma::vec& cumu_weight_all , 
                                           arma::vec& num_stocks_all ) ;

    void calculate_criterion_one_variable_APTree_TS( State& state , 
                                                     size_t var , 
                                                     std::vector<CAPTree*>& bottom_nodes_vec , 
                                                     size_t node_ind , 
                                                     std::vector<double>& output , 
                                                     arma::vec& weighted_return_all , 
                                                     arma::vec& cumu_weight_all , 
                                                     arma::vec& num_stocks_all , 
                                                     size_t var_ind ) ;

    void node_sufficient_stat( State& state , 
                               arma::umat& Xorder , 
                               arma::vec& weighted_return_all , 
                               arma::vec& cumu_weight_all , 
                               arma::vec& num_stocks_all ) ;

    void calculate_factor( CAPTree& root , 
                           arma::vec& leaf_node_index , 
                           arma::mat& all_leaf_portfolio , 
                           arma::mat& leaf_weight , 
                           arma::mat& ft , 
                           State& state ) ;

    double calculate_R2( State& state , 
                         arma::mat& ft ) ;

} ;

#endif

