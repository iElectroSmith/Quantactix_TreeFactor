


#include "model.h"
#include "APTree.h"
#include "common.h"

////////////////////////////
//
//
//      APTree global split criterion
//
//
////////////////////////////



void CAPTreeModel::check_node_splitability( CState& state ,
                                            std::vector<CAPTree*>& vec_bottom_nodes ,
                                            std::vector<bool>& vec_node_splitability )
{
   

    DEBUG_PRINT("");

    //raise(SIGTRAP) ;  // Insert this line to trigger a breakpoint    
    
    CAPTree::APTree_Pt node ;

    DEBUG_PRINT( "vec_bottom_nodes size : " << vec_bottom_nodes.size( ) );

    // check node depth and number of data observations
    for( size_t i = 0 ; i < vec_bottom_nodes.size( ) ; i++ )
    {

        node = vec_bottom_nodes[ i ] ;

        if( node->get_treeDepth( ) >= state.max_depth )  //获取节点的深度（树的层级）
        {
            
            DEBUG_PRINT( "node depth  >=  max_depth: 4 , unSplit "   );

            vec_node_splitability[ i ] = false ;

        }
        else if( node->get_numData_inNode( ) <= state.min_leaf_size ) //获取节点中的数据观测数目
        {

            DEBUG_PRINT( "node numData  <=  min_leaf_size: 10  , unSplit "   );

            vec_node_splitability[ i ] = false ;
        
        }
        else
        {

            DEBUG_PRINT( "set this BTM node Split "   );

            vec_node_splitability[ i ] = true ;
        
        }
    
    }

    return ;
}


//用于计算分裂标准值，并确定最佳分裂节点、变量和分裂点
/*
这个函数通过计算每个节点和每个候选分裂点的分裂标准值，确定最佳分裂点并返回相关信息。
该算法通过遍历所有底部节点和候选分裂点，寻找分裂标准值最小的分裂点，
以此来优化树模型的分裂决策。
*/
void CAPTreeModel::calculate_criterion( CState& state , 
                                       std::vector<CAPTree*>& vec_bottom_nodes , 
                                       std::vector<bool>& vec_node_splitability , 
                                       size_t& split_node , 
                                       size_t& split_var , 
                                       size_t& split_point , 
                                       bool& splitable , 
                                       std::vector<double>& vec_criterion_values )
{

    DEBUG_PRINT("");


    size_t num_btmNodes   = vec_bottom_nodes.size( ) ;
    //每个节点的分裂候选点数量，等于分裂点数量乘以特征数量
    size_t num_candidates = state.num_cutPoints * state.numOfCharitisc ;

    // a vector to save split criterion valuation of all nodes, all candidates
    // initialized at infinity
    // the first num_cutpoints * p is for the first node, etc
    //用于保存所有节点和所有候选点的分裂标准值，初始值为无穷大    
    vec_criterion_values.resize( num_btmNodes * num_candidates ) ;
    std::fill( vec_criterion_values.begin( ) ,  vec_criterion_values.end( ) ,  std::numeric_limits<double>::max( ) ) ;
    // std::vector<double> criterion_values(num_nodes * num_candidates, std::numeric_limits<double>::max()) ;

    // temp_vector stores criterion evaluation of ONE variable
    //临时向量，用于存储单个变量的分裂标准值    
    std::vector<double> temp_vec_criterio( state.num_cutPoints ) ;


    // three major sufficient statistics, calculate one for each month
    // sum of weighted returns, sum(w * R)
    //用于存储节点的加权收益、累计权重和股票数量的统计量，初始值为0
    arma::vec vec_weighted_return_all( state.num_months , arma::fill::zeros ) ;
    // sum of weights, sum(w)
    arma::vec vec_cumu_weight_all( state.num_months , arma::fill::zeros ) ;
    // number of stocks
    arma::vec vec_num_stocks_all( state.num_months , arma::fill::zeros ) ;

    size_t temp_split_var ;

    DEBUG_PRINT("loop over all current leaf nodes");

    // loop over all current leaf nodes
    for( size_t nodeID_index = 0 ; nodeID_index < num_btmNodes ; nodeID_index++ )
    {

        if( !vec_node_splitability[ nodeID_index ] )
        {

            DEBUG_PRINT(" this leaf node cannot split here, do nothing, the split criterion value will remain infinite ");
            // if cannot split here, do nothing, the split criterion value will remain infinite
	    //检查 vec_node_splitability，如果节点不可分裂，则跳过。
        
        }
        else
        {
            DEBUG_PRINT(" this leaf node is splitable, checkout split candidates ,calculate sufficient statistics for a node  ");

            CAPTree* tmppNode =  vec_bottom_nodes[ nodeID_index ] ;
            //printMat(  *(tmppNode->m_pMatx_Xorder)   , 0 , 0 , 20 ,  tmppNode->m_pMatx_Xorder->n_cols ) ;
            //pirntMat(  tmppNode->m   , 0 , 0 , 20 ,  tmppNode->m_pMatx_Xorder->n_cols ) ;

            // this node is splitable, checkout split candidates
            // calculate sufficient statistics for a node
            node_sufficient_stat( state , 
                                  *( vec_bottom_nodes[ nodeID_index ]->m_pMatx_Xorder ) , 
                                  vec_weighted_return_all , 
                                  vec_cumu_weight_all , 
                                  vec_num_stocks_all ) ;


            if( vec_bottom_nodes[ nodeID_index ]->get_treeDepth( ) == 1 )
            {

                DEBUG_PRINT_SPACE
                DEBUG_PRINT_SPACE
                DEBUG_PRINT("depth 1, this is the root ") ;

                //对于根节点，遍历 state.first_split_var 指定的变量，
                // 计算分裂标准值并保存到 vec_criterion_values 中。

                // depth 1, this is the root

                DEBUG_PRINT_SPACE

                DEBUG_PRINT("loop over first_split_var->n_elem  times 1 : "  <<  state.m_vec_first_split_var->n_elem ) ;
                DEBUG_PRINT(" vec_first_split_var : " ) ;
                printVec(  *state.m_vec_first_split_var , 0 , 20 );

                for( size_t var = 0 ; var < state.m_vec_first_split_var->n_elem ; var++ )
                {

                    DEBUG_PRINT("");
 

                    // loop over variables, note the constraint on variables for the root
                    temp_split_var = ( size_t ) ( *state.m_vec_first_split_var )( var ) ;

                   DEBUG_PRINT("======================================= treeDepth == 1 "  << " split_var = "  <<  temp_split_var  );

                    DEBUG_PRINT("curt split_var @1 : "  <<  temp_split_var ) ;

                    this->calculate_criterion_one_variable( state , 
                                                            temp_split_var , 
                                                            vec_bottom_nodes , 
                                                            nodeID_index , 
                                                            temp_vec_criterio , 
                                                            vec_weighted_return_all , 
                                                            vec_cumu_weight_all , 
                                                            vec_num_stocks_all ) ;



                    for( size_t ind = 0 ; ind < state.num_cutPoints ; ind++ )
                    {
                        vec_criterion_values[ num_candidates * nodeID_index + temp_split_var * state.num_cutPoints + ind ] = temp_vec_criterio[ ind ] ;

                    }

                }

            }
            else if( vec_bottom_nodes[ nodeID_index ]->get_treeDepth( ) == 2 )
            {

                 DEBUG_PRINT_SPACE
                 DEBUG_PRINT_SPACE
                 DEBUG_PRINT("depth 2  ") ;

                 DEBUG_PRINT("loop over first_split_var->n_elem  times 2 : "  <<  state.m_vec_first_split_var->n_elem ) ;
                 DEBUG_PRINT(" vec_second_split_var : " ) ;
                 printVec(  *state.m_vec_second_split_var , 0 , 20 );

                //对于深度为 2 的节点，遍历 state.second_split_var 指定的变量，
                // 计算分裂标准值并保存到 vec_criterion_values 中
                // depth 2
                for( size_t var = 0 ; var < state.m_vec_second_split_var->n_elem ; var++ )
                {

                    // loop over variables, note the constraint on variables for depth 2
                    temp_split_var = ( size_t ) ( *state.m_vec_second_split_var )( var ) ;
                    DEBUG_PRINT("curt temp_split_var 2 : "  <<  temp_split_var ) ;

                    DEBUG_PRINT("======================================= treeDepth == 2 "  << " split_var = "  <<  temp_split_var  );

                    this->calculate_criterion_one_variable( state , 
                                                            temp_split_var , 
                                                            vec_bottom_nodes , 
                                                            nodeID_index , 
                                                            temp_vec_criterio , 
                                                            vec_weighted_return_all , 
                                                            vec_cumu_weight_all , 
                                                            vec_num_stocks_all ) ;

                    for( size_t ind = 0 ; ind < state.num_cutPoints ; ind++ )
                    {
                        vec_criterion_values[ num_candidates * nodeID_index + temp_split_var * state.num_cutPoints + ind ] = temp_vec_criterio[ ind ] ;
                    }
                }

            }
            else
            {

                DEBUG_PRINT_SPACE
                DEBUG_PRINT_SPACE
                DEBUG_PRINT("Other depth node  ") ;

                DEBUG_PRINT("loop over state->numOfCharitisc  times   : "  <<  state.numOfCharitisc ) ;

                //对于其他节点，遍历所有变量，计算分裂标准值并保存到 vec_criterion_values 中
                // all other following nodes
                for( size_t var = 0 ; var < state.numOfCharitisc ; var++ )
                {

                    DEBUG_PRINT("======================================= treeDepth >= 3 "  << " split_var = "  <<  var  );

                    // loop over variables, there is no constraint, loop over all variables
                    this->calculate_criterion_one_variable( state , 
                                                            var , 
                                                            vec_bottom_nodes , 
                                                            nodeID_index , 
                                                            temp_vec_criterio , 
                                                            vec_weighted_return_all , 
                                                            vec_cumu_weight_all , 
                                                            vec_num_stocks_all ) ;

                    for( size_t ind = 0 ; ind < state.num_cutPoints ; ind++ )
                    {
                        vec_criterion_values[ num_candidates * nodeID_index + var * state.num_cutPoints + ind ] = temp_vec_criterio[ ind ] ;
                    }

                }

            }

        }

    }


    DEBUG_PRINT_SPACE
    DEBUG_PRINT("loop over vec_criterion_values , to find the lowest split criterion and its index ") ;


    //遍历 vec_criterion_values，寻找最小的分裂标准值及其索引。
    // find the lowest split criterion
    size_t lowest_index = 0 ;
    double temp_criterion_value = vec_criterion_values[ 0 ] ;
    //遍历 criterion_values，找到最小的分裂标准值 temp 及其对应的索引 lowest_index
    for( size_t i = 1 ; i < vec_criterion_values.size( ) ; i++ )
    {
        if( vec_criterion_values[ i ] <= temp_criterion_value )
        {
            temp_criterion_value = vec_criterion_values[ i ] ;
            lowest_index = i ;
        }
    }


    DEBUG_PRINT_SPACE
    //如果所有分裂点的标准值都是无穷大，设置 splitable 为 false 并返回。
    if( temp_criterion_value == std::numeric_limits<double>::max( ) )
    {

        DEBUG_PRINT("if all cutpoints have loss infinite, stop split") ;

        //如果所有分裂点的分裂标准值都为正无穷大，则无法分裂，设置 splitable 为 false。
        // if all cutpoints have loss infinite, stop split
        splitable = false ;

        return ;
    }


    DEBUG_PRINT_SPACE

    DEBUG_PRINT("restore corresponding index of node, cutpoint variable and data index") ;


    //否则，更新 state.overall_loss。
    state.m_d_overall_loss = temp_criterion_value ;
    DEBUG_PRINT( "state.m_d_overall_loss : " << temp_criterion_value ) ;


    //lowest_index 并不是用来对当前所有的叶节点进行分裂的，而是用于确定在当前所有的叶节点中，
    //哪个叶节点、在哪个分裂变量、以及在哪个分裂点上进行最优的一次分裂。

    //这段代码的核心思想是将一维索引 lowest_index 分解为三部分：
    //叶节点索引 split_node、分裂变量索引 split_var 和分裂点索引 split_point。
    //通过这种方式，算法可以从一维数组中精确定位出哪个叶节点的哪个分裂变量以及在哪个分裂点上最优化全局标准。
    //假设：
    //num_btmNodes 表示当前树中的叶节点数量。
    //numOfCharitisc 表示可用的分裂变量数量。
    //num_cutPoints 表示每个变量的分裂候选数量。

    //对于每个叶节点 i，分裂变量 var，以及分裂点 cutpoint，
    // 它在 vec_criterion_values 中的索引计算公式是：
    //index = 𝑖 × ( numOfCharitisc × num_cutPoints ) + 𝑣𝑎𝑟 × num_cutPoints + cutpoint

    // restore corresponding index of node, cutpoint variable and data index
    //恢复分裂点信息：
    size_t temp2 ;

    //lowest_index 是 criterion_values 中最小值的索引，num_candidates 
    //是每个节点的候选分裂点数量，因此通过整数除法可以得到分裂点所属的节点索引。
    split_node  = lowest_index / num_candidates ;
    DEBUG_PRINT( "split_node : " << split_node ) ;

    //通过取模运算（%）得到 lowest_index 在当前节点内的具体索引位置。
    //temp2 现在表示在该节点中，分裂变量和分裂点的综合索引
    temp2       = lowest_index % num_candidates ;
    
    //将 temp2 除以 state.num_cutPoints 得到具体的分裂变量索引。
    //每个变量有 state.num_cutPoints 个候选分裂点，
    //因此通过整数除法可以得到分裂变量的索引。
    split_var   = temp2 / state.num_cutPoints ;
    DEBUG_PRINT( "split_var : " << split_var ) ;

    //通过取模运算得到具体的分裂点索引。
    //split_point 表示在分裂变量 split_var 下，具体的分裂点位置。
    split_point = temp2 % state.num_cutPoints ;
    DEBUG_PRINT( "split_point : " << split_point ) ;


    return ;

}

//计算每个节点的分裂标准，并确定最佳分裂点
/*
这段代码的作用是遍历所有叶子节点并计算其可能的分裂点，确定最佳分裂点及其对应的
分裂变量和分裂点。它通过计算每个节点的足够统计量，并依次计算每个分裂变量的
分裂标准值，最终选出最优的分裂方案。
*/
void CAPTreeModel::calculate_criterion_APTree_TS( CState& state , 
                                                  std::vector<CAPTree*>& vec_bottom_nodes , 
                                                  std::vector<bool>& vec_node_splitability , 
                                                  size_t& split_node , 
                                                  size_t& split_var , 
                                                  size_t& split_point , 
                                                  bool& splitable )
{

    DEBUG_PRINT("");

    size_t num_nodes = vec_bottom_nodes.size( ) ;
    //每个节点的分裂候选点数量，等于分裂点数量乘以特征数量
    size_t num_candidates = state.num_cutPoints * state.numOfCharitisc ;

    // a vector to save split criterion valuation of all nodes, all candidates
    // initialized at infinity
    // the first num_cutpoints * p is for the first node, etc
    //用于保存所有节点和所有候选点的分裂标准值，初始值为无穷大
    std::vector<double> criterion_values( num_nodes * num_candidates , std::numeric_limits<double>::max( ) ) ; 
    //临时向量，用于存储单个变量的分裂标准值
    std::vector<double> temp_vector( state.num_cutPoints ) ;

    //用于存储节点的加权收益、累计权重和股票数量的统计量，初始值为0
    arma::vec weighted_return_all( state.num_months , arma::fill::zeros ) ;
    arma::vec cumu_weight_all( state.num_months , arma::fill::zeros ) ;
    arma::vec num_stocks_all( state.num_months , arma::fill::zeros ) ;

    size_t temp_index ;

    for( size_t i = 0 ; i < num_nodes ; i++ )
    {
        if( !vec_node_splitability[ i ] )
        {
            // if cannot split here, do nothing, the split criterion value will remained infinite
            //检查 vec_node_splitability，如果节点不可分裂，则跳过。
        }
        else
        {
            // this node is splitable, checkout split candidates
            // calculate sufficient statistics for a node

            node_sufficient_stat( state , 
                                  *( vec_bottom_nodes[ i ]->m_pMatx_Xorder ) , 
                                  weighted_return_all , 
                                  cumu_weight_all , 
                                  num_stocks_all ) ;

            // depth 1, root
            for( size_t var = 0 ; var < state.m_matx_first_split_mat->n_cols ; var++ )
            {

                // loop over variables
                temp_index = ( size_t ) ( *state.m_vec_first_split_var )( var ) ;
                this->calculate_criterion_one_variable_APTree_TS( state , 
                                                                    temp_index , 
                                                                    vec_bottom_nodes , 
                                                                    i , 
                                                                    temp_vector , 
                                                                    weighted_return_all , 
                                                                    cumu_weight_all , 
                                                                    num_stocks_all , 
                                                                    var ) ;
                for( size_t ind = 0 ; ind < state.num_cutPoints ; ind++ )
                {
                    criterion_values[ num_candidates * i + temp_index * state.num_cutPoints + ind ] = temp_vector[ ind ] ;
                }
            }
        }
    }

    // find the lowest split criterion
    size_t lowest_index = 0 ;
    double temp = criterion_values[ 0 ] ;
    //遍历 criterion_values，找到最小的分裂标准值 temp 及其对应的索引 lowest_index
    for( size_t i = 1 ; i < criterion_values.size( ) ; i++ )
    {
        if( criterion_values[ i ] <= temp )
        {
            temp = criterion_values[ i ] ;
            lowest_index = i ;
        }
    }

    //如果所有分裂点的标准值都是无穷大，设置 splitable 为 false 并返回。
    if( temp == std::numeric_limits<double>::max( ) )
    {
        // if all cutpoints have loss infinite, stop split
        splitable = false ;
        return ;
    }

    //否则，更新 state.overall_loss。
    state.m_d_overall_loss = temp ;

    // restore corresponding index of node, cutpoint variable and data index
    //恢复分裂点信息：
    size_t temp2 ;



    //这段代码的核心思想是将一维索引 lowest_index 分解为三部分：
    //叶节点索引 split_node、分裂变量索引 split_var 和分裂点索引 split_point。
    //通过这种方式，算法可以从一维数组中精确定位出哪个叶节点的哪个分裂变量以及在哪个分裂点上最优化全局标准。
    //假设：
    //num_btmNodes 表示当前树中的叶节点数量。
    //numOfCharitisc 表示可用的分裂变量数量。
    //num_cutPoints 表示每个变量的分裂候选数量。

    //对于每个叶节点 i，分裂变量 var，以及分裂点 cutpoint，它在 vec_criterion_values 中的索引计算公式是：
    //index = 𝑖 × ( numOfCharitisc × num_cutPoints ) + 𝑣𝑎𝑟 × num_cutPoints + cutpoint

    //num_candidates 是每个叶节点的分裂候选总数，即 numOfCharitisc* num_cutPoints。
    //这行代码通过除法操作，计算出 lowest_index 所对应的叶节点索引 split_node。
    split_node  = lowest_index / num_candidates ;

    //取模操作（% ）用来确定 lowest_index 在当前叶节点中的相对位置，
    //即它在当前叶节点的所有候选分裂点中的索引。
    temp2       = lowest_index % num_candidates ;

    //将 temp2 除以 state.num_cutPoints 得到具体的分裂变量索引。
    //每个变量有 state.num_cutPoints 个候选分裂点，
    //因此通过整数除法可以得到分裂变量的索引。
    split_var   = temp2 / state.num_cutPoints ;

    //通过取模运算得到具体的分裂点索引。
    //split_point 表示在分裂变量 split_var 下，具体的分裂点位置。
    split_point = temp2 % state.num_cutPoints ;

    return ;
}



//计算给定节点的三个足够统计量：加权收益的累计值、累计权重的累计值和每个月份的股票数量
void CAPTreeModel::node_sufficient_stat( CState& state , 
                                        arma::umat& matx_Xorder , 
                                        arma::vec& vec_weighted_return_all , 
                                        arma::vec& vec_cumu_weight_all , 
                                        arma::vec& vec_num_stocks_all )
{


    DEBUG_PRINT("");

    // This function create basis portfolio for the node
    // Use R not Y
    size_t num_obs = matx_Xorder.n_rows ; //获取数据排序顺序 Xorder 的行数，表示数据观测数目。
    size_t temp_index ;
    size_t temp_month ;
    size_t temp_month_index ;

    // three sufficient statistics
    // weighted return, w * Rt
    vec_weighted_return_all.fill( 0.0 ) ;  //加权收益的累计值
    // cumulative weight, sum of w
    vec_cumu_weight_all.fill( 0.0 ) ;  //累计权重的累计值
    // number of stocks
    vec_num_stocks_all.fill( 0.0 ) ;  //每个月份的股票数量


    DEBUG_PRINT(" cal all the DATA , loop times : " << num_obs );

    //循环遍历 Xorder 的每一行（数据观测）
    for( size_t i = 0 ; i < num_obs ; i++ )
    {

        temp_index = matx_Xorder( i , 0 ) ;
        temp_month = ( *state.m_vec_months )( temp_index ) ;
        temp_month_index = state.m_map_months_list->at( temp_month ) ;

        //std::cout << " temp_index: " << temp_index  << std::endl ;
        //std::cout << " temp_month: " << temp_month   << std::endl ;
        //std::cout << " temp_month_index: " << temp_month_index  << std::endl ;

        //累加加权收益，乘以对应的权重
        vec_weighted_return_all( temp_month_index ) += ( *state.m_vec_R_train )( temp_index ) * ( *state.m_vec_portfolio_weight )( temp_index ) ;
        
        //累加权重。
        vec_cumu_weight_all( temp_month_index )     += ( *state.m_vec_portfolio_weight )( temp_index ) ;
        
        //增加月份对应的股票数量
        vec_num_stocks_all( temp_month_index )      += 1.0 ;


    }

    DEBUG_PRINT("vec_weighted_return_all");
    //printVec( vec_weighted_return_all , 0 , 20  ) ; 
    DEBUG_PRINT("vec_cumu_weight_all");
    //printVec( vec_cumu_weight_all , 0 , 20  ) ; 
    DEBUG_PRINT("vec_num_stocks_all");
    //printVec( vec_num_stocks_all , 0 , 20  ) ; 


    /**********************************************************************
    DEBUG_PRINT("  print vec_weighted_return_all : "  );
    for( size_t i = 0 ; i < vec_weighted_return_all.size( ) ; i++ )
    {

        std::cout  << vec_weighted_return_all[i]   << std::endl ;

    }

    DEBUG_PRINT("  print vec_cumu_weight_all : "  );
    for( size_t i = 0 ; i < vec_cumu_weight_all.size( ) ; i++ )
    {

        std::cout << vec_cumu_weight_all[i]   << std::endl ;

    }

    DEBUG_PRINT("  print vec_num_stocks_all : "  );
    for( size_t i = 0 ; i < vec_num_stocks_all.size( ) ; i++ )
    {

        std::cout << vec_num_stocks_all[i]   << std::endl ;

    }
    /**********************************************************************/


    return ;

}


//calculate_criterion 函数会遍历当前所有的叶节点（即 vec_bottom_nodes 中的节点），
// 并对每个叶节点的所有可能分裂候选进行评估。
//对于每个叶节点，算法会计算在不同变量、不同分裂点上的分裂标准（例如损失函数的减少）。
//所有这些分裂标准会被存储在一个一维数组 vec_criterion_values 中。

//通过比较 vec_criterion_values 中的值，算法会找到全局最优的分裂点，
//也就是 lowest_index 所指向的分裂点。
//lowest_index 对应的是某一个特定的叶节点、变量和分裂点的组合，
//它能在全局范围内最大化预定义的性能指标（如最小化损失函数）。

//找到 lowest_index 后，算法不会同时分裂所有的叶节点，而是仅分裂那个由 lowest_index 指定的叶节点。
//这个叶节点会被分成两个新的子节点，并成为新的叶节点。
//其他叶节点保持不变，等待在后续步骤中可能被分裂。

//生成的新叶节点会被添加到叶节点列表中(vec_bottom_nodes)。
//在下一次迭代中，算法将重新评估所有当前的叶节点，并重复上述步骤。

//用于计算特定变量在特定节点上的分裂标准值。 
void CAPTreeModel::calculate_criterion_one_variable( CState& state , 
                                                     size_t var , 
                                                     std::vector<CAPTree*>& vec_bottom_nodes , 
                                                     size_t node_ind , 
                                                     std::vector<double>& vec_output_criterion , 
                                                     arma::vec& vec_weighted_return_all , 
                                                     arma::vec& vec_cumu_weight_all , 
                                                     arma::vec& vec_num_stocks_all )
{

    DEBUG_PRINT("");



    // calculate split criterion for one variable at a specific node
    //获取当前处理的节点。
    CAPTree* pCurtNode = vec_bottom_nodes[ node_ind ] ;

    //pirntMat(  *(pCurtNode->m_pMatx_Xorder)   , 0 , 0 , 20 ,  pCurtNode->m_pMatx_Xorder->n_cols ) ;

    // initialize split criterion, start from infinity
    //初始化 vec_output，用于存储分裂标准值，设置为无穷大
    std::fill( vec_output_criterion.begin( ) , vec_output_criterion.end( ) , std::numeric_limits<double>::max( ) ) ;

    // essentially, the sufficient statistics are two vectors with length num_months ;
    // first vector: weight * return
    // second vector: cumulative weight
    // the portfolio is just elementwise ratio of the two vectors
    size_t num_nodes = vec_bottom_nodes.size( ) ;
    // 获取节点的排序矩阵
    arma::umat* Xorder = pCurtNode->m_pMatx_Xorder ;

    //printMat(  *(Xorder)   , 0 , 0 , 20 ,  Xorder->n_cols ) ;

    // calculate sufficient statistics of all data here
    //temp_index，temp_month，temp_month_index 用于存储临时索引和月份信息。
    size_t temp_index = 0 ;
    size_t temp_month = 0 ;
    size_t temp_month_index = 0 ;

    // next loop over cutpoints, calculate sufficient statistics on left / right side
    // basis porfolio, use R not Y
    // 初始化左右节点的统计量 
    // 初始化左右节点的加权回报、累计权重和股票数量
    arma::vec vec_weighted_return_left( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_cumu_weight_left( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_num_stocks_left( state.num_months , arma::fill::zeros ) ;

    arma::vec vec_weighted_return_right( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_cumu_weight_right( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_num_stocks_right( state.num_months , arma::fill::zeros ) ;

    double curtCutPoint ;
    size_t loop_index = 0 ;
    arma::mat matx_mu ;
    arma::mat matx_sigma ;
    arma::mat matx_weight ;
    arma::mat matx_FactorReturns ;
    double weight_sum ;

    //初始化所有投资组合矩阵
    //初始化 all_portfolio 矩阵，用于存储所有投资组合的回报。
    arma::mat matx_all_portfolio( state.num_months , num_nodes + 1 , arma::fill::zeros ) ;
    temp_index = 2 ; // the FIRST two columns for the candidate split // 前两列用于候选分裂

    // 填充非当前节点的投资组合回报
    //将非当前节点的投资组合回报填充到 all_portfolio 矩阵中。
    for( size_t i = 0 ; i < num_nodes ; i++ )
    {
        if( i != node_ind )
        {
            // if it is not the current node
            // copy portfolio return from the leaf directly
            for( size_t ind = 0 ; ind < state.num_months ; ind++ )
            {
                matx_all_portfolio( ind , temp_index ) = ( vec_bottom_nodes[ i ]->m_vec_month_theta )[ ind ] ;
            }
            temp_index++ ;
        }
    }


    DEBUG_PRINT("matx_all_portfolio");
    //printMat( matx_all_portfolio , 0 , 0 , 20 , matx_all_portfolio.n_cols);


    // next calculate portfolio returns for current candidate
    // 遍历所有候选分裂点，计算当前节点的投资组合回报
    for( size_t indexCutPoint = 0 ; indexCutPoint < state.num_cutPoints ; indexCutPoint++ )
    {

        // reset all vectors for a new cutpoint
        vec_weighted_return_left.fill( 0.0 ) ;
        vec_weighted_return_right.fill( 0.0 ) ;
        vec_cumu_weight_left.fill( 0.0 ) ;
        vec_cumu_weight_right.fill( 0.0 ) ;
        vec_num_stocks_left.fill( 0.0 ) ;
        vec_num_stocks_right.fill( 0.0 ) ;

        // cout << "dim of all portfolio " << all_portfolio.n_rows << " " << all_portfolio.n_cols << endl ;


        // loop over candidates
        curtCutPoint = state.m_vec_split_candidates[ indexCutPoint ] ;
        DEBUG_PRINT(" m_vec_split_candidates  : "   ) ;
        printVec( state.m_vec_split_candidates ,  0 , 20 ) ;

        DEBUG_PRINT_SPACE
        DEBUG_PRINT("cutPoints_1 : " << curtCutPoint ) ;

        // while ((*state.X)((*Xorder)(loop_index, var), var) <= cutpoint)
        // {
        //     // the observation is on the left side
        //     temp_index = (*Xorder)(loop_index, var) ;              // convert from sorted index (rank) to the original index
        //     temp_month = (*state.months)(temp_index) ;             // find corresponding month
        //     temp_month_index = state.months_list->at(temp_month) ; // index of the month in the month_list
        //     // update weighted return, cumulative weight and count of stocks
        //     weighted_return_left(temp_month_index) += (*state.R)(temp_index) * (*state.weight)(temp_index) ;
        //     cumu_weight_left(temp_month_index) += (*state.weight)(temp_index) ;
        //     num_stocks_left(temp_month_index) += 1.0 ;
        //     loop_index++ ; // index of the current obs in the original Xorder matrix, will be used in the next round until it reaches total number of obs
        //     if (loop_index == (*Xorder).n_rows)
        //     {
        //         // terminating condition, avoid overflow
        //         break ;
        //     }
        // }

        // weighted_return_right = weighted_return_all - weighted_return_left ;
        // cumu_weight_right = cumu_weight_all - weighted_return_left ;
        // num_stocks_right = num_stocks_all - num_stocks_left ;

        DEBUG_PRINT(" This time , traverse the column   : " << var << " of matx_X_train ");

        DEBUG_PRINT( " loop ( *Xorder ).n_rows times : " << ( *Xorder ).n_rows ) ; 

        DEBUG_PRINT_SPACE ;
        DEBUG_PRINT("========================" << " column =  "  <<  var <<  "   cutPoints = " << curtCutPoint ); 
        DEBUG_PRINT_SPACE ;  

        // 根据候选分裂点分配数据到左右节点
        for( size_t jj = 0 ; jj < ( *Xorder ).n_rows ; jj ++ )
        {
                       
            //double  tmp = ( ( *Xorder )( jj , var ) ) ;             DEBUG_PRINT(" 1 tmp : " << tmp );
            //     tmp = ( *state.m_matx_X_train )( ( *Xorder )( jj , var ) , var ) ;       DEBUG_PRINT(" 2 tmp : " << tmp );


            if( ( *state.m_matx_X_train )( ( *Xorder )( jj , var ) , var ) <= curtCutPoint )
            {

                temp_index = ( *Xorder )( jj , var ) ;                     // convert from sorted index (rank) to the original index
                temp_month = ( *state.m_vec_months )( temp_index ) ;             // find corresponding month
                temp_month_index = state.m_map_months_list->at( temp_month ) ;   // index of the month in the month_list

                // update weighted return, cumulative weight and count of stocks
                //P11 (4)式
                vec_weighted_return_left( temp_month_index ) += ( *state.m_vec_R_train )( temp_index ) * ( *state.m_vec_portfolio_weight )( temp_index ) ;
                vec_cumu_weight_left( temp_month_index )     += ( *state.m_vec_portfolio_weight )( temp_index ) ;
                vec_num_stocks_left( temp_month_index )      += 1.0 ;
            }
            else
            {

                temp_index = ( *Xorder )( jj , var ) ;                    // convert from sorted index (rank) to the original index
                temp_month = ( *state.m_vec_months )( temp_index ) ;            // find corresponding month
                temp_month_index = state.m_map_months_list->at( temp_month ) ;  // index of the month in the month_list
                
                // update weighted return, cumulative weight and count of stocks
                //P11 (4)式
                vec_weighted_return_right( temp_month_index ) += ( *state.m_vec_R_train )( temp_index ) * ( *state.m_vec_portfolio_weight )( temp_index ) ;
                vec_cumu_weight_right( temp_month_index )     += ( *state.m_vec_portfolio_weight )( temp_index ) ;
                vec_num_stocks_right( temp_month_index )      += 1.0 ;

            }

        }

        DEBUG_PRINT_SPACE
        DEBUG_PRINT("vec_weighted_return_left");
        //printVec( vec_weighted_return_left , 0 , 20  ) ; 
        DEBUG_PRINT("vec_cumu_weight_left");
        //printVec( vec_cumu_weight_left , 0 , 20  ) ; 
        DEBUG_PRINT("vec_num_stocks_left");
        //printVec( vec_num_stocks_left , 0 , 20  ) ; 

        DEBUG_PRINT_SPACE
        DEBUG_PRINT("vec_weighted_return_right");
        //printVec( vec_weighted_return_right , 0 , 20  ) ; 
        DEBUG_PRINT("vec_cumu_weight_right");
        //printVec( vec_cumu_weight_right , 0 , 20  ) ; 
        DEBUG_PRINT("vec_num_stocks_right");
        //printVec( vec_num_stocks_right , 0 , 20  ) ; 

        DEBUG_PRINT_SPACE
        //printMat(*(Xorder) , 0 , 0 , 20 , Xorder->n_cols);

        // cout << " ---- " << endl ;
        // cout << arma::join_rows(num_stocks_all - num_stocks_left - num_stocks_right, weighted_return_all - weighted_return_left - weighted_return_right, cumu_weight_all - cumu_weight_left - cumu_weight_right) << endl ;

        // check stopping conditions such as minimal leaf size, number of stocks
        // 检查停止条件，例如最小叶子大小、股票数量等
        if(     vec_num_stocks_right.min( ) < state.min_leaf_size 
             || vec_num_stocks_left.min( )  < state.min_leaf_size 
             || arma::accu( vec_num_stocks_right ) == 0 
             || arma::accu( vec_num_stocks_left ) == 0 )
        {

            // too few data in the leaf, set criterion as infinity
            vec_output_criterion[ indexCutPoint ] = std::numeric_limits<double>::max( ) ;

        }
        else
        {

            DEBUG_PRINT(" candidate is splitable, calculate split criterion : "  );

            // if this candidate is splitable, calculate split criterion
             // 计算左右子叶的投资组合回报
            for( size_t ind = 0 ; ind < state.num_months ; ind++ )
            {

                // calculate weighted return for the candidate left / right child leaves
                // first column for the left portfolio
                matx_all_portfolio( ind , 0 ) = ( vec_num_stocks_left( ind )  == 0 ) ? 0 : vec_weighted_return_left( ind )  / vec_cumu_weight_left( ind  ) ;
		        // second column for the right portfolio
                matx_all_portfolio( ind , 1 ) = ( vec_num_stocks_right( ind ) == 0 ) ? 0 : vec_weighted_return_right( ind ) / vec_cumu_weight_right( ind ) ;
            
            }

            DEBUG_PRINT("matx_all_portfolio");
            //printMat( matx_all_portfolio , 0 , 0 , 20 , matx_all_portfolio.n_cols);

           /*******************************************************************************
            DEBUG_PRINT_SPACE;
            DEBUG_PRINT_SPACE;
            DEBUG_PRINT("  print vec_num_stocks_left : "  );
            for( size_t i = 0 ; i < vec_num_stocks_left.size( ) ; i++ )
            {
                std::cout  << vec_num_stocks_left[i]   << std::endl ;
            }

            DEBUG_PRINT("  print vec_weighted_return_left : "  );
            for( size_t i = 0 ; i < vec_weighted_return_all.size( ) ; i++ )
            {
                std::cout  << vec_weighted_return_all[i]   << std::endl ;
            }

            DEBUG_PRINT("  print vec_cumu_weight_left : "  );
            for( size_t i = 0 ; i < vec_cumu_weight_left.size( ) ; i++ )
            {
                std::cout << vec_cumu_weight_left[i]   << std::endl ;
            }

            DEBUG_PRINT("  print matx_all_portfolio : "  );
            for( size_t i = 0 ; i < matx_all_portfolio.size( ) ; i++ )
            {
                std::cout << matx_all_portfolio[i]   << std::endl ;
            }


            /*******************************************************************************
            DEBUG_PRINT_SPACE;
            DEBUG_PRINT_SPACE;

            DEBUG_PRINT("  print vec_num_stocks_right : "  );
            for( size_t i = 0 ; i < vec_num_stocks_right.size( ) ; i++ )
            {
                std::cout  << vec_num_stocks_right[i]   << std::endl ;
            }

            DEBUG_PRINT("  print vec_weighted_return_right : "  );
            for( size_t i = 0 ; i < vec_weighted_return_right.size( ) ; i++ )
            {
                std::cout  << vec_weighted_return_right[i]   << std::endl ;
            }

            DEBUG_PRINT("  print vec_cumu_weight_right : "  );
            for( size_t i = 0 ; i < vec_cumu_weight_right.size( ) ; i++ )
            {
                std::cout << vec_cumu_weight_right[i]   << std::endl ;
            }

            DEBUG_PRINT("  print matx_all_portfolio : "  );
            for( size_t i = 0 ; i < matx_all_portfolio.size( ) ; i++ )
            {
                std::cout << matx_all_portfolio[i]   << std::endl ;
            }
            DEBUG_PRINT_SPACE;
            DEBUG_PRINT_SPACE;
           /*******************************************************************************/

            // 计算均值、协方差、权重和投资组合回报  
            // matx_all_portfolio: 80R*2C
            // matx_mu: 1R*2C  ->  2R*1C
            // matx_sigma : 2R*2C
            //P11 (4)式
            matx_mu    = arma::mean( matx_all_portfolio , 0 ) ; // 计算所有投资组合的均值（列均值） // 0 for column mean
            matx_mu    = arma::trans( matx_mu ) ;    // 转置为列向量 // transpose to column vectors

            size_t n_leafs = matx_mu.n_elem ;

            // 计算所有投资组合的协方差矩阵 协方差矩阵描述了不同投资组合回报之间的关系。
            matx_sigma = arma::cov( matx_all_portfolio ) ;
            DEBUG_PRINT("matx_sigma");
            //printMat( matx_sigma , 0 , 0 , 20 , matx_sigma.n_cols);


            // mean variance efficient weight
            // 计算均值方差有效权重
            // 使用协方差矩阵 sigma 和均值 mu 计算均值方差有效权重（最优权重）
			// eye 是 "identity"（单位）的缩写，eye 函数用于生成单位矩阵。
			// arma::ones 是一个用于生成全为 1 的矩阵或向量的函数
			// arma::inv 是一个用于计算矩阵逆矩阵
            // P11底下： 收缩参数
            std::cout << "lambda_mean: " << state.m_d_lambda_mean << "\n";
            std::cout << "lambda_cov: " << state.m_d_lambda_cov << "\n";
            matx_weight =   arma::inv(  matx_sigma + state.m_d_lambda_cov  * arma::eye(  n_leafs   , n_leafs   ) ) 
                            *         ( matx_mu    + state.m_d_lambda_mean * arma::ones( matx_mu.n_rows , matx_mu.n_cols ) ) ;

            DEBUG_PRINT("matx_weight 1 ");
            //printMat( matx_weight , 0 , 0 , 20 , matx_weight.n_cols);

            // 创建等权重向量 所有权重相等
            arma::vec vec_equal_weight( n_leafs ) ;
            vec_equal_weight.fill( 1.0 / n_leafs ) ;

            DEBUG_PRINT("vec_equal_weight");
            //printMat( vec_equal_weight , 0 ,  20 );

            // 混合权重，结合均值方差有效权重和等权重
            // 将均值方差有效权重和等权重按照比例 state.eta 进行混合。
            matx_weight = matx_weight * state.m_d_eta + ( 1.0 - state.m_d_eta ) * vec_equal_weight ;

            DEBUG_PRINT("matx_weight 2 ");
            //printMat( matx_weight , 0 , 0 , 20 , matx_weight.n_cols);

            // 如果需要绝对值归一化
            if( state.m_b_abs_normalize )
            {
				//arma::abs(matx_weight)：计算 matx_weight 向量中每个元素的绝对值。
				//arma::accu(...)：计算绝对值向量的所有元素之和
                weight_sum = arma::accu( arma::abs( matx_weight ) ) ;
            }
            else
            {
                weight_sum = arma::accu( ( matx_weight ) ) ;
            }



            // 归一化权重 
			// 通过这种方式，确保了 matx_weight 向量的所有元素之和为1，即：
            matx_weight = matx_weight / weight_sum ;

            // mean variance efficient portfolio 
            // 计算均值方差有效投资组合回报
            matx_FactorReturns = matx_all_portfolio * matx_weight ;   //80R*2C *  2R*1 -> 80R*1C

            DEBUG_PRINT("matx_FactorReturns");
            //printMat( matx_FactorReturns , 0 , 0 , 20 , matx_FactorReturns.n_cols);

            DEBUG_PRINT_SPACE
            DEBUG_PRINT("calc matx_regressor = Z_{it} * ft ， traverse all of columns of  matx_Z_train ");

            //遍历所有观察值（state.num_obs_all）和所有特征列（(*state.Z).n_cols）。
            //计算交互项 Z_{it} * ft
            // P11 式(5) 
            for( size_t i = 0 ; i < state.num_obs_all ; i++ )
            {
                // 获取当前月份的索引
                for( size_t j = 0 ; j < ( *state.m_matx_Z_train ).n_cols ; j++ )
                {

                    // interaction term, Z_{it} * ft
                    //Z_{ it } 是特征矩阵 Z 的第 i 行第 j 列元素，ft 是投资组合回报。两者相乘，得到交互项。
                    temp_month_index = state.m_map_months_list->at( ( *state.m_vec_months )( i ) ) ;

                    this->m_matx_regressor( i , j ) = ( *state.m_matx_Z_train )( i , j ) * matx_FactorReturns( temp_month_index , 0 ) ;
                }

            }


            DEBUG_PRINT("m_matx_regressor");
            //printMat( m_matx_regressor , 0 , 0 , 20 , m_matx_regressor.n_cols);

            DEBUG_PRINT("*state.m_matx_Y_train( in fact the input xret ) ");
            //printMat( *state.m_matx_Y_train , 0 , 0 , 20 , ( *state.m_matx_Y_train).n_cols ) ;

            // 根据是否加权选择损失函数
            if( state.m_b_weighted_loss )
            {

                // Loss function, Use Y instead of R
                // pricing error of Y  
                // 计算带权重的 OLS 回归模型的加权残差平方和
                vec_output_criterion[ indexCutPoint ] = fastLm_weighted( ( *state.m_matx_Y_train ) , this->m_matx_regressor , ( *state.m_vec_loss_weight ) ) ;

            }
            else
            {

                // no weight on loss function, standard regression
                // 计算 OLS 回归模型的残差平方和
                vec_output_criterion[ indexCutPoint ] = fastLm( ( *state.m_matx_Y_train ) , this->m_matx_regressor ) ;

            }

            // 检查是否停止分裂  //检查分割是否改善整体损失:
            if( state.m_b_stop_no_gain )
            {

                 DEBUG_PRINT(" state.m_d_overall_loss  = " << state.m_d_overall_loss );

                // compare with overall loss, stop split if no gain
                if( vec_output_criterion[ indexCutPoint ] >= state.m_d_overall_loss )
                {

                    // if cannot improve overall pricing error, discard this split candidate
                    vec_output_criterion[ indexCutPoint ] = std::numeric_limits<double>::max( ) ;
                }

            }

        }


        DEBUG_PRINT("vec_output_criterion");
        printVec( vec_output_criterion , 0 , 20  ) ; 


        //printMat( *(Xorder) , 0 , 0 , 20 , Xorder->n_cols );

 	     //终止条件检查:
        // 如果所有数据都属于左侧，不需要遍历下一个更大的分裂点
        if( loop_index == ( *Xorder ).n_rows )
        {

            DEBUG_PRINT("loop_index : " << loop_index );

            // if loop_index = number of data, means that all observations belongs to left side
            // not necessary to loop over the next larger cutpoint
            break ;

        }

        DEBUG_PRINT("cutPoints_2 : " << curtCutPoint ) ;
        DEBUG_PRINT_SPACE

    }


    return ;


}


/*
它通过遍历切点并计算左右两侧的数据统计量来找到最优的分割点。
加权的部分则是为了提高模型的鲁棒性和准确性。
通过这种方式，函数能够在特定节点上基于一个变量来优化模型的拟合效果。

TS 可能表示 "Time Series"，即时间序列的意思。
因此，函数名的整体意思可能是基于时间序列的某种树模型（APTree）中计算特定变量的分割标准。
*/
//计算了在特定节点处基于一个变量的分割标准，其目的是找到最优的分割点以
//最小化回归模型的残差平方和（Sum of Squared Residuals, SSR）。
//它的加权部分使得不同的样本点对模型的贡献不同，从而更好地反映数据的特性和重要性。
//State& state: 包含了状态信息和参数。
//size_t var : 变量的索引。
//std::vector<CAPTree*>& vec_bottom_nodes : 包含所有叶节点的向量。
//size_t node_ind : 当前节点在 vec_bottom_nodes 中的索引。
//std::vector<double>& vec_output : 存储输出的向量。
//arma::vec& weighted_return_all : 所有数据的加权收益向量。
//arma::vec& cumu_weight_all : 所有数据的累计权重向量。
//arma::vec& num_stocks_all : 所有数据的股票数量向量。
//size_t var_ind : 变量在切点矩阵中的索引。
void CAPTreeModel::calculate_criterion_one_variable_APTree_TS( CState& state , 
                                                              size_t var , 
                                                              std::vector<CAPTree*>& vec_bottom_nodes , 
                                                              size_t node_ind , 
                                                              std::vector<double>& vec_output , 
                                                              arma::vec& vec_weighted_return_all , 
                                                              arma::vec& vec_cumu_weight_all , 
                                                              arma::vec& vec_num_stocks_all , 
                                                              size_t var_ind )
{

    DEBUG_PRINT("");

    // calculate split criterion for one variable at a specific node
    //获取当前节点指针：
    CAPTree* node = vec_bottom_nodes[ node_ind ] ;

    // initialize split criterion, start from infinity
    //初始化分割标准:
    std::fill( vec_output.begin( ) , vec_output.end( ) , std::numeric_limits<double>::max( ) ) ;

    // essentially, the sufficient statistics are two vectors with length num_months ;
    // first vector: weight * return
    // second vector: cumulative weight
    // the portfolio is just elementwise ratio of the two vectors
    size_t num_nodes = vec_bottom_nodes.size( ) ;
    arma::umat* Xorder = node->m_pMatx_Xorder ;  //获取变量排序矩阵

    // // calculate sufficient statistics of all data here
    //temp_index，temp_month，temp_month_index 用于存储临时索引和月份信息。
    size_t temp_index ;
    size_t temp_month ;
    size_t temp_month_index ;

    // next loop over cutpoints, calculate sufficient statistics on left / right side
    //初始化用于存储左、右侧加权收益、累计权重和股票数量的向量。
    arma::vec vec_weighted_return_left( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_cumu_weight_left( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_num_stocks_left( state.num_months , arma::fill::zeros ) ;

    arma::vec vec_weighted_return_right( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_cumu_weight_right( state.num_months , arma::fill::zeros ) ;
    arma::vec vec_num_stocks_right( state.num_months , arma::fill::zeros ) ;

    double cutpoint ;
    size_t loop_index = 0 ;
    arma::mat mu ;
    arma::mat sigma ;
    arma::mat weight ;
    arma::mat ft ;
    double weight_sum ;

    // for time series split, the months on the left / right sides are not the same
    size_t num_months_left = 0 ;
    size_t num_months_right = 0 ;

    // number of original data observations on the left / right side
    size_t num_obs_left = 0 ;
    size_t num_obs_right = 0 ;

    arma::mat all_portfolio( state.num_months , num_nodes + 1 , arma::fill::zeros ) ;
    temp_index = 2 ; // first two columns for the candidate split

    for( size_t i = 0 ; i < num_nodes ; i++ )
    {
        if( i != node_ind )
        {
            for( size_t ind = 0 ; ind < state.num_months ; ind++ )
            {
                all_portfolio( ind , temp_index ) = ( vec_bottom_nodes[ i ]->m_vec_month_theta )[ ind ] ;
            }

            temp_index++ ;
        }
    }

    for( size_t i = 0 ; i < state.num_cutPoints ; i++ )
    {

        //对每个切点，重置左、右侧向量。
        
        // reset all vectors for a new cutpoint
        vec_weighted_return_left.fill( 0.0 ) ;
        vec_weighted_return_right.fill( 0.0 ) ;
        vec_cumu_weight_left.fill( 0.0 ) ;
        vec_cumu_weight_right.fill( 0.0 ) ;
        vec_num_stocks_left.fill( 0.0 ) ;
        vec_num_stocks_right.fill( 0.0 ) ;

        // 获取当前切点： 
        cutpoint = ( *state.m_matx_first_split_mat )( i , var_ind ) ;

        //遍历数据点，根据切点将数据点分到左侧或右侧，并更新相应的向量
        while( ( *state.m_matx_X_train )( ( *Xorder )( loop_index , var ) , var ) <= cutpoint )
        {

            // the observation is on the left side
            temp_index = ( *Xorder )( loop_index , var ) ;
            temp_month = ( *state.m_vec_months )( temp_index ) ;
            temp_month_index = state.m_map_months_list->at( temp_month ) ;

            vec_weighted_return_left( temp_month_index ) += ( *state.m_vec_R_train )( temp_index ) * ( *state.m_vec_portfolio_weight )( temp_index ) ;
            vec_cumu_weight_left( temp_month_index ) += ( *state.m_vec_portfolio_weight )( temp_index ) ;
            vec_num_stocks_left( temp_month_index ) += 1.0 ;

            loop_index++ ;

            if( loop_index == ( *Xorder ).n_rows )
            {
                // terminating condition, avoid overflow
                break ;
            }
        }

        //计算右侧的统计量:
        vec_weighted_return_right = vec_weighted_return_all - vec_weighted_return_left ;
        vec_cumu_weight_right = vec_cumu_weight_all - vec_cumu_weight_left ;
        vec_num_stocks_right = vec_num_stocks_all - vec_num_stocks_left ;

        // cout << "some conditions " << endl ;
        // cout << (!state.flag_first_cut) << endl ;
        // cout << (num_stocks_right.min() < state.min_leaf_size) << endl ;
        // cout << (num_stocks_left.min() < state.min_leaf_size) << endl ;
        // cout << (arma::accu(num_stocks_right) == 0) << endl ;
        // cout << (arma::accu(num_stocks_left) == 0) << endl ;

        //检查叶节点的最小数据量要求:
        if( ( !state.m_b_flag_first_cut ) && ( vec_num_stocks_right.min( ) < state.min_leaf_size || vec_num_stocks_left.min( ) < state.min_leaf_size || arma::accu( vec_num_stocks_right ) == 0 || arma::accu( vec_num_stocks_left ) == 0 ) )
        {
            // too few data in the leaf
            //如果叶节点的数据量不足，则跳过该切点。
            vec_output[ i ] = std::numeric_limits<double>::max( ) ;
        }
        else
        {
            //计算左右两侧的投资组合:
            for( size_t ind = 0 ; ind < state.num_months ; ind++ )
            {
                // first column for the left portfolio
                all_portfolio( ind , 0 ) = ( vec_cumu_weight_left( ind ) == 0 ) ? 0 : vec_weighted_return_left( ind ) / vec_cumu_weight_left( ind ) ;

                // second column for the right portfolio
                all_portfolio( ind , 1 ) = ( vec_cumu_weight_right( ind ) == 0 ) ? 0 : vec_weighted_return_right( ind ) / vec_cumu_weight_right( ind ) ;
            }

            // count how many months on the left
            num_months_left = 0 ;
            for( size_t tt = 0 ; tt < vec_num_stocks_left.n_elem ; tt++ )
            {
                if( vec_num_stocks_left[ tt ] == 0 )
                {
                    num_months_left++ ;
                }
            }
            num_months_right = state.num_months - num_months_left ;

            // cout << "number of months on left " << num_months_left << " " << num_months_right << endl ;

            //计算分割后的投资组合的均值和协方差:
            mu = arma::mean( all_portfolio , 0 ) ; // 0 for column mean
            mu = arma::trans( mu ) ;              // transpose to column vectors
            sigma = arma::cov( all_portfolio ) ;

            size_t n_leafs = mu.n_elem ;

            //计算权重和调整权重:
            weight = arma::inv( sigma + state.m_d_lambda * arma::eye( n_leafs , n_leafs ) ) * mu ;

            arma::vec equal_weight( n_leafs ) ;

            equal_weight.fill( 1.0 / n_leafs ) ;

            weight = weight * state.m_d_eta + ( 1.0 - state.m_d_eta ) * equal_weight ;

            if( state.m_b_abs_normalize )
            {
                weight_sum = arma::accu( arma::abs( weight ) ) ;
            }
            else
            {
                weight_sum = arma::accu( ( weight ) ) ;
            }
            weight = weight / weight_sum ;

            ft = all_portfolio * weight ;

            num_obs_left = num_obs_right = 0 ;

            //更新回归矩阵:
            for( size_t i = 0 ; i < state.num_obs_all ; i++ )
            {
                for( size_t j = 0 ; j < ( *state.m_matx_Z_train ).n_cols ; j++ )
                {
                    temp_month_index = state.m_map_months_list->at( ( *state.m_vec_months )( i ) ) ;

                    if( vec_num_stocks_left( temp_month_index ) == 0 )
                    {
                        // no stock in the left leaf, look at the right one
                        this->m_matx_regressor( i , j ) = ( *state.m_matx_Z_train )( i , j ) * all_portfolio( temp_month_index , 1 ) ;
                        num_obs_right++ ;
                    }
                    else
                    {
                        // have stocks in the left leaf, empty at the right one
                        this->m_matx_regressor( i , j ) = ( *state.m_matx_Z_train )( i , j ) * all_portfolio( temp_month_index , 0 ) ;
                        num_obs_left++ ;
                    }
                }
            }

            size_t num_obs_left = loop_index ;
            size_t num_obs_right = state.num_obs_all - num_obs_left ;
            size_t num_regressor_cols = this->m_matx_regressor.n_cols ;

            arma::vec Y_left( num_obs_left ) ;
            arma::vec Y_right( num_obs_right ) ;
            arma::mat regressor_left( num_obs_left , num_regressor_cols ) ;
            arma::mat regressor_right( num_obs_right , num_regressor_cols ) ;

            //计算残差平方和（加权或未加权）:
            if( state.m_b_weighted_loss )
            {
                vec_output[ i ] = fastLm_weighted( ( *state.m_matx_Y_train ) , this->m_matx_regressor , ( *state.m_vec_loss_weight ) ) ;
            }
            else
            {
                vec_output[ i ] = fastLm( ( *state.m_matx_Y_train ) , this->m_matx_regressor ) ;
            }

            //检查分割是否改善整体损失:
            if( state.m_b_stop_no_gain )
            {
                // compare with overall loss, stop split if no gain
                if( vec_output[ i ] >= state.m_d_overall_loss )
                {
                    // if cannot improve overall pricing error, discard this split candidate
                    vec_output[ i ] = std::numeric_limits<double>::max( ) ;
                }
            }
        }

        //终止条件检查:
        if( loop_index == ( *Xorder ).n_rows )
        {
            // if loop_index = number of data, means that all observations belongs to left side
            // not necessary to loop over the next larger cutpoint
            break ;
        }
    }

    return ;
}



/*
 函数用于在时间序列分析（Time Series Analysis）中对一个树节点 pAPTreeNode 
 按照指定的分裂变量 split_var 和分裂点 split_point 进行分裂，生成左右子节点。
*/
void CAPTreeModel::split_node_APTree_TS( CState& state , 
                                         CAPTree* pAPTreeNode , 
                                         size_t split_var , 
                                         size_t split_point )
{

    DEBUG_PRINT("");

    // first, figure out how many are on the left side and right side
    arma::umat* Xorder = pAPTreeNode->m_pMatx_Xorder ; // 获取当前节点的数据顺序矩阵
    size_t num_obs_left = 0 ;
    size_t num_obs_right = 0 ;

    // first find the corresponding index in the first_split_var vector
    // 在 first_split_var 向量中查找 split_var 的索引
    //遍历 state.first_split_var 向量，找到 split_var 在向量中的索引 var。
    //这个步骤确保我们在分裂时使用正确的分裂点。
    size_t var = 999 ; 
    for( size_t i = 0 ; i < state.m_vec_first_split_var->n_elem ; i++ )
    {
        if( split_var == ( *state.m_vec_first_split_var )( i ) )
        {
            var = i ;
        }
    }

    // 统计左右子节点的样本数
    //遍历当前节点的数据，统计在分裂变量 split_var 上小于等于分裂点的样本数 
    //num_obs_left，以及大于分裂点的样本数 num_obs_right
    for( size_t i = 0 ; i < pAPTreeNode->get_numData_inNode( ) ; i++ )
    {
        ( ( *state.m_matx_X_train )( ( *Xorder )( i , split_var ) , split_var ) <= ( *state.m_matx_first_split_mat )( split_point , var ) ) ? num_obs_left++ : num_obs_right++ ;
    }

    double temp_split = ( *state.m_matx_first_split_mat )( split_point , var ) ;

    // 设置分裂点信息
    pAPTreeNode->set_varIndex2Split( split_var ) ;
    pAPTreeNode->set_valueIndex2Split( split_point ) ;
    pAPTreeNode->set_rawValue2Split( temp_split ) ;

    // 初始化左右子节点的 Xorder 矩阵
    arma::umat* Xorder_left = new arma::umat( num_obs_left , state.numOfCharitisc , arma::fill::zeros ) ;
    arma::umat* Xorder_right = new arma::umat( num_obs_right , state.numOfCharitisc , arma::fill::zeros ) ;

    // node->split_Xorder((*Xorder_left), (*Xorder_right), (*Xorder), split_point, split_var, state) ;

    // 创建左右子节点
    CAPTree::APTree_Pt lchild = new CAPTree( state.num_months , pAPTreeNode->get_treeDepth( ) + 1 , num_obs_left , pAPTreeNode->get_nodeID( ) * 2 , pAPTreeNode , Xorder_left ) ;
    CAPTree::APTree_Pt rchild = new CAPTree( state.num_months , pAPTreeNode->get_treeDepth( ) + 1 , num_obs_right , pAPTreeNode->get_nodeID( ) * 2 + 1 , pAPTreeNode , Xorder_right ) ;

    // 设置左右子节点
    pAPTreeNode->set_leftChild( lchild ) ;
    pAPTreeNode->set_rightChild( rchild ) ;

    // 初始化左右子节点的投资组合
    this->initialize_portfolio( state , lchild ) ;
    this->initialize_portfolio( state , rchild ) ;

    return ;
}


//用于将一个树节点 pAPTreeNode 按照指定的分裂变量 split_var 
//和分裂点 split_point 进行分裂，生成左右子节点
/*
正则化在 split_node 中的体现
在这个函数中，正则化的思想体现在数据分裂和节点初始化过程中：

数据分裂：
分裂过程中，通过选择合适的分裂变量和分裂点，可以有效地减少每个节点的数据量，
从而防止模型过拟合。

节点初始化：
在初始化左右子节点的投资组合时，可能会使用正则化方法来确保每个节点的投资组合
不受少量样本或异常值的过度影响。
*/
void CAPTreeModel::split_node( CState& state , 
                               CAPTree* pAPTreeNode , 
                               size_t split_var , 
                               size_t split_point )
{


    DEBUG_PRINT_SPACE ;
    DEBUG_PRINT("");


    DEBUG_PRINT("first, figure out how many are on the left side and right side");

    // first, figure out how many are on the left side and right side
    arma::umat* matx_pXorder = pAPTreeNode->m_pMatx_Xorder ; // 获取当前节点的数据顺序矩阵
    size_t num_obs_left = 0 ;
    size_t num_obs_right = 0 ;

 
    // 统计左右子节点的样本数
    DEBUG_PRINT(" get_numData_inNode : " << pAPTreeNode->get_numData_inNode( ) );

    printVec( state.m_vec_split_candidates , 0 , state.m_vec_split_candidates.size( ) ) ; 

    //printMat(  *matx_pXorder , 0 , 0 ,  20 , matx_pXorder->n_cols  ) ; 
    //printMat(  ( *state.m_matx_X_train ) , 0 , 0 ,  20 , (*state.m_matx_X_train).n_cols  ) ; 

    // 统计左右子节点的样本数
    //遍历当前节点的数据，统计在分裂变量 split_var 上小于等于分裂点的样本数 
    //num_obs_left，以及大于分裂点的样本数 num_obs_right
    for( size_t i = 0 ; i < pAPTreeNode->get_numData_inNode( ) ; i++ )
    {
        ( ( *state.m_matx_X_train )( ( *matx_pXorder )( i , split_var ) , split_var ) <= state.m_vec_split_candidates[ split_point ] )  ?  num_obs_left++ : num_obs_right++ ;
    }


    double temp_split = state.m_vec_split_candidates[ split_point ] ;
    DEBUG_PRINT(" temp_split : " << temp_split );

    // 设置分裂点信息
    DEBUG_PRINT(" Set split node infor " ) ;
    pAPTreeNode->set_varIndex2Split( split_var ) ;
    pAPTreeNode->set_valueIndex2Split( split_point ) ;
    pAPTreeNode->set_rawValue2Split( temp_split ) ;

    DEBUG_PRINT(" Init Xorder_left and Xorder_right ");

    // 初始化左右子节点的 Xorder 矩阵
    //Xorder_left 和 Xorder_right 分别存储左、右子节点的数据顺序索引矩阵。
    arma::umat* matx_Xorder_left  = new arma::umat( num_obs_left  , state.numOfCharitisc , arma::fill::zeros ) ;
    arma::umat* matx_Xorder_right = new arma::umat( num_obs_right , state.numOfCharitisc , arma::fill::zeros ) ;

    DEBUG_PRINT(" Split matx_pXorder into Xorder_left and Xorder_right ");
    //分裂 Xorder 矩阵
    //将 Xorder 矩阵分割成 Xorder_left 和 Xorder_right，分别对应左、右子节点的数据。
    pAPTreeNode->split_Xorder( ( *matx_Xorder_left ) , 
                               ( *matx_Xorder_right ) , 
                               ( *matx_pXorder ) , 
                               split_point , 
                               split_var , 
                               state ) ;

    DEBUG_PRINT_SPACE ;
    DEBUG_PRINT("===================================== create  pAPTree_LeftChild and pAPTree_RightChild ");
    DEBUG_PRINT_SPACE ;

    DEBUG_PRINT("create  pAPTree_LeftChild and pAPTree_RightChild ");
    // 创建左右子节点
    CAPTree::APTree_Pt pAPTree_LeftChild  = new CAPTree( state.num_months , 
                                                         pAPTreeNode->get_treeDepth( ) + 1 , 
                                                         num_obs_left , 
                                                         pAPTreeNode->get_nodeID( ) * 2 , 
                                                         pAPTreeNode , 
                                                         matx_Xorder_left ) ;
    CAPTree::APTree_Pt pAPTree_RightChild = new CAPTree( state.num_months , 
                                                         pAPTreeNode->get_treeDepth( ) + 1 , 
                                                         num_obs_right , 
                                                         pAPTreeNode->get_nodeID( ) * 2 + 1 , 
                                                         pAPTreeNode , 
                                                         matx_Xorder_right ) ;

    DEBUG_PRINT(" Add  pAPTree_LeftChild and pAPTree_RightChild ");
    // 设置左右子节点
    pAPTreeNode->set_leftChild( pAPTree_LeftChild ) ;
    pAPTreeNode->set_rightChild( pAPTree_RightChild ) ;

    // 初始化左右子节点的投资组合
    DEBUG_PRINT("initialize_portfolio   left child ");
    this->initialize_portfolio( state , pAPTree_LeftChild ) ;

    DEBUG_PRINT("initialize_portfolio   right child ");
    this->initialize_portfolio( state , pAPTree_RightChild ) ;

    return ;

}


//用于在给定节点初始化投资组合
/*
State& state: 表示当前状态，包括月数、投资组合回报率、权重等信息。
CAPTree* pAPTreeNode: 指向 CAP 树节点的指针，包含节点的排序信息和计算结果。

theta 是希腊字母 θ（θήτα）。在数学、统计学和金融领域，希腊字母经常用来表示变量、
参数和特定的数学概念。
在这个函数中，theta 被用作变量名，表示计算出的投资组合回报率或某种预测值。
使用希腊字母作为变量名是一种常见的做法，因为它们简洁且广泛认可。例如：

在统计学中，θ 通常表示模型的参数或估计值。
在金融中，θ 有时用来表示时间衰减（如期权定价模型中的“theta”）。
因此，在这个函数中，theta 作为变量名可能是为了表示其特殊意义，即每个月的投资组合回报率或预测值。
*/
void CAPTreeModel::initialize_portfolio( CState& state , CAPTree* pAPTreeNode )
{

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("");

    // initialize Rt at the given node
    // calculate equal weight / value weight portfolio return of a node
    size_t num_obs = ( *pAPTreeNode->m_pMatx_Xorder ).n_rows ;  // 当前节点的观测数。
    size_t month ;  //当前观测的月份
    size_t row_index ; //当前观测在数据中的行索引。
    size_t temp_month_index ;  //当前月份在月列表中的索引。
    std::vector<double> vec_month_weight_sum( state.num_months ) ;  //一个向量，用于存储每个月的权重和。


    arma::umat *tmppXOrder = pAPTreeNode->m_pMatx_Xorder ;
    //printMat( *tmppXOrder , 0 , 0 , 20 , tmppXOrder->n_cols);


    // 如果是等权重投资组合
    if( state.m_b_equal_weight )
    {

        DEBUG_PRINT(" equal_weight , loop times : " << num_obs);

        //如果 equal_weight 为真：遍历所有观测，将每个月的回报率相加，并计算权重和。
        for( size_t i = 0 ; i < num_obs ; i++ )
        {

            row_index = ( *pAPTreeNode->m_pMatx_Xorder )( i , 0 ) ;
            month   = ( *state.m_vec_months )[ row_index ] ;

            unsigned long long  month1 = (*state.m_vec_months)[row_index];
            //std::cout << " month: " << month1   ;

            temp_month_index = state.m_map_months_list->at( month ) ;
            //std::cout << "  temp_month_index: " << temp_month_index << std::endl  ;

            // 将当前观测的回报率累加到对应月份的theta中
            ( pAPTreeNode->m_vec_month_theta )[ temp_month_index ] += ( *state.m_vec_R_train )[ row_index ] ;
            
            // 对应月份的权重和加1
            vec_month_weight_sum[ temp_month_index ] = vec_month_weight_sum[ temp_month_index ] + 1 ;

        }

        DEBUG_PRINT_SPACE;
        std::cout << "print month weight : " << std::endl ; 
        for (size_t i = 0; i < vec_month_weight_sum.size() ; i++ )
        {

            //std::cout <<  i << " : " << vec_weight_sum[i] << std::endl;
        }
        DEBUG_PRINT_SPACE;

    }
    else  // 如果是加权投资组合
    {

        DEBUG_PRINT(" NOT equal_weight , loop times : " << num_obs );

        //如果 equal_weight 为假：遍历所有观测，将每个月的加权回报率相加，并计算加权权重和。
        for( size_t i = 0 ; i < num_obs ; i++ )
        {

            row_index = ( *pAPTreeNode->m_pMatx_Xorder )( i , 0 ) ;
            month   = ( *state.m_vec_months )[ row_index ] ;
            temp_month_index = state.m_map_months_list->at( month ) ;
            
            // 将加权回报率累加到对应月份的theta中
            ( pAPTreeNode->m_vec_month_theta )[ temp_month_index ] += ( *state.m_vec_R_train )[ row_index ] * ( *state.m_vec_portfolio_weight )[ row_index ] ;
           
            // 对应月份的权重和加上当前观测的权重
            vec_month_weight_sum[ temp_month_index ] = vec_month_weight_sum[ temp_month_index ] + ( *state.m_vec_portfolio_weight )[ row_index ] ;
        
        }

    }

    DEBUG_PRINT(" vec_month_weight_sum : "  );
    printVec( vec_month_weight_sum  , 0 , 20 ) ;

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT_SPACE;
    std::cout << "print month theta : " << std::endl;

    // 标准化每个月的回报率
    //遍历所有月份，如果权重和不为零，则将回报率除以权重和；否则，将回报率设为零。
    for( size_t i = 0 ; i < state.num_months ; i++ )
    {
        ( pAPTreeNode->m_vec_month_theta )[ i ] = ( vec_month_weight_sum[ i ] == 0 ) ? 0.0 : ( pAPTreeNode->m_vec_month_theta )[ i ] / vec_month_weight_sum[ i ] ;
    
        //std::cout << i << " : " << (pAPTreeNode->m_vec_theta)[i] << std::endl ;
    }


    printVec( pAPTreeNode->m_vec_month_theta  , 0 , 20 ) ;


    DEBUG_PRINT_SPACE;

    return ;

}



/*
用于初始化回归矩阵，目的是在模型中计算定价误差。
回归模型形式为 Yt ~ Zt * Ft + Ht，其中 Y 替代了 R 作为响应变量。
*/
void CAPTreeModel::initialize_regressor_matrix( CState& state )
{

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("");

    // initialize the regressor matrix in the model class
    // used in calculating pricing error
    // pre allocate space to save computing time
    // regress Yt ~ Zt * Ft + Ht 
    // Use Y instead of R

    size_t num_obs = state.num_obs_all ; //观测总数。
    size_t num_H   = ( *state.m_matx_H_train ).n_cols ;  //矩阵 H 的列数
    size_t num_Z   = ( *state.m_matx_Z_train ).n_cols ;  //矩阵 Z 的列数

     
    if( state.m_b_no_H ) // 如果没有矩阵H
    {
        
        DEBUG_PRINT(" state.no_H is TRUE : "  );

        //如果 state.no_H 为真，则回归矩阵的大小为 (num_obs, num_Z)。
        this->m_matx_regressor.resize( num_obs , num_Z ) ;
        this->m_matx_regressor.fill( arma::fill::zeros ) ;
    
    }
    else
    {

        DEBUG_PRINT(" state.no_H is FALSE , loop times : " << num_obs );

        printMat(   *state.m_matx_H_train , 0 , 0 , 20 , (*state.m_matx_H_train).n_cols ) ;

        //如果 state.no_H 为假，则回归矩阵的大小为(num_obs, num_H + num_Z)，
        //并将矩阵 H 的内容复制到回归矩阵的适当位置。
        this->m_matx_regressor.resize( num_obs , num_H + num_Z ) ;
        this->m_matx_regressor.fill( arma::fill::zeros ) ;

        for( size_t i = 0 ; i < num_obs ; i++ )
        {
            for( size_t j = 0 ; j < num_H ; j++ )
            {

                // first columns leave for the SDF
                // 在初始化回归矩阵时，前几列是预留给 SDF（Stochastic Discount Factor，随机折现因子）
                // 在金融和资产定价模型中，SDF 是一个重要的概念，用于将未来的现金流折现到当前值。
                // 它通常表示为一个回归模型中的解释变量，反映市场的风险和时间价值。
               
                // 将H的列存储在回归矩阵中Z之后的位置
                this->m_matx_regressor( i , j + num_Z ) = ( *state.m_matx_H_train )( i , j ) ;
            
            }
        
        }


    }


    //printMat( this->m_matx_regressor , 0 , 0 , 20 , this->m_matx_regressor.n_cols);

    return ;

}



/*
predict_AP 函数通过在给定的数据矩阵 X 上进行预测，找到每个数据点在树中的叶子节点，
并将其节点 ID 存储在 leaf_index 中。这使得我们可以知道每个数据点在决策树中的具体位置，
从而进一步进行分析和处理。
*/
void CAPTreeModel::predict_AP( arma::mat& matx_X , 
                               CAPTree& root , 
                               arma::vec& vec_months , 
                               arma::vec& vec_leaf_index )
{

    DEBUG_PRINT("");

    CAPTree* leaf ;

    // 遍历数据矩阵 X 的每一行
    for( size_t i = 0 ; i < matx_X.n_rows ; i++ )
    {
        // 找到数据点 X 第 i 行所在的叶子节点
        leaf = root.find_BtmNodeOfData( matx_X , i ) ;
        // 存储叶子节点的 ID
        vec_leaf_index( i ) = leaf->nid( ) ;
    }
    
    return ;

}


/*
用于计算给定 CAP 树模型的因子权重（leaf_weight）和因子回报（ft）。
具体来说，该函数首先确定所有叶节点的权重，然后根据这些权重计算每个月的因子回报。

输入参数：
CAPTree& root: CAP树的根节点。
arma::vec& leaf_node_index: 叶节点索引向量。
arma::mat& all_leaf_portfolio: 所有叶节点的投资组合矩阵。
arma::mat& leaf_weight: 叶节点的权重矩阵。
arma::mat& ft: 因子回报矩阵。
State& state: 包含计算所需的状态信息，如月份数、正则化参数等。

*/

void CAPTreeModel::calculate_factor( CAPTree& root , 
                                     arma::vec& vec_leaf_node_index , 
                                     arma::mat& matx_all_leaf_portfolio , 
                                     arma::mat& matx_leaf_weight , 
                                     arma::mat& matx_fator , 
                                     CState& state )
{

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("");


    // 获取所有叶节点
    std::vector<CAPTree*> vec_bottom_nodes ;
    // once fitting is done, calculate weight of all leaf nodes
    vec_bottom_nodes.resize( 0 ) ;
    root.get_vecOfBtmNodes( vec_bottom_nodes ) ;

    // 初始化叶节点索引和投资组合矩阵
    vec_leaf_node_index.resize( vec_bottom_nodes.size( ) ) ;
    vec_leaf_node_index.fill( arma::fill::zeros ) ;
    matx_all_leaf_portfolio.resize( state.num_months , vec_bottom_nodes.size( ) ) ;
    matx_all_leaf_portfolio.fill( arma::fill::zeros ) ;

    // 填充叶节点索引和投资组合矩阵
    for( size_t i = 0 ; i < vec_bottom_nodes.size( ) ; i++ )
    {
        vec_leaf_node_index( i ) = vec_bottom_nodes[ i ]->nid( ) ;

        for( size_t j = 0 ; j < state.num_months ; j++ )
        {
            matx_all_leaf_portfolio( j , i ) = ( vec_bottom_nodes[ i ]->m_vec_month_theta )[ j ] ;
        }
    }

 
    DEBUG_PRINT("vec_leaf_node_index ");
    printVec( vec_leaf_node_index   , 0 , 20  ) ; 

    DEBUG_PRINT("matx_all_leaf_portfolio ");
    printMat( matx_all_leaf_portfolio   , 0 , 0 , 20 ,  matx_all_leaf_portfolio.n_cols ) ;


    // 计算叶节点投资组合的均值向量 mu 和协方差矩阵sigma
    arma::mat mu = arma::mean( matx_all_leaf_portfolio , 0 ) ;
    mu = arma::trans( mu ) ;

    DEBUG_PRINT("mu ");
    printMat( mu   , 0 , 0 , 20 ,  mu.n_cols ) ;

    size_t n_leafs = mu.n_elem ;
    arma::mat sigma = arma::cov( matx_all_leaf_portfolio ) ;

    DEBUG_PRINT("sigma ");
    printMat( sigma   , 0 , 0 , 20 ,  sigma.n_cols ) ;

    //通过加权最小二乘法计算叶节点的权重，并结合等权重进行平滑处理。
    // 计算叶节点权重
    // P11底下： 收缩参数
    /*
    arma::inv 计算矩阵的逆。
    sigma + state.lambda_cov * arma::eye(n_leafs, n_leafs) 是带有正则化项的协方差矩阵。
    mu + state.lambda_mean * arma::ones(mu.n_rows, mu.n_cols) 是带有正则化项的均值向量。
    */
    std::cout << "lambda_mean: " << state.m_d_lambda_mean << "\n";
    std::cout << "lambda_cov: "  << state.m_d_lambda_cov << "\n";

    matx_leaf_weight = arma::inv( sigma + state.m_d_lambda_cov * arma::eye( n_leafs , n_leafs ) ) 
                       * ( mu + state.m_d_lambda_mean * arma::ones( mu.n_rows , mu.n_cols ) ) ;

    // 计算等权重
    arma::vec vec_equal_weight( n_leafs ) ;
    vec_equal_weight.fill( 1.0 / n_leafs ) ;

    // 综合考虑等权重和计算的权重
    // state.eta 是平滑参数，控制计算权重和等权重的比例。
    matx_leaf_weight = matx_leaf_weight * state.m_d_eta + ( 1.0 - state.m_d_eta ) * vec_equal_weight ;


    //对叶节点权重进行归一化处理，确保其和为 1  
    double weight_sum ;

    if( state.m_b_abs_normalize )
    {
        weight_sum = arma::accu( arma::abs( matx_leaf_weight ) ) ;
    }
    else
    {
        weight_sum = arma::accu( ( matx_leaf_weight ) ) ;
    }

    DEBUG_PRINT("matx_leaf_weight ");
    printMat( matx_leaf_weight   , 0 , 0 , 20 ,  matx_leaf_weight.n_cols ) ;


    matx_leaf_weight = matx_leaf_weight / weight_sum ;

    DEBUG_PRINT("matx_leaf_weight ");
    printMat( matx_leaf_weight   , 0 , 0 , 20 ,  matx_leaf_weight.n_cols ) ;

    // 计算因子回报
    matx_fator = matx_all_leaf_portfolio * matx_leaf_weight ;

    // 如果因子回报为负，则进行调整
    if( arma::accu( matx_fator ) < 0 )
    {
        // if the average return is negative, short it
        matx_leaf_weight = matx_leaf_weight * ( -1.0 ) ;
        matx_fator = matx_all_leaf_portfolio * matx_leaf_weight ;
    }


    DEBUG_PRINT("matx_all_leaf_portfolio ");
    printMat( matx_leaf_weight   , 0 , 0 , 20 ,  matx_leaf_weight.n_cols ) ;

    DEBUG_PRINT("matx_all_leaf_portfolio ");
    printMat( matx_fator   , 0 , 0 , 20 ,  matx_fator.n_cols ) ;

    DEBUG_PRINT_SPACE

    return ;


}



//函数用于计算模型的 R平方值
//这是一种评估模型拟合效果的指标。它通过将实际值与预测值之间的误差进行比较来衡量模型的解释力
double CAPTreeModel::calculate_R2( CState& state , 
                                   arma::mat& matx_factor )
{

    DEBUG_PRINT("");

    arma::mat regressor ;  // 用于存储回归矩阵。
    double loss = 0.0 ;  //用于存储损失值
    size_t temp_month_index ;  //临时存储月份索引。

    if( !state.m_b_no_H )
    {
        //如果 state.no_H 为 false，则 regressor 的列数为 Z 和 H 矩阵列数之和。
        regressor.resize( state.num_obs_all , ( *state.m_matx_Z_train ).n_cols + ( *state.m_matx_H_train ).n_cols ) ;
        regressor.fill( arma::fill::zeros ) ;  //初始化为全零矩阵

        // 将 H 矩阵的数据填入 regressor 中
        for( size_t i = 0 ; i < state.num_obs_all ; i++ )
        {
            for( size_t j = 0 ; j < ( *state.m_matx_H_train ).n_cols ; j++ )
            {
                regressor( i , j + ( *state.m_matx_Z_train ).n_cols ) = ( *state.m_matx_H_train )( i , j ) ;
            }
        }
    }
    else
    {
        //如果 state.no_H 为 true，则 regressor 的列数仅为 Z 矩阵的列数。
        regressor.resize( state.num_obs_all , ( *state.m_matx_Z_train ).n_cols ) ;
        regressor.fill( arma::fill::zeros ) ; 
    }


    DEBUG_PRINT("regressor 1: ");
    printMat( regressor   , 0 , 0 , 20 ,  regressor.n_cols ) ;

    //将 Z 矩阵与 ft 矩阵的乘积填入 regressor 的相应列。
    for( size_t i = 0 ; i < state.num_obs_all ; i++ )
    {
        for( size_t j = 0 ; j < ( *state.m_matx_Z_train ).n_cols ; j++ )
        {
            temp_month_index = state.m_map_months_list->at( ( *state.m_vec_months )( i ) ) ;
            regressor( i , j ) = ( *state.m_matx_Z_train )( i , j ) * matx_factor( temp_month_index , 0 ) ;
        }
    }

    DEBUG_PRINT("regressor 2 : ");
    printMat( regressor   , 0 , 0 , 20 ,  regressor.n_cols ) ;

    // 根据是否使用加权损失函数，计算损失
    if( state.m_b_weighted_loss )
    {
        //如果 state.weighted_loss 为 true，使用加权线性回归 fastLm_weighted 计算损失。
        loss = fastLm_weighted( ( *state.m_matx_Y_train ) , regressor , ( *state.m_vec_loss_weight ) ) ;
    }
    else
    {
        //使用普通线性回归 fastLm 计算损失
        loss = fastLm( ( *state.m_matx_Y_train ) , regressor ) ;
    }

    // 计算 R^2 值
    // 分子 loss 是通过线性回归计算的残差平方和，
    // 分母 arma::accu( pow( *state.Y , 2 ) ) 是实际值的平方和。
    loss = 1 - loss / arma::accu( pow( *state.m_matx_Y_train , 2 ) ) ;

    return loss ;
}