#include "APTree.h"
#include <chrono>
#include <ctime>

#include "common.h"



//--------------------------------------------------
std::ostream& operator<<( std::ostream& os , const CAPTree& APTree )
{

    CAPTree::vec_APTree_cnstPt nds ;

    APTree.get_AllNodes_const( nds ) ;
    os << nds.size( ) << std::endl ;

    // size_t theta_length = nds[0]->getthetasize() ;
    // cout << "theta length is " << theta_length << endl ;
    for( size_t i = 0 ; i < nds.size( ) ; i++ )
    {
        os << nds[ i ]->nid( ) << " " ;
        os << nds[ i ]->get_varIndex2Split( ) << " " ;
        os << nds[ i ]->get_rawValue2Split( ) << " " ;
        os << nds[ i ]->get_valueIndex2Split( ) << " " ;
        os << nds[ i ]->get_iter( ) ;
        // for (size_t j = 0 ; j < theta_length ; j++)
        // {
        //     os << " " << nds[i]->gettheta(j) ;
        //     // os << " " << nds[i]->getRt(j) ;
        // }
        os << std::endl ;
    
    }
    
    return os ;
}


std::istream& operator>>( std::istream& is , CAPTree& APTree )
{

    size_t curtID , paretID ;                        //tid: id of current node, pid: parent's id
    std::map<size_t , CAPTree::APTree_Pt> map_ID2PT ; //pointers to nodes indexed by node id
    size_t numOfNodes ;                              //number of nodes

    APTree.toNull( ) ; // obliterate old tree (if there)

    //read number of nodes----------
    is >> numOfNodes ;
    if( !is )
    {
        return is ;
    }

    // The idea is to dump string to a lot of node_info structure first, then link them as a tree, by nid

    //read in vector of node information----------
    std::vector<node_info> vec_NodeInfor( numOfNodes ) ;
    for( size_t i = 0 ; i != numOfNodes ; i++ )
    {
        is >> vec_NodeInfor[ i ].id >> vec_NodeInfor[ i ].var >> vec_NodeInfor[ i ].cutPoint >> vec_NodeInfor[ i ].vec_theta[ 0 ] ; // Only works on first theta for now, fix latex if needed
        if( !is )
        {
            return is ;
        }
    }

    //first node has to be the top one
    map_ID2PT[ 1 ] = &APTree ; //be careful! this is not the first pts, it is pointer of id 1.
    APTree.set_varIndex2Split( vec_NodeInfor[ 0 ].var ) ;
    APTree.set_rawValue2Split( vec_NodeInfor[ 0 ].cutPoint ) ;
    APTree.set_theta( vec_NodeInfor[ 0 ].vec_theta ) ;
    APTree.m_parentNode = 0 ;

    //now loop through the rest of the nodes knowing parent is already there.
    for( size_t i = 1 ; i != vec_NodeInfor.size( ) ; i++ )
    {

        CAPTree::APTree_Pt pTreeNode = new CAPTree ;

        pTreeNode->m_varIndex2Split = vec_NodeInfor[ i ].var ;
        pTreeNode->m_rawValue2Split = vec_NodeInfor[ i ].cutPoint ;
        pTreeNode->m_vec_month_theta = vec_NodeInfor[ i ].vec_theta ;
        curtID = vec_NodeInfor[ i ].id ;
        map_ID2PT[ curtID ] = pTreeNode ;
        paretID = curtID / 2 ;

        if( curtID % 2 == 0 )
        { //left child has even id
            map_ID2PT[ paretID ]->m_leftChild = pTreeNode ;
        }
        else
        {
            map_ID2PT[ paretID ]->m_rightChild = pTreeNode ;
        }

        pTreeNode->m_parentNode = map_ID2PT[ paretID ] ;

    }

    return is ;

}



size_t CAPTree::nid( ) const
{

    if( !m_parentNode )
        return 1 ; //if you don't have a parent, you are the top

    if( this == m_parentNode->m_leftChild )
        return 2 * ( m_parentNode->nid( ) ) ; //if you are a left child
    else
        return 2 * ( m_parentNode->nid( ) ) + 1 ; //else you are a right child

}

CAPTree::APTree_Pt CAPTree::get_pt2ParentByNodeID( size_t nid )
{

    if( this->nid( ) == nid )
        return this ; //found it

    if( m_leftChild == 0 )
        return 0 ; //no children, did not find it

    APTree_Pt lp = m_leftChild->get_pt2ParentByNodeID( nid ) ;
    if( lp )
        return lp ; //found on left

    APTree_Pt rp = m_rightChild->get_pt2ParentByNodeID( nid ) ;
    if( rp )
        return rp ; //found on right
    
    return 0 ;      //never found it

}

size_t CAPTree::treeSize( )
{
    if( m_leftChild == 0 )
        return 1 ; //if bottom node, tree size is 1
    else
        return ( 1 + m_leftChild->treeSize( ) + m_rightChild->treeSize( ) ) ;
}

char CAPTree::nodeType( )
{
    //t:top, b:bottom, n:no grandchildren, i:internal
    if( !m_parentNode )
        return 't' ;

    if( !m_leftChild )
        return 'b' ;

    if( !( m_leftChild->m_leftChild ) && !( m_rightChild->m_leftChild ) )
        return 'n' ;

    return 'i' ;

}



void CAPTree::printScreen( bool pc )
{

    size_t depth = this->m_treeDepth ;
    size_t id = nid( ) ;
    size_t pid ;

    if( !m_parentNode )
        pid = 0 ; //parent of top node
    else
        pid = m_parentNode->nid( ) ;

    std::string pad( 2 * depth , ' ' ) ;
    std::string strSplit( ", " ) ;

    if( pc && ( nodeType( ) == 't' ) )
        std::cout << "tree size: " << treeSize( ) << std::endl ;

    std::cout << pad << "(id,parent): " << id << strSplit << pid ;
    std::cout << strSplit << "(v,c): " << m_varIndex2Split << strSplit << m_rawValue2Split ;
    // std::cout << sp << "theta: " << theta ;
    std::cout << strSplit << "type: " << nodeType( ) ;
    std::cout << strSplit << "depth: " << this->m_treeDepth ;
    std::cout << strSplit << "pointer: " << this << std::endl ;

    if( pc )
    {

        if( m_leftChild )
        {
            m_leftChild->printScreen( pc ) ;
            m_rightChild->printScreen( pc ) ;
        }

    }

}



bool CAPTree::isNoGrandChildren( )
{

    bool isnog = true ;

    if( m_leftChild )
    {
        if( m_leftChild->m_leftChild || m_rightChild->m_leftChild )
            isnog = false ; //one of the children has children.
    }
    else
    {
        isnog = false ; //no children
    }

    return isnog ;

}


size_t CAPTree::numNoGrandChildsNodes( )
{

    if( !m_leftChild )
        return 0 ; //bottom node

    if( m_leftChild->m_leftChild || m_rightChild->m_leftChild )
    {   
        //not a nog
        return ( m_leftChild->numNoGrandChildsNodes( ) + m_rightChild->numNoGrandChildsNodes( ) ) ;
    }
    else
    { 
        //is a nog
        return 1 ;
    }

}

size_t CAPTree::numLeafNodes( )
{

    DEBUG_PRINT("");

    if( m_leftChild == 0 )
    { //if a bottom node
        return 1 ;
    }
    else
    {
        return m_leftChild->numLeafNodes( ) + m_rightChild->numLeafNodes( ) ;
    }

}

void CAPTree::get_vecOfBtmNodes( vec_APTree_Pt& vec_BtmNodes )
{

    DEBUG_PRINT("");

    if( m_leftChild )
    { 
        DEBUG_PRINT( " m_leftChild  exsit "   );

        //have children
        m_leftChild->get_vecOfBtmNodes( vec_BtmNodes ) ;
        m_rightChild->get_vecOfBtmNodes( vec_BtmNodes ) ;
    
    }
    else
    {
        DEBUG_PRINT( " no children yet :  add this CAPTree OBJECT "   );
    
        vec_BtmNodes.push_back( this ) ;
    
    }

}

void CAPTree::get_vecNoGrandChildsNodes( vec_APTree_Pt& vec_NoGrandChildsNodes )
{

    DEBUG_PRINT("");

    if( m_leftChild )
    { 

        DEBUG_PRINT( " m_leftChild  exsit :  have children "   );

        //have children
        if( ( m_leftChild->m_leftChild ) || ( m_rightChild->m_leftChild ) )
        { 
            DEBUG_PRINT( " and have grandchildren "   );

            //have grandchildren
            if( m_leftChild->m_leftChild )
                m_leftChild->get_vecNoGrandChildsNodes( vec_NoGrandChildsNodes ) ;

            if( m_rightChild->m_leftChild )
                m_rightChild->get_vecNoGrandChildsNodes( vec_NoGrandChildsNodes ) ;
        }
        else
        {

            DEBUG_PRINT( " NO grandchildren : add this CAPTree OBJECT "   ) ;

            vec_NoGrandChildsNodes.push_back( this ) ;
        }
    }

}

CAPTree::APTree_Pt CAPTree::get_pt2TopNode( )
{

    DEBUG_PRINT("");

    if( !m_parentNode )
    {
        return this ;
    }
    else
    {
        return m_parentNode->get_pt2TopNode( ) ;
    }

}

void CAPTree::get_AllNodes( vec_APTree_Pt& vec_AllNodes )
{

    DEBUG_PRINT("");

    vec_AllNodes.push_back( this ) ;

    if( m_leftChild )
    {
        m_leftChild->get_AllNodes( vec_AllNodes ) ;
        m_rightChild->get_AllNodes( vec_AllNodes ) ;
    }

}


void CAPTree::get_AllNodes_const( vec_APTree_cnstPt& vec_AllNodesConst ) const
{

    DEBUG_PRINT("");

    vec_AllNodesConst.push_back( this ) ;

    if( m_leftChild )
    {
        m_leftChild->get_AllNodes_const( vec_AllNodesConst ) ;
        m_rightChild->get_AllNodes_const( vec_AllNodesConst ) ;
    }

}



CAPTree::APTree_Pt CAPTree::find_BtmNodeOfData( arma::mat& matx_x , size_t& row_ind )
{

    DEBUG_PRINT("");

    // v is variable to split, c is raw value
    // not index in matrix<double>, so compare x[v] with c directly
    if( m_leftChild == 0 )
        return this ;

    if( matx_x( row_ind , m_varIndex2Split ) <= m_rawValue2Split )
    {
        return m_leftChild->find_BtmNodeOfData( matx_x , row_ind ) ;
    }
    else
    {
        return m_rightChild->find_BtmNodeOfData( matx_x , row_ind ) ;
    }

}

void CAPTree::toNull( )
{

    DEBUG_PRINT("");

    size_t tree_size = treeSize( ) ;

    //loop invariant: ts>=1
    while( tree_size > 1 )
    { 
        //if false ts=1
        vec_APTree_Pt vec_APTreePt ;
        get_vecNoGrandChildsNodes( vec_APTreePt ) ;

        for( size_t i = 0 ; i < vec_APTreePt.size( ) ; i++ )
        {
            delete vec_APTreePt[ i ]->m_leftChild ;
            delete vec_APTreePt[ i ]->m_rightChild ;

            vec_APTreePt[ i ]->m_leftChild = 0 ;
            vec_APTreePt[ i ]->m_rightChild = 0 ;

        }

        tree_size = treeSize( ) ; //make invariant true

    }

    m_varIndex2Split = 0 ;
    m_rawValue2Split = 0 ;

    m_parentNode = 0 ;
    m_leftChild = 0 ;
    m_rightChild = 0 ;

}


//copy tree tree o to tree n
void CAPTree::copyTree( APTree_Pt pNewTree , APTree_cnstPt pOldTree )
//assume n has no children (so we don't have to kill them)
//recursion down
// create a new copy of tree in NEW memory space
{

    DEBUG_PRINT("");


    if( pNewTree->m_leftChild )
    {
        std::cout << "cp:error node has children\n" ;
        return ;
    }

    pNewTree->m_varIndex2Split = pOldTree->m_varIndex2Split ;
    pNewTree->m_rawValue2Split = pOldTree->m_rawValue2Split ;
    pNewTree->m_vec_month_theta = pOldTree->m_vec_month_theta ;

    if( pOldTree->m_leftChild )
    { 
        //if o has children
        pNewTree->m_leftChild = new CAPTree ;
        ( pNewTree->m_leftChild )->m_parentNode = pNewTree ;
        copyTree( pNewTree->m_leftChild , pOldTree->m_leftChild ) ;

        pNewTree->m_rightChild = new CAPTree ;
        ( pNewTree->m_rightChild )->m_parentNode = pNewTree ;
        copyTree( pNewTree->m_rightChild , pOldTree->m_rightChild ) ;
    }

}

void CAPTree::copy_only_root( APTree_Pt pOldTree )
//assume n has no children (so we don't have to kill them)
//NOT LIKE cp() function
//this function pointer new root to the OLD structure
{

    DEBUG_PRINT("");

    this->m_varIndex2Split = pOldTree->m_varIndex2Split ;
    this->m_rawValue2Split = pOldTree->m_rawValue2Split ;
    this->m_vec_month_theta = pOldTree->m_vec_month_theta ;

    if( pOldTree->m_leftChild )
    {
        // keep the following structure, rather than create a new tree in memory
        this->m_leftChild = pOldTree->m_leftChild ;
        this->m_rightChild = pOldTree->m_rightChild ;

        // also update pointers to parents
        this->m_leftChild->m_parentNode = this ;
        this->m_rightChild->m_parentNode = this ;
    }
    else
    {
        this->m_leftChild = 0 ;
        this->m_rightChild = 0 ;
    }

}



//--------------------------------------------------
//operators
CAPTree& CAPTree::operator=( const CAPTree& rightHandSide )
{

    if( &rightHandSide != this )
    {
        toNull( ) ;       //kill left hand side (this)
        copyTree( this , &rightHandSide ) ; //copy right hand side to left hand side
    }

    return *this ;

}


/*
split_Xorder 函数通过将当前节点的数据顺序矩阵 Xorder 根据给定的分裂点和分裂变量
分割为左右子节点的数据顺序矩阵 Xorder_left 和 Xorder_right，实现了决策树的分裂过程。
*/
void CAPTree::split_Xorder( arma::umat& matx_Xorder_left , 
                            arma::umat& matx_Xorder_right , 
                            arma::umat& matx_Xorder , 
                            size_t split_point , 
                            size_t split_var , 
                            CState& state )
{


    DEBUG_PRINT_SPACE ;

    DEBUG_PRINT("") ;


    size_t num_obs = matx_Xorder.n_rows ;
    
    //printMat(  matx_Xorder , 0 , 0 ,  20 , matx_Xorder.n_cols  ) ; 

    //确定分裂值 cutvalue。
    double curtCutValue = state.m_vec_split_candidates[ split_point ] ;

    size_t left_index ;
    size_t right_index ;
    for( size_t i = 0 ; i < state.numOfCharitisc ; i++ )
    {
        //初始化左右子节点数据的索引 left_index 和 right_index。
        left_index = 0 ;
        right_index = 0 ;

        // loop over variables  遍历每个观察样本
        //对于每一个变量，遍历所有样本，检查其在分裂变量上的值是否小于等于 curtCutValue，
        //将样本索引分配到 Xorder_left 或 Xorder_right 中
        for( size_t j = 0 ; j < num_obs ; j++ )
        {

            // loop over observations
            if( ( *state.m_matx_X_train )( matx_Xorder( j , i ) , split_var ) <= curtCutValue )
            {
                //将符合条件的样本索引存储在左子节点的 Xorder
                // left side  
                matx_Xorder_left( left_index , i ) = matx_Xorder( j , i ) ;
                left_index++ ;

            }
            else
            {
                // 将不符合条件的样本索引存储在右子节点的 Xorder
                // right side
                matx_Xorder_right( right_index , i ) = matx_Xorder( j , i ) ;
                right_index++ ;

            }
        }

    }

    return ;

}

json CAPTree::to_json( )
{

    json j ;
    
    if( m_leftChild == 0 )
    {
        j = this->m_vec_month_theta ;
    }
    else
    {
        j[ "variable" ] = this->m_varIndex2Split ;
        j[ "cutpoint" ] = this->m_rawValue2Split ;
        j[ "cutpoint_index" ] = this->m_valueIndex2Split ;
        j[ "nodeid" ] = this->nid( ) ;
        j[ "depth" ] = this->m_treeDepth ;
        j[ "left" ] = this->m_leftChild->to_json( ) ;
        j[ "right" ] = this->m_rightChild->to_json( ) ;
    }

    return j ;

}


void CAPTree::from_json( json& j3 , size_t dim_theta )
{

    if( j3.is_array( ) )
    {

        // this is the leaf
        std::vector<double> temp ;
        j3.get_to( temp ) ;

        if( temp.size( ) > 1 )
        {
            this->m_vec_month_theta = temp ;
        }
        else
        {
            this->m_vec_month_theta[ 0 ] = temp[ 0 ] ;
        }

    }
    else
    {

        // this is an intermediate node
        j3.at( "variable" ).get_to( this->m_varIndex2Split ) ;
        j3.at( "cutpoint" ).get_to( this->m_rawValue2Split ) ;
        j3.at( "cutpoint_index" ).get_to( this->m_valueIndex2Split ) ;
        j3.at( "depth" ).get_to( this->m_treeDepth ) ;

        CAPTree* lchild = new CAPTree( dim_theta ) ;
        lchild->from_json( j3[ "left" ] , dim_theta ) ;
        CAPTree* rchild = new CAPTree( dim_theta ) ;
        rchild->from_json( j3[ "right" ] , dim_theta ) ;

        lchild->m_parentNode = this ;
        rchild->m_parentNode = this ;
        this->m_leftChild = lchild ;
        this->m_rightChild = rchild ;
    }

}

  
void CAPTree::grow( bool& b_breakFlag , 
                    CAPTreeModel& model , 
                    CState& state , 
                    size_t& iter , 
                    std::vector<double>& vec_criterion_values )
{

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("");
    DEBUG_PRINT_SPACE;

    std::vector<CAPTree*> vec_bottom_nodes ;
    std::vector<bool> vec_node_splitability ;

    size_t split_node =  0;
    size_t split_var  = 0 ;
    size_t split_point = 0  ;
    bool splitable = true ;

    // grow a tree by iteration instead of recursion

    cout << "first, find all leaves"   << endl ;

    // first, find all leaves
    vec_bottom_nodes.resize( 0 ) ;
    this->get_vecOfBtmNodes( vec_bottom_nodes ) ;

    cout << "second, check splitability "   << endl ;

    // second, check splitability, 1 for splitable, 0 for terminated
    vec_node_splitability.resize( vec_bottom_nodes.size( ) ) ;
    model.check_node_splitability( state , vec_bottom_nodes , vec_node_splitability ) ;

    if( sum( vec_node_splitability ) )
    {

        cout << "exist at least one node for split : " <<  sum( vec_node_splitability ) << endl ;
        cout << "third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point "   << endl ;

        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion( state , 
                                   vec_bottom_nodes , 
                                   vec_node_splitability , 
                                   split_node , 
                                   split_var , 
                                   split_point , 
                                   splitable , 
                                   vec_criterion_values ) ;
        // split the selected node

        if( splitable )
        {

            arma::umat *tmppXOrder = vec_bottom_nodes[ split_node ]->m_pMatx_Xorder ;
            //printMat( *tmppXOrder , 0 , 0 , 20 , tmppXOrder->n_cols);

            vec_bottom_nodes[ split_node ]->set_iter( iter ) ;

            //printMat( *tmppXOrder , 0 , 0 , 20 , tmppXOrder->n_cols);

            model.split_node( state , vec_bottom_nodes[ split_node ] , split_var , split_point ) ;
        
        }
        else
        {

            cout << "break of no GOOD candidate 1 " << endl ;
        
            b_breakFlag = true ;
        
        }

    }
    else
    {

        cout << "break of no NODE splitable 2 " << endl ;

        b_breakFlag = true ;
    
    }

    return ;
}

void CAPTree::grow_APTree_TS( bool& b_breakFlag , CAPTreeModel& model , CState& state )
{

    DEBUG_PRINT("");

    std::vector<CAPTree*> vec_bottom_nodes ;
    std::vector<bool> vec_node_splitability ;

    size_t split_node ;
    size_t split_var ;
    size_t split_point ;
    bool bSplitable = true ;

    // grow a tree by iteration instead of recursion
    // first, find all leaves
    vec_bottom_nodes.resize( 0 ) ;
    this->get_vecOfBtmNodes( vec_bottom_nodes ) ;

    // second, check splitability, 1 for splitable, 0 for terminated
    vec_node_splitability.resize( vec_bottom_nodes.size( ) ) ;
    model.check_node_splitability( state , 
                                   vec_bottom_nodes , 
                                   vec_node_splitability ) ;

    if( sum( vec_node_splitability ) )
    {
        
        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion_APTree_TS(  state , 
                                              vec_bottom_nodes , 
                                              vec_node_splitability , 
                                              split_node , 
                                              split_var , 
                                              split_point , 
                                              bSplitable ) ;
        // split the selected node
        if( bSplitable )
        {
            model.split_node_APTree_TS( state , 
                                        vec_bottom_nodes[ split_node ] , 
                                        split_var , 
                                        split_point ) ;
        }
        else
        {
            cout << "break of no good candidate" << endl ;
            b_breakFlag = true ;
        }
        
    }
    else
    {
        cout << "break of no node splitable" << endl ;
        b_breakFlag = true ;
    }

    return ;
}
