#include "APTree.h"
#include <chrono>
#include <ctime>

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

    if( m_leftChild == 0 )
    { //if a bottom node
        return 1 ;
    }
    else
    {
        return m_leftChild->numLeafNodes( ) + m_rightChild->numLeafNodes( ) ;
    }

}

void CAPTree::get_vecOfBtmNodes( vec_APTree_Pt& bv )
{

    if( m_leftChild )
    { //have children
        m_leftChild->get_vecOfBtmNodes( bv ) ;
        m_rightChild->get_vecOfBtmNodes( bv ) ;
    }
    else
    {
        bv.push_back( this ) ;
    }

}

void CAPTree::getNoGrandChildsNodes( vec_APTree_Pt& nv )
{

    if( m_leftChild )
    { 
        //have children
        if( ( m_leftChild->m_leftChild ) || ( m_rightChild->m_leftChild ) )
        { 
            //have grandchildren
            if( m_leftChild->m_leftChild )
                m_leftChild->getNoGrandChildsNodes( nv ) ;

            if( m_rightChild->m_leftChild )
                m_rightChild->getNoGrandChildsNodes( nv ) ;
        }
        else
        {
            nv.push_back( this ) ;
        }
    }

}

CAPTree::APTree_Pt CAPTree::get_pt2TopNode( )
{

    if( !m_parentNode )
    {
        return this ;
    }
    else
    {
        return m_parentNode->get_pt2TopNode( ) ;
    }

}

void CAPTree::get_AllNodes( vec_APTree_Pt& v )
{

    v.push_back( this ) ;

    if( m_leftChild )
    {
        m_leftChild->get_AllNodes( v ) ;
        m_rightChild->get_AllNodes( v ) ;
    }

}


void CAPTree::get_AllNodes_const( vec_APTree_cnstPt& v ) const
{

    v.push_back( this ) ;

    if( m_leftChild )
    {
        m_leftChild->get_AllNodes_const( v ) ;
        m_rightChild->get_AllNodes_const( v ) ;
    }

}



CAPTree::APTree_Pt CAPTree::findBtmNodeOfData( arma::mat& x , size_t& row_ind )
{

    // v is variable to split, c is raw value
    // not index in matrix<double>, so compare x[v] with c directly
    if( m_leftChild == 0 )
        return this ;

    if( x( row_ind , m_varIndex2Split ) <= m_rawValue2Split )
    {
        return m_leftChild->findBtmNodeOfData( x , row_ind ) ;
    }
    else
    {
        return m_rightChild->findBtmNodeOfData( x , row_ind ) ;
    }

}

void CAPTree::toNull( )
{

    size_t ts = treeSize( ) ;

    //loop invariant: ts>=1
    while( ts > 1 )
    { 
        //if false ts=1
        vec_APTree_Pt nv ;
        getNoGrandChildsNodes( nv ) ;

        for( size_t i = 0 ; i < nv.size( ) ; i++ )
        {
            delete nv[ i ]->m_leftChild ;
            delete nv[ i ]->m_rightChild ;
            nv[ i ]->m_leftChild = 0 ;
            nv[ i ]->m_rightChild = 0 ;
        }

        ts = treeSize( ) ; //make invariant true

    }

    m_varIndex2Split = 0 ;
    m_rawValue2Split = 0 ;

    m_parentNode = 0 ;
    m_leftChild = 0 ;
    m_rightChild = 0 ;

}


//copy tree tree o to tree n
void CAPTree::copyTree( APTree_Pt n , APTree_cnstPt o )
//assume n has no children (so we don't have to kill them)
//recursion down
// create a new copy of tree in NEW memory space
{

    if( n->m_leftChild )
    {
        std::cout << "cp:error node has children\n" ;
        return ;
    }

    n->m_varIndex2Split = o->m_varIndex2Split ;
    n->m_rawValue2Split = o->m_rawValue2Split ;
    n->theta = o->theta ;

    if( o->m_leftChild )
    { 
        //if o has children
        n->m_leftChild = new CAPTree ;
        ( n->m_leftChild )->m_parentNode = n ;
        copyTree( n->m_leftChild , o->m_leftChild ) ;

        n->m_rightChild = new CAPTree ;
        ( n->m_rightChild )->m_parentNode = n ;
        copyTree( n->m_rightChild , o->m_rightChild ) ;
    }

}

void CAPTree::copy_only_root( APTree_Pt o )
//assume n has no children (so we don't have to kill them)
//NOT LIKE cp() function
//this function pointer new root to the OLD structure
{

    this->m_varIndex2Split = o->m_varIndex2Split ;
    this->m_rawValue2Split = o->m_rawValue2Split ;
    this->theta = o->theta ;

    if( o->m_leftChild )
    {
        // keep the following structure, rather than create a new tree in memory
        this->m_leftChild = o->m_leftChild ;
        this->m_rightChild = o->m_rightChild ;
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
CAPTree& CAPTree::operator=( const CAPTree& rhs )
{
    if( &rhs != this )
    {
        toNull( ) ;       //kill left hand side (this)
        copyTree( this , &rhs ) ; //copy right hand side to left hand side
    }
    return *this ;
}

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


std::istream& operator>>( std::istream& is , CAPTree& AP_tree )
{

    size_t tid , pid ;                        //tid: id of current node, pid: parent's id
    std::map<size_t , CAPTree::APTree_Pt> pts ; //pointers to nodes indexed by node id
    size_t numOfNodes ;                              //number of nodes

    AP_tree.toNull( ) ; // obliterate old tree (if there)

    //read number of nodes----------
    is >> numOfNodes ;
    if( !is )
    {
        return is ;
    }

    // The idea is to dump string to a lot of node_info structure first, then link them as a tree, by nid

    //read in vector of node information----------
    std::vector<node_info> nv( numOfNodes ) ;
    for( size_t i = 0 ; i != numOfNodes ; i++ )
    {
        is >> nv[ i ].id >> nv[ i ].v >> nv[ i ].c >> nv[ i ].theta[ 0 ] ; // Only works on first theta for now, fix latex if needed
        if( !is )
        {
            return is ;
        }
    }

    //first node has to be the top one
    pts[ 1 ] = &AP_tree ; //be careful! this is not the first pts, it is pointer of id 1.
    AP_tree.set_varIndex2Split( nv[ 0 ].v ) ;
    AP_tree.set_rawValue2Split( nv[ 0 ].c ) ;
    AP_tree.set_theta( nv[ 0 ].theta ) ;
    AP_tree.m_parentNode = 0 ;

    //now loop through the rest of the nodes knowing parent is already there.
    for( size_t i = 1 ; i != nv.size( ) ; i++ )
    {
        CAPTree::APTree_Pt np = new CAPTree ;
        np->m_varIndex2Split = nv[ i ].v ;
        np->m_rawValue2Split = nv[ i ].c ;
        np->theta = nv[ i ].theta ;
        tid = nv[ i ].id ;
        pts[ tid ] = np ;
        pid = tid / 2 ;

        if( tid % 2 == 0 )
        { //left child has even id
            pts[ pid ]->m_leftChild = np ;
        }
        else
        {
            pts[ pid ]->m_rightChild = np ;
        }

        np->m_parentNode = pts[ pid ] ;

    }

    return is ;

}

void CAPTree::split_Xorder( arma::umat& Xorder_left , 
                            arma::umat& Xorder_right , 
                            arma::umat& Xorder , 
                            size_t split_point , 
                            size_t split_var , 
                            State& state )
{

    size_t num_obs = Xorder.n_rows ;

    double cutvalue = state.split_candidates[ split_point ] ;

    size_t left_index ;
    size_t right_index ;
    for( size_t i = 0 ; i < state.numOfCharitisc ; i++ )
    {
        left_index = 0 ;
        right_index = 0 ;

        // loop over variables
        for( size_t j = 0 ; j < num_obs ; j++ )
        {

            // loop over observations
            if( ( *state.X )( Xorder( j , i ) , split_var ) <= cutvalue )
            {
                // left side
                Xorder_left( left_index , i ) = Xorder( j , i ) ;
                left_index++ ;

            }
            else
            {
                // right side
                Xorder_right( right_index , i ) = Xorder( j , i ) ;
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
        j = this->theta ;
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
            this->theta = temp ;
        }
        else
        {
            this->theta[ 0 ] = temp[ 0 ] ;
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

  
void CAPTree::grow( bool& break_flag , 
                    CAPTreeModel& model , 
                    State& state , 
                    size_t& iter , 
                    std::vector<double>& criterion_values )
{

    std::vector<CAPTree*> bottom_nodes_vec ;
    std::vector<bool> node_splitability ;

    size_t split_node ;
    size_t split_var ;
    size_t split_point ;
    bool splitable = true ;

    // grow a tree by iteration instead of recursion
    // first, find all leaves
    bottom_nodes_vec.resize( 0 ) ;
    this->get_vecOfBtmNodes( bottom_nodes_vec ) ;

    // second, check splitability, 1 for splitable, 0 for terminated
    node_splitability.resize( bottom_nodes_vec.size( ) ) ;
    model.check_node_splitability( state , bottom_nodes_vec , node_splitability ) ;

    if( sum( node_splitability ) )
    {
        
        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion( state , 
                                   bottom_nodes_vec , 
                                   node_splitability , 
                                   split_node , 
                                   split_var , 
                                   split_point , 
                                   splitable , 
                                   criterion_values ) ;
        // split the selected node

        if( splitable )
        {
            bottom_nodes_vec[ split_node ]->set_iter( iter ) ;
            model.split_node( state , bottom_nodes_vec[ split_node ] , split_var , split_point ) ;
        }
        else
        {
            cout << "break of no good candidate" << endl ;
            break_flag = true ;
        }
    }
    else
    {
        cout << "break of no node splitable" << endl ;
        break_flag = true ;
    }

    return ;
}

void CAPTree::grow_APTree_TS( bool& break_flag , CAPTreeModel& model , State& state )
{

    std::vector<CAPTree*> bottom_nodes_vec ;
    std::vector<bool> node_splitability ;

    size_t split_node ;
    size_t split_var ;
    size_t split_point ;
    bool splitable = true ;

    // grow a tree by iteration instead of recursion
    // first, find all leaves
    bottom_nodes_vec.resize( 0 ) ;
    this->get_vecOfBtmNodes( bottom_nodes_vec ) ;

    // second, check splitability, 1 for splitable, 0 for terminated
    node_splitability.resize( bottom_nodes_vec.size( ) ) ;
    model.check_node_splitability( state , 
                                   bottom_nodes_vec , 
                                   node_splitability ) ;

    if( sum( node_splitability ) )
    {
        
        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion_APTree_TS(  state , 
                                              bottom_nodes_vec , 
                                              node_splitability , 
                                              split_node , 
                                              split_var , 
                                              split_point , 
                                              splitable ) ;
        // split the selected node
        if( splitable )
        {
            model.split_node_APTree_TS( state , 
                                        bottom_nodes_vec[ split_node ] , 
                                        split_var , 
                                        split_point ) ;
        }
        else
        {
            cout << "break of no good candidate" << endl ;
            break_flag = true ;
        }
        
    }
    else
    {
        cout << "break of no node splitable" << endl ;
        break_flag = true ;
    }

    return ;
}
