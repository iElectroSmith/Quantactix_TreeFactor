#include "APTree.h"
#include <chrono>
#include <ctime>

size_t CAPTree::nid( ) const
{

    if( !parentNode )
        return 1; //if you don't have a parent, you are the top

    if( this == parentNode->leftChild )
        return 2 * ( parentNode->nid( ) ); //if you are a left child
    else
        return 2 * ( parentNode->nid( ) ) + 1; //else you are a right child

}

CAPTree::APTree_Pt CAPTree::getParent( size_t nid )
{

    if( this->nid( ) == nid )
        return this; //found it

    if( leftChild == 0 )
        return 0; //no children, did not find it

    APTree_Pt lp = leftChild->getParent( nid );
    if( lp )
        return lp; //found on left

    APTree_Pt rp = rightChild->getParent( nid );
    if( rp )
        return rp; //found on right
    
    return 0;      //never found it

}

size_t CAPTree::treeSize( )
{
    if( leftChild == 0 )
        return 1; //if bottom node, tree size is 1
    else
        return ( 1 + leftChild->treeSize( ) + rightChild->treeSize( ) );
}

char CAPTree::ntype( )
{
    //t:top, b:bottom, n:no grandchildren, i:internal
    if( !parentNode )
        return 't';

    if( !leftChild )
        return 'b';

    if( !( leftChild->leftChild ) && !( rightChild->leftChild ) )
        return 'n';

    return 'i';

}



void CAPTree::printScreen( bool pc )
{

    size_t depth = this->depth ;
    size_t id = nid( );
    size_t pid;

    if( !parentNode )
        pid = 0; //parent of top node
    else
        pid = parentNode->nid( );

    std::string pad( 2 * depth , ' ' );
    std::string strSplit( ", " );

    if( pc && ( ntype( ) == 't' ) )
        std::cout << "tree size: " << treeSize( ) << std::endl;

    std::cout << pad << "(id,parent): " << id << strSplit << pid;
    std::cout << strSplit << "(v,c): " << v << strSplit << c;
    // std::cout << sp << "theta: " << theta;
    std::cout << strSplit << "type: " << ntype( );
    std::cout << strSplit << "depth: " << this->depth;
    std::cout << strSplit << "pointer: " << this << std::endl;

    if( pc )
    {

        if( leftChild )
        {
            leftChild->printScreen( pc );
            rightChild->printScreen( pc );
        }

    }

}



bool CAPTree::isnog( )
{

    bool isnog = true;

    if( leftChild )
    {
        if( leftChild->leftChild || rightChild->leftChild )
            isnog = false; //one of the children has children.
    }
    else
    {
        isnog = false; //no children
    }

    return isnog;

}


size_t CAPTree::numNoGrandChildsNodes( )
{

    if( !leftChild )
        return 0; //bottom node

    if( leftChild->leftChild || rightChild->leftChild )
    { //not a nog
        return ( leftChild->numNoGrandChildsNodes( ) + rightChild->numNoGrandChildsNodes( ) );
    }
    else
    { //is a nog
        return 1;
    }

}

size_t CAPTree::numLeafNodes( )
{

    if( leftChild == 0 )
    { //if a bottom node
        return 1;
    }
    else
    {
        return leftChild->numLeafNodes( ) + rightChild->numLeafNodes( );
    }

}

void CAPTree::getbots( vec_APTree_Pt& bv )
{

    if( leftChild )
    { //have children
        leftChild->getbots( bv );
        rightChild->getbots( bv );
    }
    else
    {
        bv.push_back( this );
    }

}

void CAPTree::getNoGrandChildsNodes( vec_APTree_Pt& nv )
{

    if( leftChild )
    { 
        //have children
        if( ( leftChild->leftChild ) || ( rightChild->leftChild ) )
        { 
            //have grandchildren
            if( leftChild->leftChild )
                leftChild->getNoGrandChildsNodes( nv );

            if( rightChild->leftChild )
                rightChild->getNoGrandChildsNodes( nv );
        }
        else
        {
            nv.push_back( this );
        }
    }

}

CAPTree::APTree_Pt CAPTree::gettop( )
{

    if( !parentNode )
    {
        return this;
    }
    else
    {
        return parentNode->gettop( );
    }

}

void CAPTree::getnodes( vec_APTree_Pt& v )
{

    v.push_back( this );

    if( leftChild )
    {
        leftChild->getnodes( v );
        rightChild->getnodes( v );
    }

}


void CAPTree::getnodes( vec_APTree_cnstPt& v ) const
{

    v.push_back( this );

    if( leftChild )
    {
        leftChild->getnodes( v );
        rightChild->getnodes( v );
    }

}



CAPTree::APTree_Pt CAPTree::bn( arma::mat& x , size_t& row_ind )
{

    // v is variable to split, c is raw value
    // not index in matrix<double>, so compare x[v] with c directly
    if( leftChild == 0 )
        return this;

    if( x( row_ind , v ) <= c )
    {
        return leftChild->bn( x , row_ind );
    }
    else
    {
        return rightChild->bn( x , row_ind );
    }

}

void CAPTree::toNull( )
{

    size_t ts = treeSize( );

    //loop invariant: ts>=1
    while( ts > 1 )
    { 
        //if false ts=1
        vec_APTree_Pt nv;
        getNoGrandChildsNodes( nv );
        for( size_t i = 0; i < nv.size( ); i++ )
        {
            delete nv[ i ]->leftChild;
            delete nv[ i ]->rightChild;
            nv[ i ]->leftChild = 0;
            nv[ i ]->rightChild = 0;
        }

        ts = treeSize( ); //make invariant true

    }

    v = 0;
    c = 0;

    parentNode = 0;
    leftChild = 0;
    rightChild = 0;

}


//copy tree tree o to tree n
void CAPTree::copyTree( APTree_Pt n , APTree_cnstPt o )
//assume n has no children (so we don't have to kill them)
//recursion down
// create a new copy of tree in NEW memory space
{

    if( n->leftChild )
    {
        std::cout << "cp:error node has children\n";
        return;
    }

    n->v = o->v;
    n->c = o->c;
    n->theta = o->theta;

    if( o->leftChild )
    { 
        //if o has children
        n->leftChild = new CAPTree;
        ( n->leftChild )->parentNode = n;
        copyTree( n->leftChild , o->leftChild );

        n->rightChild = new CAPTree;
        ( n->rightChild )->parentNode = n;
        copyTree( n->rightChild , o->rightChild );
    }

}

void CAPTree::copy_only_root( APTree_Pt o )
//assume n has no children (so we don't have to kill them)
//NOT LIKE cp() function
//this function pointer new root to the OLD structure
{

    this->v = o->v;
    this->c = o->c;
    this->theta = o->theta;

    if( o->leftChild )
    {
        // keep the following structure, rather than create a new tree in memory
        this->leftChild = o->leftChild;
        this->rightChild = o->rightChild;
        // also update pointers to parents
        this->leftChild->parentNode = this;
        this->rightChild->parentNode = this;
    }
    else
    {
        this->leftChild = 0;
        this->rightChild = 0;
    }

}



//--------------------------------------------------
//operators
CAPTree& CAPTree::operator=( const CAPTree& rhs )
{
    if( &rhs != this )
    {
        toNull( );       //kill left hand side (this)
        copyTree( this , &rhs ); //copy right hand side to left hand side
    }
    return *this;
}

//--------------------------------------------------
std::ostream& operator<<( std::ostream& os , const CAPTree& APTree )
{

    CAPTree::vec_APTree_cnstPt nds;

    APTree.getnodes( nds );
    os << nds.size( ) << std::endl;

    // size_t theta_length = nds[0]->getthetasize();
    // cout << "theta length is " << theta_length << endl;
    for( size_t i = 0; i < nds.size( ); i++ )
    {
        os << nds[ i ]->nid( ) << " ";
        os << nds[ i ]->getv( ) << " ";
        os << nds[ i ]->getc( ) << " ";
        os << nds[ i ]->getc_index( ) << " ";
        os << nds[ i ]->getiter( );
        // for (size_t j = 0; j < theta_length; j++)
        // {
        //     os << " " << nds[i]->gettheta(j);
        //     // os << " " << nds[i]->getRt(j);
        // }
        os << std::endl;
    
    }
    
    return os;
}


std::istream& operator>>( std::istream& is , CAPTree& t )
{

    size_t tid , pid;                        //tid: id of current node, pid: parent's id
    std::map<size_t , CAPTree::APTree_Pt> pts; //pointers to nodes indexed by node id
    size_t nn;                              //number of nodes

    t.toNull( ); // obliterate old tree (if there)

    //read number of nodes----------
    is >> nn;
    if( !is )
    {
        return is;
    }

    // The idea is to dump string to a lot of node_info structure first, then link them as a tree, by nid

    //read in vector of node information----------
    std::vector<node_info> nv( nn );
    for( size_t i = 0; i != nn; i++ )
    {
        is >> nv[ i ].id >> nv[ i ].v >> nv[ i ].c >> nv[ i ].theta[ 0 ]; // Only works on first theta for now, fix latex if needed
        if( !is )
        {
            return is;
        }
    }

    //first node has to be the top one
    pts[ 1 ] = &t; //be careful! this is not the first pts, it is pointer of id 1.
    t.setv( nv[ 0 ].v );
    t.setc( nv[ 0 ].c );
    t.settheta( nv[ 0 ].theta );
    t.parentNode = 0;

    //now loop through the rest of the nodes knowing parent is already there.
    for( size_t i = 1; i != nv.size( ); i++ )
    {
        CAPTree::APTree_Pt np = new CAPTree;
        np->v = nv[ i ].v;
        np->c = nv[ i ].c;
        np->theta = nv[ i ].theta;
        tid = nv[ i ].id;
        pts[ tid ] = np;
        pid = tid / 2;

        if( tid % 2 == 0 )
        { //left child has even id
            pts[ pid ]->leftChild = np;
        }
        else
        {
            pts[ pid ]->rightChild = np;
        }

        np->parentNode = pts[ pid ];

    }

    return is;

}

void CAPTree::split_Xorder( arma::umat& Xorder_left , 
                            arma::umat& Xorder_right , 
                            arma::umat& Xorder , 
                            size_t split_point , 
                            size_t split_var , 
                            State& state )
{

    size_t num_obs = Xorder.n_rows;

    double cutvalue = state.split_candidates[ split_point ];

    size_t left_index;
    size_t right_index;
    for( size_t i = 0; i < state.p; i++ )
    {
        left_index = 0;
        right_index = 0;

        // loop over variables
        for( size_t j = 0; j < num_obs; j++ )
        {

            // loop over observations
            if( ( *state.X )( Xorder( j , i ) , split_var ) <= cutvalue )
            {
                // left side
                Xorder_left( left_index , i ) = Xorder( j , i );
                left_index++;

            }
            else
            {
                // right side
                Xorder_right( right_index , i ) = Xorder( j , i );
                right_index++;

            }
        }

    }

    return;

}

json CAPTree::to_json( )
{

    json j;
    
    if( leftChild == 0 )
    {
        j = this->theta;
    }
    else
    {
        j[ "variable" ] = this->v;
        j[ "cutpoint" ] = this->c;
        j[ "cutpoint_index" ] = this->c_index;
        j[ "nodeid" ] = this->nid( );
        j[ "depth" ] = this->depth;
        j[ "left" ] = this->leftChild->to_json( );
        j[ "right" ] = this->rightChild->to_json( );
    }

    return j;

}


void CAPTree::from_json( json& j3 , size_t dim_theta )
{

    if( j3.is_array( ) )
    {

        // this is the leaf
        std::vector<double> temp;
        j3.get_to( temp );

        if( temp.size( ) > 1 )
        {
            this->theta = temp;
        }
        else
        {
            this->theta[ 0 ] = temp[ 0 ];
        }

    }
    else
    {

        // this is an intermediate node
        j3.at( "variable" ).get_to( this->v );
        j3.at( "cutpoint" ).get_to( this->c );
        j3.at( "cutpoint_index" ).get_to( this->c_index );
        j3.at( "depth" ).get_to( this->depth );

        CAPTree* lchild = new CAPTree( dim_theta );
        lchild->from_json( j3[ "left" ] , dim_theta );
        CAPTree* rchild = new CAPTree( dim_theta );
        rchild->from_json( j3[ "right" ] , dim_theta );

        lchild->parentNode = this;
        rchild->parentNode = this;
        this->leftChild = lchild;
        this->rightChild = rchild;
    }

}

  
void CAPTree::grow( bool& break_flag , 
                    CAPTreeModel& model , 
                    State& state , 
                    size_t& iter , 
                    std::vector<double>& criterion_values )
{

    std::vector<CAPTree*> bottom_nodes_vec;
    std::vector<bool> node_splitability;

    size_t split_node;
    size_t split_var;
    size_t split_point;
    bool splitable = true;

    // grow a tree by iteration instead of recursion
    // first, find all leaves
    bottom_nodes_vec.resize( 0 );
    this->getbots( bottom_nodes_vec );

    // second, check splitability, 1 for splitable, 0 for terminated
    node_splitability.resize( bottom_nodes_vec.size( ) );
    model.check_node_splitability( state , bottom_nodes_vec , node_splitability );

    if( sum( node_splitability ) )
    {
        
        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion( state , bottom_nodes_vec , node_splitability , 
                                   split_node , split_var , split_point , splitable , criterion_values );
        // split the selected node

        if( splitable )
        {
            bottom_nodes_vec[ split_node ]->setiter( iter );
            model.split_node( state , bottom_nodes_vec[ split_node ] , split_var , split_point );
        }
        else
        {
            cout << "break of no good candidate" << endl;
            break_flag = true;
        }
    }
    else
    {
        cout << "break of no node splitable" << endl;
        break_flag = true;
    }

    return;
}

void CAPTree::grow_APTree_TS( bool& break_flag , CAPTreeModel& model , State& state )
{
    std::vector<CAPTree*> bottom_nodes_vec;
    std::vector<bool> node_splitability;

    size_t split_node;
    size_t split_var;
    size_t split_point;
    bool splitable = true;

    // grow a tree by iteration instead of recursion
    // first, find all leaves
    bottom_nodes_vec.resize( 0 );
    this->getbots( bottom_nodes_vec );

    // second, check splitability, 1 for splitable, 0 for terminated
    node_splitability.resize( bottom_nodes_vec.size( ) );
    model.check_node_splitability( state , bottom_nodes_vec , node_splitability );

    if( sum( node_splitability ) )
    {
        
        // if there exist at least one node for split
        // third, loop  over those splitabiliable nodes, calculate split criterion, figure out split node, var and point
        model.calculate_criterion_APTree_TS(  state , bottom_nodes_vec , 
                                              node_splitability , split_node , 
                                              split_var , split_point , splitable );
        // split the selected node
        if( splitable )
        {
            model.split_node_APTree_TS( state , bottom_nodes_vec[ split_node ] , split_var , split_point );
        }
        else
        {
            cout << "break of no good candidate" << endl;
            break_flag = true;
        }
        
    }
    else
    {
        cout << "break of no node splitable" << endl;
        break_flag = true;
    }

    return;
}
