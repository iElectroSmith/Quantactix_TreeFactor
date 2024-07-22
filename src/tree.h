#ifndef GUARD_tree_h
#define GUARD_tree_h

#include "common.h"
#include "state.h"
#include "model.h"
#include "json.h"

using json = nlohmann::json ;

class CTree
{

public:
    // define types
    typedef CTree* tree_p ;
    typedef const CTree* tree_cp ;
    typedef std::vector<tree_p> npv ;
    typedef std::vector<tree_cp> cnpv ;

    //leaf parameters and sufficient statistics
    std::vector<double> theta ;
    arma::umat* Xorder ;
    arma::mat beta ;

    // constructors
    CTree( ) : 
        theta( 1 , 0.0 ) , 
        Xorder( 0 ) , 
        num_inNode( 0 ) , 
        node_ID( 1 ) , 
        var_Index( 0 ) , 
        c_index( 0 ) , 
        rawValue( 0.0 ) , 
        depth( 0 ) , 
        parentNode( 0 ) , 
        leftChild( 0 ) , 
        rightChild( 0 ) , 
        beta( 0 ) { }

    CTree( size_t dim_theta ) : 
        theta( dim_theta , 0.0 ) , 
        Xorder( 0 ) , 
        num_inNode( 0 ) , 
        node_ID( 1 ) , 
        var_Index( 0 ) , 
        c_index( 0 ) , 
        rawValue( 0.0 ) , 
        depth( 0 ) , 
        parentNode( 0 ) , 
        leftChild( 0 ) , 
        rightChild( 0 ) { ( this->beta ).resize( dim_theta , 1 ) ; }

    CTree( size_t dim_theta , arma::umat* Xordermat ) : 
        theta( dim_theta , 0.0 ) , 
        Xorder( Xordermat ) , 
        num_inNode( 0 ) , 
        node_ID( 1 ) , 
        var_Index( 0 ) , 
        c_index( 0 ) , 
        rawValue( 0.0 ) , 
        depth( 0 ) , 
        parentNode( 0 ) , 
        leftChild( 0 ) , 
        rightChild( 0 ) { }

    CTree( size_t dim_theta , size_t depth , size_t N , size_t ID , tree_p p , arma::umat* Xordermat ) : 
        theta( dim_theta , 0.0 ) , 
        Xorder( Xordermat ) , 
        num_inNode( N ) , 
        node_ID( ID ) , 
        var_Index( 0 ) , 
        c_index( 0 ) , 
        rawValue( 0.0 ) , 
        depth( depth ) , 
        parentNode( p ) , 
        leftChild( 0 ) , 
        rightChild( 0 ) 
    { ( this->beta ).resize( dim_theta , 1 ) ; }

    // functions
    void settheta( std::vector<double>& theta ) { this->theta = theta ; }
    void setv( size_t v ) { this->var_Index = v ; }
    void setc( double c ) { this->rawValue = c ; }
    void setc_index( size_t c_index ) { this->c_index = c_index ; }
    void setN( size_t N ) { this->num_inNode = N ; }
    void setID( size_t ID ) { this->node_ID = ID ; }
    void setdepth( size_t depth ) { this->depth = depth ; }

    size_t getv( ) const { return var_Index ; }
    size_t getdepth( ) const { return depth ; }
    std::vector<double> gettheta( ) { return theta ; }
    double gettheta( size_t ind ) const { return theta[ ind ] ; }
    size_t getthetasize( ) const { return theta.size( ) ; }
    size_t getc_index( ) const { return c_index ; }
    double getc( ) const { return rawValue ; }
    size_t getN( ) const { return num_inNode ; }
    tree_p getParent( ) { return parentNode ; }
    tree_p getRightChild( ) { return rightChild ; }
    tree_p getLeftChild( ) { return leftChild ; }
    size_t getID( ) { return node_ID ; }

    CTree& operator=( const CTree& rhs ) ;

    // tree operation functions
    void toNull( ) ;                              // delete the tree
    tree_p getptr( size_t nid ) ;                  // get node pointer from node ID, 0 if not there ;
    void pr( bool pc = true ) ;                    // to screen, pc is "print child"
    size_t treesize( ) ;                          // return number of nodes in the tree
    size_t nnogs( ) ;                             // number of nog (no grandchildren nodes)
    size_t nbots( ) ;                             // number of leaf nodes
    void getbots( npv& v ) ;                       // return a vector of bottom nodes
    void getnogs( npv& v ) ;                       // get nog nodes
    void getnodes( npv& v ) ;                      // get ALL nodes
    void getnodes( cnpv& v ) const ;               // get ALL nodes (const)
    tree_p gettop( ) ;                            // get pointer to the top node (root node) of the tree
    tree_p bn( arma::mat& x , size_t& row_index ) ; // search tree, find bottom node of the data
    size_t nid( ) const ;                         // nid of the node
    char ntype( ) ;                               // node type, t:top, b:bot, n:no grandchildren, i:interior (t can be b) ;
    bool isnog( ) ;
    void cp( tree_p n , tree_cp o ) ;  // copy tree from o to n
    void copy_only_root( tree_p o ) ; // copy tree, point new root to old structure
    friend std::istream& operator>>( std::istream& , CTree& ) ;

    // growing functions
    void grow( State& state , CModel& model , arma::umat& Xorder ) ;
    void split_Xorder( arma::umat& Xorder_left , arma::umat& Xorder_right , arma::umat& Xorder , size_t split_point , size_t split_var , State& state , CModel& model ) ;
    void predict( arma::mat X , arma::vec months , arma::vec& output ) ;

    // input and output to json
    json to_json( ) ;
    void from_json( json& j3 , size_t dim_theta ) ;

private:
    size_t num_inNode ;       // number of training data observation in this node
    size_t node_ID ;      // ID of the node
    size_t var_Index ;       // index of the variable to split
    size_t c_index ; // index of the value to split (index in the Xorder matrix)
    double rawValue ;       // raw value to split
    size_t depth ;   // depth of the tree

    tree_p parentNode ; // pointer to the parent node
    tree_p leftChild  ; // pointer to left child
    tree_p rightChild ; // pointer to right child

} ;


// io functions
std::istream& operator>>( std::istream& , CTree& ) ;
std::ostream& operator<<( std::ostream& , const CTree& ) ;

#endif


