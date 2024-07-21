#ifndef GUARD_APTree_h
#define GUARD_APTree_h

#include "common.h"
#include "state.h"
#include "model.h"
#include "json.h"

using json = nlohmann::json;

class State;
class APTreeModel;

class CAPTree
{

public:
    // define types
    typedef CAPTree* APTree_Pt;
    typedef const CAPTree* APTree_cnstPt;
    typedef std::vector<APTree_Pt> vec_APTree_Pt;
    typedef std::vector<APTree_cnstPt> vec_APTree_cnstPt;

    //leaf parameters and sufficient statistics
    std::vector<double> theta;
    arma::umat* Xorder;

    // constructors
    CAPTree( ) : theta( 1 , 0.0 ) , 
                 Xorder( 0 ) , 
                 N( 0 ) , 
                 ID( 1 ) , 
                 v( 0 ) , 
                 c_index( 0 ) , 
                 c( 0.0 ) , 
                 depth( 0 ) , 
                 parentNode( 0 ) , 
                 leftChild( 0 ) , 
                 rightChild( 0 ) , 
                 iter( 0 ) { }

    CAPTree( size_t dim_theta ) : theta( dim_theta , 0.0 ) , 
                                  Xorder( 0 ) , 
                                  N( 0 ) , 
                                  ID( 1 ) , 
                                  v( 0 ) , 
                                  c_index( 0 ) , 
                                  c( 0.0 ) , 
                                  depth( 0 ) , 
                                  parentNode( 0 ) , 
                                  leftChild( 0 ) , 
                                  rightChild( 0 ) , 
                                  iter( 0 ) { }
    
    CAPTree( size_t dim_theta , arma::umat* Xordermat ) : theta( dim_theta , 0.0 ) , 
                                                          Xorder( Xordermat ) , 
                                                          N( 0 ) , 
                                                          ID( 1 ) , 
                                                          v( 0 ) , 
                                                          c_index( 0 ) , 
                                                          c( 0.0 ) , 
                                                          depth( 0 ) , 
                                                          parentNode( 0 ) , 
                                                          leftChild( 0 ) , 
                                                          rightChild( 0 ) , 
                                                          iter( 0 ) { }
    CAPTree( size_t dim_theta ,   //state.num_months
             size_t depth ,       // 1 , 
             size_t N ,           //state.num_obs_all  
             size_t ID ,          //1 ,
             APTree_Pt p ,         //0 ,
             arma::umat* Xordermat ) : theta( dim_theta , 0.0 ) , 
                                       Xorder( Xordermat ) , 
                                       N( N ) , 
                                       ID( ID ) , 
                                       v( 0 ) , 
                                       c_index( 0 ) , 
                                       c( 0.0 ) , 
                                       depth( depth ) , 
                                       parentNode( p ) , 
                                       leftChild( 0 ) , 
                                       rightChild( 0 ) , 
                                       iter( 0 ) { }


    // functions
    void settheta( std::vector<double>& theta ) { this->theta = theta; }
    void setv( size_t v ) { this->v = v; }
    void setc( double c ) { this->c = c; }
    void setc_index( size_t c_index ) { this->c_index = c_index; }
    void setN( size_t N ) { this->N = N; }
    void setID( size_t ID ) { this->ID = ID; }
    void setdepth( size_t depth ) { this->depth = depth; }
    void setl( APTree_Pt l ) { this->leftChild = l; }
    void setr( APTree_Pt r ) { this->rightChild = r; }
    void setp( APTree_Pt p ) { this->parentNode = p; }
    void setiter( size_t iter ) { this->iter = iter; }

    size_t getv( ) const { return v; }
    size_t getdepth( ) const { return depth; }
    std::vector<double> gettheta( ) { return theta; }
    double gettheta( size_t ind ) const { return theta[ ind ]; }
    size_t getthetasize( ) const { return theta.size( ); }
    size_t getc_index( ) const { return c_index; }
    double getc( ) const { return c; }
    size_t getN( ) const { return N; }
    APTree_Pt getp( ) { return parentNode; }
    APTree_Pt getr( ) { return rightChild; }
    APTree_Pt getl( ) { return leftChild; }
    size_t getID( ) { return ID; }
    size_t getiter( ) const { return iter; }

    CAPTree& operator=( const CAPTree& rhs );

    // tree operation functions
    void toNull( );                                // delete the tree
    APTree_Pt getParent( size_t nid );              // get node pointer from node ID, 0 if not there;
    void printScreen( bool pc = true );            // to screen, pc is "print child"
    size_t treeSize( );                            // return number of nodes in the tree
    size_t numNoGrandChildsNodes( );                               // number of nog (no grandchildren nodes)
    size_t numLeafNodes( );                               // number of leaf nodes
    void getbots( vec_APTree_Pt& v );                        // return a vector of bottom nodes

    void getNoGrandChildsNodes( vec_APTree_Pt& v );                         // get nog nodes
    void getnodes( vec_APTree_Pt& v );                        // get ALL nodes
    void getnodes( vec_APTree_cnstPt& v ) const;                 // get ALL nodes (const)
    APTree_Pt gettop( );                              // get pointer to the top node (root node) of the tree
    APTree_Pt bn( arma::mat& x , size_t& row_index ); // search tree, find bottom node of the data
    
    size_t nid( ) const;                           // nid of the node
    char ntype( );                                 // node type, t:top, b:bot, n:no grandchildren, i:interior (t can be b);
    bool isnog( );
    void copyTree( APTree_Pt n , APTree_cnstPt o ); // copy tree from o to n
    void copy_only_root( APTree_Pt o );  // copy tree, point new root to old structure
    friend std::istream& operator>>( std::istream& , CAPTree& );

    void split_Xorder( arma::umat& Xorder_left , 
                       arma::umat& Xorder_right , 
                       arma::umat& Xorder , 
                       size_t split_point , 
                       size_t split_var , 
                       State& state );
    void predict( arma::mat X , arma::vec months , arma::vec& output );

    void grow( bool& break_flag , 
               CAPTreeModel& model , 
               State& state , 
               size_t& iter , 
               std::vector<double>& criterion_values );

    void grow_APTree_TS( bool& break_flag , CAPTreeModel& model , State& state );
 
    // input and output to json
    json to_json( );
    void from_json( json& j3 , size_t dim_theta );

private:
    size_t N;       // number of training data observation in this node
    size_t ID;      // ID of the node
    size_t v;       // index of the variable to split
    size_t c_index; // index of the value to split (index in the Xorder matrix)
    double c;       // raw value to split
    size_t depth;   // depth of the tree
    size_t iter;    // which iteration this node splits

    APTree_Pt parentNode; // pointer to the parent node
    APTree_Pt leftChild; // pointer to left child
    APTree_Pt rightChild; // pointer to right child

};

// io functions
std::istream& operator>>( std::istream& , CAPTree& );
std::ostream& operator<<( std::ostream& , const CAPTree& );

#endif