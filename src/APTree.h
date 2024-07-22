#ifndef GUARD_APTree_h
#define GUARD_APTree_h

#include "common.h"
#include "state.h"
#include "model.h"
#include "json.h"

using json = nlohmann::json ;

class State ;
class APTreeModel ;

class CAPTree
{

public:
    // define types
    typedef CAPTree* APTree_Pt ;
    typedef const CAPTree* APTree_cnstPt ;
    typedef std::vector<APTree_Pt> vec_APTree_Pt ;
    typedef std::vector<APTree_cnstPt> vec_APTree_cnstPt ;

    //leaf parameters and sufficient statistics
    std::vector<double> theta ;
    arma::umat* Xorder ;

    // constructors
    CAPTree( ) : theta( 1 , 0.0 ) , 
                 Xorder( 0 ) , 
                 m_numData_inNode( 0 ) , 
                 m_nodeID( 1 ) , 
                 m_varIndex2Split( 0 ) , 
                 m_valueIndex2Split( 0 ) , 
                 m_rawValue2Split( 0.0 ) , 
                 m_treeDepth( 0 ) , 
                 m_parentNode( 0 ) , 
                 m_leftChild( 0 ) , 
                 m_rightChild( 0 ) , 
                 iter( 0 ) { }

    CAPTree( size_t dim_theta ) : theta( dim_theta , 0.0 ) , 
                                  Xorder( 0 ) , 
                                  m_numData_inNode( 0 ) , 
                                  m_nodeID( 1 ) , 
                                  m_varIndex2Split( 0 ) , 
                                  m_valueIndex2Split( 0 ) , 
                                  m_rawValue2Split( 0.0 ) , 
                                  m_treeDepth( 0 ) , 
                                  m_parentNode( 0 ) , 
                                  m_leftChild( 0 ) , 
                                  m_rightChild( 0 ) , 
                                  iter( 0 ) { }
    
    CAPTree( size_t dim_theta , arma::umat* Xordermat ) : theta( dim_theta , 0.0 ) , 
                                                          Xorder( Xordermat ) , 
                                                          m_numData_inNode( 0 ) , 
                                                          m_nodeID( 1 ) , 
                                                          m_varIndex2Split( 0 ) , 
                                                          m_valueIndex2Split( 0 ) , 
                                                          m_rawValue2Split( 0.0 ) , 
                                                          m_treeDepth( 0 ) , 
                                                          m_parentNode( 0 ) , 
                                                          m_leftChild( 0 ) , 
                                                          m_rightChild( 0 ) , 
                                                          iter( 0 ) { }
    CAPTree( size_t dim_theta ,   //state.num_months
             size_t depth ,       // 1 , 
             size_t N ,           //state.num_obs_all  
             size_t ID ,          //1 ,
             APTree_Pt p ,         //0 ,
             arma::umat* Xordermat ) : theta( dim_theta , 0.0 ) , 
                                       Xorder( Xordermat ) , 
                                       m_numData_inNode( N ) , 
                                       m_nodeID( ID ) , 
                                       m_varIndex2Split( 0 ) , 
                                       m_valueIndex2Split( 0 ) , 
                                       m_rawValue2Split( 0.0 ) , 
                                       m_treeDepth( depth ) , 
                                       m_parentNode( p ) , 
                                       m_leftChild( 0 ) , 
                                       m_rightChild( 0 ) , 
                                       iter( 0 ) { }


    // functions
    void set_theta( std::vector<double>& theta ) { this->theta = theta ; }
    void set_varIndex2Split( size_t v ) { this->m_varIndex2Split = v ; }
    void set_rawValue2Split( double c ) { this->m_rawValue2Split = c ; }
    void set_valueIndex2Split( size_t c_index ) { this->m_valueIndex2Split = c_index ; }
    void set_numData_inNode( size_t N ) { this->m_numData_inNode = N ; }
    void set_nodeID( size_t ID ) { this->m_nodeID = ID ; }
    void set_treeDepth( size_t depth ) { this->m_treeDepth = depth ; }
    void set_leftChild( APTree_Pt l ) { this->m_leftChild = l ; }
    void set_rightChild( APTree_Pt r ) { this->m_rightChild = r ; }
    void set_parentNode( APTree_Pt p ) { this->m_parentNode = p ; }
    void set_iter( size_t iter ) { this->iter = iter ; }

    size_t get_varIndex2Split( ) const { return m_varIndex2Split ; }
    size_t get_treeDepth( ) const { return m_treeDepth ; }
    std::vector<double> get_theta( ) { return theta ; }
    double get_theta( size_t ind ) const { return theta[ ind ] ; }
    size_t get_thetaSize( ) const { return theta.size( ) ; }
    size_t get_valueIndex2Split( ) const { return m_valueIndex2Split ; }
    double get_rawValue2Split( ) const { return m_rawValue2Split ; }
    size_t get_numData_inNode( ) const { return m_numData_inNode ; }
    APTree_Pt get_parentNode( ) { return m_parentNode ; }
    APTree_Pt get_rightChild( ) { return m_rightChild ; }
    APTree_Pt get_leftChild( ) { return m_leftChild ; }
    size_t get_nodeID( ) { return m_nodeID ; }
    size_t get_iter( ) const { return iter ; }

    CAPTree& operator=( const CAPTree& rhs ) ;

    // tree operation functions
    void toNull( ) ;                                // delete the tree
    APTree_Pt get_pt2ParentByNodeID( size_t nid ) ;              // get node pointer from node ID, 0 if not there ;
    void printScreen( bool pc = true ) ;            // to screen, pc is "print child"
    size_t treeSize( ) ;                            // return number of nodes in the tree
    size_t numNoGrandChildsNodes( ) ;                               // number of nog (no grandchildren nodes)
    size_t numLeafNodes( ) ;                               // number of leaf nodes
    void get_vecOfBtmNodes( vec_APTree_Pt& v ) ;                        // return a vector of bottom nodes

    void getNoGrandChildsNodes( vec_APTree_Pt& v ) ;           // get nog nodes
    void get_AllNodes( vec_APTree_Pt& v ) ;                        // get ALL nodes
    void get_AllNodes_const( vec_APTree_cnstPt& v ) const ;              // get ALL nodes (const)
    APTree_Pt get_pt2TopNode( ) ;                              // get pointer to the top node (root node) of the tree
    APTree_Pt findBtmNodeOfData( arma::mat& x , size_t& row_index ) ; // search tree, find bottom node of the data
    
    size_t nid( ) const ;                           // nid of the node
    char nodeType( ) ;                                 // node type, t:top, b:bot, n:no grandchildren, i:interior (t can be b) ;
    bool isNoGrandChildren( ) ;
    void copyTree( APTree_Pt n , APTree_cnstPt o ) ; // copy tree from o to n
    void copy_only_root( APTree_Pt o ) ;  // copy tree, point new root to old structure
    friend std::istream& operator>>( std::istream& , CAPTree& ) ;

    void split_Xorder( arma::umat& Xorder_left , 
                       arma::umat& Xorder_right , 
                       arma::umat& Xorder , 
                       size_t split_point , 
                       size_t split_var , 
                       State& state ) ;
    void predict( arma::mat X , arma::vec months , arma::vec& output ) ;

    void grow( bool& break_flag , 
               CAPTreeModel& model , 
               State& state , 
               size_t& iter , 
               std::vector<double>& criterion_values ) ;

    void grow_APTree_TS( bool& break_flag , CAPTreeModel& model , State& state ) ;
 
    // input and output to json
    json to_json( ) ;
    void from_json( json& j3 , size_t dim_theta ) ;

private:
    size_t m_numData_inNode ;       // number of training data observation in this node
    size_t m_nodeID ;               // ID of the node
    size_t m_varIndex2Split ;       // index of the variable to split
    size_t m_valueIndex2Split ; // index of the value to split (index in the Xorder matrix)
    double m_rawValue2Split ;       // raw value to split
    size_t m_treeDepth ;   // depth of the tree
    size_t iter ;    // which iteration this node splits

    APTree_Pt m_parentNode ; // pointer to the parent node
    APTree_Pt m_leftChild ; // pointer to left child
    APTree_Pt m_rightChild ; // pointer to right child

} ;

// io functions
std::istream& operator>>( std::istream& , CAPTree& ) ;
std::ostream& operator<<( std::ostream& , const CAPTree& ) ;

#endif