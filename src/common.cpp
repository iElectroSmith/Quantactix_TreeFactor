#include "common.h"

// overload to print vectors and vector<vector>

std::ostream& operator<<( std::ostream& out , const std::vector<double>& v )
{

    size_t last = v.size( ) - 1 ;

    for( size_t i = 0 ; i < v.size( ) ; ++i )
    {
        out << v[ i ] ;

        if( i != last )
            out << ", " ;

    }

    return out ;

}

std::ostream& operator<<( std::ostream& out , const std::vector<bool>& v )
{
    size_t last = v.size( ) - 1 ;

    for( size_t i = 0 ; i < v.size( ) ; ++i )
    {
        out << v[ i ] ;

        if( i != last )
            out << ", " ;
    }

    return out ;

}

std::ostream& operator<<( std::ostream& out , const std::vector<size_t>& v )
{

    size_t last = v.size( ) - 1 ;

    for( size_t i = 0 ; i < v.size( ) ; ++i )
    {
        out << v[ i ] ;

        if( i != last )
            out << ", " ;

    }

    return out ;

}

std::ostream& operator<<( std::ostream& out , const std::vector<std::vector<double>>& v )
{

    // size_t last = v.size() - 1 ;
    for( size_t i = 0 ; i < v.size( ) ; ++i )
    {
        out << v[ i ] << endl ;
    }

    return out ;

}

std::ostream& operator<<( std::ostream& out , const std::vector<std::vector<size_t>>& v )
{
    // size_t last = v.size() - 1 ;
    for( size_t i = 0 ; i < v.size( ) ; ++i )
    {
        out << v[ i ] << endl ;
    }

    return out ;

}

//函数 fastLm 接受两个参数 y 和 X，分别表示响应变量向量和预测变量矩阵。
//函数返回一个 double 类型的值，即残差平方和
/*
该函数 fastLm 主要用于计算 OLS 回归模型的残差平方和。
它通过求解线性方程组得到回归系数，计算残差，并最终返回残差的平方和。以下是该函数的步骤总结：

计算回归系数 coef。
计算残差 resid。
计算残差方差 sig2 和回归系数的标准误差 stderrest。
计算并返回残差的平方和 output。
*/
//fastLm: Fast Linear Model（快速线性模型计算）
double fastLm( const arma::vec& y , const arma::mat& X )
{


    DEBUG_PRINT_SPACE;
    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("  "  );

    // this function calculate sum of residual squares for OLS
    //获取样本数和预测变量的数量：n 表示样本数，k 表示预测变量的数量。
    size_t n = X.n_rows ;
    size_t k = X.n_cols ;

    //使用 Armadillo 的 solve 函数计算回归系数 coef，相当于 coef = (X^T * X)^(-1) * X^T * y。
	// 用于解线性方程组。它可以用来求解形如 AX=B 的方程，
	// 其中 A 是系数矩阵，X 是未知数向量（或矩阵），B 是已知数向量（或矩阵）
    arma::colvec coef = arma::solve( X , y ) ;

    //计算残差 resid，即实际值 y 与预测值 X * coef 之间的差。
    arma::colvec resid = y - X * coef ;

    //计算残差方差 sig2，即残差的平方和除以自由度 n - k。
    double sig2 = arma::as_scalar( arma::trans( resid ) * resid / ( n - k ) ) ;
    
    //计算回归系数的标准误差 stderrest。
    //arma::diagvec 函数提取协方差矩阵的对角线元素，arma::inv 函数计算矩阵的逆。
    arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv( arma::trans( X ) * X ) ) ) ;

    //计算残差的平方 temp，并使用 arma::accu 函数求和得到 output。
    arma::colvec temp = arma::pow( resid , 2 ) ;
    double output = arma::accu( temp ) ;

    return output ;

}

/*
该函数 fastLm_weighted 主要用于计算带权重的 OLS 回归模型的加权残差平方和。
通过求解线性方程组得到回归系数，计算残差，并对残差的平方与权重向量逐元素相乘，
最终返回加权残差平方和。以下是该函数的步骤总结：

计算回归系数 coef。
计算残差 resid。
计算残差方差 sig2 和回归系数的标准误差 stderrest。
计算加权残差的平方 temp。
返回加权残差平方和 output。
*/
//三个参数 y、X 和 weight，分别表示响应变量向量、预测变量矩阵和权重向量。
//函数返回一个 double 类型的值，即加权残差平方和
//fastLm_weighted: Fast Weighted Linear Model（快速加权线性模型计算）
double fastLm_weighted( const arma::vec& y , const arma::mat& X , const arma::vec& weight )
{

    DEBUG_PRINT_SPACE;
    DEBUG_PRINT_SPACE;
    DEBUG_PRINT("  "  );

    // this function calculate sum of residual squares for OLS
    //获取样本数和预测变量的数量：n 表示样本数，k 表示预测变量的数量。
    size_t n = X.n_rows ;
    size_t k = X.n_cols ;

    //使用 Armadillo 的 solve 函数计算回归系数 coef，相当于 coef = (X^T * X)^(-1) * X^T * y
    //计算残差 resid，即实际值 y 与预测值 X * coef 之间的差。
    arma::colvec coef = arma::solve( X , y ) ;
    arma::colvec resid = y - X * coef ;

    //计算残差方差 sig2，即残差的平方和除以自由度 n - k。
    double sig2 = arma::as_scalar( arma::trans( resid ) * resid / ( n - k ) ) ;
    //计算回归系数的标准误差 stderrest。
    //arma::diagvec 函数提取协方差矩阵的对角线元素，arma::inv 函数计算矩阵的逆。
    arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv( arma::trans( X ) * X ) ) ) ;

    //计算残差的平方 temp，并使用元素逐一相乘操作 % 与权重向量 weight 相乘。
    // 与前者的差别： 对残差的平方进行加权求和：arma::colvec temp = arma::pow(resid, 2) % weight，
    arma::colvec temp = arma::pow( resid , 2 ) % weight ;
    //计算并返回加权残差平方和：
    double output = arma::accu( temp ) ;

    return output ;

}

bool sum( std::vector<bool>& v )
{

    bool output = false ;
    for( size_t i = 0 ; i < v.size( ) ; i++ )
    {
        output = output + v[ i ] ;
    }

    return output ;

}

double log_normal_density( arma::vec& R , arma::mat& cov )
{
    double output = 0.0 ;

    return output ;
}

double soft_c( double a , double lambda )
{

    // soft threshold
    if( a > lambda )
    {
        return ( a - lambda ) ;
    }
    else if( a < lambda )
    {
        return ( a + lambda ) ;
    }
    else
    {
        return 0.0 ;
    }

}

double lasso_loss( const arma::mat& X , 
                    const arma::mat& Y , 
                    const arma::vec& beta , 
                    double lambda )
{

    size_t n = X.n_rows ;
    size_t p = X.n_cols ;
    double output = accu( square( Y - X * beta ) / ( 2 * n ) ) + lambda * accu( abs( beta ) ) ;

    return output ;
}

arma::vec lasso_fit_standardized( const arma::mat& X , 
                                    const arma::mat& Y , 
                                    double lambda ,
                                    const arma::vec& beta_ini , 
                                    double eps = 0.0001 )
{

    // solve Lasso by coordinate descent method
    // not the closed form solution
    size_t n = X.n_rows ;
    size_t p = X.n_cols ;
    arma::vec beta_last = beta_ini ;
    arma::vec beta_new = beta_ini ;

    double loss_diff = 100.00 ;

    arma::vec r = Y - X * beta_ini ;

    double loss_old ;

    while( loss_diff >= eps )
    {

        beta_last = beta_new ;

        loss_old = lasso_loss( X , Y , beta_last , lambda ) ;

        for( size_t i = 0 ; i < p ; i++ )
        {
            beta_new( i ) = soft_c( arma::as_scalar( beta_last( i ) + X.col( i ).t( ) * r / n ) , lambda ) ;

            r = r + X.col( i ) * ( beta_last( i ) - beta_new( i ) ) ;
        }

        double loss_new = lasso_loss( X , Y , beta_new , lambda ) ;

        loss_diff = loss_old - loss_new ;

    }

    return beta_new ;

}

void int_to_bin( size_t num , std::vector<size_t>& s )
{

    size_t p = s.size( ) ;

    // i has to be int type here, cannot be size_t
    // otherwise it overflow to large positive number, the loop cannot stop
    for( int i = ( p - 1 ) ; i >= 0 ; i-- )
    {

        if( num & 1 )
        {
            s.at( i ) = 1 ;
        }
        else
        {
            s.at( i ) = 0 ;
        }

        num = num >> 1 ;

    }

    return ;

}


#include <sstream>


void printMat(const arma::mat& m , size_t start_row  , size_t start_col , 
                                   size_t num_rows  , size_t num_cols   ) 
{

    DEBUG_PRINT_SPACE;
	
    std::ostringstream oss;
	
 // 设置固定的浮点数格式和精度
    oss << std::fixed << std::setprecision( 6 );
	
    // 打印矩阵的内存地址
    oss << "Matrix address: " << &m << "\n";
	
	
  // 检查行列范围是否有效
    if (start_row >= m.n_rows || start_col >= m.n_cols) 
	{
        oss << "Invalid start index.";
        return   ;
    }

	  // 如果 num_rows 或 num_cols 为 0 或超过矩阵范围，则打印到矩阵末尾
    if (num_rows == 0 || start_row + num_rows > m.n_rows) 
	{
        num_rows = m.n_rows - start_row;
    }
    if (num_cols == 0 || start_col + num_cols > m.n_cols) 
	{
        num_cols = m.n_cols - start_col;
    }

  // 打印指定范围的元素
    for (size_t i = start_row; i < start_row + num_rows; ++i) 
	{
        for (size_t j = start_col; j < start_col + num_cols; ++j) 
		{
            oss << std::setw(15) << m(i, j) << "\t";
        }
        oss << "\n";
    }
	
    std::string matStr = oss.str();
	
	std::cout << matStr << std::endl;	
	
}

void printMat(const arma::umat& m , size_t start_row  , size_t start_col , 
                                   size_t num_rows  , size_t num_cols   ) 							   
{
	
    DEBUG_PRINT_SPACE;

	 // 将 umat 转换为 mat
    arma::mat M = arma::conv_to<arma::mat>::from( m );
	
	printMat( M ,   start_row  ,   start_col ,  num_rows  ,   num_cols ) ;
	
}

void  printVec(const arma::vec& v , size_t start , size_t length  ) 
{
    
    DEBUG_PRINT_SPACE;

    std::ostringstream oss;
	
 // 设置固定的浮点数格式和精度
    oss << std::fixed << std::setprecision( 6 );
	

    // 打印向量的内存地址
    oss << "Vector address: " << &v << "\n";
	

    // 如果 length 为 0 或超过向量长度，则使用向量的剩余长度
    if (length == 0 || start + length > v.n_elem) {
        length = v.n_elem - start;
    }	
	
    // 检查起始索引是否有效
    if (start >= v.n_elem) {
        oss << "Invalid start index.";
        return  ;
    } 
	
	for (size_t i = start; i < start + length; ++i) 
	{	
        oss << std::setw(15)<<  v(i) << "\n";
	}
    
	std::string vecStr = oss.str();
	
	std::cout << vecStr << std::endl;	
	
}