
# The data-generating-process is simple.
# We assume the returns are governed by the market factor $mkt_t$ and two characteristics $c1_{i,t}, c2_{i,t}$.
# $ r_{i,t} = \beta_i mkt_t + b1 * c1_{i,t} + b2_i*c2_{i,t} + \epsilon_{i,t} $.
# The characteristics follow Uniform[-1,1] with a normal(0,0.1) fluctuation.
# The market beta is a function of characteristics $c2, c3$ and macroeconomic indicator $m2$.
# The characteristics are independent.
# The return and volatility of the market factor is conditional on a macroeconomic indicator $m1$.
# Other noise variables include $c5, m2$

# parameters

set.seed( 20220215 )
# set.seed(89)

N <- 1000  # number of asset
T <- 100   # number of time period
b1 <- 1/100
b2 <- 2/100

# functions

## function that generates characteristics
#这个函数char生成一个T行N列的矩阵c_matrix，其每个元素是由以下部分组成：
#一个均匀分布在[-1, 1]之间的随机数c（在每行中重复出现）。
#一个均值为0、标准差为0.1的正态分布随机数（独立生成的噪声）。
char <- function( N ,T )
{

  #均匀分布随机数
  #runif(n,min=0,max=1) n表示生成的随机数数量，min表示均匀分布的下限，max表示均匀分布的上限
  c <- runif( n = N , min = -1 , max = 1 )    # 1*1000  # 生成N个均匀分布在[-1, 1]之间的随机数

  x1<- rep( c , T )     # 1*1000*100 ，赋值100次 # 重复向量c，重复T次，得到长度为N*T的向量
  x2<- matrix( x1 ,  nrow = N )   #1000* 100 ，C行变列，变成100列，每一列都一样 ## 将向量x1重组为一个N行的矩阵，每列都相同
  x3<- t( x2  )   # 转置 ， 变成100*1000 ，每一行都一样   # 转置矩阵x2，得到一个T行N列的矩阵，每行都相同
  
  x4<- rnorm( n = N*T , mean=0 , sd=0.1 )  # 1000*100 ，rnorm 随机正态分布，mean 是平均数， sd 是标准差 ，然后随机抽样 或者取值 n 次，
  x5<- matrix( x4 , nrow = T )  #变成 100* 1000矩阵 # 将这些随机数重组为一个T行N列的矩阵
 
  
  #叠加
  c_matrix <- t( matrix( rep( c , T ) ,  nrow = N ) )  + matrix( rnorm( n = N*T , mean=0 , sd=0.1 ) , nrow = T )

  return( c_matrix )

}

## function that generate market beta
#这个函数通过参数c2和c3的线性组合以及参数m1的条件判断来计算市场beta。
#如果m1大于0，计算出的beta值将减少1；否则，beta值将增加1。
#该函数生成的beta值以1为基础，然后加上0.5 * (c2 + c3)的结果，
#再加上一个根据m1值的条件调整。
betaF <- function( c1 , c2 , c3 , c4 , c5 , m1 , m2 )
{

    #(m1 > 0)返回一个布尔值，如果m1大于0，则返回TRUE（对应数值1），否则返回FALSE（对应数值0）。
    #(m1 > 0) * 2将布尔值转换为数值并乘以2。
    #(1 - (m1 > 0) * 2)最终结果为：如果m1 > 0，则结果为1 - 2 = -1；如果m1 <= 0，则结果为1 - 0 = 1。
    beta <- 1 + 0.5 * ( c2 + c3 + ( 1 - ( m1 > 0 ) * 2 ) )

    return( beta )

}



## function that reshape matrix into a stacked list
# 它按照列的顺序依次将矩阵的每一列元素连接起来，形成一个长向量
m2l <- function( m )
{

  d <- dim( m )  # 获取矩阵的维度
  t <- d[1]  # nrow
  k <- d[2]  # ncol
  l <- rep( 0 , t*k )  # 初始化长度为 t*k 的零向量
  
  for ( i in c( 1:k ) )   #遍历每一列
  {
  
    left <- t*(i-1) + 1  # 计算当前列在向量中的起始位置
    right <- t*i         # 计算当前列在向量中的结束位置
  
    # print(left, right)
    l[ left:right ] <- m[ , i ] # 将当前列的元素赋值给向量对应的位置
  
  }
  
  return( l )

}

## function that reshape matrix into a stacked list

# simulation

## simulate the characteristics of assets

c1 <- char(N,T)
c2 <- char(N,T)
c3 <- char(N,T)
c4 <- char(N,T)
c5 <- char(N,T)


# simulate macroeconomic indicator

m1 <- rnorm( n=T, mean=0, sd=1/100) #从正态分布中生成随机数
m2 <- rnorm( n=T, mean=0, sd=1/100)


## simualte market factor
#一个长度为T的随机数向量，其中每个随机数都来自一个均值为1/100，标准差为2/100的正态分布
mkt <- rnorm( n=T, mean=1/100 , sd=2/100 )
#将向量mkt转换为一个T行1列的矩阵。因为nrow = T，
#所以生成的矩阵的每一行对应向量mkt的一个元素。
#注意，这里并不是转置操作，只是将向量转换为一个列矩阵。
#如果想要明确地转置矩阵，可以使用t函数。
mkt <- matrix( mkt , nrow = T )  #转置了


## simulate the market beta of assets

beta <- betaF( c1 , c2 , c3 , c4 , c5 , m1 , m2 )


## simulate asset returns

r <- matrix( 0 , T , N ) # 100 * 1000 全0矩阵 
dim_xy <-dim( r )  #，r矩阵有100行和10列。通过调用dim(r)， 确认矩阵的维度是否符合预期。
print(dim_xy)
cx <- c(1:N)  # [1 : 1000]
for ( i in c(1:N) )
{
    
    x0 <- i

    x11 <- beta[ , i ]
    x12 <- beta[ , i ] *mkt  #表示矩阵beta的第i列与向量mkt的按元素乘积。
    
    x13 <- c1[ , i ]
    x14 <- c1[ , i ]*b1 #表示矩阵c1的第i列与向量b1的按元素乘积。
    
    x15 <- c2[ , i ]
    x16 <- c2[ , i ]*b2 #表示矩阵c2的第i列与向量b2的按元素乘积

    x17<- rnorm(n=T, mean=0, sd=10/100) #生成一个正态分布的随机噪声向量，长度为T，均值为0，标准差为0.1。
    
    #表示矩阵r的第i列。
    r[,i] <- beta[ , i ]*mkt + c1[ , i ]*b1 + c2[ , i ]*b2 + rnorm(n=T, mean=0, sd=10/100)
}

# save
## stack data and save as rda

da <- data.frame(

      xret = m2l(r),  # 将矩阵 r 转换为向量
      
      #ceiling函数的作用是将小数部分舍去，并将整数部分加1，除非它本身是整数
      # 这生成了一个长度为N*T的向量，其中每个整数重复N次，从1到T
      # ceiling(1:(N*T)/N)生成了一个id列，它的作用是为每一组N个观测值分配一个ID，使得数据框中每N个观测值共享一个ID
      id   = ceiling(1:(N*T)/N), #生成一个从1到T重复N次的整数序列，
      
      date = rep(c(1:T),N), #创建一个日期变量，通过rep(c(1:T), N)生成，它的长度也是N*T，表示每个观测值的日期。
      mkt  = rep(mkt,N), #创建一个市场变量，通过rep(mkt, N)生成，它的长度为N*T
      m1   = rep(m1,N), #长度为N*T
      m2   = rep(m2,N), #长度为N*T
      
      c1   = m2l(c1),
      c2   = m2l(c2),
      c3   = m2l(c3),
      c4   = m2l(c4),
      c5   = m2l(c5)
      
)

f  = as.matrix( cbind( mkt ) )
xt = as.matrix( cbind( m1 , m2 ) )
save(da, f, xt, beta, file = "simu_data.rda")


