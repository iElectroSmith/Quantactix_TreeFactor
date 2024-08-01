
library(TreeFactor)
library(rpart)
library(ranger)


#定义了一个名为tf_residual的函数，用于计算残差。
#这个函数的主要目的是在使用树模型因子的回归分析中，计算实际值与预测值之间的差异。
#fit：一个包含树模型因子等信息的对象。
#Y：响应变量（目标变量）。
#Z：回归变量矩阵。
#H：额外的回归变量矩阵。
#months：一个索引向量，用于选择fit$ft中的月份。
#no_H：一个布尔值，决定是否使用H矩阵。
tf_residual = function(fit,Y,Z,H,months,no_H)
{

  # Tree Factor Models
  regressor = Z   #将Z赋值给regressor，初始化回归变量矩阵
  #对于Z的每一列，将其乘以fit$ft中的对应月份因子。
  #这一部分代码的目的是对Z中的每个回归变量进行缩放处理。
  for( j in 1:dim(Z)[2] )
  {
      
    #从fit对象中提取名为ft的元素，这可能是一个数值向量或列表  
    regressor[ , j ] = Z[ , j ] * fit$ft[ months + 1 ]
    
  }

  #如果no_H为假（即no_H为FALSE），则将H矩阵合并到regressor中，作为额外的回归变量。
  if(!no_H)
  {
    regressor = cbind(regressor, H)
  }

  # print(fit$R2*100)
  x <- as.matrix(regressor)
  y <- Y
  
  #使用矩阵运算计算回归系数b_tf。计算的是线性回归模型中回归系数的最小二乘估计。
  # t(x) 表示矩阵 𝑋的转置。
  # %*% 表示矩阵乘法。
  # solve 被用来计算矩阵的逆
  b_tf = solve( t(x)%*%x )%*%t(x)%*%y 
  
  #这里 x 是自变量矩阵，b_tf 是之前通过最小二乘法估计得到的回归系数向量。
  #矩阵乘法 x %*% b_tf 计算的是预测值向量, 即每个样本的预测值
  #x %*% b_tf 的结果是一个矩阵，但通常它会是一个 n×1 的矩阵（即一列）。
  #[,1] 表示提取该矩阵的第一列。这样做是为了将结果转换为一个向量，
  #因为 R 中即使是一个列向量，进行矩阵乘法的结果仍然会是矩阵形式。
  # haty 包含了根据线性回归模型对每个样本进行预测的所有预测值。
  haty <- (x%*%b_tf)[,1]
  print(b_tf)
  
  #返回残差，即实际值Y减去预测值haty。
  return( Y - haty )

}

###### parameters #####

start = 1
split = 80
end   = 100

case='demo' 
max_depth=4   # 决策树的最大深度。限制树的深度可以防止过拟合。
min_leaf_size = 10 #决策树叶节点的最小样本数。控制树的复杂度和防止过拟合。
max_depth_boosting = 3 #提升模型（如梯度提升树）的最大深度。
num_iter = 1000  # 迭代次数，通常用于提升方法中的迭代轮数
num_cutpoints = 4 #每个特征的切分点数目。可能用于离散化连续特征
equal_weight = TRUE #是否对样本使用相等的权重
no_H = TRUE  #可能表示某个特定的约束或配置，具体取决于模型实现。
abs_normalize = TRUE  #是否对特征进行绝对值归一化。
weighted_loss = FALSE #是否使用加权损失函数。
stop_no_gain = FALSE  #是否在没有收益时停止训练。
nu = 1 #学习率或步长，通常用于梯度提升方法中

# this tiny regularization ensures the matrix inversion
# penalty for the sigma (sigma + lambda I)^{-1} * mu
lambda = 1e-4
eta = 1


##### load data #####

#load("../../data/simu_data.rda")
#load("E:/GitHub/Quantactix_TreeFactor/data/simu_data.rda")
load("D:/PKU_work/Quantactix_TreeFactor/data/simu_data.rda")

#print(names(da))打印数据框da的列名，以便查看其结构。
print( names( da ) )

data <- da # 将da数据框赋值给data，即data是da的副本。
data['lag_me'] = 1 #在data中添加一列名为lag_me的列，所有行的值都设置为1。

#这行代码从data中选择特定的列并创建一个新的数据框tmp。所选的列包括：
#'id'：可能是标识符列。
#'date'：日期列。
#'xret'：可能是某种回报率列。
#'lag_me'：前面添加的列，值为1。
#'c1'至'c5'：五个列，可能代表某些类别或特征。
#'m1'和'm2'：两个列，可能代表某些度量或指标。
#'mkt'：可能代表市场相关的一个列。
tmp = data[, c('id', 'date','xret','lag_me', 
              # 5
              'c1', 'c2', 'c3', 'c4', 'c5', # 9
              'm1','m2', # 11
              'mkt' # 12
              ) ]
data = tmp  #将tmp赋值回data，即data现在只包含所选的列。
rm( tmp ) #删除tmp变量，释放内存
rm( da )

# chars

all_chars <- names(data)[c(5:9)]  #提取data数据框中第5到第9列的列名，
top5chars <- c(1:5) #生成一个向量，包含从1到5的数字。
instruments = all_chars[top5chars]  #提取all_chars向量中第1到第5个元素
splitting_chars <- all_chars

first_split_var  = c(1:5)-1  #生成一个向量，包含从0到4的数字，
second_split_var = c(1:5)-1  #生成一个向量，包含从0到4的数字，

first_split_var_boosting = c(1:5)-1 #生成一个向量，包含从0到4的数字
second_split_var_boosting = c(1:5)-1

##### train-test split #####
#将原始数据框 data 按照 date 列拆分成两个子数据框 data1 和 data2，
#分别包含 start 到 split 和 split 到 end 之间的记录。

#data[,c('date')]: 提取数据框 data 中 date 列的所有值
#逻辑表达式，筛选 date 列中的值，使其在 start 和 split 之间（包括 start 和 split）
#data 使用逻辑表达式作为索引，筛选出符合条件的行
#data1 包含 data 中所有 date 列的值在 start 和 split 之间（包括 start 和 split）的行
data1 <- data[ (data[ , c('date')]>=start ) & ( data[ , c('date')]<=split) ,  ]

#data2 包含 data 中所有 date 列的值在 split 和 end 之间（包括 end 但不包括 split）的行
data2 <- data[ (data[ , c('date')]>split ) & ( data[ , c('date')]<=end ) , ]


# rm(data)

###### train data for all boosting steps #####
#段代码从数据框 data1 中提取和转换所需的特征变量、响应变量、
#时间和股票ID变量、协变量、组合权重和损失权重，
#并计算出唯一月份和股票的数量，为提升步骤的训练准备数据。

#从数据框 data1 中提取 splitting_chars 中指定的列，作为训练数据的特征变量 X_train
X_train = data1[,splitting_chars] 
#从数据框 data1 中提取 xret 列，作为训练数据的响应变量 R_train。
R_train = data1[,c("xret")]

# 将 data1 中的 date 列转换为因子（factor）
#将因子转换为数值型变量。结果是每个日期对应一个唯一的整数值
months_train = as.numeric(as.factor(data1[,c("date")]))
#将所有整数值减1，使月份从0开始。
months_train = months_train - 1 # start from 0

#将 data1 中的 id 列转换为因子。
#因子转换为数值型变量。结果是每个股票ID对应一个唯一的整数值
#将所有整数值减1，使股票ID从0开始
stocks_train = as.numeric(as.factor(data1[,c("id")])) - 1

#从数据框 data1 中提取 instruments 中指定的列，作为训练数据的协变量 Z_train。
Z_train = data1[, instruments]
# 在 Z_train 的第一列添加常数项1，用于线性模型的截距项
Z_train = cbind(1, Z_train)

#数据框 data1 中提取 lag_me 列，作为训练数据的组合权重 portfolio_weight_train
portfolio_weight_train = data1[,c("lag_me")]

#从数据框 data1 中提取 lag_me 列，作为训练数据的损失权重 loss_weight_train。
loss_weight_train = data1[,c("lag_me")]

#提取 months_train 向量中的唯一值
#计算唯一值的数量，即月份的数量 num_months。
num_months = length(unique(months_train))

#提取 stocks_train 向量中的唯一值。
#计算唯一值的数量，即股票的数量 num_stocks
num_stocks = length(unique(stocks_train))


###### train data 1 #####

# the first H is the mkt 表示第一个变量 H 是市场变量（mkt）。

#数据框 data1 中提取 xret 列，作为第一阶段训练数据的响应变量 Y_train1。
Y_train1 = data1[,c("xret")]
#从数据框 data1 中提取 mkt 列，作为第一阶段训练数据的市场变量 H_train1。
H_train1 = data1[,c("mkt")]

#对市场变量 H_train1 和协变量 Z_train 进行逐元素相乘。
#这里，H_train1 是一个向量，Z_train 是一个矩阵。
#相乘的结果是一个矩阵，每个元素是 H_train1 中对应元素
#和 Z_train 中对应元素的乘积。
#这个操作通常用在模型中，用于引入交互效应或构建新特征。
H_train1 = H_train1 * Z_train

# train 1 
t = proc.time()

#R_train：从 data1 中提取的 xret 列，作为训练数据的响应变量。
#Y_train1：从 data1 中提取的 xret 列，作为第一阶段训练数据的响应变量。
#X_train：从 data1 中提取的 splitting_chars 列，作为特征变量。
#Z_train：包含协变量的矩阵，其中第一列是常数1，其余列是从 instruments 提取的变量。
#H_train1：市场变量 mkt 与 Z_train 逐元素相乘的结果矩阵。
#portfolio_weight_train：从 data1 中提取的 lag_me 列，作为组合权重。
#loss_weight_train：从 data1 中提取的 lag_me 列，作为损失权重。
#stocks_train：将 data1 中 id 列转换为因子并减去1的结果，用于表示股票的索引。
#months_train：将 data1 中 date 列转换为因子并减去1的结果，用于表示月份的索引。
#first_split_var：定义第一个分裂变量的索引向量。
#second_split_var：定义第二个分裂变量的索引向量。
#num_stocks：唯一股票的数量。
#num_months：唯一月份的数量。
#min_leaf_size：树的最小叶节点大小。
#max_depth：树的最大深度。
#num_iter：迭代次数。
#num_cutpoints：用于分裂的切点数量。
#lambda：正则化参数。
#eta：学习率。
#equal_weight：布尔值，是否对每个样本赋予相同权重。
#no_H：布尔值，是否在模型中使用 H 矩阵。
#abs_normalize：布尔值，是否对数据进行绝对值归一化。
#weighted_loss：布尔值，是否使用加权损失。
#stop_no_gain：布尔值，是否在无增益时停止训练。


# 创建参数列表
params <- list(
    R_train = R_train,
    Y_train1 = Y_train1,
    X_train = X_train,
    Z_train = Z_train,
    H_train1 = H_train1,
    portfolio_weight_train = portfolio_weight_train,
    loss_weight_train = loss_weight_train,
    stocks_train = stocks_train,
    months_train = months_train,
    first_split_var = first_split_var,
    second_split_var = second_split_var,
    num_stocks = num_stocks,
    num_months = num_months,
    min_leaf_size = min_leaf_size,
    max_depth = max_depth,
    num_iter = num_iter,
    num_cutpoints = num_cutpoints,
    lambda = lambda,
    eta = eta,
    equal_weight = equal_weight,
    no_H = no_H,
    abs_normalize = abs_normalize,
    weighted_loss = weighted_loss,
    stop_no_gain = stop_no_gain
)

# 将参数列表写入文件
saveRDS(params, "./params.rds")


dir.create("params", showWarnings = FALSE)
write_vec <- function(vec, filename) {
    filepath <- file.path("params", filename)
    write.table(vec, file = filepath, row.names = FALSE, col.names = FALSE)
}

write_mat <- function(mat, filename) {
    filepath <- file.path("params", filename)
    write.table(mat, file = filepath, row.names = FALSE, col.names = FALSE)
}

# 写入向量和矩阵到文件
write_vec(R_train, "R_train.txt")
write_vec(Y_train1, "Y_train1.txt")
write_mat(X_train, "X_train.txt")
write_mat(Z_train, "Z_train.txt")
write_mat(H_train1, "H_train1.txt")
write_vec(portfolio_weight_train, "portfolio_weight_train.txt")
write_vec(loss_weight_train, "loss_weight_train.txt")
write_vec(stocks_train, "stocks_train.txt")
write_vec(months_train, "months_train.txt")
write_vec(unique(months_train), "unique_months_train.txt")
write_vec(first_split_var, "first_split_var.txt")
write_vec(second_split_var, "second_split_var.txt")

# 写入单一数值到文件
write.table(num_stocks, file.path("params", "num_stocks.txt"), row.names = FALSE, col.names = FALSE)
write.table(num_months, file.path("params", "num_months.txt"), row.names = FALSE, col.names = FALSE)
write.table(min_leaf_size, file.path("params", "min_leaf_size.txt"), row.names = FALSE, col.names = FALSE)
write.table(max_depth, file.path("params", "max_depth.txt"), row.names = FALSE, col.names = FALSE)
write.table(num_iter, file.path("params", "num_iter.txt"), row.names = FALSE, col.names = FALSE)
write.table(num_cutpoints, file.path("params", "num_cutpoints.txt"), row.names = FALSE, col.names = FALSE)
write.table(eta, file.path("params", "eta.txt"), row.names = FALSE, col.names = FALSE)
write.table(equal_weight, file.path("params", "equal_weight.txt"), row.names = FALSE, col.names = FALSE)
write.table(no_H, file.path("params", "no_H.txt"), row.names = FALSE, col.names = FALSE)
write.table(abs_normalize, file.path("params", "abs_normalize.txt"), row.names = FALSE, col.names = FALSE)
write.table(weighted_loss, file.path("params", "weighted_loss.txt"), row.names = FALSE, col.names = FALSE)
write.table(stop_no_gain, file.path("params", "stop_no_gain.txt"), row.names = FALSE, col.names = FALSE)
write.table(lambda_mean, file.path("params", "lambda_mean.txt"), row.names = FALSE, col.names = FALSE)
write.table(lambda_cov, file.path("params", "lambda_cov.txt"), row.names = FALSE, col.names = FALSE)


fit1 = TreeFactor_APTree( R_train, 
                          Y_train1, 
                          X_train, 
                          Z_train, 
                          H_train1, 
                          portfolio_weight_train, 
                          loss_weight_train, 
                          stocks_train, 
                          months_train, 
                          first_split_var, 
                          second_split_var, 
                          num_stocks, 
                          num_months, 
                          min_leaf_size, 
                          max_depth, 
                          num_iter, 
                          num_cutpoints, 
                          lambda, 
                          eta, 
                          equal_weight, 
                          no_H, 
                          abs_normalize, 
                          weighted_loss, 
                          stop_no_gain )
                          
t = proc.time() - t
print(t)

# in sample check
#使用模型 fit1 进行预测，
#传入训练数据的特征变量 X_train、响应变量 R_train、月份索引 months_train 
#和组合权重 portfolio_weight_train。
insPred1 = predict(fit1, X_train, R_train, months_train, portfolio_weight_train)
#insPred1$ft: 预测结果中的 ft 分量
#fit1$ft: 原模型中的 ft 分量
#计算预测结果与原模型结果之间的平方和差异。这是一个简单的误差衡量指标。
sum((insPred1$ft - fit1$ft)^2)
#打印模型 fit1 的 R² 值，用于衡量模型的拟合优度
print(fit1$R2)

# residual 1
#计算模型 fit1 的残差。传入的参数包括响应变量 Y_train1、协变量 Z_train、
#市场变量 H_train1、月份索引 months_train 和布尔值 no_H。
res1 = tf_residual(fit1,Y_train1,Z_train,H_train1,months_train,no_H)

# beta
#fit1$ft[months_train + 1]: 根据月份索引 months_train 获取 fit1 的 ft 分量，
#并加1以确保索引从1开始。
#将 fit1 的 ft 分量与协变量 Z_train 逐元素相乘。
#将结果转换为矩阵，存储在变量 x 中
x <- as.matrix(fit1$ft[months_train+1]*Z_train)
#将响应变量转换为矩阵，存储在变量 y 中
y <- as.matrix(Y_train1)
#计算 beta 系数的估计值 beta_bf1，这是普通最小二乘法（OLS）估计的一部分。
beta_bf1 = solve(t(x)%*%x)%*%t(x)%*%y
# print(beta_bf1)

###### train data 2 #####
# the first H is the mkt
Y_train2 = res1
no_H2 = TRUE

# train
t = proc.time()
fit2 = TreeFactor_APTree(R_train, Y_train2, X_train, Z_train, H_train1, portfolio_weight_train, 
                        loss_weight_train, stocks_train, months_train, first_split_var_boosting, second_split_var_boosting, num_stocks, 
                        num_months, min_leaf_size, max_depth_boosting, num_iter, num_cutpoints, lambda, eta, equal_weight, 
                        no_H2, abs_normalize, weighted_loss, stop_no_gain)
t = proc.time() - t
print(t)
print(fit2$R2)

# in sample check
insPred2 = predict(fit2, X_train, R_train, months_train, portfolio_weight_train)
sum((insPred2$ft - fit2$ft)^2)

# residual
res2 = tf_residual(fit2, Y_train2 ,Z_train,H_train1,months_train,no_H2)

# beta
x <- as.matrix(fit2$ft[months_train+1]*Z_train)
y <- as.matrix(Y_train2)
beta_bf2 = solve(t(x)%*%x)%*%t(x)%*%y
# print(beta_bf2)

############# Train Period tf2 #############

###### test data #####

X_test = data2[,splitting_chars]
R_test = data2[,c("xret")]
months_test = as.numeric(as.factor(data2[,c("date")]))
months_test = months_test - 1 # start from 0
stocks_test = as.numeric(as.factor(data2[,c("id")])) - 1
Z_test = data2[,instruments]
Z_test = cbind(1, Z_test)
H_test = data2[,c("mkt")]
H_test = H_test * Z_test
portfolio_weight_test = data2[,c("lag_me")]
loss_weight_test = data2[,c("lag_me")]

rm(data2)

############# Test Period Factors #############

pred1 = predict(fit1, X_test, R_test, months_test, portfolio_weight_test)
pred2 = predict(fit2, X_test, R_test, months_test, portfolio_weight_test)

tf_test <- cbind(pred1$ft, pred2$ft)
avg_tf_test <- colMeans(tf_test)

print("####### 1/N ########")
ewsdf = rowMeans(pred1$portfolio)
wsdf = pred1$portfolio %*% fit1$leaf_weight
print("1/N sdf sr: \n")
print(mean(ewsdf)/sd(ewsdf)*sqrt(length(ewsdf)))
print("weighted sdf sr: \n")
print(mean(wsdf)/sd(wsdf)*sqrt(length(wsdf)))
print("####### 1/N ########")


# ############# START OUTPUT #############

# tf_train <- cbind(fit1$ft, fit2$ft)
# # output factors
# write.csv(data.frame(tf_train), paste0(case,"_bf_train",".csv"))
# write.csv(data.frame(tf_test), paste0(case,"_bf_test",".csv"))

# write.csv(data.frame(insPred1$leaf_index), paste0(case,"_leaf_index",".csv"))
# write.csv(data.frame(months_train), paste0(case,"_months_train",".csv"))
# write.csv(data.frame(R_train), paste0(case,"_R_train",".csv"))
# write.csv(data.frame(portfolio_weight_train), paste0(case,"_weight_train",".csv"))

# write.csv(data.frame(fit1$leaf_weight), paste0(case,"_leaf_weight",".csv"))

# ### output test,for plot mve ###

# write.csv(data.frame(pred1$leaf_index), paste0(case,"_leaf_index_test",".csv"))
# write.csv(data.frame(months_test), paste0(case,"_months_test",".csv"))

# write.csv(data.frame(R_test), paste0(case,"_R_test",".csv"))
# write.csv(data.frame(portfolio_weight_test), paste0(case,"_weight_test",".csv"))

# ### output basis portfolio returns

# write.csv(data.frame(fit1$portfolio), paste0(case,"_portfolio_fit1",".csv"))
# write.csv(data.frame(pred1$portfolio), paste0(case,"_portfolio_pred1",".csv"))

# write.csv(data.frame(fit2$portfolio), paste0(case,"_portfolio_fit2",".csv"))
# write.csv(data.frame(pred2$portfolio), paste0(case,"_portfolio_pred2",".csv"))

# ############# END OUTPUT #############