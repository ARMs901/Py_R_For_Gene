# 定义核心算法函数
CoreAlg <- function(X, y){
  # 设置迭代次数
  svn_itor <- 3
  
  # 定义一个内部函数，用于为给定的 i 返回相应的 nu 值
  res <- function(i){
    nus <- 0 # 为 nus 设置一个默认值
    if(i == 1) {nus <- 0.25}
    if(i == 2) {nus <- 0.5}
    if(i == 3) {nus <- 0.75}
    # 使用给定的 nu 值训练 SVM 模型
    model <- svm(X, y, type="nu-regression", kernel="linear", nu=nus, scale=F)
    return(model)
  }
  
  # 根据操作系统决定如何并行地运行 res 函数
  if(Sys.info()['sysname'] == 'Windows') {
    out <- mclapply(1:svn_itor, res, mc.cores=1) 
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  }
  
  # 初始化存储结果的向量
  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)
  t <- 1
  while(t <= svn_itor) {
    # 计算权重
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    
    # 计算 k
    u <- sweep(X, MARGIN=2, w, '*')
    k <- apply(u, 1, sum)
    
    # 计算均方根误差和相关性
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  w <- (q / sum(q))
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  return(list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r))
}

# 定义置换测试函数
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  while(itor <= perm){
    # 对 Y 进行随机排序
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    yr <- (yr - mean(yr)) / sd(yr)
    
    # 调用核心算法函数
    result <- CoreAlg(X, yr)
    mix_r <- result$mix_r
    
    if(itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }
    itor <- itor + 1
  }
  return(list("dist" = dist))
}

# 定义 CIBERSORT 函数
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  # 加载必要的包
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  # 读取输入文件
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)
  
  # 将 X 和 Y 转换为数据矩阵并按行名称排序
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  # 如果 Y 的最大值小于 50，则对 Y 进行 2 的对数变换
  if(max(Y) < 50) {Y <- 2^Y}
  
  # 如果 QN 为 TRUE，则对 Y 进行定量归一化
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  # 确保 X 和 Y 有相同的行名
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  # 对 X 进行标准化
  X <- (X - mean(X)) / sd(as.vector(X))
  
  # 如果指定了置换次数，则执行置换测试
  if(perm > 0) {
    nulldist <- sort(doPerm(perm, X, Y)$dist)
  }
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  while(itor <= mixtures){
    # 提取 Y 的当前列并进行标准化
    y <- Y[,itor]
    y <- (y - mean(y)) / sd(y)
    
    # 调用核心算法函数
    result <- CoreAlg(X, y)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    # 如果执行了置换测试，则计算 p 值
    if(perm > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }
    
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }
    itor <- itor + 1
  }
  
  # 将结果写入文件
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  return(obj)
}

# 设置工作目录
setwd("C:\\Users\\MOMO\\OneDrive\\桌面\\机器学习数据重做\\UC数据\\015.CIBERSORT")

# 调用 CIBERSORT 函数
# 请确保 inputFile 已被定义
inputFile = "normalize.txt"

outTab <- CIBERSORT("ref.txt", inputFile, perm=1000, QN=TRUE)

# 保留 p 值小于 0.05 的结果
outTab <- outTab[outTab[,"P-value"] < 0.05,]

# 删除最后三列并将结果写入文件
outTab <- as.matrix(outTab[, 1:(ncol(outTab) - 3)])
outTab <- rbind(colnames(outTab), outTab) # 使用列名代替未定义的 id
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)
