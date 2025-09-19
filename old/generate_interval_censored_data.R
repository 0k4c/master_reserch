source("C:/Users/81801/Documents/R/code/senior_thesis/hazard_functions.R")
source("C:/Users/81801/Documents/R/code/senior_thesis/window_mean.R")
# source("C:/Users/81801/Documents/R/code/senior_thesis/analyze_wmst1.R")
#分布関数の種類を追加
generate_interval_censored_data <- function(n=100,K=5,p_dropout="Medium",p_exact=0.75,d=1){
	#初期化
	t_end <- 1
	epsilon <- 1e-6

	#エラーと警告を表示
	if (!is.element(n,c(100,200,400))) stop("only (100,200,400) is permited for n")
	if (!is.element(K,c(3,5,10,20))) stop("only (3,5,10,20) is permited for K")
	if (!is.element(p_dropout,c("None","Low","Medium","High"))) stop('only ("None","Low","Medium","High") is permited for p_dropout')
	if (!is.element(p_exact,c(0,0.25,0.5,0.75,1.0))) stop("only (0,0.25,0.5,1.0) is permited for p_exact")
	if (!is.element(d,1:11)) stop("only 1:11 is permited for d")

	#区間幅Lを求める
	L <- t_end/(K+1)

	#グリッドの長さを決定
	grid_length <- K+3

	#n個のグリッドを初期化
	grids <- array(numeric(2*(K+3)*n), dim=c(K+3,2,n))

	#n個の観測開始時点を一様分布から生成
	tau_0s <- runif(n, 0, L)

	#観測時点行列を作成
	mat1 <- array(numeric((K+3)*n), dim=c(n,K+3))
	mat1[,3:(K+2)] <- L
	mat2 <- array(numeric((K+3)*n), dim=c(n,K+3))
	mat2[,2:(K+2)] <- tau_0s
	mat2[,K+3] <- Inf
	#観測時点配列をグリッドに代入
	grids[,1,] <- apply(mat1,1,cumsum)+aperm(mat2,perm=c(2,1))

	#観測欠落配列を生成
	if (p_dropout=="None") p_drops <- numeric(grid_length)
	else if (p_dropout=="Low") p_drops <- c(numeric(2),rep(0.1,K-1),0.2, 0)
	else if (p_dropout=="Medium") p_drops <- c(numeric(2), rep(0.2, K-1), 0.4, 0)
	else if (p_dropout=="High") p_drops <- c(numeric(2), rep(0.3, K-1), 0.6, 0)
	else print("Error")
	#観測欠落配列をグリッドに代入(1:欠落,0:正常)
	grids[,2,] <- rbinom(grid_length*n, size=1, prob=p_drops)

	#真の観測データを生成(nに拡張)
	if (d==1) T_vals <- rweibull(n,scale=1,shape=1)
	else if (d==2) T_vals <- rweibull(n,scale=1,shape=0.5)
	else if (d==3) T_vals <- rweibull(n,scale=1,shape=2)
	else if (d==4) T_vals <- rweibull(n,scale=0.5,shape=1)
	else if (d==5) T_vals <- rweibull(n,scale=2,shape=1)
	else if (d==6) T_vals <- rf6(n)
	else if (d==7) T_vals <- rf7(n)
	else if (d==8) T_vals <- rf8(n)
	else if (d==9) T_vals <- rf9(n)
	else if (d==10) T_vals <- rf10(n)
	else if (d==11) T_vals <- rf11(n)
	else print("Error")

	#正確な観測かどうかをベルヌーイ分布からサンプリング(1:正確,0:区間打ち切り)
	xis <- rbinom(n,size=1,prob=p_exact)
  # print(xis)
	#区間打ち切りデータ生成
	interval_datas <- array(numeric(2*n),dim=c(n,2))

	#正確な観測データを代入
	nums1 <- which(xis==1)
	
	interval_datas[nums1,] <- c(T_vals[xis==1],T_vals[xis==1])
	interval_datas[nums1,2] <- interval_datas[nums1,1]+epsilon

	#T_val>1の区間打ち切りデータを代入
	nums2 <- c(which(xis==0&T_vals>t_end))
  # print(nums2)
	#このdrop_maskは欠落をFALSEとする
	drop_mask <- as.logical(grids[,2,nums2])
	#ここ
	# print(which(drop_mask|grids[,1,nums2]==Inf))
	# print(grids[,1,nums2])
	# print(grids[,1,nums2,drop=F])
	interval_datas[nums2,1] <- apply(replace(grids[,1,nums2,drop=F],c(which(drop_mask|grids[,1,nums2]==Inf)),0),2,max)
	interval_datas[nums2,2] <- Inf

	#T_val<=1の区間打ち切りデータを代入
	nums3 <- c(which(xis==0&T_vals<=t_end))
	l_logic_mask <- aperm(aperm(grids[,1,nums3], perm=c(2,1))-T_vals[nums3],perm=c(2,1))<0
	r_logic_mask <- !l_logic_mask
	#このdrop_maskは欠落をTUREとする
	drop_mask <- as.logical(grids[,2,nums3])
	interval_datas[nums3,1] <- apply(replace(grids[,1,nums3],c(which(r_logic_mask | drop_mask)),0),2,max)
	interval_datas[nums3,2] <- apply(replace(grids[,1,nums3],c(which(l_logic_mask | drop_mask)),Inf),2,min)

	#1列目が区間の左端、2列目が区間の右端、行はi番目の被験者
	return(invisible(interval_datas))
}

# #中点代入して生存時間データを生成(t,delta)delta==0が右側打ち切り(中点代入法)
# midpoint_assignment <- function(matrix){
# 	#結果を格納するresultを初期化
# 	result <- array(numeric(nrow(matrix)*2),dim=c(nrow(matrix),2))
# 	#区間打ち切りデータを中点代入処理(delta=1:イベント発生)
# 	result[,1][matrix[,2]<=1] <- (matrix[,1][matrix[,2]<=1]+matrix[,2][matrix[,2]<=1])/2
# 	result[,2][matrix[,2]<=1] <- 1
# 	#観察中途打ち切り(delta=0:打ち切り処理)
# 	result[,1][matrix[,2]>1&matrix[,1]<=1] <- matrix[,1][matrix[,2]>1&matrix[,1]<=1]
# 	result[,2][matrix[,2]>1&matrix[,1]<=1] <- 0
# 	#観察完了打ち切り(delta=0:打ち切り処理)
# 	result[,1][matrix[,1]>1] <- 1
# 	result[,2][matrix[,1]>1] <- 0
# 	#結果をまとめて返す
# 	result_df <- data.frame(time=result[,1],cens=result[,2])
# 	return(result_df)
# }
# 
# #右端法
# rightpoint_assignment <- function(matrix){
# 	#結果を格納するresultを初期化
# 	result <- array(numeric(nrow(matrix)*2),dim=c(nrow(matrix),2))
# 	#区間打ち切りデータを右点代入処理(delta=1:イベント発生)
# 	result[,1][matrix[,2]<=1] <- matrix[,2][matrix[,2]<=1]
# 	result[,2][matrix[,2]<=1] <- 1
# 	#観察中途打ち切り(delta=0:打ち切り処理)
# 	result[,1][matrix[,2]>1&matrix[,1]<=1] <- matrix[,1][matrix[,2]>1&matrix[,1]<=1]
# 	result[,2][matrix[,2]>1&matrix[,1]<=1] <- 0
# 	#観察完了打ち切り(delta=0:打ち切り処理)
# 	result[,1][matrix[,1]>1] <- 1
# 	result[,2][matrix[,1]>1] <- 0
# 	#結果をまとめて返す
# 	result_df <- data.frame(time=result[,1],cens=result[,2])
# 	return(result_df)
# }

# 台形の面積を計算
calc_trapezoid <- function(base, top, height){
  area <- (base + top) * height / 2
  return(area)
}

#ターンブル法
turnbull <- function(mtx){
  # View(df)
  mtx_L <- cbind(mtx[,1], 0)
  mtx_R <- cbind(mtx[,2], 1)
  # 1列目は左端or右端の数値、2列目は左端か右端かの識別子0:左、1:右
  new_mtx <- rbind(mtx_L, mtx_R)
  rm(mtx_L,mtx_R)
  # View(new_df)
  # 1列目の値で行列を並び替え
  order_mtx <- new_mtx[order(new_mtx[,1]),]
  # 3列目に一つ後のid、4列目に一つ前のid
  order_mtx2 <- cbind(order_mtx, append(order_mtx[-1,2],NA),append(NA,order_mtx[-nrow(order_mtx),2]))
  # 左-右の連続する組を抽出
  equal_set <- order_mtx2[(order_mtx2[,2]==0&order_mtx2[,3]==1)|(order_mtx2[,2]==1&order_mtx2[,4]==0),]
  # View(equal_set)
  # 同等集合を作成。1列目が左端、2列目は右端
  equal_set_wide <- matrix(equal_set[,1],nrow(equal_set)/2,2,byrow=T)#cbind(equal_set[equal_set[,2]==0,1], equal_set[equal_set[,2]==1,1])
  # 3列目に右端と左端の差を追加
  equal_set_wide <- cbind(equal_set_wide,equal_set_wide[,2] - equal_set_wide[,1])
  
  # for (j in 1:length(equal_set_wide$L)){
  #   assign(paste0("s", j), function(x){x})
  #   equal_set_wide$Survdiff[j] <- paste0("s", j)
  # }
  # View(equal_set_wide)
  # 同等集合の要素数
  func_length <- nrow(equal_set_wide)
  # 被験者数
  n <- nrow(mtx)
  
  # for (i in 1:length(df$L)){
  #   for (j in 1:length(equal_set_wide$L)){
  #     if (equal_set_wide$L[j] >= df$L[i] & equal_set_wide$R[j] <= df$R[i]){
  #       equal_set_wide$sumidx[j] <- 1
  #     }else equal_set_wide$sumidx[j] <- 0
  #   }
  # }
  
  make_likelihood <- function(surv_length, func_length) {
    likelihood <- function(params){
      # partial_likelihood = 1
      # for (i in 1:surv_length){
      #   sum_val <- 0
      #   for (j in 1:func_length){
      #     if (equal_set_wide$L[j] >= df$L[i] & equal_set_wide$R[j] <= df$R[i]){
      #       equal_set_wide$sumidx[j] <- 1
      #     }else equal_set_wide$sumidx[j] <- 0
      #   }
      #   for (j in 1:func_length){
      #     sum_val <- sum_val + params[j]*equal_set_wide$sumidx[j]
      #   }
      #   partial_likelihood <- partial_likelihood * sum_val
      #   #print(partial_likelihood)
      # 
      # }
      
      # 最初にすべて0の行列を作っとく
      mtx_ZeroOne <- matrix(rep(0, surv_length*func_length), func_length, surv_length)
      for (i in 1:func_length){
        mtx_ZeroOne[i,] <- ifelse(mtx[,1] <= equal_set_wide[i,1] & mtx[,2] >= equal_set_wide[i,2], 1, 0)
      }
      # 中の足し算
      sums_vec <- t(params) %*% mtx_ZeroOne
      # 外の掛け算
      prod <- prod(sums_vec)
      return(-prod)
    }
    return(likelihood)
  }
  likelihood <- make_likelihood(n, func_length)
  # initial_params <- rep(1/func_length, func_length)
  # ui <- rbind(diag(func_length), rep(1, func_length))
  # ci <- c(rep(0, func_length), 1-epsilon)
  # check_constraints <- ui %*% initial_params - ci

  # result <- constrOptim(theta = initial_params, f = likelihood, grad = NULL, ui = ui, ci=ci, control = list(reltol = -1e-200, trace=T, REPORT=500, maxit= 100000))
  
  constraint <- function(params){
    return(sum(params)-1)
  }
  result <- nloptr(
    x0 = rep(1/func_length, func_length), 
    eval_f = likelihood,  # 目的関数
    lb = rep(0, func_length),
    eval_g_eq = constraint,  # 不等式制約（非線形制約）
    opts = list(
      print_level = 0,  # 進行状況表示
      algorithm = "NLOPT_LN_COBYLA",  # 制約付きの最適化アルゴリズム
      xtol_rel = 1e-4,  # 相対収束基準
      xtol_abs = 1e-4,
      maxeval = 1  # 最大反復回数
    )
  )
  # 4列目に生存関数の差(params)を代入
  equal_set_wide <- cbind(equal_set_wide,result$solution)
  surv_turnbull <- 1-cumsum(result$solution)
  # 5列目に生存関数が1段下がった後の生存割合を代入
  equal_set_wide <- cbind(equal_set_wide, surv_turnbull)
  surv_turnbull <- c(1, surv_turnbull)
  surv_turnbull <- surv_turnbull[-length(surv_turnbull)]
  # 6列目に生存関数が1段下がる前の生存割合を代入
  equal_set_wide <- cbind(equal_set_wide,surv_turnbull)
  equal_set_wide_df <- data.frame(equal_set_wide)
  colnames(equal_set_wide_df) <- c("L","R","LRdiff","surv_turnbull_diff","surv_turnbull_next","surv_turnbull_lag")
  return(equal_set_wide_df)
}

wmst_turnbull <- function(df, tau_0, tau_1=NULL){
  if (is.null(tau_1) == T & tail(df$R, 1)==Inf){
    tau_1 = tail(df$L, 1)
  }else if (is.null(tau_1) == T & tail(df$R, 1)!=Inf){
    tau_1 = tail(df$R, 1)
  }
  # print(paste0("tau_1:",tau_1))

  calc_df <- df %>% filter(R > tau_0 & L <= tau_1)
  # View(calc_df)
  #間の長方形を埋める
  calc_df2 <- calc_df[1, ]  # 初期値として1行目
  
  for (i in 1:(nrow(calc_df) - 1)) {
    # 新しい行を挿入する
    new_row <- data.frame(
      L = calc_df$R[i],  # Lは1行目のR
      R = calc_df$L[i + 1],  # Rは次の行のL
      LRdiff = calc_df$L[i + 1] - calc_df$R[i],  # LRdiffは0
      surv_turnbull_diff = 0,  # surv_turnbull_diffは0
      surv_turnbull_next = calc_df$surv_turnbull_next[i],  # surv_turnbull_nextは1行目のsurv_turnbull_next
      surv_turnbull_lag = calc_df$surv_turnbull_next[i]   # surv_turnbull_lagは1行目のsurv_turnbull_next
    )
    
    # 1行目と2行目の間に新しい行を追加
    calc_df2 <- rbind(calc_df2, new_row, calc_df[i + 1, ])
  }
  row.names(calc_df2) <- NULL
  # View(calc_df2)
  # tau_0によって最初の長方形を作ったり、台形を削ったり
  if (tau_0 < head(calc_df2$L, 1)){
    new_raw <- data.frame(L=tau_0, R=head(calc_df2$L, 1), LRdiff=head(calc_df2$L, 1)-tau_0, surv_turnbull_diff=0, surv_turnbull_lag=head(calc_df2$surv_turnbull_lag, 1), surv_turnbull_next=head(calc_df2$surv_turnbull_lag, 1))
    calc_df2 <- rbind(new_raw, calc_df2)
  }else if(tau_0 > head(calc_df2$L, 1)){
    calc_df2$surv_turnbull_lag[1] <- calc_df2$surv_turnbull_lag[1] - calc_df2$surv_turnbull_diff[1] * ((calc_df2$R[1]-tau_0) / (calc_df2$R[1]-calc_df2$L[1]))
    calc_df2$surv_turnbull_diff[1] <- calc_df2$surv_turnbull_lag[1] - calc_df2$surv_turnbull_next[1]
    calc_df2$L[1] <- tau_0
    calc_df2$LRdiff[1] <- calc_df2$R[1] - calc_df2$L[1]
  }
  
  #最後が右側打ち切りなら削除
  if (tail(calc_df2$R, 1) == Inf){
    calc_df2 <- calc_df2[-nrow(calc_df2), ]
  }
  # 各台形の面積を計算
  calc_df2$area <- calc_trapezoid(base=calc_df2$surv_turnbull_lag, top=calc_df2$surv_turnbull_next, height=calc_df2$LRdiff)
  # View(calc_df2)
  wmst_est <- sum(calc_df2$area)
  return(wmst_est)
}

plot_turnbull <- function(df){
  
}