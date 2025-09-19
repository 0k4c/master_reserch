# ===============================================================================
# 柔軟な区間打ち切りデータ生成関数（制限緩和版）
#
# 元の関数の制限を緩和し、テスト用の小さいサンプルサイズにも対応
# ===============================================================================

# 元の関数を読み込み
source("interval_censord_data_function.R")
source("distributions_2arm.R")

# 制限を緩和した区間打ち切りデータ生成関数
generate_interval_censored_data_flexible <- function(n=100, K=5, p_dropout="Medium", p_exact=0, d=1){
  # パラメータ設定
  t_end <- 1  # 観察終了時点
  epsilon <- 1e-6

  # 緩和されたエラーチェック
  if (n < 10 || n > 1000)
    stop("n must be between 10 and 1000")
  if (!is.element(K, c(3,5,10,20)))
    stop('only (3,5,10,20) is permitted for K')
  if (!is.element(p_dropout, c("None","Low","Medium","High")))
    stop('only ("None","Low","Medium","High") is permitted for p_dropout')
  if (!is.element(p_exact, c(0,0.2,0.5,1.0)))
    stop("only (0,0.2,0.5,1.0) is permitted for p_exact")
  if (!is.element(d, 1:24))
    stop("only 1:24 is permitted for d")

  # 観測間隔の長さを決定
  L <- t_end/(K+1)

  # グリッドの長さ（観測時点数）
  grid_length <- K+3

  # n人分のグリッド配列
  grids <- array(numeric(2*(K+3)*n), dim=c(K+3,2,n))

  # n人の観測開始時点を一様分布から生成
  tau_0s <- 0
    #runif(n, 0, L)

  # 観測時点の作成
  mat1 <- array(numeric((K+3)*n), dim=c(n,K+3))
  mat1[,3:(K+2)] <- L
  mat2 <- array(numeric((K+3)*n), dim=c(n,K+3))
  mat2[,2:(K+2)] <- tau_0s
  mat2[,K+3] <- Inf

  # 観測時点をグリッドに代入
  grids[,1,] <- apply(mat1, 1, cumsum) + aperm(mat2, perm=c(2,1))

  # 観測脱落確率の設定
  if (p_dropout=="None") {
    p_drops <- numeric(grid_length)
  } else if (p_dropout=="Low") {
    p_drops <- c(numeric(2), rep(0.1, K-1), 0.2, 0)
  } else if (p_dropout=="Medium") {
    p_drops <- c(numeric(2), rep(0.2, K-1), 0.4, 0)
  } else if (p_dropout=="High") {
    p_drops <- c(numeric(2), rep(0.3, K-1), 0.6, 0)
  } else {
    print("Error")
  }

  # 観測脱落をグリッドに代入（1:脱落, 0:観測継続）
  grids[,2,] <- rbinom(grid_length*n, size=1, prob=p_drops)

  # 真の生存時間データを生成（分布タイプdに基づく）
  if (d==1) T_vals <- rweibull(n,scale=1,shape=1)
  else if (d==2) T_vals <- rweibull(n,scale=1,shape=0.5)
  else if (d==3) T_vals <- rweibull(n,scale=1,shape=2)
  else if (d==4) T_vals <- rT4(n)
  else if (d==5) T_vals <- rT5(n)
  else if (d==6) T_vals <- rweibull(n,scale=0.5,shape=1)
  else if (d==7) T_vals <- rweibull(n,scale=0.5,shape=0.5)
  else if (d==8) T_vals <- rweibull(n,scale=0.75,shape=2)
  else if (d==9) T_vals <- rT9(n)
  else if (d==10) T_vals <- rT10(n)
  else if (d==11) T_vals <- rT11(n)
  else if (d==12) T_vals <- rT12(n)
  else if (d==13) T_vals <- rT13(n)
  else if (d==14) T_vals <- rT14(n)
  else if (d==15) T_vals <- rT15(n)
  else if (d==16) T_vals <- rT16(n)
  else if (d==17) T_vals <- rT17(n)
  else if (d==18) T_vals <- rT18(n)
  else if (d==19) T_vals <- rT19(n)
  else if (d==20) T_vals <- rT20(n)
  else if (d==21) T_vals <- rT21(n)
  else if (d==22) T_vals <- rT22(n)
  else if (d==23) T_vals <- rT23(n)
  else if (d==24) T_vals <- rT24(n)
  else print("Error")

  # 正確な観測の有無をベルヌーイサンプリング（1:正確, 0:区間打ち切り）
  xis <- rbinom(n, size=1, prob=p_exact)

  # 区間打ち切りデータの配列を初期化
  interval_datas <- array(numeric(2*n), dim=c(n,2))

  # 正確な観測データの処理
  nums1 <- which(xis==1)
  if(length(nums1) > 0) {
    interval_datas[nums1,] <- cbind(T_vals[nums1], T_vals[nums1])
  }

  # T_val > 1の区間打ち切りデータの処理
  nums2 <- c(which(xis==0 & T_vals>t_end))
  if(length(nums2)>0){
    drop_mask <- as.logical(grids[,2,nums2])
    cens_point_candidate <- replace(grids[,1,nums2],
                                    c(which(drop_mask | grids[,1,nums2]==Inf)), 0)
    if (is.null(dim(cens_point_candidate))){
      cens_point_candidate <- array(cens_point_candidate,
                                    dim=c(length(cens_point_candidate),1))
    }
    interval_datas[nums2,1] <- apply(cens_point_candidate, 2, max)
    interval_datas[nums2,2] <- Inf
  }

  # T_val <= 1の区間打ち切りデータの処理
  nums3 <- c(which(xis==0 & T_vals<=t_end))
  if(length(nums3)>0){
    small_grid <- grids[,1,nums3]
    if (is.null(dim(small_grid))){
      small_grid <- array(small_grid, dim=c(length(small_grid),1))
    }
    l_logic_mask <- aperm(aperm(small_grid, perm=c(2,1))-T_vals[nums3], perm=c(2,1)) < 0
    r_logic_mask <- !l_logic_mask
    drop_mask <- as.logical(grids[,2,nums3])

    # 左端点の計算
    l_candidate <- replace(grids[,1,nums3], c(which(r_logic_mask | drop_mask)), 0)
    if (is.null(dim(l_candidate))){
      l_candidate <- array(l_candidate, dim=c(length(l_candidate),1))
    }
    interval_datas[nums3,1] <- apply(l_candidate, 2, max)

    # 右端点の計算
    r_candidate <- replace(grids[,1,nums3], c(which(l_logic_mask | drop_mask)), Inf)
    if (is.null(dim(r_candidate))){
      r_candidate <- array(r_candidate, dim=c(length(r_candidate),1))
    }
    interval_datas[nums3,2] <- apply(r_candidate, 2, min)
  }

  return(invisible(interval_datas))
}

# 既存の関数をオーバーライド（テスト用）
generate_interval_censored_data_2arm <- generate_interval_censored_data_flexible

cat("柔軟なデータ生成関数を読み込みました。\n")
cat("サンプルサイズの制限が 10-1000 に緩和されました。\n")