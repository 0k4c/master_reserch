source("C:/Users/81801/Documents/R/code/senior_thesis/distributions_2arm.R")

#z֐̎ނǉ
generate_interval_censored_data_2arm <- function(n=100,K=5,p_dropout="Medium",p_exact=0,d=1){
  #
  t_end <- 1
  epsilon <- 1e-6
  
  #G[ƌx\
  if (!is.element(n,c(100,200,400))) stop("only (100,200,400) is permited for n")
  if (!is.element(K,c(3,5,10,20))) stop("only (3,5,10,20) is permited for K")
  if (!is.element(p_dropout,c("None","Low","Medium","High"))) stop('only ("None","Low","Medium","High") is permited for p_dropout')
  if (!is.element(p_exact,c(0,0.2,0.5,1.0))) stop("only (0,0.2,0.5,1.0) is permited for p_exact")
  if (!is.element(d,1:24)) stop("only 1:24 is permited for d")
  
  #ԕL����߂
  L <- t_end/(K+1)
  
  #Obh̒����
  grid_length <- K+3
  
  #n̃Obh����
  grids <- array(numeric(2*(K+3)*n), dim=c(K+3,2,n))
  
  #n̊ϑJn_lz琶
  tau_0s <- runif(n, 0, L)
  
  #ϑ_s쐬
  mat1 <- array(numeric((K+3)*n), dim=c(n,K+3))
  mat1[,3:(K+2)] <- L
  mat2 <- array(numeric((K+3)*n), dim=c(n,K+3))
  mat2[,2:(K+2)] <- tau_0s
  mat2[,K+3] <- Inf
  #ϑ_zObhɑ
  grids[,1,] <- apply(mat1,1,cumsum)+aperm(mat2,perm=c(2,1))
  
  #ϑz𐶐
  if (p_dropout=="None") p_drops <- numeric(grid_length)
  else if (p_dropout=="Low") p_drops <- c(numeric(2),rep(0.1,K-1),0.2, 0)
  else if (p_dropout=="Medium") p_drops <- c(numeric(2), rep(0.2, K-1), 0.4, 0)
  else if (p_dropout=="High") p_drops <- c(numeric(2), rep(0.3, K-1), 0.6, 0)
  else print("Error")
  #ϑzObhɑ(1:,0:)
  grids[,2,] <- rbinom(grid_length*n, size=1, prob=p_drops)
  #print("grids")
  #print(grids)
  
  #^̊ϑf[^𐶐(nɊg)
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
  
  #mȊϑǂxk[CzTvO(1:m,0:ԑł؂)
  xis <- rbinom(n,size=1,prob=p_exact)
  #print("xis")
  #print(xis)
  
  #ԑł؂f[^
  #-------------------------------------------------------------
  interval_datas <- array(numeric(2*n),dim=c(n,2))
  #print("interval_datas")
  #print(interval_datas)
  
  #mȊϑf[^
  nums1 <- which(xis==1)
  #print("nums1")
  #print(nums1)
  interval_datas[nums1,] <- c(T_vals[xis==1],T_vals[xis==1])
  #print("interval_datas")
  #print(interval_datas)
  
  #T_val>1̋ԑł؂f[^
  nums2 <- c(which(xis==0&T_vals>t_end))
  #print("nums2")
  #print(nums2)
  #if(length(nums2)>0){
  #drop_mask͌FALSEƂ
  drop_mask <- as.logical(grids[,2,nums2])
  #print("drop_mask")
  #print(drop_mask)
  cens_point_candidate <- replace(grids[,1,nums2],c(which(drop_mask|grids[,1,nums2]==Inf)),0)
  #print("cens_point_carndidate")
  #print(array(cens_point_candidate))
  #OȂapplyKp邽߂matrixϊ
  #cens_point_candidate_mat = matrix(cens_point_candidate)
  if (is.null(dim(cens_point_candidate))){
    #arrayłȂappplyȂ̂applyKp邽߂arrayɕϊ
    cens_point_candidate <- array(cens_point_candidate,dim=c(length(cens_point_candidate),1))
  }
  interval_datas[nums2,1] <- apply(cens_point_candidate,2,max)
  interval_datas[nums2,2] <- Inf
  #print("interval_datas")
  #print(interval_datas)
  #}
  
  #T_val<=1̋ԑł؂f[^
  nums3 <- c(which(xis==0&T_vals<=t_end))
  #if(length(nums3)>0){
  #OȂapermg߂ɕϊ
  small_grid <- grids[,1,nums3]
  if (is.null(dim(small_grid))){
    small_grid <- array(small_grid,dim=c(length(small_grid),1))
  }
  l_logic_mask <- aperm(aperm(small_grid, perm=c(2,1))-T_vals[nums3],perm=c(2,1))<0
  #print("l_logic_mask")
  #print(l_logic_mask)
  r_logic_mask <- !l_logic_mask
  #print("r_logic_mask")
  #print(r_logic_mask)
  #drop_mask͌TUREƂ
  drop_mask <- as.logical(grids[,2,nums3])
  l_candidate <- replace(grids[,1,nums3],c(which(r_logic_mask | drop_mask)),0)
  #OȂapplyKp邽߂matrixϊ
  l_candidate_mat <- matrix(l_candidate)
  if (is.null(dim(l_candidate))){
    l_candidate = array(l_candidate,dim=c(length(l_candidate),1))
  }
  #print(l_candidate)
  #print(l_candidate_mat)
  interval_datas[nums3,1] <- apply(l_candidate,2,max)
  r_candidate <- replace(grids[,1,nums3],c(which(l_logic_mask | drop_mask)),Inf)
  #OȂapplyKp邽߂matrixϊ
  r_candidate_mat <- matrix(r_candidate)
  if (is.null(dim(r_candidate))){
    r_candidate = array(r_candidate,dim=c(length(r_candidate),1))
  }
  interval_datas[nums3,2] <- apply(r_candidate,2,min)
  #}
  return(invisible(interval_datas))
}

#_Đԃf[^𐶐(t,delta)delta==0Eł؂(_@)
midpoint_assignment <- function(matrix){
  #ʂi[result����
  result <- array(numeric(nrow(matrix)*2),dim=c(nrow(matrix),2))
  #ԑł؂f[^𒆓_(delta=1:Cxg)
  result[,1][matrix[,2]<=1] <- (matrix[,1][matrix[,2]<=1]+matrix[,2][matrix[,2]<=1])/2
  result[,2][matrix[,2]<=1] <- 1
  #ώ@rł؂(delta=0:ł؂菈)
  result[,1][matrix[,2]>1&matrix[,1]<=1] <- matrix[,1][matrix[,2]>1&matrix[,1]<=1]
  result[,2][matrix[,2]>1&matrix[,1]<=1] <- 0
  #ώ@ł؂(delta=0:ł؂菈)
  result[,1][matrix[,1]>1] <- 1
  result[,2][matrix[,1]>1] <- 0
  #ʂ܂Ƃ߂ĕԂ
  result_df <- data.frame(time=result[,1],cens=result[,2])
  return(result_df)
}

#E[@
rightpoint_assignment <- function(matrix){
  #ʂi[result����
  result <- array(numeric(nrow(matrix)*2),dim=c(nrow(matrix),2))
  #ԑł؂f[^E_(delta=1:Cxg)
  result[,1][matrix[,2]<=1] <- matrix[,2][matrix[,2]<=1]
  result[,2][matrix[,2]<=1] <- 1
  #ώ@rł؂(delta=0:ł؂菈)
  result[,1][matrix[,2]>1&matrix[,1]<=1] <- matrix[,1][matrix[,2]>1&matrix[,1]<=1]
  result[,2][matrix[,2]>1&matrix[,1]<=1] <- 0
  #ώ@ł؂(delta=0:ł؂菈)
  result[,1][matrix[,1]>1] <- 1
  result[,2][matrix[,1]>1] <- 0
  #ʂ܂Ƃ߂ĕԂ1],cens=result[,2])
  return(result_df)
}
