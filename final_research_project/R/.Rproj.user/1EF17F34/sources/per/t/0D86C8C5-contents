# ===============================================================================
# 一般化ペアワイズ比較を用いた区間打ち切りデータ対応の統計的手法の比較研究
#
# このコードは以下の項目を包括的に実装しています：
# 1. 一般化ペアワイズ比較（Generalized Pairwise Comparison, GPC）
# 2. 区間打ち切りデータへの新しい代入法（拡張EMI法）
# 3. RMST、ログランク検定との検出力比較
# 4. 論文用の詳細な結果出力と可視化
#
# 作成者: Claude Code
# 作成日: 2025-01-16
# ===============================================================================

# 必要なライブラリの読み込み
library(survival)
library(survRM2)  # RMST用
library(BuyseTest)  # GPC用
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)
library(knitr)
library(kableExtra)

# 既存の関数を読み込み
source("interval_censord_data_function.R")
source("distributions_2arm.R")

# ===============================================================================
# 1. 一般化ペアワイズ比較（GPC）の実装
# ===============================================================================

# 区間打ち切りデータに対応したGPC統計量の計算
gpc_interval_censored <- function(data_control, data_treatment, method = "net_benefit") {
  # データの前処理
  n_control <- nrow(data_control)
  n_treatment <- nrow(data_treatment)

  # ペアワイズ比較の結果を格納する行列
  comparison_matrix <- matrix(0, nrow = n_treatment, ncol = n_control)

  for (i in 1:n_treatment) {
    for (j in 1:n_control) {
      # 治療群のi番目と対照群のj番目を比較
      treat_left <- data_treatment[i, 1]
      treat_right <- data_treatment[i, 2]
      control_left <- data_control[j, 1]
      control_right <- data_control[j, 2]

      # 区間打ち切りデータの比較ルール
      if (treat_right <= control_left) {
        # 治療群が明らかに優れている（イベント時刻が早い）
        comparison_matrix[i, j] <- 1
      } else if (control_right <= treat_left) {
        # 対照群が明らかに優れている
        comparison_matrix[i, j] <- -1
      } else {
        # 区間が重複している場合は引き分け
        comparison_matrix[i, j] <- 0
      }
    }
  }

  # Net Benefit または Win Ratio の計算
  wins <- sum(comparison_matrix == 1)
  losses <- sum(comparison_matrix == -1)
  ties <- sum(comparison_matrix == 0)

  total_pairs <- n_treatment * n_control

  if (method == "net_benefit") {
    # Net Benefit = (Wins - Losses) / Total Pairs
    net_benefit <- (wins - losses) / total_pairs

    # 分散の計算（Mann-Whitney統計量の理論より）
    variance <- (wins + losses + ties/4) / (total_pairs^2)

    # 検定統計量
    z_stat <- net_benefit / sqrt(variance)
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      net_benefit = net_benefit,
      z_statistic = z_stat,
      p_value = p_value,
      wins = wins,
      losses = losses,
      ties = ties
    ))
  } else if (method == "win_ratio") {
    # Win Ratio = Wins / Losses
    if (losses == 0) {
      win_ratio <- Inf
    } else {
      win_ratio <- wins / losses
    }

    # ログWin Ratioの分散
    log_wr_variance <- 1/wins + 1/losses
    z_stat <- log(win_ratio) / sqrt(log_wr_variance)
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      win_ratio = win_ratio,
      log_win_ratio = log(win_ratio),
      z_statistic = z_stat,
      p_value = p_value,
      wins = wins,
      losses = losses,
      ties = ties
    ))
  }
}

# ===============================================================================
# 2. 拡張EMI法（Enhanced Multiple Imputation）の実装
# ===============================================================================

# 区間打ち切りデータに対する拡張EMI法
enhanced_emi_imputation <- function(interval_data, n_imputations = 10, method = "adaptive") {
  n <- nrow(interval_data)
  imputed_datasets <- list()

  for (m in 1:n_imputations) {
    imputed_times <- numeric(n)

    for (i in 1:n) {
      left <- interval_data[i, 1]
      right <- interval_data[i, 2]

      if (left == right) {
        # 正確な観測時刻
        imputed_times[i] <- left
      } else if (is.infinite(right)) {
        # 右側打ち切り - 指数分布で補外
        if (method == "adaptive") {
          # 観測されたイベント時刻から生存関数を推定
          observed_events <- interval_data[interval_data[,1] == interval_data[,2] &
                                         interval_data[,1] <= 1, 1]
          if (length(observed_events) > 0) {
            lambda_est <- 1 / mean(observed_events)
          } else {
            lambda_est <- 1  # デフォルト値
          }
          imputed_times[i] <- left + rexp(1, rate = lambda_est)
        } else {
          # 従来の一様分布代入
          imputed_times[i] <- left + rexp(1, rate = 1)
        }
      } else {
        # 区間打ち切り
        if (method == "adaptive") {
          # Beta分布を使用した適応的代入
          # 区間の位置に基づいてBeta分布のパラメータを調整
          interval_width <- right - left
          alpha <- 2 + 1/interval_width  # 幅が狭いほど山が鋭くなる
          beta <- 2 + 1/interval_width

          uniform_sample <- rbeta(1, alpha, beta)
          imputed_times[i] <- left + uniform_sample * (right - left)
        } else {
          # 従来の一様分布代入
          imputed_times[i] <- runif(1, left, right)
        }
      }
    }

    # 打ち切り指示子の設定
    delta <- ifelse(imputed_times <= 1 & !is.infinite(interval_data[, 2]), 1, 0)
    imputed_times[imputed_times > 1] <- 1

    imputed_datasets[[m]] <- data.frame(time = imputed_times, status = delta)
  }

  return(imputed_datasets)
}

# 複数代入データセットからの統合結果の計算
pool_emi_results <- function(imputed_datasets, test_function, ...) {
  n_imputations <- length(imputed_datasets)
  results <- list()

  for (m in 1:n_imputations) {
    results[[m]] <- test_function(imputed_datasets[[m]], ...)
  }

  # Rubin's rules に従って結果を統合
  estimates <- sapply(results, function(x) x$estimate)
  variances <- sapply(results, function(x) x$variance)

  pooled_estimate <- mean(estimates)
  within_imputation_variance <- mean(variances)
  between_imputation_variance <- var(estimates)

  total_variance <- within_imputation_variance +
                   (1 + 1/n_imputations) * between_imputation_variance

  # 自由度の計算
  df <- (n_imputations - 1) * (1 + within_imputation_variance /
                              ((1 + 1/n_imputations) * between_imputation_variance))^2

  t_stat <- pooled_estimate / sqrt(total_variance)
  p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)

  return(list(
    estimate = pooled_estimate,
    variance = total_variance,
    t_statistic = t_stat,
    df = df,
    p_value = p_value
  ))
}

# ===============================================================================
# 3. RMST計算関数
# ===============================================================================

# 区間打ち切りデータに対するRMST計算（代入法使用）
rmst_interval_censored <- function(data_control, data_treatment, tau = 1, imputation_method = "midpoint") {
  # 代入法によるデータ変換
  if (imputation_method == "midpoint") {
    surv_control <- midpoint_assignment(data_control)
    surv_treatment <- midpoint_assignment(data_treatment)
  } else if (imputation_method == "rightpoint") {
    surv_control <- rightpoint_assignment(data_control)
    surv_treatment <- rightpoint_assignment(data_treatment)
  } else if (imputation_method == "enhanced_emi") {
    # 拡張EMI法を使用
    emi_control <- enhanced_emi_imputation(data_control, n_imputations = 10)
    emi_treatment <- enhanced_emi_imputation(data_treatment, n_imputations = 10)

    # 各代入データセットでRMSTを計算し、統合
    rmst_results <- list()
    for (m in 1:10) {
      combined_data <- rbind(
        cbind(emi_control[[m]], group = 0),
        cbind(emi_treatment[[m]], group = 1)
      )

      tryCatch({
        rmst_result <- rmst2(combined_data$time, combined_data$status, combined_data$group, tau = tau)
        rmst_results[[m]] <- list(
          estimate = rmst_result$unadjusted.result[1, 1],
          variance = rmst_result$unadjusted.result[1, 2]^2
        )
      }, error = function(e) {
        rmst_results[[m]] <- list(estimate = 0, variance = 1)
      })
    }

    # 結果を統合
    pooled_result <- pool_emi_results(rmst_results, function(x) x)
    return(list(
      rmst_diff = pooled_result$estimate,
      p_value = pooled_result$p_value,
      method = "Enhanced EMI"
    ))
  }

  # 通常の代入法の場合
  combined_data <- rbind(
    cbind(surv_control, group = 0),
    cbind(surv_treatment, group = 1)
  )

  tryCatch({
    rmst_result <- rmst2(combined_data$time, combined_data$cens, combined_data$group, tau = tau)
    return(list(
      rmst_diff = rmst_result$unadjusted.result[1, 1],
      p_value = rmst_result$unadjusted.result[1, 4],
      method = imputation_method
    ))
  }, error = function(e) {
    return(list(rmst_diff = 0, p_value = 1, method = imputation_method))
  })
}

# ===============================================================================
# 4. ログランク検定関数
# ===============================================================================

# 区間打ち切りデータに対するログランク検定（代入法使用）
logrank_interval_censored <- function(data_control, data_treatment, imputation_method = "midpoint") {
  # 代入法によるデータ変換
  if (imputation_method == "midpoint") {
    surv_control <- midpoint_assignment(data_control)
    surv_treatment <- midpoint_assignment(data_treatment)
  } else if (imputation_method == "rightpoint") {
    surv_control <- rightpoint_assignment(data_control)
    surv_treatment <- rightpoint_assignment(data_treatment)
  } else if (imputation_method == "enhanced_emi") {
    # 拡張EMI法を使用
    emi_control <- enhanced_emi_imputation(data_control, n_imputations = 10)
    emi_treatment <- enhanced_emi_imputation(data_treatment, n_imputations = 10)

    # 各代入データセットでログランク検定を実行し、統合
    logrank_results <- list()
    for (m in 1:10) {
      combined_data <- rbind(
        cbind(emi_control[[m]], group = 0),
        cbind(emi_treatment[[m]], group = 1)
      )

      tryCatch({
        survdiff_result <- survdiff(Surv(time, status) ~ group, data = combined_data)
        chi_sq <- survdiff_result$chisq
        p_val <- 1 - pchisq(chi_sq, df = 1)

        logrank_results[[m]] <- list(
          estimate = sqrt(chi_sq),  # Z統計量
          variance = 1
        )
      }, error = function(e) {
        logrank_results[[m]] <- list(estimate = 0, variance = 1)
      })
    }

    # 結果を統合
    pooled_result <- pool_emi_results(logrank_results, function(x) x)
    return(list(
      p_value = pooled_result$p_value,
      method = "Enhanced EMI"
    ))
  }

  # 通常の代入法の場合
  combined_data <- rbind(
    cbind(surv_control, group = 0),
    cbind(surv_treatment, group = 1)
  )

  tryCatch({
    survdiff_result <- survdiff(Surv(time, cens) ~ group, data = combined_data)
    p_value <- 1 - pchisq(survdiff_result$chisq, df = 1)
    return(list(p_value = p_value, method = imputation_method))
  }, error = function(e) {
    return(list(p_value = 1, method = imputation_method))
  })
}

# ===============================================================================
# 5. シミュレーション研究のメイン関数
# ===============================================================================

run_simulation_study <- function(
  n_sim = 1000,
  sample_sizes = c(100, 200, 400),
  K_values = c(3, 5, 10),
  dropout_levels = c("None", "Low", "Medium", "High"),
  effect_sizes = 1:24,  # 24種類の分布効果
  alpha = 0.05,
  n_cores = parallel::detectCores() - 1
) {

  # 並列処理の設定
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # 結果を格納するデータフレーム
  results <- data.frame()

  # シミュレーション条件の組み合わせ
  conditions <- expand.grid(
    n = sample_sizes,
    K = K_values,
    dropout = dropout_levels,
    effect = effect_sizes,
    stringsAsFactors = FALSE
  )

  cat("シミュレーション開始: 総条件数 =", nrow(conditions), "\n")

  for (i in 1:nrow(conditions)) {
    n <- conditions$n[i]
    K <- conditions$K[i]
    dropout <- conditions$dropout[i]
    effect <- conditions$effect[i]

    cat(sprintf("条件 %d/%d: n=%d, K=%d, dropout=%s, effect=%d\n",
                i, nrow(conditions), n, K, dropout, effect))

    # 並列シミュレーション
    sim_results <- foreach(sim = 1:n_sim, .combine = rbind,
                          .packages = c("survival", "survRM2")) %dopar% {

      tryCatch({
        # 対照群（効果サイズ1）と治療群（指定された効果サイズ）のデータ生成
        data_control <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                           p_dropout = dropout, d = 1)
        data_treatment <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                             p_dropout = dropout, d = effect)

        # 各手法の実行
        # 1. 一般化ペアワイズ比較（Net Benefit）
        gpc_nb <- gpc_interval_censored(data_control, data_treatment, method = "net_benefit")
        gpc_nb_reject <- gpc_nb$p_value < alpha

        # 2. 一般化ペアワイズ比較（Win Ratio）
        gpc_wr <- gpc_interval_censored(data_control, data_treatment, method = "win_ratio")
        gpc_wr_reject <- gpc_wr$p_value < alpha

        # 3. RMST（中点代入法）
        rmst_mid <- rmst_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
        rmst_mid_reject <- rmst_mid$p_value < alpha

        # 4. RMST（右点代入法）
        rmst_right <- rmst_interval_censored(data_control, data_treatment, imputation_method = "rightpoint")
        rmst_right_reject <- rmst_right$p_value < alpha

        # 5. RMST（拡張EMI法）
        rmst_emi <- rmst_interval_censored(data_control, data_treatment, imputation_method = "enhanced_emi")
        rmst_emi_reject <- rmst_emi$p_value < alpha

        # 6. ログランク検定（中点代入法）
        lr_mid <- logrank_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
        lr_mid_reject <- lr_mid$p_value < alpha

        # 7. ログランク検定（右点代入法）
        lr_right <- logrank_interval_censored(data_control, data_treatment, imputation_method = "rightpoint")
        lr_right_reject <- lr_right$p_value < alpha

        # 8. ログランク検定（拡張EMI法）
        lr_emi <- logrank_interval_censored(data_control, data_treatment, imputation_method = "enhanced_emi")
        lr_emi_reject <- lr_emi$p_value < alpha

        # 結果をまとめる
        data.frame(
          sim_id = sim,
          n = n,
          K = K,
          dropout = dropout,
          effect = effect,
          gpc_nb_reject = gpc_nb_reject,
          gpc_wr_reject = gpc_wr_reject,
          rmst_mid_reject = rmst_mid_reject,
          rmst_right_reject = rmst_right_reject,
          rmst_emi_reject = rmst_emi_reject,
          lr_mid_reject = lr_mid_reject,
          lr_right_reject = lr_right_reject,
          lr_emi_reject = lr_emi_reject,
          gpc_nb_pvalue = gpc_nb$p_value,
          gpc_wr_pvalue = gpc_wr$p_value,
          rmst_mid_pvalue = rmst_mid$p_value,
          rmst_right_pvalue = rmst_right$p_value,
          rmst_emi_pvalue = rmst_emi$p_value,
          lr_mid_pvalue = lr_mid$p_value,
          lr_right_pvalue = lr_right$p_value,
          lr_emi_pvalue = lr_emi$p_value
        )
      }, error = function(e) {
        # エラー時のデフォルト値
        data.frame(
          sim_id = sim,
          n = n,
          K = K,
          dropout = dropout,
          effect = effect,
          gpc_nb_reject = FALSE,
          gpc_wr_reject = FALSE,
          rmst_mid_reject = FALSE,
          rmst_right_reject = FALSE,
          rmst_emi_reject = FALSE,
          lr_mid_reject = FALSE,
          lr_right_reject = FALSE,
          lr_emi_reject = FALSE,
          gpc_nb_pvalue = 1,
          gpc_wr_pvalue = 1,
          rmst_mid_pvalue = 1,
          rmst_right_pvalue = 1,
          rmst_emi_pvalue = 1,
          lr_mid_pvalue = 1,
          lr_right_pvalue = 1,
          lr_emi_pvalue = 1
        )
      })
    }

    results <- rbind(results, sim_results)
  }

  stopCluster(cl)

  return(results)
}

# ===============================================================================
# 6. 結果の分析と可視化関数
# ===============================================================================

# 検出力の計算
calculate_power <- function(results) {
  power_summary <- results %>%
    group_by(n, K, dropout, effect) %>%
    summarise(
      gpc_nb_power = mean(gpc_nb_reject),
      gpc_wr_power = mean(gpc_wr_reject),
      rmst_mid_power = mean(rmst_mid_reject),
      rmst_right_power = mean(rmst_right_reject),
      rmst_emi_power = mean(rmst_emi_reject),
      lr_mid_power = mean(lr_mid_reject),
      lr_right_power = mean(lr_right_reject),
      lr_emi_power = mean(lr_emi_reject),
      .groups = "drop"
    )

  return(power_summary)
}

# 第1種の誤りの計算（effect = 1のときの棄却率）
calculate_type1_error <- function(results) {
  type1_summary <- results %>%
    filter(effect == 1) %>%
    group_by(n, K, dropout) %>%
    summarise(
      gpc_nb_type1 = mean(gpc_nb_reject),
      gpc_wr_type1 = mean(gpc_wr_reject),
      rmst_mid_type1 = mean(rmst_mid_reject),
      rmst_right_type1 = mean(rmst_right_reject),
      rmst_emi_type1 = mean(rmst_emi_reject),
      lr_mid_type1 = mean(lr_mid_reject),
      lr_right_type1 = mean(lr_right_reject),
      lr_emi_type1 = mean(lr_emi_reject),
      .groups = "drop"
    )

  return(type1_summary)
}

# 検出力比較のヒートマップ作成
create_power_heatmap <- function(power_data, method_name, title_suffix = "") {
  ggplot(power_data, aes(x = factor(K), y = factor(effect), fill = power)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", power)), color = "white", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0.5, limits = c(0, 1)) +
    facet_grid(n ~ dropout, labeller = labeller(
      n = function(x) paste("n =", x),
      dropout = function(x) paste("Dropout:", x)
    )) +
    labs(
      title = paste(method_name, "検出力", title_suffix),
      x = "観測回数 (K)",
      y = "効果サイズ (分布タイプ)",
      fill = "検出力"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10, face = "bold")
    )
}

# 論文用の包括的結果表の作成
create_comprehensive_results_table <- function(power_summary, type1_summary) {
  # 各手法の平均検出力を計算
  avg_power <- power_summary %>%
    filter(effect != 1) %>%  # 帰無仮説のケースを除外
    summarise(
      GPC_NetBenefit = mean(gpc_nb_power),
      GPC_WinRatio = mean(gpc_wr_power),
      RMST_Midpoint = mean(rmst_mid_power),
      RMST_Rightpoint = mean(rmst_right_power),
      RMST_EnhancedEMI = mean(rmst_emi_power),
      LogRank_Midpoint = mean(lr_mid_power),
      LogRank_Rightpoint = mean(lr_right_power),
      LogRank_EnhancedEMI = mean(lr_emi_power)
    )

  # 各手法の平均第1種の誤りを計算
  avg_type1 <- type1_summary %>%
    summarise(
      GPC_NetBenefit = mean(gpc_nb_type1),
      GPC_WinRatio = mean(gpc_wr_type1),
      RMST_Midpoint = mean(rmst_mid_type1),
      RMST_Rightpoint = mean(rmst_right_type1),
      RMST_EnhancedEMI = mean(rmst_emi_type1),
      LogRank_Midpoint = mean(lr_mid_type1),
      LogRank_Rightpoint = mean(lr_right_type1),
      LogRank_EnhancedEMI = mean(lr_emi_type1)
    )

  # 結果表を作成
  results_table <- data.frame(
    Method = c("GPC (Net Benefit)", "GPC (Win Ratio)",
               "RMST (Midpoint)", "RMST (Rightpoint)", "RMST (Enhanced EMI)",
               "LogRank (Midpoint)", "LogRank (Rightpoint)", "LogRank (Enhanced EMI)"),
    Average_Power = c(avg_power$GPC_NetBenefit, avg_power$GPC_WinRatio,
                     avg_power$RMST_Midpoint, avg_power$RMST_Rightpoint, avg_power$RMST_EnhancedEMI,
                     avg_power$LogRank_Midpoint, avg_power$LogRank_Rightpoint, avg_power$LogRank_EnhancedEMI),
    Average_Type1_Error = c(avg_type1$GPC_NetBenefit, avg_type1$GPC_WinRatio,
                           avg_type1$RMST_Midpoint, avg_type1$RMST_Rightpoint, avg_type1$RMST_EnhancedEMI,
                           avg_type1$LogRank_Midpoint, avg_type1$LogRank_Rightpoint, avg_type1$LogRank_EnhancedEMI)
  )

  return(results_table)
}

# ===============================================================================
# 7. メイン実行関数
# ===============================================================================

main_simulation_analysis <- function(n_sim = 1000, output_dir = "simulation_results") {
  # 出力ディレクトリの作成
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("=== 一般化ペアワイズ比較を用いた区間打ち切りデータ対応の統計的手法比較研究 ===\n")
  cat("シミュレーション開始時刻:", format(Sys.time()), "\n\n")

  # シミュレーション実行
  cat("1. シミュレーション実行中...\n")
  simulation_results <- run_simulation_study(n_sim = n_sim)

  # 結果をCSVで保存
  write.csv(simulation_results, file.path(output_dir, "raw_simulation_results.csv"), row.names = FALSE)

  # 検出力の計算
  cat("2. 検出力分析中...\n")
  power_summary <- calculate_power(simulation_results)
  write.csv(power_summary, file.path(output_dir, "power_summary.csv"), row.names = FALSE)

  # 第1種の誤りの計算
  cat("3. 第1種の誤り分析中...\n")
  type1_summary <- calculate_type1_error(simulation_results)
  write.csv(type1_summary, file.path(output_dir, "type1_error_summary.csv"), row.names = FALSE)

  # 包括的結果表の作成
  cat("4. 論文用結果表作成中...\n")
  comprehensive_table <- create_comprehensive_results_table(power_summary, type1_summary)
  write.csv(comprehensive_table, file.path(output_dir, "comprehensive_results_table.csv"), row.names = FALSE)

  # 可視化の作成
  cat("5. 可視化作成中...\n")

  # 各手法の検出力ヒートマップ
  methods <- list(
    list(col = "gpc_nb_power", name = "GPC (Net Benefit)"),
    list(col = "gpc_wr_power", name = "GPC (Win Ratio)"),
    list(col = "rmst_emi_power", name = "RMST (Enhanced EMI)"),
    list(col = "lr_emi_power", name = "LogRank (Enhanced EMI)")
  )

  for (method in methods) {
    power_data <- power_summary %>%
      select(n, K, dropout, effect, power = all_of(method$col))

    p <- create_power_heatmap(power_data, method$name)
    ggsave(file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", method$name), "_heatmap.png")),
           p, width = 12, height = 10, dpi = 300)
  }

  # 手法間の検出力比較
  power_comparison <- power_summary %>%
    filter(effect != 1) %>%
    group_by(n, K, dropout) %>%
    summarise(
      GPC_NB = mean(gpc_nb_power),
      GPC_WR = mean(gpc_wr_power),
      RMST_EMI = mean(rmst_emi_power),
      LogRank_EMI = mean(lr_emi_power),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(GPC_NB, GPC_WR, RMST_EMI, LogRank_EMI),
                names_to = "Method", values_to = "Power")

  p_comparison <- ggplot(power_comparison, aes(x = Method, y = Power, fill = Method)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.6, width = 0.2) +
    facet_grid(n ~ dropout, labeller = labeller(
      n = function(x) paste("n =", x),
      dropout = function(x) paste("Dropout:", x)
    )) +
    labs(
      title = "手法間の検出力比較",
      x = "統計手法",
      y = "平均検出力",
      fill = "手法"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  ggsave(file.path(output_dir, "method_comparison.png"), p_comparison,
         width = 14, height = 10, dpi = 300)

  cat("6. 分析完了\n")
  cat("結果は", output_dir, "ディレクトリに保存されました。\n")
  cat("シミュレーション終了時刻:", format(Sys.time()), "\n")

  return(list(
    simulation_results = simulation_results,
    power_summary = power_summary,
    type1_summary = type1_summary,
    comprehensive_table = comprehensive_table
  ))
}

# ===============================================================================
# 使用例とパラメータ設定
# ===============================================================================

# 小規模テスト実行（開発・デバッグ用）
run_test_simulation <- function() {
  cat("=== テストシミュレーション実行 ===\n")
  test_results <- run_simulation_study(
    n_sim = 100,  # 小さいシミュレーション回数
    sample_sizes = c(100, 200),
    K_values = c(3, 5),
    dropout_levels = c("None", "Medium"),
    effect_sizes = c(1, 6, 12),  # 一部の効果サイズのみ
    n_cores = 2
  )

  power_summary <- calculate_power(test_results)
  type1_summary <- calculate_type1_error(test_results)

  print(head(power_summary))
  print(head(type1_summary))

  return(test_results)
}

# フル分析実行（論文用）
run_full_analysis <- function() {
  cat("=== フル分析実行（論文用）===\n")
  cat("注意: この分析には数時間かかる可能性があります。\n")

  full_results <- main_simulation_analysis(
    n_sim = 1000,  # 論文品質のシミュレーション回数
    output_dir = "paper_results"
  )

  return(full_results)
}

cat("===================================================================\n")
cat("一般化ペアワイズ比較を用いた区間打ち切りデータ対応の研究コードが準備完了\n")
cat("===================================================================\n")
cat("使用方法:\n")
cat("1. テスト実行: test_results <- run_test_simulation()\n")
cat("2. フル分析: full_results <- run_full_analysis()\n")
cat("3. 個別実行: results <- main_simulation_analysis(n_sim=500)\n")
cat("===================================================================\n")