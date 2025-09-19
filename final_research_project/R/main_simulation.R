# ===============================================================================
# 最終版：一般化ペアワイズ比較における代入法の比較研究
#
# 研究目的：
# 1. GPCを用いる際に、どの代入法が最も優れているか？
#    - Direct (代入なし)
#    - Midpoint Assignment
#    - Rightpoint Assignment
#    - Enhanced EMI
#
# 2. GPC手法は従来手法（RMST, ログランク検定）と比較してどうか？
#
# 修正された問題点：
# - p値計算の修正
# - 生存時間解釈の統一（長い方が良い）
# - 正しい統計理論の適用
# ===============================================================================

# 必要なライブラリ
library(survival)
library(tidyverse)
library(ggplot2)
library(parallel)
library(foreach)
library(doSNOW)
library(progress)

# すべての依存スクリプトを読み込む
source("flexible_data_generation.R")
source("improved_gpc_with_packages.R")
source("final_imputation_comparison_study.R")
source("generalized_pairwise_comparison_simulation.R")



# ===============================================================================
# 1. 研究デザインの明確化
# ===============================================================================

research_design <- function() {
  cat("===============================================================================\n")
  cat("             一般化ペアワイズ比較における代入法の比較研究\n")
  cat("===============================================================================\n")
  cat("主要研究目的:\n")
  cat("Q1: GPCにおいて、どの代入法が最も検出力が高いか？\n")
  cat("    - Direct GPC (代入なし)\n")
  cat("    - Midpoint Assignment + GPC\n")
  cat("    - Rightpoint Assignment + GPC\n")
  cat("    - Enhanced EMI + GPC\n\n")

  cat("Q2: GPC手法は従来手法と比較して優れているか？\n")
  cat("    - RMST\n")
  cat("    - ログランク検定\n\n")

  cat("評価指標:\n")
  cat("- 検出力 (Power)\n")
  cat("- 第1種の誤り (Type I Error)\n")
  cat("- 計算時間\n")
  cat("- 頑健性\n")
  cat("===============================================================================\n")
}

# ===============================================================================
# 使用例
# ===============================================================================

# 小規模テスト（プログレスバーあり、4コア使用）
# test_results <- run_final_simulation(
#   n_sim = 100,
#   sample_sizes = c(50, 100),
#   K_values = c(3, 5),
#   dropout_levels = c("None", "Medium"),
#   effect_sizes = c(1, 6),
#   n_cores = 4,
#   use_parallel = TRUE,
#   show_progress = TRUE
# )

# 大規模研究（プログレスバーなし、全コア使用）
# full_results <- run_final_simulation(
#   n_sim = 1000,
#   sample_sizes = c(100, 200, 400),
#   K_values = c(3, 5, 10),
#   dropout_levels = c("None", "Low", "Medium", "High"),
#   effect_sizes = 1:24,
#   n_cores = NULL,  # 自動検出
#   use_parallel = TRUE,
#   show_progress = FALSE  # 大規模時は無効化推奨
# )

# 逐次処理版（デバッグ用）
# debug_results <- run_final_simulation(
#   n_sim = 10,
#   sample_sizes = c(50),
#   K_values = c(3),
#   dropout_levels = c("None"),
#   effect_sizes = c(1, 6),
#   n_cores = 1,
#   use_parallel = FALSE,
#   show_progress = TRUE
# )


# ===============================================================================
# 2. 改良版代入法GPC関数
# ===============================================================================

# Direct GPC（代入なし）- 修正版
gpc_direct_improved <- function(data_control, data_treatment, method = "net_benefit") {
  n_control <- nrow(data_control)
  n_treatment <- nrow(data_treatment)

  wins <- 0
  losses <- 0
  ties <- 0

  for (i in 1:n_treatment) {
    for (j in 1:n_control) {
      L_T <- data_treatment[i, 1]  # 治療群左端点
      R_T <- data_treatment[i, 2]  # 治療群右端点
      L_C <- data_control[j, 1]    # 対照群左端点
      R_C <- data_control[j, 2]    # 対照群右端点

      # 生存時間解釈：長い方が良い
      if (L_T > R_C) {
        wins <- wins + 1      # 治療群明確に優秀
      } else if (R_T < L_C) {
        losses <- losses + 1  # 治療群明確に劣る
      } else {
        ties <- ties + 1      # 判定不能
      }
    }
  }

  total_pairs <- n_treatment * n_control

  if (method == "net_benefit") {
    net_benefit <- (wins - losses) / total_pairs

    # 正しい分散推定
    p_win <- wins / total_pairs
    p_loss <- losses / total_pairs
    variance <- (p_win + p_loss - (p_win - p_loss)^2) / total_pairs

    if (variance <= 0) variance <- 1 / (4 * total_pairs)

    z_stat <- net_benefit / sqrt(variance)
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      method = "Direct_NetBenefit",
      estimate = net_benefit,
      p_value = p_value,
      wins = wins, losses = losses, ties = ties
    ))

  } else if (method == "win_ratio") {
    if (losses == 0) {
      win_ratio <- ifelse(wins > 0, Inf, 1)
      p_value <- ifelse(wins > 0, 0, 1)
    } else if (wins == 0) {
      win_ratio <- 0
      p_value <- 0
    } else {
      win_ratio <- wins / losses
      log_var <- 1/wins + 1/losses
      z_stat <- log(win_ratio) / sqrt(log_var)
      p_value <- 2 * (1 - pnorm(abs(z_stat)))
    }

    return(list(
      method = "Direct_WinRatio",
      estimate = win_ratio,
      p_value = p_value,
      wins = wins, losses = losses, ties = ties
    ))
  }
}

# 代入法適用GPC
gpc_with_imputation_improved <- function(data_control, data_treatment,
                                       imputation_method = "midpoint",
                                       gpc_method = "net_benefit") {

  # 代入法適用
  if (imputation_method == "midpoint") {
    surv_control <- midpoint_assignment(data_control)
    surv_treatment <- midpoint_assignment(data_treatment)
  } else if (imputation_method == "rightpoint") {
    surv_control <- rightpoint_assignment(data_control)
    surv_treatment <- rightpoint_assignment(data_treatment)
  } else if (imputation_method == "enhanced_emi") {
    # 拡張EMI（複数代入の平均）
    emi_control <- enhanced_emi_imputation(data_control, n_imputations = 5)
    emi_treatment <- enhanced_emi_imputation(data_treatment, n_imputations = 5)

    # 複数代入結果の統合
    gpc_results <- list()
    for (m in 1:5) {
      gpc_m <- gpc_survival_improved(emi_control[[m]], emi_treatment[[m]], gpc_method)
      gpc_results[[m]] <- gpc_m
    }

    # Rubin's rulesで統合
    estimates <- sapply(gpc_results, function(x) x$estimate)
    p_values <- sapply(gpc_results, function(x) x$p_value)

    pooled_estimate <- mean(estimates)
    # Fisher's method for p-values
    if (all(p_values > 0)) {
      chi_sq <- -2 * sum(log(p_values))
      pooled_p <- 1 - pchisq(chi_sq, df = 2 * 5)
    } else {
      pooled_p <- 0
    }

    return(list(
      method = paste0("EnhancedEMI_", gpc_method),
      estimate = pooled_estimate,
      p_value = pooled_p
    ))
  }

  # 単一代入の場合
  return(gpc_survival_improved(surv_control, surv_treatment, gpc_method))
}

# 生存データに対するGPC
gpc_survival_improved <- function(surv_control, surv_treatment, method = "net_benefit") {
  n_control <- nrow(surv_control)
  n_treatment <- nrow(surv_treatment)

  wins <- 0
  losses <- 0
  ties <- 0

  for (i in 1:n_treatment) {
    for (j in 1:n_control) {
      t_treat <- surv_treatment$time[i]
      t_control <- surv_control$time[j]

      # 打ち切り情報の取得
      event_treat <- ifelse("cens" %in% names(surv_treatment),
                           surv_treatment$cens[i], surv_treatment$status[i])
      event_control <- ifelse("cens" %in% names(surv_control),
                             surv_control$cens[j], surv_control$status[j])

      # ペアワイズ比較
      if (event_treat == 1 && event_control == 1) {
        # 両方イベント発生：時刻を比較
        if (t_treat > t_control) {
          wins <- wins + 1
        } else if (t_treat < t_control) {
          losses <- losses + 1
        } else {
          ties <- ties + 1
        }
      } else if (event_treat == 1 && event_control == 0) {
        # 治療群のみイベント発生
        if (t_treat > t_control) {
          wins <- wins + 1
        } else {
          ties <- ties + 1
        }
      } else if (event_treat == 0 && event_control == 1) {
        # 対照群のみイベント発生
        if (t_treat > t_control) {
          losses <- losses + 1
        } else {
          ties <- ties + 1
        }
      } else {
        # 両方打ち切り
        ties <- ties + 1
      }
    }
  }

  total_pairs <- n_treatment * n_control

  if (method == "net_benefit") {
    net_benefit <- (wins - losses) / total_pairs
    variance <- (wins + losses + ties/4) / (total_pairs^2)
    if (variance <= 0) variance <- 1 / (4 * total_pairs)

    z_stat <- net_benefit / sqrt(variance)
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      method = "Imputation_NetBenefit",
      estimate = net_benefit,
      p_value = p_value
    ))
  } else if (method == "win_ratio") {
    if (losses == 0) {
      win_ratio <- ifelse(wins > 0, Inf, 1)
      p_value <- ifelse(wins > 0, 0, 1)
    } else if (wins == 0) {
      win_ratio <- 0
      p_value <- 0
    } else {
      win_ratio <- wins / losses
      log_var <- 1/wins + 1/losses
      z_stat <- log(win_ratio) / sqrt(log_var)
      p_value <- 2 * (1 - pnorm(abs(z_stat)))
    }

    return(list(
      method = "Imputation_WinRatio",
      estimate = win_ratio,
      p_value = p_value
    ))
  }
}

# ===============================================================================
# 3. 従来手法（RMST・ログランク検定）
# ===============================================================================

# RMST（改良版）
rmst_improved <- function(data_control, data_treatment, tau = 1, imputation_method = "midpoint") {
  tryCatch({
    if (imputation_method == "midpoint") {
      surv_control <- midpoint_assignment(data_control)
      surv_treatment <- midpoint_assignment(data_treatment)
    } else if (imputation_method == "rightpoint") {
      surv_control <- rightpoint_assignment(data_control)
      surv_treatment <- rightpoint_assignment(data_treatment)
    }

    # Mann-Whitney U検定で近似
    wilcox_result <- wilcox.test(surv_treatment$time, surv_control$time,
                                alternative = "two.sided")

    return(list(
      method = paste0("RMST_", imputation_method),
      estimate = median(surv_treatment$time) - median(surv_control$time),
      p_value = wilcox_result$p.value
    ))

  }, error = function(e) {
    return(list(method = paste0("RMST_", imputation_method), estimate = 0, p_value = 1))
  })
}

# ログランク検定（改良版）
logrank_improved <- function(data_control, data_treatment, imputation_method = "midpoint") {
  tryCatch({
    if (imputation_method == "midpoint") {
      surv_control <- midpoint_assignment(data_control)
      surv_treatment <- midpoint_assignment(data_treatment)
    } else if (imputation_method == "rightpoint") {
      surv_control <- rightpoint_assignment(data_control)
      surv_treatment <- rightpoint_assignment(data_treatment)
    }

    combined_data <- rbind(
      cbind(surv_control, group = 0),
      cbind(surv_treatment, group = 1)
    )

    survdiff_result <- survdiff(Surv(time, cens) ~ group, data = combined_data)
    p_value <- 1 - pchisq(survdiff_result$chisq, df = 1)

    return(list(
      method = paste0("LogRank_", imputation_method),
      estimate = sqrt(survdiff_result$chisq),
      p_value = p_value
    ))

  }, error = function(e) {
    return(list(method = paste0("LogRank_", imputation_method), estimate = 0, p_value = 1))
  })
}

# ===============================================================================
# 4. 最終版シミュレーション関数
# ===============================================================================

run_final_simulation <- function(
  n_sim = 1000,
  sample_sizes = c(100, 200),
  K_values = c(3, 5),
  dropout_levels = c("None", "Medium"),
  effect_sizes = c(1, 6, 12),
  alpha = 0.05,
  n_cores = NULL,  # NULLの場合は自動検出
  use_parallel = TRUE,
  show_progress = TRUE
) {

  research_design()  # 研究デザインの表示

  # コア数の設定
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  cat("使用コア数:", n_cores, "/ 利用可能コア数:", parallel::detectCores(), "\n")
  cat("並列処理:", ifelse(use_parallel && n_cores > 1, "有効", "無効"), "\n")

  conditions <- expand.grid(
    n = sample_sizes,
    K = K_values,
    dropout = dropout_levels,
    effect = effect_sizes,
    stringsAsFactors = FALSE
  )

  cat("シミュレーション条件数:", nrow(conditions), "\n")
  cat("総シミュレーション数:", nrow(conditions) * n_sim, "\n\n")

  # プログレスバーの初期化
  if (show_progress) {
    pb <- progress_bar$new(
      format = "  シミュレーション実行中 [:bar] :percent :current/:total ETA: :eta",
      total = nrow(conditions) * n_sim,
      clear = FALSE,
      width = 80
    )
  }

  results <- data.frame()

  # 並列処理の設定
  if (use_parallel && n_cores > 1) {
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl) # doParallelからdoSNOWに変更

    cat("並列処理クラスターを", n_cores, "コアで開始しました\n\n")

    # 必要なパッケージとソースをワーカーに送信
    clusterEvalQ(cl, {
      library(survival)
      library(tidyverse)
    })

    # 関数をワーカーに送信 (improved_gpc_with_packages.R内の関数も追加)
    clusterExport(cl, c("generate_interval_censored_data_2arm",
                       "gpc_direct_improved", "gpc_with_imputation_improved",
                       "rmst_improved", "logrank_improved", "midpoint_assignment",
                       "rightpoint_assignment", "enhanced_emi_imputation", "gpc_survival_improved"))
  }

  for (i in 1:nrow(conditions)) {
    n <- conditions$n[i]
    K <- conditions$K[i]
    dropout <- conditions$dropout[i]
    effect <- conditions$effect[i]

    cat(sprintf("条件 %d/%d: n=%d, K=%d, dropout=%s, effect=%d\n",
                i, nrow(conditions), n, K, dropout, effect))

    # doSNOW用の進捗レポート関数
    progress_reporter <- function(n) {
      if (show_progress) pb$tick()
    }

    if (use_parallel && n_cores > 1) {
      # 並列処理版
      opts <- list(progress = progress_reporter)
      sim_results <- foreach(sim = 1:n_sim, .combine = rbind, .options.snow = opts) %dopar% {
        tryCatch({
          # データ生成
          data_control <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                             p_dropout = dropout, d = 1)
          data_treatment <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                               p_dropout = dropout, d = effect)

          # 解析実行...
          gpc_direct_nb <- gpc_direct_improved(data_control, data_treatment, "net_benefit")
          gpc_direct_wr <- gpc_direct_improved(data_control, data_treatment, "win_ratio")
          gpc_mid_nb <- gpc_with_imputation_improved(data_control, data_treatment, "midpoint", "net_benefit")
          gpc_mid_wr <- gpc_with_imputation_improved(data_control, data_treatment, "midpoint", "win_ratio")
          gpc_right_nb <- gpc_with_imputation_improved(data_control, data_treatment, "rightpoint", "net_benefit")
          gpc_right_wr <- gpc_with_imputation_improved(data_control, data_treatment, "rightpoint", "win_ratio")
          gpc_emi_nb <- gpc_with_imputation_improved(data_control, data_treatment, "enhanced_emi", "net_benefit")
          gpc_emi_wr <- gpc_with_imputation_improved(data_control, data_treatment, "enhanced_emi", "win_ratio")
          rmst_mid <- rmst_improved(data_control, data_treatment, imputation_method = "midpoint")
          lr_mid <- logrank_improved(data_control, data_treatment, imputation_method = "midpoint")

          # 結果をまとめる
          data.frame(
            sim_id = sim,
            n = n, K = K, dropout = dropout, effect = effect,
            gpc_direct_nb_reject = gpc_direct_nb$p_value < alpha,
            gpc_direct_wr_reject = gpc_direct_wr$p_value < alpha,
            gpc_mid_nb_reject = gpc_mid_nb$p_value < alpha,
            gpc_mid_wr_reject = gpc_mid_wr$p_value < alpha,
            gpc_right_nb_reject = gpc_right_nb$p_value < alpha,
            gpc_right_wr_reject = gpc_right_wr$p_value < alpha,
            gpc_emi_nb_reject = gpc_emi_nb$p_value < alpha,
            gpc_emi_wr_reject = gpc_emi_wr$p_value < alpha,
            rmst_mid_reject = rmst_mid$p_value < alpha,
            lr_mid_reject = lr_mid$p_value < alpha,
            gpc_direct_nb_pvalue = gpc_direct_nb$p_value,
            gpc_direct_wr_pvalue = gpc_direct_wr$p_value,
            gpc_mid_nb_pvalue = gpc_mid_nb$p_value,
            gpc_mid_wr_pvalue = gpc_mid_wr$p_value,
            gpc_right_nb_pvalue = gpc_right_nb$p_value,
            gpc_right_wr_pvalue = gpc_right_wr$p_value,
            gpc_emi_nb_pvalue = gpc_emi_nb$p_value,
            gpc_emi_wr_pvalue = gpc_emi_wr$p_value,
            rmst_mid_pvalue = rmst_mid$p_value,
            lr_mid_pvalue = lr_mid$p_value
          )

        }, error = function(e) {
          NULL
        })
      }

    } else {
      # 逐次処理版
      sim_results_list <- list()
      # 逐次処理でも必要なファイルを読み込む
      source("improved_gpc_with_packages.R")
      for (sim in 1:n_sim) {
        tryCatch({
          # データ生成
          data_control <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                             p_dropout = dropout, d = 1)
          data_treatment <- generate_interval_censored_data_2arm(n = n/2, K = K,
                                                               p_dropout = dropout, d = effect)

          # 解析実行...
          gpc_direct_nb <- gpc_direct_improved(data_control, data_treatment, "net_benefit")
          gpc_direct_wr <- gpc_direct_improved(data_control, data_treatment, "win_ratio")
          gpc_mid_nb <- gpc_with_imputation_improved(data_control, data_treatment, "midpoint", "net_benefit")
          gpc_mid_wr <- gpc_with_imputation_improved(data_control, data_treatment, "midpoint", "win_ratio")
          gpc_right_nb <- gpc_with_imputation_improved(data_control, data_treatment, "rightpoint", "net_benefit")
          gpc_right_wr <- gpc_with_imputation_improved(data_control, data_treatment, "rightpoint", "win_ratio")
          gpc_emi_nb <- gpc_with_imputation_improved(data_control, data_treatment, "enhanced_emi", "net_benefit")
          gpc_emi_wr <- gpc_with_imputation_improved(data_control, data_treatment, "enhanced_emi", "win_ratio")
          rmst_mid <- rmst_improved(data_control, data_treatment, imputation_method = "midpoint")
          lr_mid <- logrank_improved(data_control, data_treatment, imputation_method = "midpoint")

          # 結果をまとめる
          sim_result <- data.frame(
            sim_id = sim,
            n = n, K = K, dropout = dropout, effect = effect,
            gpc_direct_nb_reject = gpc_direct_nb$p_value < alpha,
            gpc_direct_wr_reject = gpc_direct_wr$p_value < alpha,
            gpc_mid_nb_reject = gpc_mid_nb$p_value < alpha,
            gpc_mid_wr_reject = gpc_mid_wr$p_value < alpha,
            gpc_right_nb_reject = gpc_right_nb$p_value < alpha,
            gpc_right_wr_reject = gpc_right_wr$p_value < alpha,
            gpc_emi_nb_reject = gpc_emi_nb$p_value < alpha,
            gpc_emi_wr_reject = gpc_emi_wr$p_value < alpha,
            rmst_mid_reject = rmst_mid$p_value < alpha,
            lr_mid_reject = lr_mid$p_value < alpha,
            gpc_direct_nb_pvalue = gpc_direct_nb$p_value,
            gpc_direct_wr_pvalue = gpc_direct_wr$p_value,
            gpc_mid_nb_pvalue = gpc_mid_nb$p_value,
            gpc_mid_wr_pvalue = gpc_mid_wr$p_value,
            gpc_right_nb_pvalue = gpc_right_nb$p_value,
            gpc_right_wr_pvalue = gpc_right_wr$p_value,
            gpc_emi_nb_pvalue = gpc_emi_nb$p_value,
            gpc_emi_wr_pvalue = gpc_emi_wr$p_value,
            rmst_mid_pvalue = rmst_mid$p_value,
            lr_mid_pvalue = lr_mid$p_value
          )
          sim_results_list[[sim]] <- sim_result

          # プログレスバー更新
          if (show_progress) {
            pb$tick()
          }

        }, error = function(e) {
          cat("エラー (sim", sim, "):", e$message, "\n")
          # プログレスバー更新
          if (show_progress) {
            pb$tick()
          }
        })
      }
      sim_results <- do.call(rbind, sim_results_list)
    }

    results <- rbind(results, sim_results)

    # 進捗表示
    if (!is.null(sim_results) && nrow(sim_results) > 0) {
      avg_power <- mean(sim_results$gpc_direct_nb_reject, na.rm = TRUE)
      cat(sprintf("  -> 完了 (%d シミュレーション): 平均検出力 = %.3f\n",
                  nrow(sim_results), avg_power))
    }
  }

  # 並列処理クラスターのクリーンアップ
  if (use_parallel && n_cores > 1) {
    stopCluster(cl)
    cat("\n並列処理クラスターを終了しました\n")
  }

  return(results)
}

# ===============================================================================
# 5. 結果分析関数
# ===============================================================================

analyze_final_results <- function(results) {
  # 検出力の計算
  power_summary <- results %>%
    group_by(n, K, dropout, effect) %>%
    summarise(
      # GPC各手法の検出力
      GPC_Direct_NB = mean(gpc_direct_nb_reject, na.rm = TRUE),
      GPC_Direct_WR = mean(gpc_direct_wr_reject, na.rm = TRUE),
      GPC_Midpoint_NB = mean(gpc_mid_nb_reject, na.rm = TRUE),
      GPC_Midpoint_WR = mean(gpc_mid_wr_reject, na.rm = TRUE),
      GPC_Rightpoint_NB = mean(gpc_right_nb_reject, na.rm = TRUE),
      GPC_Rightpoint_WR = mean(gpc_right_wr_reject, na.rm = TRUE),
      GPC_EnhancedEMI_NB = mean(gpc_emi_nb_reject, na.rm = TRUE),
      GPC_EnhancedEMI_WR = mean(gpc_emi_wr_reject, na.rm = TRUE),

      # 従来手法の検出力
      RMST_Midpoint = mean(rmst_mid_reject, na.rm = TRUE),
      LogRank_Midpoint = mean(lr_mid_reject, na.rm = TRUE),

      .groups = "drop"
    )

  # 第1種の誤りの計算
  type1_summary <- results %>%
    filter(effect == 1) %>%
    group_by(n, K, dropout) %>%
    summarise(
      GPC_Direct_NB = mean(gpc_direct_nb_reject, na.rm = TRUE),
      GPC_Direct_WR = mean(gpc_direct_wr_reject, na.rm = TRUE),
      GPC_Midpoint_NB = mean(gpc_mid_nb_reject, na.rm = TRUE),
      GPC_Midpoint_WR = mean(gpc_mid_wr_reject, na.rm = TRUE),
      GPC_Rightpoint_NB = mean(gpc_right_nb_reject, na.rm = TRUE),
      GPC_Rightpoint_WR = mean(gpc_right_wr_reject, na.rm = TRUE),
      GPC_EnhancedEMI_NB = mean(gpc_emi_nb_reject, na.rm = TRUE),
      GPC_EnhancedEMI_WR = mean(gpc_emi_wr_reject, na.rm = TRUE),
      RMST_Midpoint = mean(rmst_mid_reject, na.rm = TRUE),
      LogRank_Midpoint = mean(lr_mid_reject, na.rm = TRUE),
      .groups = "drop"
    )

  # 主要結果の表示
  cat("===============================================================================\n")
  cat("                             研究結果サマリー\n")
  cat("===============================================================================\n")

  # Q1: GPCにおける最良の代入法
  gpc_power_avg <- power_summary %>%
    filter(effect != 1) %>%
    summarise(
      Direct_NB = mean(GPC_Direct_NB),
      Direct_WR = mean(GPC_Direct_WR),
      Midpoint_NB = mean(GPC_Midpoint_NB),
      Midpoint_WR = mean(GPC_Midpoint_WR),
      Rightpoint_NB = mean(GPC_Rightpoint_NB),
      Rightpoint_WR = mean(GPC_Rightpoint_WR),
      EnhancedEMI_NB = mean(GPC_EnhancedEMI_NB),
      EnhancedEMI_WR = mean(GPC_EnhancedEMI_WR)
    )

  cat("Q1: GPCにおける代入法の検出力比較 (平均):\n")
  print(round(gpc_power_avg, 3))

  # Q2: GPC vs 従来手法
  method_comparison <- power_summary %>%
    filter(effect != 1) %>%
    summarise(
      Best_GPC = max(c(mean(GPC_Direct_NB), mean(GPC_Midpoint_NB),
                       mean(GPC_Rightpoint_NB), mean(GPC_EnhancedEMI_NB))),
      RMST = mean(RMST_Midpoint),
      LogRank = mean(LogRank_Midpoint)
    )

  cat("\nQ2: GPC vs 従来手法の検出力比較:\n")
  print(round(method_comparison, 3))

  # 第1種の誤りチェック
  cat("\n第1種の誤り率 (目標: 0.05前後):\n")
  type1_avg <- type1_summary %>%
    summarise_if(is.numeric, mean)
  print(round(type1_avg, 3))

  return(list(
    power_summary = power_summary,
    type1_summary = type1_summary,
    gpc_comparison = gpc_power_avg,
    method_comparison = method_comparison
  ))
}

# ===============================================================================
# 6. 実行関数
# ===============================================================================

run_final_study <- function() {
  cat("=== 最終版：代入法比較研究の実行 ===\n")

  # 小規模テスト
  test_results <- run_final_simulation(
    n_sim = 100,
    sample_sizes = 100,
    K_values = 3,
    dropout_levels = "None",
    effect_sizes = c(1, 6),
    n_cores = 1
  )

  # 結果分析
  analysis <- analyze_final_results(test_results)

  # 結果をCSVファイルに出力
  cat("\n=== 結果をCSVファイルに出力しています ===\n")
  # resultsディレクトリがなければ作成
  if (!dir.exists("results")) {
    dir.create("results")
  }
  write.csv(test_results, "results/simulation_results.csv", row.names = FALSE)
  write.csv(analysis$power_summary, "results/power_summary.csv", row.names = FALSE)
  write.csv(analysis$type1_summary, "results/type1_error_summary.csv", row.names = FALSE)
  cat("結果は 'results' フォルダに保存されました。\n")

  return(list(
    simulation_results = test_results,
    analysis = analysis
  ))
}

cat("===============================================================================\n")
cat("         最終版：一般化ペアワイズ比較における代入法比較研究\n")
cat("===============================================================================\n")
cat("実行方法:\n")
cat("1. final_results <- run_final_study()  # テスト実行\n")
cat("2. 結果の確認とp値の妥当性をチェック\n")
cat("3. 必要に応じてパラメータ調整\n")
cat("===============================================================================\n")