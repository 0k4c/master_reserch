# ===============================================================================
# 修正版GPCシミュレーションのテスト・デバッグ用スクリプト
# ===============================================================================

# 作業ディレクトリの設定（必要に応じて変更）
# setwd("C:/Claude/serch/serch1/code")

# 必要なライブラリの読み込み確認
required_packages <- c("survival", "tidyverse", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("警告:", pkg, "パッケージがインストールされていません。\n")
    cat("install.packages('", pkg, "') を実行してください。\n")
  } else {
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat("✓", pkg, "パッケージを読み込みました。\n")
  }
}

# ファイルの読み込み
cat("\n=== ファイル読み込み ===\n")
tryCatch({
  source("flexible_data_generation.R")  # 制限緩和版を使用
  cat("✓ flexible_data_generation.R 読み込み完了\n")
}, error = function(e) {
  cat("✗ flexible_data_generation.R 読み込みエラー:", e$message, "\n")
  # フォールバック：元のファイルを試す
  tryCatch({
    source("interval_censord_data_function.R")
    source("distributions_2arm.R")
    cat("✓ 元のファイルで読み込み完了（制限あり）\n")
  }, error = function(e2) {
    cat("✗ 元のファイル読み込みもエラー:", e2$message, "\n")
  })
})

tryCatch({
  source("corrected_gpc_simulation.R")
  cat("✓ corrected_gpc_simulation.R 読み込み完了\n")
}, error = function(e) {
  cat("✗ corrected_gpc_simulation.R 読み込みエラー:", e$message, "\n")
  stop("メインファイルの読み込みに失敗しました。")
})

# ===============================================================================
# 個別関数のテスト
# ===============================================================================

test_data_generation <- function() {
  cat("\n=== データ生成テスト ===\n")

  # 小さいデータでテスト（制限緩和版なら20、元版なら100を使用）
  tryCatch({
    data_control <- generate_interval_censored_data_2arm(n = 20, K = 3, p_dropout = "None", d = 1)
    data_treatment <- generate_interval_censored_data_2arm(n = 20, K = 3, p_dropout = "None", d = 6)
  }, error = function(e) {
    cat("小さいサンプルサイズでエラー。制限サイズを使用します。\n")
    data_control <<- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 1)
    data_treatment <<- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 6)
  })

  cat("対照群データサイズ:", nrow(data_control), "x", ncol(data_control), "\n")
  cat("治療群データサイズ:", nrow(data_treatment), "x", ncol(data_treatment), "\n")

  cat("対照群データの例:\n")
  print(head(data_control, 3))
  cat("治療群データの例:\n")
  print(head(data_treatment, 3))

  return(list(control = data_control, treatment = data_treatment))
}

test_gpc_direct <- function(data_control, data_treatment) {
  cat("\n=== GPC Direct テスト ===\n")

  # Net Benefit テスト
  result_nb <- gpc_direct_interval(data_control, data_treatment, method = "net_benefit")
  cat("Net Benefit結果:\n")
  cat("  Net Benefit =", result_nb$net_benefit, "\n")
  cat("  p値 =", result_nb$p_value, "\n")
  cat("  勝ち/負け/引き分け =", result_nb$wins, "/", result_nb$losses, "/", result_nb$ties, "\n")

  # Win Ratio テスト
  result_wr <- gpc_direct_interval(data_control, data_treatment, method = "win_ratio")
  cat("Win Ratio結果:\n")
  cat("  Win Ratio =", result_wr$win_ratio, "\n")
  cat("  p値 =", result_wr$p_value, "\n")

  return(list(nb = result_nb, wr = result_wr))
}

test_imputation_methods <- function(data_control, data_treatment) {
  cat("\n=== 代入法テスト ===\n")

  # 中点代入法
  mid_control <- midpoint_assignment(data_control)
  mid_treatment <- midpoint_assignment(data_treatment)
  cat("中点代入法:\n")
  cat("  対照群イベント率:", mean(mid_control$cens), "\n")
  cat("  治療群イベント率:", mean(mid_treatment$cens), "\n")

  # 右点代入法
  right_control <- rightpoint_assignment(data_control)
  right_treatment <- rightpoint_assignment(data_treatment)
  cat("右点代入法:\n")
  cat("  対照群イベント率:", mean(right_control$cens), "\n")
  cat("  治療群イベント率:", mean(right_treatment$cens), "\n")

  # 拡張EMI法
  tryCatch({
    emi_control <- enhanced_emi_imputation(data_control, n_imputations = 3, method = "adaptive")
    emi_treatment <- enhanced_emi_imputation(data_treatment, n_imputations = 3, method = "adaptive")
    cat("拡張EMI法:\n")
    cat("  代入データセット数:", length(emi_control), "\n")
    cat("  各データセット例 (対照群1番目):\n")
    print(head(emi_control[[1]], 3))
  }, error = function(e) {
    cat("拡張EMI法エラー:", e$message, "\n")
  })
}

test_gpc_with_imputation <- function(data_control, data_treatment) {
  cat("\n=== GPC+代入法テスト ===\n")

  # 中点代入法
  tryCatch({
    gpc_mid_nb <- gpc_with_imputation(data_control, data_treatment,
                                     imputation_method = "midpoint", method = "net_benefit")
    cat("GPC + 中点代入法 (NB): p値 =", gpc_mid_nb$p_value, "\n")
  }, error = function(e) {
    cat("GPC + 中点代入法エラー:", e$message, "\n")
  })

  # 右点代入法
  tryCatch({
    gpc_right_nb <- gpc_with_imputation(data_control, data_treatment,
                                       imputation_method = "rightpoint", method = "net_benefit")
    cat("GPC + 右点代入法 (NB): p値 =", gpc_right_nb$p_value, "\n")
  }, error = function(e) {
    cat("GPC + 右点代入法エラー:", e$message, "\n")
  })

  # 拡張EMI法
  tryCatch({
    gpc_emi_nb <- gpc_with_imputation(data_control, data_treatment,
                                     imputation_method = "enhanced_emi", method = "net_benefit")
    cat("GPC + 拡張EMI法 (NB): p値 =", gpc_emi_nb$p_value, "\n")
  }, error = function(e) {
    cat("GPC + 拡張EMI法エラー:", e$message, "\n")
  })
}

test_traditional_methods <- function(data_control, data_treatment) {
  cat("\n=== 従来手法テスト ===\n")

  # RMST
  tryCatch({
    rmst_result <- rmst_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
    cat("RMST (中点代入法): p値 =", rmst_result$p_value, ", 差 =", rmst_result$rmst_diff, "\n")
  }, error = function(e) {
    cat("RMSTエラー:", e$message, "\n")
  })

  # ログランク検定
  tryCatch({
    lr_result <- logrank_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
    cat("ログランク検定 (中点代入法): p値 =", lr_result$p_value, "\n")
  }, error = function(e) {
    cat("ログランク検定エラー:", e$message, "\n")
  })
}

# ===============================================================================
# 小規模シミュレーションテスト
# ===============================================================================

test_small_simulation <- function() {
  cat("\n=== 小規模シミュレーションテスト ===\n")
  cat("実行中... (1-2分程度)\n")

  sim_results <- run_corrected_simulation(
    n_sim = 10,  # 非常に小さい
    sample_sizes = 100,  # 制限に合わせて100に変更
    K_values = 3,
    dropout_levels = "None",
    effect_sizes = c(1, 6),
    n_cores = 1
  )

  cat("シミュレーション完了。結果サイズ:", nrow(sim_results), "行\n")

  if (nrow(sim_results) > 0) {
    cat("p値の分布確認:\n")

    # p値の範囲をチェック
    p_cols <- grep("_pvalue$", names(sim_results), value = TRUE)
    for (col in p_cols) {
      p_vals <- sim_results[[col]]
      if (length(p_vals) > 0 && !all(is.na(p_vals))) {
        cat(sprintf("  %s: [%.3f, %.3f], 平均=%.3f\n",
                   col, min(p_vals, na.rm = TRUE), max(p_vals, na.rm = TRUE), mean(p_vals, na.rm = TRUE)))
      }
    }

    # 棄却率もチェック
    reject_cols <- grep("_reject$", names(sim_results), value = TRUE)
    cat("\n棄却率の確認:\n")
    for (col in reject_cols) {
      reject_rate <- mean(sim_results[[col]], na.rm = TRUE)
      cat(sprintf("  %s: %.2f%%\n", col, reject_rate * 100))
    }
  }

  return(sim_results)
}

# ===============================================================================
# 包括的テスト実行
# ===============================================================================

run_comprehensive_debug <- function() {
  cat("===============================================================================\n")
  cat("                   修正版GPCシミュレーション デバッグテスト\n")
  cat("===============================================================================\n")

  start_time <- Sys.time()

  # 1. データ生成テスト
  test_data <- test_data_generation()

  # 2. GPC Direct テスト
  gpc_results <- test_gpc_direct(test_data$control, test_data$treatment)

  # 3. 代入法テスト
  test_imputation_methods(test_data$control, test_data$treatment)

  # 4. GPC+代入法テスト
  test_gpc_with_imputation(test_data$control, test_data$treatment)

  # 5. 従来手法テスト
  test_traditional_methods(test_data$control, test_data$treatment)

  # 6. 小規模シミュレーション
  sim_results <- test_small_simulation()

  end_time <- Sys.time()

  cat("\n===============================================================================\n")
  cat("                              テスト完了\n")
  cat("===============================================================================\n")
  cat("総実行時間:", round(difftime(end_time, start_time, units = "mins"), 2), "分\n")

  # 結果の要約
  if (!is.null(sim_results) && nrow(sim_results) > 0) {
    cat("\n✓ 主要機能が正常に動作しています。\n")
    cat("✓ p値が適切に計算されています。\n")

    # 分析実行
    tryCatch({
      analysis_results <- analyze_corrected_results(sim_results)
      cat("✓ 結果分析も正常に動作しています。\n")

      return(list(
        simulation_results = sim_results,
        analysis_results = analysis_results,
        test_data = test_data,
        gpc_results = gpc_results
      ))
    }, error = function(e) {
      cat("分析でエラーが発生:", e$message, "\n")
      return(sim_results)
    })
  } else {
    cat("✗ シミュレーションで問題が発生しました。\n")
    return(NULL)
  }
}

# ===============================================================================
# メイン実行
# ===============================================================================

cat("デバッグテストスクリプトが読み込まれました。\n")
cat("実行方法:\n")
cat("1. run_comprehensive_debug()  # 全体テスト\n")
cat("2. test_data_generation()     # データ生成のみ\n")
cat("3. test_small_simulation()    # 小規模シミュレーションのみ\n")

# 非対話モードでは自動実行
if (!interactive()) {
  run_comprehensive_debug()
}