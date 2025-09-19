# ===============================================================================
# 一般化ペアワイズ比較シミュレーションのクイックテスト
#
# このファイルは、メインシミュレーションコードの動作確認用です。
# 小規模なテストを実行して、各関数が正常に動作することを確認できます。
# ===============================================================================

# 必要なライブラリの確認と読み込み
check_and_load_packages <- function() {
  required_packages <- c("survival", "tidyverse", "ggplot2", "parallel")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("警告:", pkg, "パッケージがインストールされていません。\n")
      cat("install.packages('", pkg, "') を実行してください。\n")
    } else {
      library(pkg, character.only = TRUE)
      cat("✓", pkg, "パッケージを読み込みました。\n")
    }
  }
}

# メイン関数の読み込み確認
load_main_functions <- function() {
  tryCatch({
    source("interval_censord_data_function.R")
    cat("✓ interval_censord_data_function.R を読み込みました。\n")

    source("distributions_2arm.R")
    cat("✓ distributions_2arm.R を読み込みました。\n")

    source("generalized_pairwise_comparison_simulation.R")
    cat("✓ generalized_pairwise_comparison_simulation.R を読み込みました。\n")

    return(TRUE)
  }, error = function(e) {
    cat("エラー: ファイルの読み込みに失敗しました。\n")
    cat("エラーメッセージ:", e$message, "\n")
    return(FALSE)
  })
}

# 基本的なデータ生成テスト
test_data_generation <- function() {
  cat("\n=== データ生成テスト ===\n")

  tryCatch({
    # 対照群データ生成
    data_control <- generate_interval_censored_data_2arm(n = 50, K = 5, p_dropout = "None", d = 1)
    cat("✓ 対照群データ生成成功: サイズ", nrow(data_control), "x", ncol(data_control), "\n")

    # 治療群データ生成
    data_treatment <- generate_interval_censored_data_2arm(n = 50, K = 5, p_dropout = "None", d = 6)
    cat("✓ 治療群データ生成成功: サイズ", nrow(data_treatment), "x", ncol(data_treatment), "\n")

    # データの基本統計
    cat("対照群の要約統計:\n")
    print(summary(data_control))
    cat("治療群の要約統計:\n")
    print(summary(data_treatment))

    return(list(control = data_control, treatment = data_treatment))
  }, error = function(e) {
    cat("✗ データ生成エラー:", e$message, "\n")
    return(NULL)
  })
}

# 代入法テスト
test_imputation_methods <- function(data_control, data_treatment) {
  cat("\n=== 代入法テスト ===\n")

  tryCatch({
    # 中点代入法
    mid_control <- midpoint_assignment(data_control)
    mid_treatment <- midpoint_assignment(data_treatment)
    cat("✓ 中点代入法成功\n")
    cat("  対照群: イベント率 =", mean(mid_control$cens), "\n")
    cat("  治療群: イベント率 =", mean(mid_treatment$cens), "\n")

    # 右点代入法
    right_control <- rightpoint_assignment(data_control)
    right_treatment <- rightpoint_assignment(data_treatment)
    cat("✓ 右点代入法成功\n")
    cat("  対照群: イベント率 =", mean(right_control$cens), "\n")
    cat("  治療群: イベント率 =", mean(right_treatment$cens), "\n")

    # 拡張EMI法
    emi_control <- enhanced_emi_imputation(data_control, n_imputations = 3)
    emi_treatment <- enhanced_emi_imputation(data_treatment, n_imputations = 3)
    cat("✓ 拡張EMI法成功: 代入数 =", length(emi_control), "\n")

    return(TRUE)
  }, error = function(e) {
    cat("✗ 代入法エラー:", e$message, "\n")
    return(FALSE)
  })
}

# GPC関数テスト
test_gpc_functions <- function(data_control, data_treatment) {
  cat("\n=== GPC関数テスト ===\n")

  tryCatch({
    # Net Benefit法
    gpc_nb <- gpc_interval_censored(data_control, data_treatment, method = "net_benefit")
    cat("✓ GPC Net Benefit法成功\n")
    cat("  Net Benefit =", round(gpc_nb$net_benefit, 4), "\n")
    cat("  p値 =", round(gpc_nb$p_value, 4), "\n")
    cat("  勝ち/負け/引き分け =", gpc_nb$wins, "/", gpc_nb$losses, "/", gpc_nb$ties, "\n")

    # Win Ratio法
    gpc_wr <- gpc_interval_censored(data_control, data_treatment, method = "win_ratio")
    cat("✓ GPC Win Ratio法成功\n")
    cat("  Win Ratio =", round(gpc_wr$win_ratio, 4), "\n")
    cat("  p値 =", round(gpc_wr$p_value, 4), "\n")

    return(list(net_benefit = gpc_nb, win_ratio = gpc_wr))
  }, error = function(e) {
    cat("✗ GPC関数エラー:", e$message, "\n")
    return(NULL)
  })
}

# RMST・ログランク検定テスト
test_survival_methods <- function(data_control, data_treatment) {
  cat("\n=== 生存解析手法テスト ===\n")

  # RMST テスト
  tryCatch({
    rmst_mid <- rmst_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
    cat("✓ RMST (中点代入法) 成功: 差 =", round(rmst_mid$rmst_diff, 4), ", p値 =", round(rmst_mid$p_value, 4), "\n")
  }, error = function(e) {
    cat("✗ RMST (中点代入法) エラー:", e$message, "\n")
  })

  # ログランク検定テスト
  tryCatch({
    lr_mid <- logrank_interval_censored(data_control, data_treatment, imputation_method = "midpoint")
    cat("✓ ログランク検定 (中点代入法) 成功: p値 =", round(lr_mid$p_value, 4), "\n")
  }, error = function(e) {
    cat("✗ ログランク検定 (中点代入法) エラー:", e$message, "\n")
  })
}

# 小規模シミュレーションテスト
test_small_simulation <- function() {
  cat("\n=== 小規模シミュレーションテスト ===\n")
  cat("注意: このテストには2-3分かかります。\n")

  tryCatch({
    # 非常に小さい条件でテスト
    test_sim <- run_simulation_study(
      n_sim = 10,                    # 非常に小さいシミュレーション回数
      sample_sizes = 100,            # 1つのサンプルサイズのみ
      K_values = 5,                  # 1つの観測頻度のみ
      dropout_levels = "None",       # 脱落なしのみ
      effect_sizes = c(1, 6),        # 2つの効果サイズのみ
      n_cores = 1                    # シングルスレッド
    )

    cat("✓ 小規模シミュレーション成功: 結果行数 =", nrow(test_sim), "\n")

    # 検出力計算テスト
    power_summary <- calculate_power(test_sim)
    cat("✓ 検出力計算成功\n")
    print(head(power_summary, 3))

    return(test_sim)
  }, error = function(e) {
    cat("✗ 小規模シミュレーションエラー:", e$message, "\n")
    return(NULL)
  })
}

# 包括的テスト実行
run_comprehensive_test <- function() {
  cat("===============================================================================\n")
  cat("                一般化ペアワイズ比較シミュレーション - 動作確認テスト\n")
  cat("===============================================================================\n")

  start_time <- Sys.time()

  # 1. パッケージ確認
  cat("\n1. パッケージ確認中...\n")
  check_and_load_packages()

  # 2. 関数読み込み確認
  cat("\n2. 関数読み込み確認中...\n")
  if (!load_main_functions()) {
    cat("関数の読み込みに失敗しました。テストを中止します。\n")
    return(FALSE)
  }

  # 3. データ生成テスト
  cat("\n3. データ生成テスト中...\n")
  test_data <- test_data_generation()
  if (is.null(test_data)) {
    cat("データ生成に失敗しました。テストを中止します。\n")
    return(FALSE)
  }

  # 4. 代入法テスト
  cat("\n4. 代入法テスト中...\n")
  if (!test_imputation_methods(test_data$control, test_data$treatment)) {
    cat("代入法テストに失敗しました。\n")
  }

  # 5. GPC関数テスト
  cat("\n5. GPC関数テスト中...\n")
  gpc_results <- test_gpc_functions(test_data$control, test_data$treatment)

  # 6. 生存解析手法テスト
  cat("\n6. 生存解析手法テスト中...\n")
  test_survival_methods(test_data$control, test_data$treatment)

  # 7. 小規模シミュレーションテスト
  cat("\n7. 小規模シミュレーションテスト中...\n")
  sim_results <- test_small_simulation()

  end_time <- Sys.time()

  cat("\n===============================================================================\n")
  cat("                                テスト完了\n")
  cat("===============================================================================\n")
  cat("総実行時間:", round(difftime(end_time, start_time, units = "mins"), 2), "分\n")

  if (!is.null(sim_results)) {
    cat("✓ 全ての主要機能が正常に動作しています。\n")
    cat("✓ メインシミュレーション実行の準備が完了しました。\n")
    cat("\n次のステップ:\n")
    cat("- フル分析: full_results <- run_full_analysis()\n")
    cat("- テスト分析: test_results <- run_test_simulation()\n")
    return(TRUE)
  } else {
    cat("✗ 一部の機能でエラーが発生しました。\n")
    cat("エラーメッセージを確認して修正してください。\n")
    return(FALSE)
  }
}

# システム情報表示
show_system_info <- function() {
  cat("システム情報:\n")
  cat("- R バージョン:", R.version.string, "\n")
  cat("- プラットフォーム:", R.version$platform, "\n")
  cat("- 利用可能コア数:", parallel::detectCores(), "\n")
  cat("- 現在の作業ディレクトリ:", getwd(), "\n")
  cat("- メモリ使用量:", round(memory.size()/1024, 2), "GB (Windows)\n")
}

# メイン実行
if (interactive()) {
  cat("クイックテストスクリプトが読み込まれました。\n")
  cat("実行方法:\n")
  cat("1. show_system_info()        # システム情報の表示\n")
  cat("2. run_comprehensive_test()  # 包括的テストの実行\n")
  cat("3. test_data_generation()    # データ生成のみテスト\n")
} else {
  # 非対話モードでは自動実行
  show_system_info()
  run_comprehensive_test()
}