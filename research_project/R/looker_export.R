# ===============================================================================
# Google Looker Studio用データ出力機能
#
# このファイルは、シミュレーション結果をGoogle Looker Studioで
# 可視化できる形式のCSVファイルとして出力する機能を提供します。
# ===============================================================================

library(tidyverse)

# ===============================================================================
# 1. Looker Studio用データ変換関数
# ===============================================================================

# メインのデータ変換関数
transform_for_looker_studio <- function(simulation_results, analysis_results = NULL) {
  cat("=== Google Looker Studio用データ変換 ===\n")

  # 基本的な結果データの準備
  main_data <- simulation_results %>%
    # long format に変換（各手法を行として展開）
    pivot_longer(
      cols = ends_with("_reject"),
      names_to = "method_type",
      values_to = "rejected"
    ) %>%
    # p値も同様に変換
    pivot_longer(
      cols = ends_with("_pvalue"),
      names_to = "method_pvalue_type",
      values_to = "p_value",
      names_prefix = "",
      values_drop_na = TRUE
    ) %>%
    # メソッド名のマッチング
    filter(
      str_remove(method_type, "_reject") == str_remove(method_pvalue_type, "_pvalue")
    ) %>%
    # 不要な列を削除
    select(-method_pvalue_type) %>%
    # メソッド名をクリーンアップ
    mutate(
      method = str_remove(method_type, "_reject"),
      method = str_replace_all(method, "_", " "),
      method = str_to_title(method)
    ) %>%
    # カテゴリ分類
    mutate(
      method_category = case_when(
        str_detect(method, "Gpc.*Nb") ~ "GPC_NetBenefit",
        str_detect(method, "Gpc.*Wr") ~ "GPC_WinRatio",
        str_detect(method, "Rmst") ~ "RMST",
        str_detect(method, "Lr") ~ "LogRank",
        TRUE ~ "Other"
      ),
      imputation_method = case_when(
        str_detect(method, "Direct") ~ "Direct",
        str_detect(method, "Mid") ~ "Midpoint",
        str_detect(method, "Right") ~ "Rightpoint",
        str_detect(method, "Emi") ~ "Enhanced_EMI",
        TRUE ~ "Traditional"
      ),
      is_gpc_method = str_detect(method_category, "GPC"),
      effect_type = case_when(
        effect == 1 ~ "Null_Effect",
        effect <= 5 ~ "Small_Effect",
        effect <= 12 ~ "Medium_Effect",
        TRUE ~ "Large_Effect"
      ),
      sample_size_category = case_when(
        n <= 100 ~ "Small",
        n <= 200 ~ "Medium",
        TRUE ~ "Large"
      )
    )

  return(main_data)
}

# 検出力サマリーデータの作成
create_power_summary_for_looker <- function(simulation_results) {
  power_data <- simulation_results %>%
    group_by(n, K, dropout, effect) %>%
    summarise(
      # GPC手法
      GPC_Direct_NB = mean(gpc_direct_nb_reject, na.rm = TRUE),
      GPC_Direct_WR = mean(gpc_direct_wr_reject, na.rm = TRUE),
      GPC_Midpoint_NB = mean(gpc_mid_nb_reject, na.rm = TRUE),
      GPC_Midpoint_WR = mean(gpc_mid_wr_reject, na.rm = TRUE),
      GPC_Rightpoint_NB = mean(gpc_right_nb_reject, na.rm = TRUE),
      GPC_Rightpoint_WR = mean(gpc_right_wr_reject, na.rm = TRUE),
      GPC_EnhancedEMI_NB = mean(gpc_emi_nb_reject, na.rm = TRUE),
      GPC_EnhancedEMI_WR = mean(gpc_emi_wr_reject, na.rm = TRUE),

      # 従来手法
      RMST_Midpoint = mean(rmst_mid_reject, na.rm = TRUE),
      LogRank_Midpoint = mean(lr_mid_reject, na.rm = TRUE),

      # 実行されたシミュレーション数
      n_simulations = n(),
      .groups = "drop"
    ) %>%
    # long format に変換
    pivot_longer(
      cols = c(GPC_Direct_NB:LogRank_Midpoint),
      names_to = "method",
      values_to = "power"
    ) %>%
    # メタデータの追加
    mutate(
      method_category = case_when(
        str_detect(method, "GPC.*NB") ~ "GPC_NetBenefit",
        str_detect(method, "GPC.*WR") ~ "GPC_WinRatio",
        str_detect(method, "RMST") ~ "RMST",
        str_detect(method, "LogRank") ~ "LogRank",
        TRUE ~ "Other"
      ),
      imputation_type = case_when(
        str_detect(method, "Direct") ~ "Direct",
        str_detect(method, "Midpoint") ~ "Midpoint",
        str_detect(method, "Rightpoint") ~ "Rightpoint",
        str_detect(method, "EnhancedEMI") ~ "Enhanced_EMI",
        TRUE ~ "Not_Applicable"
      ),
      effect_category = case_when(
        effect == 1 ~ "Null",
        effect <= 5 ~ "Small",
        effect <= 12 ~ "Medium",
        TRUE ~ "Large"
      ),
      sample_size_group = paste0("n_", n),
      observation_frequency = paste0("K_", K),
      dropout_level = dropout
    )

  return(power_data)
}

# 比較分析用データの作成
create_comparison_data_for_looker <- function(simulation_results) {
  # 手法間の直接比較
  comparison_data <- simulation_results %>%
    group_by(sim_id, n, K, dropout, effect) %>%
    summarise(
      # 最良GPC手法の特定
      best_gpc_power = max(c(
        gpc_direct_nb_reject, gpc_mid_nb_reject,
        gpc_right_nb_reject, gpc_emi_nb_reject
      ), na.rm = TRUE),

      # 従来手法
      rmst_power = rmst_mid_reject,
      logrank_power = lr_mid_reject,

      # 優劣判定
      gpc_beats_rmst = best_gpc_power > rmst_power,
      gpc_beats_logrank = best_gpc_power > logrank_power,

      .groups = "drop"
    ) %>%
    group_by(n, K, dropout, effect) %>%
    summarise(
      # 勝率計算
      gpc_vs_rmst_win_rate = mean(gpc_beats_rmst, na.rm = TRUE),
      gpc_vs_logrank_win_rate = mean(gpc_beats_logrank, na.rm = TRUE),

      # 平均検出力差
      mean_power_diff_rmst = mean(best_gpc_power - rmst_power, na.rm = TRUE),
      mean_power_diff_logrank = mean(best_gpc_power - logrank_power, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    # メタデータ追加
    mutate(
      comparison_context = paste0(
        "n", n, "_K", K, "_", dropout, "_effect", effect
      ),
      overall_gpc_advantage = (gpc_vs_rmst_win_rate + gpc_vs_logrank_win_rate) / 2
    )

  return(comparison_data)
}

# ===============================================================================
# 2. CSV出力機能
# ===============================================================================

# メインの出力関数
export_for_looker_studio <- function(simulation_results, analysis_results = NULL,
                                   output_dir = "looker_studio_data") {
  # 出力ディレクトリの作成
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("=== Google Looker Studio用データ出力開始 ===\n")

  # 1. メインの変換データ
  cat("1. メインデータの変換中...\n")
  main_data <- transform_for_looker_studio(simulation_results, analysis_results)
  write.csv(main_data, file.path(output_dir, "main_simulation_data.csv"), row.names = FALSE)
  cat("   ✓ main_simulation_data.csv 出力完了\n")

  # 2. 検出力サマリーデータ
  cat("2. 検出力サマリーデータの作成中...\n")
  power_data <- create_power_summary_for_looker(simulation_results)
  write.csv(power_data, file.path(output_dir, "power_analysis_summary.csv"), row.names = FALSE)
  cat("   ✓ power_analysis_summary.csv 出力完了\n")

  # 3. 手法比較データ
  cat("3. 手法比較データの作成中...\n")
  comparison_data <- create_comparison_data_for_looker(simulation_results)
  write.csv(comparison_data, file.path(output_dir, "method_comparison_analysis.csv"), row.names = FALSE)
  cat("   ✓ method_comparison_analysis.csv 出力完了\n")

  # 4. 第1種誤り分析データ
  cat("4. 第1種誤り分析データの作成中...\n")
  type1_data <- simulation_results %>%
    filter(effect == 1) %>%
    select(n, K, dropout, contains("_reject")) %>%
    pivot_longer(
      cols = ends_with("_reject"),
      names_to = "method",
      values_to = "type1_error_rate"
    ) %>%
    group_by(n, K, dropout, method) %>%
    summarise(
      type1_error_rate = mean(type1_error_rate, na.rm = TRUE),
      n_observations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      method = str_remove(method, "_reject"),
      method_category = case_when(
        str_detect(method, "gpc.*nb") ~ "GPC_NetBenefit",
        str_detect(method, "gpc.*wr") ~ "GPC_WinRatio",
        str_detect(method, "rmst") ~ "RMST",
        str_detect(method, "lr") ~ "LogRank",
        TRUE ~ "Other"
      ),
      is_within_nominal = abs(type1_error_rate - 0.05) <= 0.02
    )

  write.csv(type1_data, file.path(output_dir, "type1_error_analysis.csv"), row.names = FALSE)
  cat("   ✓ type1_error_analysis.csv 出力完了\n")

  # 5. ダッシュボード用統合データ
  cat("5. ダッシュボード用統合データの作成中...\n")
  dashboard_data <- power_data %>%
    left_join(
      type1_data %>%
        select(n, K, dropout, method, type1_error_rate) %>%
        rename(method_for_join = method),
      by = c("n", "K", "dropout")
    ) %>%
    # 研究の主要指標を統合
    mutate(
      research_question = case_when(
        str_detect(method, "GPC") ~ "Q1: GPC代入法比較",
        TRUE ~ "Q2: 従来手法比較"
      ),
      overall_performance = power * (1 - abs(type1_error_rate - 0.05) * 2),
      performance_rank = rank(-overall_performance, ties.method = "min")
    )

  write.csv(dashboard_data, file.path(output_dir, "looker_studio_dashboard.csv"), row.names = FALSE)
  cat("   ✓ looker_studio_dashboard.csv 出力完了\n")

  # 6. データディクショナリの作成
  cat("6. データディクショナリの作成中...\n")
  create_data_dictionary(output_dir)

  cat("=== 出力完了 ===\n")
  cat("出力先:", normalizePath(output_dir), "\n")
  cat("ファイル一覧:\n")
  files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
  for (file in files) {
    cat(sprintf("  - %s\n", file))
  }

  # 出力サマリーの返却
  return(list(
    output_directory = normalizePath(output_dir),
    files_created = files,
    main_data_rows = nrow(main_data),
    power_data_rows = nrow(power_data),
    comparison_data_rows = nrow(comparison_data),
    dashboard_data_rows = nrow(dashboard_data)
  ))
}

# データディクショナリの作成
create_data_dictionary <- function(output_dir) {
  dictionary <- data.frame(
    File = c(
      "main_simulation_data.csv",
      "power_analysis_summary.csv",
      "method_comparison_analysis.csv",
      "type1_error_analysis.csv",
      "looker_studio_dashboard.csv"
    ),
    Description = c(
      "全シミュレーション結果の詳細データ（long format）",
      "手法別・条件別の検出力サマリー",
      "GPC vs 従来手法の比較分析結果",
      "第1種誤り率の分析結果",
      "Looker Studioダッシュボード用統合データ"
    ),
    Key_Columns = c(
      "method, method_category, imputation_method, n, K, dropout, effect, rejected, p_value",
      "method, power, method_category, imputation_type, n, K, dropout, effect",
      "n, K, dropout, effect, gpc_vs_rmst_win_rate, mean_power_diff_rmst",
      "method, type1_error_rate, method_category, n, K, dropout",
      "method, power, type1_error_rate, overall_performance, research_question"
    ),
    Use_Case = c(
      "詳細分析、フィルタリング、ドリルダウン",
      "検出力比較、トレンド分析",
      "手法優劣の可視化",
      "統計的妥当性の確認",
      "エグゼクティブサマリー、KPIダッシュボード"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(dictionary, file.path(output_dir, "data_dictionary.csv"), row.names = FALSE)
  cat("   ✓ data_dictionary.csv 出力完了\n")
}

# ===============================================================================
# 3. 実行例
# ===============================================================================

# クイック出力関数
quick_export_for_looker <- function() {
  cat("=== クイック出力テスト ===\n")

  # サンプルデータの作成（実際の結果がない場合用）
  if (!exists("simulation_results") || is.null(simulation_results)) {
    cat("シミュレーション結果が見つかりません。サンプルデータを作成します。\n")

    # ダミーデータの作成
    simulation_results <- expand.grid(
      sim_id = 1:10,
      n = c(100, 200),
      K = c(3, 5),
      dropout = c("None", "Medium"),
      effect = c(1, 6),
      stringsAsFactors = FALSE
    ) %>%
    mutate(
      # ダミーの棄却結果
      gpc_direct_nb_reject = runif(n()) < 0.1,
      gpc_mid_nb_reject = runif(n()) < 0.15,
      gpc_right_nb_reject = runif(n()) < 0.12,
      gpc_emi_nb_reject = runif(n()) < 0.18,
      gpc_direct_wr_reject = runif(n()) < 0.08,
      gpc_mid_wr_reject = runif(n()) < 0.13,
      gpc_right_wr_reject = runif(n()) < 0.10,
      gpc_emi_wr_reject = runif(n()) < 0.16,
      rmst_mid_reject = runif(n()) < 0.09,
      lr_mid_reject = runif(n()) < 0.11,

      # ダミーのp値
      gpc_direct_nb_pvalue = runif(n()),
      gpc_mid_nb_pvalue = runif(n()),
      gpc_right_nb_pvalue = runif(n()),
      gpc_emi_nb_pvalue = runif(n()),
      gpc_direct_wr_pvalue = runif(n()),
      gpc_mid_wr_pvalue = runif(n()),
      gpc_right_wr_pvalue = runif(n()),
      gpc_emi_wr_pvalue = runif(n()),
      rmst_mid_pvalue = runif(n()),
      lr_mid_pvalue = runif(n())
    )
  }

  # 出力実行
  export_result <- export_for_looker_studio(simulation_results)

  cat("\n=== Looker Studio での使用方法 ===\n")
  cat("1. Google Looker Studio (https://lookerstudio.google.com/) にアクセス\n")
  cat("2. '空のレポートを作成' を選択\n")
  cat("3. データソースとして 'ファイルのアップロード' を選択\n")
  cat("4. 作成されたCSVファイルをアップロード\n")
  cat("5. 推奨ファイル: looker_studio_dashboard.csv （最も包括的）\n")

  return(export_result)
}

cat("===============================================================================\n")
cat("            Google Looker Studio用データ出力機能が準備完了\n")
cat("===============================================================================\n")
cat("使用方法:\n")
cat("1. export_for_looker_studio(simulation_results)  # 本格出力\n")
cat("2. quick_export_for_looker()                     # テスト出力\n")
cat("===============================================================================\n")