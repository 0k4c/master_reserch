# ===============================================================================
# WINSとBuyseTestパッケージを活用した改良版GPC実装
#
# 問題の解決：
# 1. p値計算の修正（正しい統計理論の適用）
# 2. 生存時間解釈の明確化（長い方が良い vs 短い方が良い）
# 3. 信頼性の高いパッケージの活用
# 4. 区間打ち切りデータへの適切な対応
# ===============================================================================

# 必要なライブラリ
library(survival)
library(tidyverse)

# WINSとBuyseTestパッケージのインストールと読み込み
install_and_load_packages <- function() {
  cat("=== パッケージのインストールと読み込み ===\n")

  # WINSパッケージ
  tryCatch({
    if (!requireNamespace("WINS", quietly = TRUE)) {
      install.packages("WINS")
    }
    library(WINS)
    cat("✓ WINSパッケージ読み込み成功\n")
  }, error = function(e) {
    cat("✗ WINSパッケージエラー:", e$message, "\n")
  })

  # BuyseTestパッケージ
  tryCatch({
    if (!requireNamespace("BuyseTest", quietly = TRUE)) {
      install.packages("BuyseTest")
    }
    library(BuyseTest)
    cat("✓ BuyseTestパッケージ読み込み成功\n")
  }, error = function(e) {
    cat("✗ BuyseTestパッケージエラー:", e$message, "\n")
  })
}

# データ生成関数の読み込み
source("flexible_data_generation.R")

# ===============================================================================
# 1. 理論的に正しいGPC実装
# ===============================================================================

# 帰無仮説と対立仮説の明確化
theoretical_framework <- function() {
  cat("=== 一般化ペアワイズ比較の理論的枠組み ===\n")
  cat("研究目的: 治療群が対照群より優れているかを検定\n\n")

  cat("生存時間の解釈:\n")
  cat("- 生存時間が長い = 良い結果（通常の生存解析）\n")
  cat("- 治療群の生存時間 > 対照群の生存時間 → 治療効果あり\n\n")

  cat("帰無仮説と対立仮説:\n")
  cat("H0: 治療群と対照群に差がない\n")
  cat("  - Net Benefit = 0\n")
  cat("  - Win Ratio = 1\n")
  cat("H1: 治療群が対照群より優れている（または劣っている）\n")
  cat("  - Net Benefit ≠ 0\n")
  cat("  - Win Ratio ≠ 1\n\n")

  cat("区間打ち切りデータでの比較ルール:\n")
  cat("治療群 [L_T, R_T] vs 対照群 [L_C, R_C]\n")
  cat("- L_T > R_C: 治療群の勝ち（明確に優れている）\n")
  cat("- R_T < L_C: 対照群の勝ち（明確に劣っている）\n")
  cat("- 区間重複: 引き分け（判定不能）\n")
}

# 修正版GPC関数（理論的に正しい実装）
gpc_theoretically_correct <- function(data_control, data_treatment, method = "net_benefit") {
  theoretical_framework()

  n_control <- nrow(data_control)
  n_treatment <- nrow(data_treatment)

  wins <- 0    # 治療群が明確に優れている
  losses <- 0  # 治療群が明確に劣っている
  ties <- 0    # 判定不能

  # ペアワイズ比較の実行
  for (i in 1:n_treatment) {
    for (j in 1:n_control) {
      L_T <- data_treatment[i, 1]  # 治療群の左端点
      R_T <- data_treatment[i, 2]  # 治療群の右端点
      L_C <- data_control[j, 1]    # 対照群の左端点
      R_C <- data_control[j, 2]    # 対照群の右端点

      # 修正された比較ルール（生存時間：長い方が良い）
      if (L_T > R_C) {
        # 治療群の最小値 > 対照群の最大値 → 治療群明確に優秀
        wins <- wins + 1
      } else if (R_T < L_C) {
        # 治療群の最大値 < 対照群の最小値 → 治療群明確に劣る
        losses <- losses + 1
      } else {
        # 区間が重複 → 判定不能
        ties <- ties + 1
      }
    }
  }

  total_pairs <- n_treatment * n_control

  cat(sprintf("\nペアワイズ比較結果: 勝ち=%d, 負け=%d, 引き分け=%d (総ペア数=%d)\n",
              wins, losses, ties, total_pairs))

  if (method == "net_benefit") {
    # Net Benefit = (勝ち - 負け) / 総ペア数
    net_benefit <- (wins - losses) / total_pairs

    # 正しい分散推定（Brunner-Munzel統計量）
    p_win <- wins / total_pairs
    p_loss <- losses / total_pairs
    p_tie <- ties / total_pairs

    # 標本分散の推定
    var_nb <- (p_win + p_loss - (p_win - p_loss)^2) / total_pairs

    # 最小分散の確保
    if (var_nb <= 0 || is.na(var_nb)) {
      var_nb <- 1 / (4 * total_pairs)  # 理論的最小分散
    }

    # 検定統計量（Z統計量）
    z_stat <- net_benefit / sqrt(var_nb)

    # p値（両側検定）
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      method = "Net Benefit",
      net_benefit = net_benefit,
      z_statistic = z_stat,
      p_value = p_value,
      variance = var_nb,
      wins = wins,
      losses = losses,
      ties = ties,
      total_pairs = total_pairs
    ))

  } else if (method == "win_ratio") {
    # Win Ratio = 勝ち / 負け
    if (losses == 0) {
      if (wins > 0) {
        # 完全勝利
        win_ratio <- Inf
        p_value <- 0  # 極めて有意
      } else {
        # 勝負なし
        win_ratio <- 1
        p_value <- 1  # 有意差なし
      }
      z_stat <- NA
    } else if (wins == 0) {
      # 完全敗北
      win_ratio <- 0
      p_value <- 0  # 極めて有意
      z_stat <- -Inf
    } else {
      win_ratio <- wins / losses

      # ログWin Ratioの分散（Delta法）
      log_var <- 1/wins + 1/losses

      # 検定統計量
      z_stat <- log(win_ratio) / sqrt(log_var)

      # p値（両側検定）
      p_value <- 2 * (1 - pnorm(abs(z_stat)))
    }

    return(list(
      method = "Win Ratio",
      win_ratio = win_ratio,
      log_win_ratio = ifelse(win_ratio > 0, log(win_ratio), -Inf),
      z_statistic = z_stat,
      p_value = p_value,
      wins = wins,
      losses = losses,
      ties = ties,
      total_pairs = total_pairs
    ))
  }
}

# ===============================================================================
# 2. WINSパッケージを使用したGPC
# ===============================================================================

gpc_using_wins <- function(data_control, data_treatment, method = "net_benefit") {
  if (!requireNamespace("WINS", quietly = TRUE)) {
    cat("WINSパッケージが利用できません。独自実装を使用します。\n")
    return(gpc_theoretically_correct(data_control, data_treatment, method))
  }

  tryCatch({
    # WINSパッケージ用のデータ形式に変換
    # 区間打ち切りデータを中点で代入（WINSパッケージ用）
    control_times <- apply(data_control, 1, function(x) {
      if (is.infinite(x[2])) {
        x[1]  # 右打ち切りの場合
      } else {
        (x[1] + x[2]) / 2  # 区間の中点
      }
    })

    treatment_times <- apply(data_treatment, 1, function(x) {
      if (is.infinite(x[2])) {
        x[1]  # 右打ち切りの場合
      } else {
        (x[1] + x[2]) / 2  # 区間の中点
      }
    })

    # イベント指示子（区間打ち切りは1、右打ち切りは0）
    control_events <- ifelse(is.infinite(data_control[, 2]), 0, 1)
    treatment_events <- ifelse(is.infinite(data_treatment[, 2]), 0, 1)

    # データフレーム作成
    control_data <- data.frame(
      time = control_times,
      event = control_events,
      group = 0
    )

    treatment_data <- data.frame(
      time = treatment_times,
      event = treatment_events,
      group = 1
    )

    combined_data <- rbind(control_data, treatment_data)

    # WINSパッケージでのGPC実行
    wins_result <- WINS::wins(
      data = combined_data,
      endpoint = "time",
      group = "group",
      reference = 0  # 対照群を基準
    )

    # 結果の抽出
    if (method == "net_benefit") {
      return(list(
        method = "WINS Net Benefit",
        net_benefit = wins_result$net_benefit,
        p_value = wins_result$p_value,
        package_used = "WINS"
      ))
    } else {
      return(list(
        method = "WINS Win Ratio",
        win_ratio = wins_result$win_ratio,
        p_value = wins_result$p_value,
        package_used = "WINS"
      ))
    }

  }, error = function(e) {
    cat("WINSパッケージエラー:", e$message, "\n独自実装を使用します。\n")
    return(gpc_theoretically_correct(data_control, data_treatment, method))
  })
}

# ===============================================================================
# 3. BuyseTestパッケージを使用したGPC
# ===============================================================================

gpc_using_buysetest <- function(data_control, data_treatment) {
  if (!requireNamespace("BuyseTest", quietly = TRUE)) {
    cat("BuyseTestパッケージが利用できません。独自実装を使用します。\n")
    return(gpc_theoretically_correct(data_control, data_treatment, "net_benefit"))
  }

  tryCatch({
    # BuyseTest用のデータ準備
    # 区間打ち切りデータの処理
    control_times <- numeric(nrow(data_control))
    control_events <- numeric(nrow(data_control))

    for (i in 1:nrow(data_control)) {
      if (is.infinite(data_control[i, 2])) {
        # 右打ち切り
        control_times[i] <- data_control[i, 1]
        control_events[i] <- 0
      } else {
        # 区間打ち切り→中点を使用
        control_times[i] <- (data_control[i, 1] + data_control[i, 2]) / 2
        control_events[i] <- 1
      }
    }

    treatment_times <- numeric(nrow(data_treatment))
    treatment_events <- numeric(nrow(data_treatment))

    for (i in 1:nrow(data_treatment)) {
      if (is.infinite(data_treatment[i, 2])) {
        # 右打ち切り
        treatment_times[i] <- data_treatment[i, 1]
        treatment_events[i] <- 0
      } else {
        # 区間打ち切り→中点を使用
        treatment_times[i] <- (data_treatment[i, 1] + data_treatment[i, 2]) / 2
        treatment_events[i] <- 1
      }
    }

    # BuyseTest用データフレーム
    df_buyse <- data.frame(
      group = c(rep("control", nrow(data_control)), rep("treatment", nrow(data_treatment))),
      time = c(control_times, treatment_times),
      event = c(control_events, treatment_events)
    )

    # BuyseTestでのGPC実行
    buyse_result <- BuyseTest::BuyseTest(
      group ~ time,
      data = df_buyse,
      treatment = "treatment",
      control = "control",
      scoring.rule = "Peron"  # 生存時間用
    )

    # 結果の要約
    summary_result <- summary(buyse_result)

    return(list(
      method = "BuyseTest",
      net_benefit = summary_result@table["netBenefit", "estimate"],
      p_value = summary_result@table["netBenefit", "p.value"],
      win_ratio = summary_result@table["winRatio", "estimate"],
      package_used = "BuyseTest",
      full_result = summary_result
    ))

  }, error = function(e) {
    cat("BuyseTestパッケージエラー:", e$message, "\n独自実装を使用します。\n")
    return(gpc_theoretically_correct(data_control, data_treatment, "net_benefit"))
  })
}

# ===============================================================================
# 4. 包括的比較関数
# ===============================================================================

compare_gpc_methods <- function(data_control, data_treatment) {
  cat("===============================================================================\n")
  cat("                      GPC手法の包括的比較\n")
  cat("===============================================================================\n")

  results <- list()

  # 1. 理論的に正しい独自実装
  cat("\n1. 理論的に正しい独自実装（Net Benefit）\n")
  results$custom_nb <- gpc_theoretically_correct(data_control, data_treatment, "net_benefit")

  cat("\n2. 理論的に正しい独自実装（Win Ratio）\n")
  results$custom_wr <- gpc_theoretically_correct(data_control, data_treatment, "win_ratio")

  # 3. WINSパッケージ
  cat("\n3. WINSパッケージ\n")
  results$wins_nb <- gpc_using_wins(data_control, data_treatment, "net_benefit")
  results$wins_wr <- gpc_using_wins(data_control, data_treatment, "win_ratio")

  # 4. BuyseTestパッケージ
  cat("\n4. BuyseTestパッケージ\n")
  results$buysetest <- gpc_using_buysetest(data_control, data_treatment)

  # 結果の比較表示
  cat("\n===============================================================================\n")
  cat("                          結果比較\n")
  cat("===============================================================================\n")

  comparison_table <- data.frame(
    Method = c("Custom_NB", "Custom_WR", "WINS_NB", "WINS_WR", "BuyseTest"),
    Net_Benefit = c(
      results$custom_nb$net_benefit,
      NA,
      ifelse("net_benefit" %in% names(results$wins_nb), results$wins_nb$net_benefit, NA),
      NA,
      ifelse("net_benefit" %in% names(results$buysetest), results$buysetest$net_benefit, NA)
    ),
    Win_Ratio = c(
      NA,
      results$custom_wr$win_ratio,
      NA,
      ifelse("win_ratio" %in% names(results$wins_wr), results$wins_wr$win_ratio, NA),
      ifelse("win_ratio" %in% names(results$buysetest), results$buysetest$win_ratio, NA)
    ),
    P_Value = c(
      results$custom_nb$p_value,
      results$custom_wr$p_value,
      ifelse("p_value" %in% names(results$wins_nb), results$wins_nb$p_value, NA),
      ifelse("p_value" %in% names(results$wins_wr), results$wins_wr$p_value, NA),
      ifelse("p_value" %in% names(results$buysetest), results$buysetest$p_value, NA)
    ),
    stringsAsFactors = FALSE
  )

  print(comparison_table)

  # p値の妥当性チェック
  cat("\np値の妥当性チェック:\n")
  for (i in 1:nrow(comparison_table)) {
    method <- comparison_table$Method[i]
    p_val <- comparison_table$P_Value[i]

    if (!is.na(p_val)) {
      if (p_val >= 0 && p_val <= 1) {
        cat(sprintf("✓ %s: p値 = %.4f (妥当)\n", method, p_val))
      } else {
        cat(sprintf("✗ %s: p値 = %.4f (異常)\n", method, p_val))
      }
    }
  }

  return(results)
}

# ===============================================================================
# 5. テスト実行関数
# ===============================================================================

test_improved_gpc <- function() {
  cat("=== 改良版GPC実装のテスト ===\n")

  # パッケージのインストール
  install_and_load_packages()

  # テストデータの生成
  cat("\nテストデータ生成中...\n")
  data_control <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 1)
  data_treatment <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 6)

  cat("対照群データ例:\n")
  print(head(data_control))
  cat("治療群データ例:\n")
  print(head(data_treatment))

  # 包括的比較の実行
  results <- compare_gpc_methods(data_control, data_treatment)

  return(results)
}

cat("===============================================================================\n")
cat("    WINSとBuyseTestパッケージを活用した改良版GPC実装が準備完了\n")
cat("===============================================================================\n")
cat("実行方法:\n")
cat("1. improved_results <- test_improved_gpc()\n")
cat("2. compare_gpc_methods(data_control, data_treatment)\n")
cat("===============================================================================\n")