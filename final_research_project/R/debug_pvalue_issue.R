# ===============================================================================
# p値計算問題のデバッグ用スクリプト
# ===============================================================================

library(survival)
library(tidyverse)

# ファイル読み込み
source("flexible_data_generation.R")

# ===============================================================================
# 1. 現在のGPC関数の問題点を特定
# ===============================================================================

# 簡単なテストデータを作成
debug_gpc_calculation <- function() {
  cat("=== GPC計算のデバッグ ===\n")

  # 小さなテストデータを手動作成
  # 対照群：区間 [0.2, 0.5], [0.3, 0.6], [0.1, 0.4]
  data_control <- matrix(c(0.2, 0.5,
                          0.3, 0.6,
                          0.1, 0.4), nrow = 3, byrow = TRUE)

  # 治療群：区間 [0.7, 0.9], [0.8, 1.0], [0.6, 0.8] (より良い結果を期待)
  data_treatment <- matrix(c(0.7, 0.9,
                            0.8, 1.0,
                            0.6, 0.8), nrow = 3, byrow = TRUE)

  cat("対照群データ:\n")
  print(data_control)
  cat("治療群データ:\n")
  print(data_treatment)

  # 手動でペアワイズ比較を計算
  cat("\n=== 手動ペアワイズ比較 ===\n")
  wins <- 0
  losses <- 0
  ties <- 0

  for (i in 1:3) {
    for (j in 1:3) {
      treat_left <- data_treatment[i, 1]
      treat_right <- data_treatment[i, 2]
      control_left <- data_control[j, 1]
      control_right <- data_control[j, 2]

      cat(sprintf("治療群%d [%.1f, %.1f] vs 対照群%d [%.1f, %.1f]: ",
                  i, treat_left, treat_right, j, control_left, control_right))

      # 現在のルール（問題がある可能性）
      if (treat_right < control_left) {
        cat("治療群勝ち (現在のルール)\n")
        wins <- wins + 1
      } else if (control_right < treat_left) {
        cat("対照群勝ち (現在のルール)\n")
        losses <- losses + 1
      } else {
        cat("引き分け\n")
        ties <- ties + 1
      }
    }
  }

  cat(sprintf("\n現在のルールでの結果: 勝ち=%d, 負け=%d, 引き分け=%d\n", wins, losses, ties)

  # 生存時間の正しい解釈（長い方が良い）で再計算
  cat("\n=== 修正されたペアワイズ比較（長い生存時間が良い）===\n")
  wins_correct <- 0
  losses_correct <- 0
  ties_correct <- 0

  for (i in 1:3) {
    for (j in 1:3) {
      treat_left <- data_treatment[i, 1]
      treat_right <- data_treatment[i, 2]
      control_left <- data_control[j, 1]
      control_right <- data_control[j, 2]

      cat(sprintf("治療群%d [%.1f, %.1f] vs 対照群%d [%.1f, %.1f]: ",
                  i, treat_left, treat_right, j, control_left, control_right))

      # 修正されたルール（長い生存時間が良い）
      if (treat_left > control_right) {
        cat("治療群勝ち (修正ルール)\n")
        wins_correct <- wins_correct + 1
      } else if (control_left > treat_right) {
        cat("対照群勝ち (修正ルール)\n")
        losses_correct <- losses_correct + 1
      } else {
        cat("引き分け\n")
        ties_correct <- ties_correct + 1
      }
    }
  }

  cat(sprintf("\n修正ルールでの結果: 勝ち=%d, 負け=%d, 引き分け=%d\n",
              wins_correct, losses_correct, ties_correct))

  return(list(
    original = list(wins = wins, losses = losses, ties = ties),
    corrected = list(wins = wins_correct, losses = losses_correct, ties = ties_correct)
  ))
}

# ===============================================================================
# 2. 統計的検定の理論確認
# ===============================================================================

test_statistical_theory <- function() {
  cat("\n=== 統計的検定の理論確認 ===\n")

  # 帰無仮説と対立仮説
  cat("帰無仮説 H0: 治療群と対照群に差がない\n")
  cat("  - Net Benefit の場合: E[NB] = 0\n")
  cat("  - Win Ratio の場合: E[WR] = 1 (log(WR) = 0)\n")

  cat("\n対立仮説 H1: 治療群と対照群に差がある\n")
  cat("  - Net Benefit の場合: E[NB] ≠ 0\n")
  cat("  - Win Ratio の場合: E[WR] ≠ 1 (log(WR) ≠ 0)\n")

  # 検定統計量の計算方法
  cat("\n検定統計量の計算:\n")
  cat("1. Net Benefit: Z = NB / sqrt(Var(NB))\n")
  cat("2. Win Ratio: Z = log(WR) / sqrt(Var(log(WR)))\n")
  cat("3. p値: P = 2 * (1 - Φ(|Z|)) (両側検定)\n")
}

# ===============================================================================
# 3. 正しいGPC実装を作成
# ===============================================================================

gpc_corrected <- function(data_control, data_treatment, method = "net_benefit",
                         survival_interpretation = TRUE) {
  n_control <- nrow(data_control)
  n_treatment <- nrow(data_treatment)

  wins <- 0
  losses <- 0
  ties <- 0

  for (i in 1:n_treatment) {
    for (j in 1:n_control) {
      treat_left <- data_treatment[i, 1]
      treat_right <- data_treatment[i, 2]
      control_left <- data_control[j, 1]
      control_right <- data_control[j, 2]

      if (survival_interpretation) {
        # 生存時間解釈：長い方が良い
        if (treat_left > control_right) {
          # 治療群の最小値 > 対照群の最大値 → 治療群明らかに良い
          wins <- wins + 1
        } else if (control_left > treat_right) {
          # 対照群の最小値 > 治療群の最大値 → 対照群明らかに良い
          losses <- losses + 1
        } else {
          # 区間が重複
          ties <- ties + 1
        }
      } else {
        # イベント時間解釈：短い方が良い（元のコード）
        if (treat_right < control_left) {
          wins <- wins + 1
        } else if (control_right < treat_left) {
          losses <- losses + 1
        } else {
          ties <- ties + 1
        }
      }
    }
  }

  total_pairs <- n_treatment * n_control

  if (method == "net_benefit") {
    net_benefit <- (wins - losses) / total_pairs

    # 修正された分散計算
    p_win <- wins / total_pairs
    p_loss <- losses / total_pairs

    # Mann-Whitney統計量の分散公式
    variance <- (p_win * (1 - p_win) + p_loss * (1 - p_loss) + 2 * p_win * p_loss) / total_pairs

    if (variance <= 0) {
      variance <- 1 / total_pairs
    }

    z_stat <- net_benefit / sqrt(variance)
    p_value <- 2 * (1 - pnorm(abs(z_stat)))

    return(list(
      net_benefit = net_benefit,
      z_statistic = z_stat,
      p_value = p_value,
      wins = wins,
      losses = losses,
      ties = ties,
      variance = variance
    ))

  } else if (method == "win_ratio") {
    if (losses == 0) {
      if (wins > 0) {
        # 完全に治療群が優勢
        return(list(
          win_ratio = Inf,
          p_value = 0,
          wins = wins,
          losses = losses,
          ties = ties
        ))
      } else {
        # 勝負がついていない
        return(list(
          win_ratio = 1,
          p_value = 1,
          wins = wins,
          losses = losses,
          ties = ties
        ))
      }
    } else if (wins == 0) {
      # 完全に対照群が優勢
      return(list(
        win_ratio = 0,
        p_value = 0,
        wins = wins,
        losses = losses,
        ties = ties
      ))
    } else {
      win_ratio <- wins / losses
      log_variance <- 1/wins + 1/losses
      z_stat <- log(win_ratio) / sqrt(log_variance)
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
}

# ===============================================================================
# 4. WINSパッケージの調査
# ===============================================================================

investigate_wins_package <- function() {
  cat("\n=== WINSパッケージの調査 ===\n")

  # WINSパッケージのインストール試行
  tryCatch({
    if (!requireNamespace("WINS", quietly = TRUE)) {
      cat("WINSパッケージをインストール中...\n")
      # CRANからのインストールを試行
      install.packages("WINS")
    }

    library(WINS)
    cat("✓ WINSパッケージの読み込み成功\n")

    # パッケージの関数を確認
    cat("利用可能な関数:\n")
    print(ls("package:WINS"))

    return(TRUE)
  }, error = function(e) {
    cat("✗ WINSパッケージの読み込み失敗:", e$message, "\n")

    # GitHubからのインストールを試行
    tryCatch({
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_github("cran/WINS")
      library(WINS)
      cat("✓ WINSパッケージ（GitHub版）の読み込み成功\n")
      return(TRUE)
    }, error = function(e2) {
      cat("✗ GitHubからのインストールも失敗:", e2$message, "\n")
      return(FALSE)
    })
  })
}

# ===============================================================================
# 5. 包括的テスト実行
# ===============================================================================

run_debug_analysis <- function() {
  cat("===============================================================================\n")
  cat("                          p値計算問題のデバッグ分析\n")
  cat("===============================================================================\n")

  # 1. 理論確認
  test_statistical_theory()

  # 2. 手動計算でのデバッグ
  manual_results <- debug_gpc_calculation()

  # 3. 修正版GPC関数のテスト
  cat("\n=== 修正版GPC関数のテスト ===\n")

  # 実データでテスト
  data_control <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 1)
  data_treatment <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 6)

  # 元の解釈（短い方が良い）
  result_original <- gpc_corrected(data_control, data_treatment,
                                  method = "net_benefit", survival_interpretation = FALSE)

  # 修正された解釈（長い方が良い）
  result_corrected <- gpc_corrected(data_control, data_treatment,
                                   method = "net_benefit", survival_interpretation = TRUE)

  cat("元の解釈（短い時間が良い）:\n")
  cat("  Net Benefit =", result_original$net_benefit, "\n")
  cat("  p値 =", result_original$p_value, "\n")
  cat("  勝ち/負け/引き分け =", result_original$wins, "/",
      result_original$losses, "/", result_original$ties, "\n")

  cat("\n修正された解釈（長い時間が良い）:\n")
  cat("  Net Benefit =", result_corrected$net_benefit, "\n")
  cat("  p値 =", result_corrected$p_value, "\n")
  cat("  勝ち/負け/引き分け =", result_corrected$wins, "/",
      result_corrected$losses, "/", result_corrected$ties, "\n")

  # 4. WINSパッケージの調査
  wins_available <- investigate_wins_package()

  # 5. 結論
  cat("\n===============================================================================\n")
  cat("                                  分析結果\n")
  cat("===============================================================================\n")

  cat("問題の特定:\n")
  if (result_original$p_value >= 0 && result_original$p_value <= 1) {
    cat("✓ p値は適切な範囲内です（元の解釈）\n")
  } else {
    cat("✗ p値が適切な範囲外です（元の解釈）\n")
  }

  if (result_corrected$p_value >= 0 && result_corrected$p_value <= 1) {
    cat("✓ p値は適切な範囲内です（修正された解釈）\n")
  } else {
    cat("✗ p値が適切な範囲外です（修正された解釈）\n")
  }

  cat("\n推奨事項:\n")
  cat("1. 生存時間の解釈を明確にする（長い=良い vs 短い=良い）\n")
  cat("2. 分散計算式を見直す\n")
  if (wins_available) {
    cat("3. WINSパッケージの利用を検討する\n")
  } else {
    cat("3. WINSパッケージは利用できないため、独自実装を改善する\n")
  }

  return(list(
    original_result = result_original,
    corrected_result = result_corrected,
    wins_available = wins_available,
    manual_results = manual_results
  ))
}

# 実行
cat("デバッグスクリプトが準備完了しました。\n")
cat("実行方法: debug_results <- run_debug_analysis()\n")