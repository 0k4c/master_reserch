# ===============================================================================
# 簡単な動作確認テスト
# ===============================================================================

cat("=== 簡単な動作確認テスト ===\n")

# 必要なライブラリの読み込み
library(survival)
library(tidyverse)

# ファイルの読み込み
cat("1. ファイル読み込み中...\n")
source("flexible_data_generation.R")
source("corrected_gpc_simulation.R")

# データ生成テスト
cat("2. データ生成テスト中...\n")
data_control <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 1)
data_treatment <- generate_interval_censored_data_2arm(n = 100, K = 3, p_dropout = "None", d = 6)

cat("データ生成成功:\n")
cat("  対照群:", nrow(data_control), "x", ncol(data_control), "\n")
cat("  治療群:", nrow(data_treatment), "x", ncol(data_treatment), "\n")

# データの内容確認
cat("\n対照群データの例:\n")
print(head(data_control))
cat("\n治療群データの例:\n")
print(head(data_treatment))

# GPC Direct テスト
cat("\n3. GPC Direct テスト中...\n")
gpc_nb <- gpc_direct_interval(data_control, data_treatment, method = "net_benefit")
gpc_wr <- gpc_direct_interval(data_control, data_treatment, method = "win_ratio")

cat("GPC Net Benefit:\n")
cat("  Net Benefit =", gpc_nb$net_benefit, "\n")
cat("  p値 =", gpc_nb$p_value, "\n")
cat("  勝ち/負け/引き分け =", gpc_nb$wins, "/", gpc_nb$losses, "/", gpc_nb$ties, "\n")

cat("\nGPC Win Ratio:\n")
cat("  Win Ratio =", gpc_wr$win_ratio, "\n")
cat("  p値 =", gpc_wr$p_value, "\n")

# 代入法テスト
cat("\n4. 代入法テスト中...\n")
mid_control <- midpoint_assignment(data_control)
mid_treatment <- midpoint_assignment(data_treatment)

cat("中点代入法結果:\n")
cat("  対照群イベント率:", mean(mid_control$cens), "\n")
cat("  治療群イベント率:", mean(mid_treatment$cens), "\n")

# GPC + 代入法テスト
cat("\n5. GPC + 代入法テスト中...\n")
gpc_mid_nb <- gpc_with_imputation(data_control, data_treatment,
                                 imputation_method = "midpoint", method = "net_benefit")

cat("GPC + 中点代入法 (Net Benefit):\n")
cat("  p値 =", gpc_mid_nb$p_value, "\n")

# 結果まとめ
cat("\n=== テスト結果まとめ ===\n")
cat("✓ データ生成: 正常\n")
cat("✓ GPC Direct: p値 =", gpc_nb$p_value, "(NB),", gpc_wr$p_value, "(WR)\n")
cat("✓ 代入法: 正常\n")
cat("✓ GPC + 代入法: p値 =", gpc_mid_nb$p_value, "\n")

if (gpc_nb$p_value < 1 && gpc_nb$p_value > 0) {
  cat("✓ p値が適切に計算されています！\n")
} else {
  cat("✗ p値に問題があります。\n")
}

cat("\n次のステップ: run_corrected_test() でフルテストを実行\n")