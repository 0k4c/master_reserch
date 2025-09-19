# 🚀 実行ガイド

このプロジェクトのシミュレーションを実行するためのクイックスタートガイドです。

## ⚡ クイックスタート

**作業ディレクトリは `final_research_project/R` を想定しています。**

### 1. 小規模研究実行（動作確認）

まず、小規模なテストシミュレーションを実行して、環境に問題がないかを確認します。

```r
# main_simulation.R を読み込み、テスト関数を実行
source("main_simulation.R")
test_results <- run_final_study()
```

実行後、`results`フォルダに3つのCSVファイルが出力されていれば成功です。

### 2. 本格研究実行（パラメータ指定）

本格的なシミュレーションは、`run_final_simulation()` 関数にパラメータを指定して実行します。

```r
# main_simulation.R を読み込む
source("main_simulation.R")

# フル条件でのシミュレーション（並列処理・プログレスバー付き）
full_results <- run_final_simulation(
  n_sim = 1000,                              # シミュレーション回数
  sample_sizes = c(100, 200, 400),           # サンプルサイズ
  K_values = c(3, 5, 10),                    # 観測回数
  dropout_levels = c("None", "Low", "Medium", "High"),  # 脱落率
  effect_sizes = 1:24,                       # 24種類の効果サイズ
  n_cores = NULL                             # 使用するCPUコア数 (NULLで自動検出)
)

# 結果の分析
analysis <- analyze_final_results(full_results)

# 結果の保存
if (!dir.exists("results")) dir.create("results")
write.csv(full_results, "results/full_simulation_results.csv", row.names = FALSE)
write.csv(analysis$power_summary, "results/full_power_summary.csv", row.names = FALSE)
write.csv(analysis$type1_summary, "results/full_type1_error_summary.csv", row.names = FALSE)
```

