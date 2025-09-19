# 🚀 実行ガイド

## クイックスタート

### 1. 動作確認（必須）
```r
source("R/verification_test.R")
verification_results <- run_final_verification()
```

### 2. 小規模研究実行
```r
source("R/main_simulation.R")
test_results <- run_final_study()
```

### 3. 本格研究実行
```r
# フル条件でのシミュレーション（並列処理・プログレスバー付き）
full_results <- run_final_simulation(
  n_sim = 1000,
  sample_sizes = c(100, 200, 400),
  K_values = c(3, 5, 10),
  dropout_levels = c("None", "Low", "Medium", "High"),
  effect_sizes = 1:24,
  n_cores = 8,          # コア数指定（NULLで自動検出）
  use_parallel = TRUE,  # 並列処理有効
  show_progress = TRUE  # プログレスバー表示
)

# 結果分析
analysis <- analyze_final_results(full_results)
```

### 3.1. 詳細オプション設定
```r
# 小規模テスト
test_results <- run_final_simulation(
  n_sim = 100,
  sample_sizes = c(50, 100),
  n_cores = 4,          # 4コア使用
  show_progress = TRUE   # プログレスバー表示
)

# デバッグ用（逐次処理）
debug_results <- run_final_simulation(
  n_sim = 10,
  n_cores = 1,          # 1コア（逐次処理）
  use_parallel = FALSE, # 並列処理無効
  show_progress = TRUE
)
```

### 4. Looker Studio用データ出力
```r
source("R/looker_export.R")
export_for_looker_studio(full_results, analysis)
```

## ファイル構成

- `R/main_simulation.R`: メイン研究コード
- `R/verification_test.R`: 動作確認テスト
- `R/data_generation.R`: データ生成関数
- `R/looker_export.R`: CSV出力機能
- `presentation/seminar_slides.tex`: Overleaf用スライド

## トラブルシューティング

問題が発生した場合は、まず verification_test.R を実行してください。