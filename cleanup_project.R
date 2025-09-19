# ===============================================================================
# プロジェクトクリーンアップスクリプト
#
# 不要なファイルを整理し、最終的なプロジェクト構造を作成
# ===============================================================================

library(fs)

cleanup_project <- function() {
  cat("=== プロジェクトクリーンアップ開始 ===\n")

  # 現在のディレクトリ
  current_dir <- getwd()
  cat("作業ディレクトリ:", current_dir, "\n")

  # 最終版プロジェクト用ディレクトリの作成
  final_dir <- "final_research_project"
  if (dir.exists(final_dir)) {
    unlink(final_dir, recursive = TRUE)
  }
  dir.create(final_dir)
  dir.create(file.path(final_dir, "R"))
  dir.create(file.path(final_dir, "docs"))
  dir.create(file.path(final_dir, "results"))
  dir.create(file.path(final_dir, "presentation"))

  # 必須ファイルのコピー
  essential_files <- list(
    # メインのRファイル
    list(src = "final_imputation_comparison_study.R", dst = "R/main_simulation.R"),
    list(src = "final_verification_test.R", dst = "R/verification_test.R"),
    list(src = "flexible_data_generation.R", dst = "R/data_generation.R"),
    list(src = "looker_studio_export.R", dst = "R/looker_export.R"),

    # 既存の必要なファイル
    list(src = "interval_censord_data_function.R", dst = "R/interval_censord_data_function.R"),
    list(src = "distributions_2arm.R", dst = "R/distributions_2arm.R"),

    # ドキュメント
    list(src = "README.md", dst = "README.md"),

    # プレゼンテーション
    list(src = "seminar_presentation.tex", dst = "presentation/seminar_slides.tex")
  )

  cat("\n必須ファイルのコピー中...\n")
  for (file_info in essential_files) {
    if (file.exists(file_info$src)) {
      file.copy(file_info$src, file.path(final_dir, file_info$dst))
      cat(sprintf("✓ %s -> %s\n", file_info$src, file_info$dst))
    } else {
      cat(sprintf("✗ %s が見つかりません\n", file_info$src))
    }
  }

  # 不要ファイルのリスト
  obsolete_files <- c(
    "generalized_pairwise_comparison_simulation.R",  # 古いバージョン
    "corrected_gpc_simulation.R",                    # 中間バージョン
    "improved_gpc_with_packages.R",                  # 部分的実装
    "debug_pvalue_issue.R",                          # デバッグ用
    "run_test_debug.R",                              # デバッグ用
    "quick_test_example.R",                          # テスト用
    "simple_test.R",                                 # テスト用
    "README_simulation_guide.md"                     # 古いREADME
  )

  # アーカイブディレクトリの作成
  archive_dir <- "archived_development_files"
  if (!dir.exists(archive_dir)) {
    dir.create(archive_dir)
  }

  cat("\n開発ファイルのアーカイブ中...\n")
  for (file in obsolete_files) {
    if (file.exists(file)) {
      file.copy(file, file.path(archive_dir, file))
      cat(sprintf("📦 %s -> archived/\n", file))
    }
  }

  # プロジェクト設定ファイルの作成
  create_project_config(final_dir)

  # 実行ガイドの作成
  create_execution_guide(final_dir)

  cat("\n=== クリーンアップ完了 ===\n")
  cat("最終プロジェクト:", file.path(current_dir, final_dir), "\n")
  cat("アーカイブ:", file.path(current_dir, archive_dir), "\n")

  return(list(
    final_directory = normalizePath(file.path(current_dir, final_dir)),
    archive_directory = normalizePath(file.path(current_dir, archive_dir))
  ))
}

# プロジェクト設定ファイルの作成
create_project_config <- function(final_dir) {
  config_content <- '
# ===============================================================================
# 一般化ペアワイズ比較における代入法の比較研究
# プロジェクト設定ファイル
# ===============================================================================

project:
  name: "GPC Imputation Method Comparison Study"
  version: "2.0.0"
  authors: ["研究者名"]
  description: "区間打ち切りデータに対する一般化ペアワイズ比較の代入法比較"

structure:
  R/: "主要なRスクリプト"
  docs/: "ドキュメント"
  results/: "結果出力ディレクトリ（実行後生成）"
  presentation/: "発表資料"

main_files:
  simulation: "R/main_simulation.R"
  verification: "R/verification_test.R"
  data_generation: "R/data_generation.R"
  export: "R/looker_export.R"

dependencies:
  required: ["survival", "tidyverse", "ggplot2", "parallel"]
  optional: ["WINS", "BuyseTest"]

execution_order:
  1: "R/verification_test.R"
  2: "R/main_simulation.R"
  3: "R/looker_export.R"
'

  writeLines(config_content, file.path(final_dir, "project_config.yaml"))
}

# 実行ガイドの作成
create_execution_guide <- function(final_dir) {
  guide_content <- '# 🚀 実行ガイド

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
# フル条件でのシミュレーション
full_results <- run_final_simulation(
  n_sim = 1000,
  sample_sizes = c(100, 200, 400),
  K_values = c(3, 5, 10),
  dropout_levels = c("None", "Low", "Medium", "High"),
  effect_sizes = 1:24
)

# 結果分析
analysis <- analyze_final_results(full_results)
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
'

  writeLines(guide_content, file.path(final_dir, "QUICKSTART.md"))
}

# 実行
cat("===============================================================================\n")
cat("                    プロジェクトクリーンアップツール\n")
cat("===============================================================================\n")
cat("実行方法: cleanup_results <- cleanup_project()\n")
cat("===============================================================================\n")

if (!interactive()) {
  cleanup_results <- cleanup_project()
}