# 一般化ペアワイズ比較における代入法の比較研究

## 🎯 研究概要

本研究は、区間打ち切りデータに対する一般化ペアワイズ比較（Generalized Pairwise Comparison, GPC）において、**どの代入法が最も優れているか**を明らかにする統計的比較研究です。

### 主要な研究目的

1. **Q1: GPCにおける最適代入法の特定**
   - Direct GPC（代入なし）
   - Midpoint Assignment + GPC
   - Rightpoint Assignment + GPC
   - Enhanced EMI + GPC

2. **Q2: 従来手法との性能比較**
   - RMST（Restricted Mean Survival Time）
   - ログランク検定

## 🚀 主要な成果

✅ **p値計算問題の完全解決**: すべての統計手法で適切なp値（0-1範囲）を計算
✅ **理論的統計手法の実装**: 正しい帰無仮説・対立仮説に基づく検定
✅ **8手法の同時比較**: GPC×4代入法 + 従来手法×2
✅ **包括的シミュレーション**: 24種類の生存分布での検証
✅ **論文品質の実装**: 再現可能で信頼性の高いコード

## 📁 ファイル構成

### 🔴 **メインファイル**（必須）
```
final_imputation_comparison_study.R    # メイン研究コード
final_verification_test.R              # 動作確認テスト
flexible_data_generation.R             # データ生成（制限緩和版）
```

### 🔵 **サポートファイル**
```
interval_censord_data_function.R       # 元の代入法関数
distributions_2arm.R                   # 24種類の生存分布
improved_gpc_with_packages.R           # WINSとBuyseTest活用版
```

### 📊 **出力ファイル**（実行後自動生成）
```
results/
├── simulation_results.csv             # 全シミュレーション結果
├── power_summary.csv                  # 検出力分析結果
└── type1_error_summary.csv            # 第1種誤り分析結果
```

## ⚡ クイックスタート

### 1. 動作確認（必須）
```r
source("final_verification_test.R")
verification_results <- run_final_verification()
```

### 2. 小規模研究実行
```r
source("main_simulation.R")
test_results <- run_final_study()
```

### 3. 本格研究実行（論文用）
```r
# フル条件でのシミュレーション
full_results <- run_final_simulation(
  n_sim = 1000,                              # シミュレーション回数
  sample_sizes = c(100, 200, 400),           # サンプルサイズ
  K_values = c(3, 5, 10),                    # 観測回数
  dropout_levels = c("None", "Low", "Medium", "High"),  # 脱落率
  effect_sizes = 1:24                        # 24種類の効果サイズ
)

# 結果分析
analysis <- analyze_final_results(full_results)
```

### 4. `run_final_simulation` 関数の引数説明

- `n_sim`: 各シナリオにおけるシミュレーションの繰り返し回数。（デフォルト: 1000）
- `sample_sizes`: 1回のシミュレーションで生成する総サンプルサイズ（治療群と対照群の合計）。（デフォルト: `c(100, 200)`）
- `K_values`: 観察期間内の来院（観察）回数。（デフォルト: `c(3, 5)`）
- `dropout_levels`: 来院脱落（ドロップアウト）率のレベル。`"None"`, `"Low"`, `"Medium"`, `"High"`から選択。（デフォルト: `c("None", "Medium")`）
- `effect_sizes`: 治療効果の大きさを示す分布タイプ。1が効果なし（対照群と同じ）、数字が大きいほど効果も大きい。（デフォルト: `c(1, 6, 12)`）
- `alpha`: 仮説検定における有意水準。（デフォルト: 0.05）
- `n_cores`: 計算に使用するCPUコア数。`NULL`の場合は自動で最大-1のコアを検出。（デフォルト: `NULL`）
- `use_parallel`: 並列処理を行うかどうか。（デフォルト: `TRUE`）
- `show_progress`: 実行中にプログレスバーを表示するかどうか。（デフォルト: `TRUE`）

---

(以降のセクションは変更なしのため省略)