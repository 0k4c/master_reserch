# 一般化ペアワイズ比較における代入法の比較研究

##  研究概要

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

##  主要な成果

✅ **p値計算問題の完全解決**: すべての統計手法で適切なp値（0-1範囲）を計算
✅ **理論的統計手法の実装**: 正しい帰無仮説・対立仮説に基づく検定
✅ **8手法の同時比較**: GPC×4代入法 + 従来手法×2
✅ **包括的シミュレーション**: 24種類の生存分布での検証
✅ **論文品質の実装**: 再現可能で信頼性の高いコード

##  ファイル構成

###  **メインファイル**（必須）
```
final_imputation_comparison_study.R    # メイン研究コード
final_verification_test.R              # 動作確認テスト
flexible_data_generation.R             # データ生成（制限緩和版）
```

###  **サポートファイル**
```
interval_censord_data_function.R       # 元の代入法関数
distributions_2arm.R                   # 24種類の生存分布
improved_gpc_with_packages.R           # WINSとBuyseTest活用版
```

###  **出力ファイル**（実行後自動生成）
```
results/
├── simulation_results.csv             # 全シミュレーション結果
├── power_summary.csv                  # 検出力分析結果
└── type1_error_summary.csv            # 第1種誤り分析結果
```

##  クイックスタート

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

##  研究手法詳細

### 統計手法（8手法を同時比較）

| カテゴリ | 手法 | 代入法 | 統計量 |
|----------|------|--------|--------|
| **GPC** | Direct | なし | Net Benefit, Win Ratio |
| **GPC** | Midpoint | 区間中点 | Net Benefit, Win Ratio |
| **GPC** | Rightpoint | 区間右端点 | Net Benefit, Win Ratio |
| **GPC** | Enhanced EMI | 均等間隔代入（論文手法） | Net Benefit, Win Ratio |
| **従来法** | RMST | 中点代入 | 平均差 |
| **従来法** | ログランク | 中点代入 | カイ二乗統計量 |

### 評価指標

- **検出力（Power）**: 真の効果を正しく検出する確率
- **第1種誤り（Type I Error）**: 効果がないのに有意と判定する確率
- **計算効率**: 実行時間とメモリ使用量
- **頑健性**: 様々な条件下での安定性

##  理論的背景

### 生存時間解釈
- **前提**: 生存時間が長い = 良い結果
- **治療効果**: 治療群の生存時間 > 対照群の生存時間

### 区間打ち切りデータの比較ルール
治療群 [L_T, R_T] vs 対照群 [L_C, R_C] において：
- **L_T > R_C**: 治療群の勝ち（明確に優れている）
- **R_T < L_C**: 対照群の勝ち（明確に劣っている）
- **区間重複**: 引き分け（判定不能）

### 統計的仮説
- **H0**: Net Benefit = 0, Win Ratio = 1（差がない）
- **H1**: Net Benefit ≠ 0, Win Ratio ≠ 1（差がある）
- **検定**: 両側検定、正規近似によるp値計算

##  トラブルシューティング

### よくある問題と解決法

1. **パッケージエラー**
```r
# 必要パッケージの一括インストール
install.packages(c("survival", "tidyverse", "ggplot2", "parallel", "WINS", "BuyseTest"))
```

2. **メモリ不足**
```r
# シミュレーション回数を減らす
run_final_simulation(n_sim = 100)  # 1000 → 100
```

3. **実行時間が長い**
```r
# 条件を絞る
run_final_simulation(
  sample_sizes = 100,        # 1つのサイズのみ
  effect_sizes = c(1, 6, 12) # 3つの効果サイズのみ
)
```

### サポート情報
- **GitHub Issues**: バグ報告・機能要求
- **研究相談**: 統計手法に関する質問
- **コード貢献**: プルリクエスト歓迎

##  参考文献

1. Buyse, M. (2010). Generalized pairwise comparisons of prioritized outcomes. *Statistics in Medicine*, 29(30), 3245-3257.
2. Peron, J. et al. (2018). The net chance of a longer survival as a patient-oriented measure of treatment benefit. *Statistics in Medicine*, 37(16), 2343-2365.
3. Sun, J. (2006). *The Statistical Analysis of Interval-censored Failure Time Data*. Springer.
4. Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.

---

**Last Updated**: 2025-01-16
**Version**: 2.0.0
**Status**:  Production Ready

---

## GitとGitHubの使い方 (簡易ガイド)

### 1. GitとGitHubとは？

- **Git**: あなたのPC上で動く「**コードのバージョン管理システム**」です。ファイルの変更履歴を「スナップショット」として記録し、いつでも過去の状態に戻したり、変更内容を確認したりできます。
- **GitHub**: Gitで作成した変更履歴を保存しておくための**Web上の保管場所**です。PCが故障してもコードが安全に保たれるバックアップとしての役割や、他の人とコードを共有・共同編集する場にもなります。

### 2. 基本的な作業サイクル

今後、コードを編集してその変更をGitHubに保存したいときは、基本的に以下の3ステップを繰り返します。

#### ステップ①: `git add` (変更をカゴに入れる)

ファイルに加えた変更を、「次のスナップショットに含めます」とGitに教える作業です。

```shell
# すべての変更を一度に追加する場合（こちらが便利です）
git add .
```

#### ステップ②: `git commit` (スナップショットを撮る)

`add`したファイル群のスナップショットを撮影し、変更履歴として記録します。このとき、`-m`に**何を変更したか**のメモ（コミットメッセージ）を必ず残します。

```shell
git commit -m "〇〇の機能を追加"
```

#### ステップ③: `git push` (GitHubにアップロード)

PC上に記録した変更履歴を、Web上のGitHubリポジトリにアップロードします。これで初めて、オンラインでのバックアップが完了します。

```shell
git push
```

この**「add → commit → push」**が、日々の基本的な作業サイクルになります。

### 3. 便利なコマンド

- **`git status`**: どのファイルが変更されたか、どのファイルが`add`されているかなど、現在の状態を確認できます。作業に迷ったらまずこれを実行してみてください。
- **`git log`**: 今までのコミット履歴（変更の記録）を一覧で表示します。
