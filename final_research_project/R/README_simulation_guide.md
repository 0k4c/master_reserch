# 一般化ペアワイズ比較を用いた区間打ち切りデータ対応の統計的手法比較研究

## 概要

本研究は、区間打ち切りデータ（interval-censored data）に対する統計的手法の比較を行うものです。特に、一般化ペアワイズ比較（Generalized Pairwise Comparison, GPC）手法を中心に、従来のRMST（Restricted Mean Survival Time）やログランク検定との検出力比較を実施します。

## 研究の背景と意義

### 区間打ち切りデータの課題
- 医学研究や信頼性工学において、正確なイベント発生時刻が観測できず、区間内での発生のみが分かる場合が多い
- 従来の生存解析手法は右側打ち切りデータを前提としており、区間打ち切りデータには適用が困難
- 代入法（imputation）による近似的手法が一般的だが、バイアスや検出力の低下が懸念される

### 本研究の新規性
1. **一般化ペアワイズ比較の区間打ち切りデータへの適用**
   - Net BenefitとWin Ratioの2つのアプローチ
   - 代入法を必要としない直接的な比較手法

2. **拡張EMI法（Enhanced Multiple Imputation）の提案**
   - 従来のEMI法の改良版
   - Beta分布を用いた適応的代入
   - 区間の幅と位置を考慮した柔軟な代入戦略

3. **包括的な検出力比較**
   - 8つの統計手法の同時比較
   - 24種類の生存分布での評価
   - 様々なサンプルサイズ、観測頻度、脱落率での検証

## ファイル構成

```
C:/Claude/serch/serch1/code/
├── generalized_pairwise_comparison_simulation.R  # メインシミュレーションコード
├── interval_censord_data_function.R              # 区間打ち切りデータ生成・代入法
├── distributions_2arm.R                          # 24種類の生存分布定義
├── README_simulation_guide.md                    # 本ファイル（研究ガイド）
└── simulation_results/                           # 結果出力ディレクトリ（実行後作成）
    ├── raw_simulation_results.csv
    ├── power_summary.csv
    ├── type1_error_summary.csv
    ├── comprehensive_results_table.csv
    └── [各種可視化ファイル].png
```

## 実装された統計手法

### 1. 一般化ペアワイズ比較（GPC）
- **Net Benefit法**: `(勝数 - 負数) / 全ペア数`
- **Win Ratio法**: `勝数 / 負数`
- 区間打ち切りデータに対する直接比較ルール：
  - 治療群の右端点 ≤ 対照群の左端点 → 治療群の勝ち
  - 対照群の右端点 ≤ 治療群の左端点 → 対照群の勝ち
  - 区間が重複 → 引き分け

### 2. RMST（制限平均生存時間）
- **中点代入法**: 区間の中点を代表値として使用
- **右点代入法**: 区間の右端点を代表値として使用
- **拡張EMI法**: Beta分布による適応的代入

### 3. ログランク検定
- **中点代入法**: 生存曲線の比較
- **右点代入法**: 保守的な生存曲線推定
- **拡張EMI法**: 複数代入による頑健な推定

## 拡張EMI法の技術的詳細

### アルゴリズム概要
```r
# 各区間 [L, R] に対して：
if (L == R) {
    # 正確な観測
    T = L
} else if (R == Inf) {
    # 右側打ち切り
    T = L + Exp(λ_estimated)  # 指数分布で補外
} else {
    # 区間打ち切り
    α = 2 + 1/(R-L)  # 幅が狭いほど集中
    β = 2 + 1/(R-L)
    U ~ Beta(α, β)
    T = L + U*(R-L)
}
```

### 従来法との違い
- **従来EMI**: 一様分布のみを使用
- **拡張EMI**: Beta分布により区間幅に応じた適応的代入
- **理論的根拠**: 区間が狭い場合はより確定的、広い場合はより不確実

## シミュレーション設定

### パラメータ設定
- **サンプルサイズ**: n = 100, 200, 400
- **観測頻度**: K = 3, 5, 10, 20回
- **脱落率**: None, Low (10%), Medium (20%), High (30%)
- **効果サイズ**: 24種類の生存分布パターン
- **シミュレーション回数**: 1000回（論文用）

### 生存分布パターン
1. **基本分布** (d=1-5): Weibull、Exponential、区分的生存関数
2. **比例ハザード** (d=6-10): 対照群に対する治療効果
3. **非比例ハザード** (d=11-24): 複雑な時間依存効果

## 実行方法

### 1. 環境準備
```r
# 必要なパッケージのインストール
install.packages(c("survival", "survRM2", "BuyseTest", "tidyverse",
                   "ggplot2", "gridExtra", "parallel", "foreach",
                   "doParallel", "knitr", "kableExtra"))

# ファイルの読み込み
source("generalized_pairwise_comparison_simulation.R")
```

### 2. テスト実行（開発・デバッグ用）
```r
# 小規模テスト（約10-15分）
test_results <- run_test_simulation()
```

### 3. フル分析実行（論文用）
```r
# 完全なシミュレーション（約3-6時間）
full_results <- run_full_analysis()
```

### 4. カスタム実行
```r
# パラメータを調整した実行
custom_results <- main_simulation_analysis(
    n_sim = 500,              # シミュレーション回数
    output_dir = "custom_results"  # 出力ディレクトリ
)
```

## 出力結果の解釈

### 1. 検出力比較（power_summary.csv）
- 各条件での各手法の検出力
- 効果サイズ別、サンプルサイズ別の比較
- **目標**: 拡張EMI法が従来法を上回ることの検証

### 2. 第1種の誤り（type1_error_summary.csv）
- 帰無仮説下（effect=1）での棄却率
- **目標**: 全手法で α=0.05 付近の値

### 3. 包括的結果表（comprehensive_results_table.csv）
- 論文用の要約統計
- 手法間の平均検出力と第1種の誤り率
- **主要評価指標**

### 4. 可視化ファイル
- **ヒートマップ**: 条件別検出力の視覚化
- **ボックスプロット**: 手法間の検出力分布比較
- **論文用図表**: 高解像度PNG形式

## 論文化に向けた分析のポイント

### 主要な検証項目
1. **GPC手法の有効性**
   - 区間打ち切りデータでの直接比較の優位性
   - Net BenefitとWin Ratioの性能差

2. **拡張EMI法の改善効果**
   - 従来の中点・右点代入法との比較
   - Beta分布代入の理論的・実証的優位性

3. **手法間の検出力ランキング**
   - 条件別の最適手法の特定
   - サンプルサイズ・脱落率・効果サイズとの関係

4. **実用的な推奨事項**
   - 実際の研究での手法選択指針
   - 計算コストと性能のトレードオフ

### 統計的検定
```r
# 手法間の検出力差の有意性検定
wilcox.test(power_gpc, power_rmst, paired = TRUE)

# 条件別の分散分析
aov(Power ~ Method + SampleSize + DropoutRate, data = results)
```

## 期待される研究成果

### 理論的貢献
1. 区間打ち切りデータに対するGPC手法の理論的基盤
2. 拡張EMI法の統計的性質の解明
3. 代入法の選択がバイアス・分散に与える影響の定量化

### 実用的貢献
1. 医学研究・信頼性工学での手法選択指針
2. 既存ソフトウェアへの実装可能なアルゴリズム
3. サンプルサイズ設計への応用

### 発表・論文化計画
1. **国内学会発表**: 日本統計学会、応用統計学会
2. **国際論文投稿**: Statistics in Medicine, Biometrics
3. **Rパッケージ化**: CRAN登録による一般化

## 技術的注意事項

### 計算効率化
- 並列処理による高速化（detectCores()で最適化）
- メモリ効率的なデータ構造
- 大規模シミュレーション対応

### 再現可能性
- set.seed()による結果の再現性確保
- 全パラメータのログ記録
- バージョン管理対応

### エラーハンドリング
- 数値的不安定性への対応
- 収束しない場合のフォールバック
- 欠損データの適切な処理

## トラブルシューティング

### よくある問題
1. **メモリ不足**: `n_sim`を小さくするか、`n_cores`を調整
2. **収束エラー**: 一部の分布で数値的問題が発生する場合あり
3. **計算時間**: 完全実行は数時間必要

### デバッグ方法
```r
# 単一条件でのテスト
debug_result <- run_simulation_study(
    n_sim = 10,
    sample_sizes = 100,
    K_values = 5,
    dropout_levels = "None",
    effect_sizes = 1
)
```

## 参考文献

1. Buyse, M. (2010). Generalized pairwise comparisons of prioritized outcomes. Statistics in Medicine.
2. Sun, J. (2006). The Statistical Analysis of Interval-censored Failure Time Data.
3. Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys.
4. Wei, L.J. et al. (2017). Alternative approaches for comparing unadjusted analyses.

---

**連絡先**: Claude Code Research Team
**更新日**: 2025-01-16
**バージョン**: 1.0