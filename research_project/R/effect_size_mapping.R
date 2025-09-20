# Effect Size マッピング関数
# 効果サイズ番号から治療群・対照群の分布情報を取得

get_effect_size_info <- function(effect_size) {
  # 効果サイズ番号に対応する治療群・対照群の分布情報
  effect_mapping <- data.frame(
    effect_size = 1:24,

    # 対照群の分布
    control_distribution = c(
      "T1 (Weibull, scale=1, shape=1)",    # 1:
      "T2 (Weibull, scale=1, shape=0.5)",  # 2:
      "T3 (Weibull, scale=1, shape=2)",    # 3:
      "T4 (Piecewise exponential)",        # 4:
      "T5 (Piecewise mixed)",              # 5:
      "T1 (Weibull, scale=1, shape=1)",    # 6: Control=T1
      "T2 (Weibull, scale=1, shape=0.5)",  # 7: Control=T2
      "T3 (Weibull, scale=1, shape=2)",    # 8: Control=T3
      "T4 (Piecewise exponential)",        # 9: Control=T4
      "T5 (Piecewise mixed)",              # 10: Control=T5
      "T11 (Piecewise)",                   # 11:
      "T13 (Piecewise [0,0.2])",           # 12:
      "T15 (Complex exponential)",         # 13:
      "T17 (Complex exponential)",         # 14:
      "T19 (Piecewise [0,0.2])",           # 15:
      "T21 (Piecewise [0,0.25])",          # 16:
      "T23 (Piecewise [0,0.8])",           # 17:
      "T1 (Weibull, scale=1, shape=1)",    # 18:
      "T2 (Weibull, scale=1, shape=0.5)",  # 19:
      "T3 (Weibull, scale=1, shape=2)",    # 20:
      "T4 (Piecewise exponential)",        # 21:
      "T5 (Piecewise mixed)",              # 22:
      "T11 (Piecewise)",                   # 23:
      "T13 (Piecewise [0,0.2])"            # 24:
    ),

    # 治療群の分布
    treatment_distribution = c(
      "Same as Control (T1)",              # 1: 効果なし
      "Same as Control (T2)",              # 2: 効果なし
      "Same as Control (T3)",              # 3: 効果なし
      "Same as Control (T4)",              # 4: 効果なし
      "Same as Control (T5)",              # 5: 効果なし
      "T6 (Weibull, scale=0.5, shape=1)", # 6: Treatment=T6
      "T7 (Modified Weibull)",             # 7: Treatment=T7
      "T8 (Modified Weibull)",             # 8: Treatment=T8
      "T9 (Modified piecewise)",           # 9: Treatment=T9
      "T10 (Modified piecewise)",          # 10: Treatment=T10
      "T12 (Modified piecewise)",          # 11:
      "T14 (Modified piecewise)",          # 12:
      "T16 (Complex exponential)",         # 13:
      "T18 (Modified complex)",            # 14:
      "T20 (Modified piecewise)",          # 15:
      "T22 (Weibull, scale=0.5, shape=1)", # 16:
      "T24 (Modified piecewise)",          # 17:
      "T6 (Weibull, scale=0.5, shape=1)", # 18:
      "T7 (Modified Weibull)",             # 19:
      "T8 (Modified Weibull)",             # 20:
      "T9 (Modified piecewise)",           # 21:
      "T10 (Modified piecewise)",          # 22:
      "T12 (Modified piecewise)",          # 23:
      "T14 (Modified piecewise)"           # 24:
    ),

    # 効果の大きさ
    effect_magnitude = c(
      "None (Null hypothesis)",            # 1:
      "None (Null hypothesis)",            # 2:
      "None (Null hypothesis)",            # 3:
      "None (Null hypothesis)",            # 4:
      "None (Null hypothesis)",            # 5:
      "Small-Medium",                      # 6:
      "Small-Medium",                      # 7:
      "Small-Medium",                      # 8:
      "Small-Medium",                      # 9:
      "Small-Medium",                      # 10:
      "Medium",                            # 11:
      "Medium",                            # 12:
      "Medium-Large",                      # 13:
      "Medium-Large",                      # 14:
      "Medium-Large",                      # 15:
      "Large",                             # 16:
      "Large",                             # 17:
      "Small-Medium (Alt)",                # 18:
      "Small-Medium (Alt)",                # 19:
      "Small-Medium (Alt)",                # 20:
      "Small-Medium (Alt)",                # 21:
      "Small-Medium (Alt)",                # 22:
      "Medium (Alt)",                      # 23:
      "Medium (Alt)"                       # 24:
    ),

    # 群の設定番号（分かりやすい表示）
    group_setting = c(
      "Control=1, Treatment=1 (No effect)", # 1:
      "Control=2, Treatment=2 (No effect)", # 2:
      "Control=3, Treatment=3 (No effect)", # 3:
      "Control=4, Treatment=4 (No effect)", # 4:
      "Control=5, Treatment=5 (No effect)", # 5:
      "Control=1, Treatment=6",             # 6:
      "Control=2, Treatment=7",             # 7:
      "Control=3, Treatment=8",             # 8:
      "Control=4, Treatment=9",             # 9:
      "Control=5, Treatment=10",            # 10:
      "Control=11, Treatment=12",           # 11:
      "Control=13, Treatment=14",           # 12:
      "Control=15, Treatment=16",           # 13:
      "Control=17, Treatment=18",           # 14:
      "Control=19, Treatment=20",           # 15:
      "Control=21, Treatment=22",           # 16:
      "Control=23, Treatment=24",           # 17:
      "Control=1, Treatment=6 (Alt)",       # 18:
      "Control=2, Treatment=7 (Alt)",       # 19:
      "Control=3, Treatment=8 (Alt)",       # 20:
      "Control=4, Treatment=9 (Alt)",       # 21:
      "Control=5, Treatment=10 (Alt)",      # 22:
      "Control=11, Treatment=12 (Alt)",     # 23:
      "Control=13, Treatment=14 (Alt)"      # 24:
    ),

    stringsAsFactors = FALSE
  )

  # 指定された効果サイズの情報を返す
  if (effect_size %in% 1:24) {
    return(effect_mapping[effect_size, ])
  } else {
    return(data.frame(
      effect_size = effect_size,
      control_distribution = "Unknown",
      treatment_distribution = "Unknown",
      effect_magnitude = "Unknown",
      group_setting = paste("Unknown effect size:", effect_size),
      stringsAsFactors = FALSE
    ))
  }
}

# 効果サイズ情報の表示関数
display_effect_size_info <- function(effect_sizes) {
  cat("=== Effect Size Information ===\n")
  for (effect in effect_sizes) {
    info <- get_effect_size_info(effect)
    cat(sprintf("Effect Size %d: %s\n", effect, info$group_setting))
    cat(sprintf("  Control Group: %s\n", info$control_distribution))
    cat(sprintf("  Treatment Group: %s\n", info$treatment_distribution))
    cat(sprintf("  Effect Magnitude: %s\n", info$effect_magnitude))
    cat("\n")
  }
}

# 結果データフレームに効果サイズ情報を追加する関数
add_effect_size_info_to_results <- function(results_df) {
  # 結果データフレームに効果サイズ情報を追加
  if ("effect" %in% names(results_df)) {
    effect_info <- lapply(results_df$effect, get_effect_size_info)

    # 新しい列を追加
    results_df$control_group <- sapply(effect_info, function(x) x$control_distribution)
    results_df$treatment_group <- sapply(effect_info, function(x) x$treatment_distribution)
    results_df$effect_magnitude <- sapply(effect_info, function(x) x$effect_magnitude)
    results_df$group_setting <- sapply(effect_info, function(x) x$group_setting)

    # 列の順序を調整（効果サイズ関連情報を効果サイズの後に配置）
    col_order <- names(results_df)
    effect_idx <- which(col_order == "effect")

    if (length(effect_idx) > 0) {
      new_order <- c(
        col_order[1:effect_idx],
        "group_setting", "control_group", "treatment_group", "effect_magnitude",
        col_order[(effect_idx+1):length(col_order)]
      )
      # 重複を除去
      new_order <- new_order[!duplicated(new_order)]
      results_df <- results_df[, new_order]
    }
  }

  return(results_df)
}