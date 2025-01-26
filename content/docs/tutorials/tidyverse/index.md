---
title: "Tidyverse"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# Tidyverse

## Tidyverseってなんぞや
tidyverseは，データ分析と可視化を効率的に行うために設計されたR言語のパッケージ群です．

## 練習問題

### 問題1

Tidyverseを用いてcsvファイルを読み込み，変数名dfに代入せよ

{{% details title="Answer" open=true %}}

### 解答1
```R
df <- read_csv("src/sample.csv")
```

{{% /details %}}
