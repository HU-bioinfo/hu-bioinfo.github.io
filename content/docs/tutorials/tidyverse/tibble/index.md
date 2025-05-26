---
title: "tibble"
weight: 6
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# tibble

tibbleって知ってるか．使ってみろよ．便利だぞ．
何が便利か知りてえか．まあ，僕もよく分かってねえよ．

Rの基本文法にあったように，tibbleなるデータ型があります．これはExcelみたいな表形式で表されるデータ型です．似たようなデータ型にDataframeがあったと思います．こいつとの違いは追々話すとして，大体高速だし，みんな使ってるし，データ解析するときデータは大抵表形式だし，何か表形式のデータを弄くり回す時は，データをtibble形式にして使います．

ここでは実際にcsvファイル(Excelファイルの廉価版みたいなやつ)を用いてtibbleの使い方を実践的に学んでいきましょう．

この[リンク](student_sample.csv)からファイルをダウンロードするか，あるいは下の内容をコピペして`student_sample.csv`として保存してください．
```csv
student_id,name,gender,math_score,english_score,science_score,passed
101,John Doe,M,78,85,69,TRUE
102,Jane Smith,F,92,89,95,TRUE
103,Bob Lee,M,56,72,60,FALSE
104,Alice Wong,F,88,91,84,TRUE
105,Chris Johnson,M,45,50,48,FALSE
106,Maria Garcia,F,74,80,70,TRUE
107,David Kim,M,66,59,73,TRUE
108,Nancy Allen,F,81,77,88,TRUE
109,Tom Brown,M,39,42,35,FALSE
110,Lisa White,F,90,94,92,TRUE
```
このcsvファイルは学生の成績を表形式にまとめたデータです．

```R
library(tidyverse)

df <- readr::read_csv("data/student_sample.csv")
print(df)
```
以上のコードでは，tibbleを含めtidyverse全体を読み込み，`readr`パッケージの`read_csv`関数を用いてcsvファイルを読み込みつつtibbleに変換したものを`df`に代入しました．

出力された結果を見てみると，この表形式のtibbleは10行×7列のデータであり，`student_id`列には`<dbl>`型(数値型を意味します)．
