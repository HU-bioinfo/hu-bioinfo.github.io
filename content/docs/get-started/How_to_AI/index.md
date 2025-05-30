---
title: "AI学習の進め方"
description: "AI学習の進め方"
weight: 5
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# AI学習の進め方

## 0. 課題を解く前に

- [Cursor](https://www.cursor.com/)のようなAIエディタは非常に効率的にAIに質問したりコードを書かせることができます。
- 初学者でも簡単なタスクであればＡＩによって簡単にコードを書くことができます。
- しかし、初学者の学習にあたって重要なことは**ＡＩが出力するコードを見て学ぶ**ことです。
- 以下の流れで学習を進めていきましょう。

1. **まずは課題をそのままAIに聞いてみる。**
2. **分からない部分をさらに深堀して聞く。**
3. **疑問が解消できた場所から実行してみる。**
4. **日常のタスクに使ってみる。**

- この学習方法の具体例を以下に示します。

---

## 例題

- 大腸がんと正常大腸粘膜の間でDEG(Differential Expression Gene)解析を行う。
- 以下のデータを使用する
- DEseq2というパッケージを使用する。

{{% highlight r %}}
# データの確認
counts_matrix[1:5, 1:5]
{{% /highlight %}}

{{% highlight r %}}
        TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
A1BG               1            4            0            5            5
A1CF             786         2748         1989          646          657
A2M             2342         6609        14075         4413        38490
A2ML1              1            9           17            7            8
A3GALT2            0            1            1            0            0
{{% /highlight %}}

{{% highlight r %}}
head(sample_type)
{{% /highlight %}}

{{% highlight r %}}
# A tibble: 6 × 2
  case_id           sample_type
  <chr>             <chr>      
1 TCGA-D5-5540      tumor      
2 TCGA-EI-6509      tumor      
3 TCGA-A6-6137      tumor      
4 TCGA-QG-A5Z2      tumor      
5 TCGA-AA-3489      tumor      
6 HCM-CSHL-0160-C18 tumor      
{{% /highlight %}}

{{% highlight r %}}
table(sample_type$sample_type)
{{% /highlight %}}

{{% highlight r %}}
normal  tumor 
    50     50 
{{% /highlight %}}

---

## 1.まずは課題をそのままAIに聞いてみる。

{{% hint info %}}

- `Manual`モード、モデルは`gemini-2.5-Pro`を使用します。
- CursorやAIチャットの使い方は[Cursorの使い方]({{% ref "/docs/get-started/Cursor/index.md" %}})を参考にしてください。

{{% /hint %}}

- AIに課題をそのまま聞いてみます。初学者向けの解説も聞いてみます。
- ただし、使うデータも教えてあげる必要があるので、今回提示されたデータをそのままコピペしましょう。

```
・大腸がんと正常大腸粘膜の間でDEG解析を行って。
・DEseq2というパッケージを使用して。
・counts_matrixとsample_typeを使用して。データの一部を以下に示す。
・初学者にもわかりやすく解説して。

counts_matrix[1:5, 1:5]
        TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
A1BG               1            4            0            5            5
A1CF             786         2748         1989          646          657
A2M             2342         6609        14075         4413        38490
A2ML1              1            9           17            7            8
A3GALT2            0            1            1            0            0
head(sample_type)
# A tibble: 6 × 2
  case_id           sample_type
  <chr>             <chr>      
1 TCGA-D5-5540      tumor      
2 TCGA-EI-6509      tumor      
3 TCGA-A6-6137      tumor      
4 TCGA-QG-A5Z2      tumor      
5 TCGA-AA-3489      tumor      
6 HCM-CSHL-0160-C18 tumor      
table(sample_type$sample_type)
normal  tumor 
    50     50 
```

- 実際に返ってきた結果が以下です。長いのでたためるようにしておきます。
- このように解説付きでコードを出力してくれます。
- 始めは長文が返ってくると読み飛ばしたくなりますが、**必ず1行ずつ読みましょう**。
- （今回は例なので一度読み飛ばして構いません。）

{{% details "AIの回答" %}}

```

# 0. 必要なライブラリを読み込みます
# DESeq2: RNA-Seqデータの発現変動解析を行うための主要なパッケージです
library(DESeq2)
# dplyr: データの操作や整形を容易にするためのパッケージです（sample_typeがtibbleの場合に便利です）
library(dplyr)

# 1. データの準備
# ここでは、ユーザーが提供してくれたcounts_matrixとsample_typeが既にRの環境に読み込まれていると仮定します。
# もしファイルから読み込む場合は、以下のようなコードになります。
# counts_matrix <- read.csv("path/to/your/counts_matrix.csv", row.names = 1)
# sample_type <- read.csv("path/to/your/sample_type.csv")

# 簡単なダミーデータを作成します（実際のデータがある場合はこの部分は不要です）
# 遺伝子名を模倣
gene_names <- c("A1BG", "A1CF", "A2M", "A2ML1", "A3GALT2", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10")
# サンプル名を模倣
sample_names_tumor <- paste0("TCGA-T", 1:5) # 5つの腫瘍サンプル
sample_names_normal <- paste0("TCGA-N", 1:5) # 5つの正常サンプル
all_sample_names <- c(sample_names_tumor, sample_names_normal)

# ランダムなカウントデータを生成 (整数値)
set.seed(123) # 再現性のために乱数シードを設定
counts_matrix <- matrix(
  rnbinom(n = length(gene_names) * length(all_sample_names), mu = 100, size = 10),
  nrow = length(gene_names),
  ncol = length(all_sample_names)
)
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- all_sample_names

# サンプル情報データフレームを作成
sample_type <- data.frame(
  case_id = all_sample_names,
  sample_type = factor(rep(c("tumor", "normal"), each = 5)) # 最初の5つがtumor, 次の5つがnormal
)

# ユーザー提供データ形式に合わせる (もしtibbleなら)
# library(tibble)
# sample_type <- as_tibble(sample_type)


# 2. DESeqDataSetオブジェクトの作成
# DESeq2で解析を行うためには、まずカウントデータとサンプル情報を行列形式でまとめ、
# DESeqDataSetオブジェクトという専用の形式に変換する必要があります。

# sample_typeのcase_id列がcounts_matrixの列名と一致しているか、
# また、その順序が同じであることを確認します。これは非常に重要です。
# もし異なっている場合は、並び替えが必要です。
# 今回のダミーデータでは作成時に一致させていますが、実際のデータでは注意してください。
if (!all(colnames(counts_matrix) == sample_type$case_id)) {
  # sample_typeの行の順番をcounts_matrixの列の順番に合わせる
  sample_type <- sample_type[match(colnames(counts_matrix), sample_type$case_id), ]
}

# sample_typeのsample_type列を因子型(factor)に変換し、比較の基準レベルを指定します。
# ここでは"normal"を基準として、"tumor"が"normal"に対してどう変動したかを見ます。
sample_type$sample_type <- factor(sample_type$sample_type, levels = c("normal", "tumor"))

# DESeqDataSetオブジェクトを作成します。
# countData: 遺伝子発現量のマトリックス (行が遺伝子、列がサンプル)
# colData: サンプルの情報 (行がサンプル、列がサンプルの属性)
# design: 解析のデザイン式。ここではsample_typeの違いによって発現量が変動するかを見たいので、"~ sample_type"とします。
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sample_type,
                              design = ~ sample_type)

# ddsオブジェクトが正しく作成されたか確認 (任意)
# print(dds)

# 3. 低発現遺伝子のフィルタリング (推奨)
# 発現量が極端に低い遺伝子は、統計的な検出力が低く、解析のノイズになることがあります。
# 例えば、「少なくともX個のサンプルでYカウント以上の遺伝子のみを残す」といった基準でフィルタリングします。
# ここでは、少なくとも3サンプルで10カウント以上の遺伝子を残す例を示します。
# (ユーザーのデータサイズや特性に応じて調整してください)
keep <- rowSums(counts(dds) >= 10) >= 3 # 今回は5サンプルずつなので、3は妥当な値
dds <- dds[keep, ]

# 4. DESeq関数の実行
# この関数一つで、以下の主要なステップが実行されます。
#   1. 正規化 (Estimate Size Factors): ライブラリサイズの違いを補正します。
#   2. 分散の推定 (Estimate Dispersions): 遺伝子ごとの発現量のばらつきを推定します。
#   3. モデルフィッティングと統計検定 (Negative Binomial GLM):
#      負の二項分布を用いた一般化線形モデルを各遺伝子に適用し、
#      サンプルタイプ間の発現差の統計的有意性を検定します。
dds <- DESeq(dds)

# 5. 解析結果の取得
# results関数で、指定した比較グループ間の発現変動解析の結果を取得します。
# contrast引数で比較する条件を指定します。
# c("要因名", "比較したい群", "基準となる群") の順で指定します。
# この場合、"sample_type"という要因について、"tumor"が"normal"に対してどう変動したかを見ます。
res <- results(dds, contrast = c("sample_type", "tumor", "normal"))

# 結果の概要を表示 (任意)
# summary(res)
# print(head(res))

# 6. 結果の解釈と出力
# resオブジェクトには、各遺伝子について以下の様な情報が含まれます。
#   - baseMean: 全サンプルでの正規化後カウント数の平均値。発現量の目安。
#   - log2FoldChange: 発現量の変化率の対数（底2）。
#                     正の値ならtumorで発現量が増加、負の値ならtumorで発現量が減少。
#                     例えば、log2FoldChange = 1 は2倍の変化、-1は1/2倍の変化を意味します。
#   - lfcSE: log2FoldChangeの標準誤差。
#   - stat: Wald検定統計量。
#   - pvalue: p値。観測された発現差が偶然生じる確率。
#   - padj: 調整済みp値 (FDR: False Discovery Rate)。多重検定の問題を補正したp値。
#           一般的に、このpadjを使って有意な遺伝子を判断します (例: padj < 0.05)。

# 有意な発現変動遺伝子を抽出する例
# padj < 0.05 かつ |log2FoldChange| > 1 (つまり2倍以上の変動) を基準とする場合
significant_degs <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
# significant_degs <- as.data.frame(significant_degs) # データフレームとして扱いやすくする
# significant_degs <- significant_degs[order(significant_degs$padj), ] # padjでソート

# 結果をソートして表示 (log2FoldChangeの絶対値が大きい順、またはpadjが小さい順など)
res_ordered <- res[order(res$padj), ] # padjが小さい順にソート
print("解析結果 (padjでソート):")
print(head(res_ordered))

# 有意な遺伝子のみを抽出 (padj < 0.05 の例)
res_significant <- subset(res_ordered, padj < 0.05)
print("有意な発現変動遺伝子 (padj < 0.05):")
print(head(res_significant))

# 結果をCSVファイルとして保存する場合
# write.csv(as.data.frame(res_ordered), file = "DEG_results_all.csv")
# write.csv(as.data.frame(res_significant), file = "DEG_results_significant.csv")

cat("\nDEG解析が完了しました。\n")
cat("結果は 'res' オブジェクトに格納されています。\n")
cat("有意な遺伝子 (padj < 0.05) は 'res_significant' オブジェクトにあります。\n")

```

{{% /details %}}

## 2. 分からない部分をさらに深堀して聞く

- AIはかなり細かいところまで解説を書いてくれるため、順番に読んでいけばだんだんと理解が深まっていきます。
- 例えば#1. データの準備 の項目について。今回はcounts_matrixやsample_typeは自分で読み込んで用意しているので、#1の部分は丸々消してしまっても問題ないことが分かります。
- しかし、**一度読んだだけではよくわからないコードや解説**も必ずあると思います。
- そういう時はわからない部分をコピペして**分からないところだけAIに聞いてみましょう**。

```
# sample_typeのcase_id列がcounts_matrixの列名と一致しているか、
# また、その順序が同じであることを確認します。これは非常に重要です。
# もし異なっている場合は、並び替えが必要です。
# 今回のダミーデータでは作成時に一致させていますが、実際のデータでは注意してください。
if (!all(colnames(counts_matrix) == sample_type$case_id)) {
  # sample_typeの行の順番をcounts_matrixの列の順番に合わせる
  sample_type <- sample_type[match(colnames(counts_matrix), sample_type$case_id), ]
}

上記の部分で具体的に並びを確認したり並べ替えたりしているのかが分からないのでわかりやすく解説して。
```

{{% details "AIの回答" %}}

```
この部分のコードは、DESeq2の解析において非常に重要な前処理ステップです。なぜなら、DESeqDataSetFromMatrix関数は、counts_matrix（遺伝子のカウントデータ）の各列と、colData（サンプル情報、ここではsample_type）の各行が、同じサンプルに対応していると仮定して処理を進めるからです。この対応付けは、それぞれの順番に依存します。
もし、counts_matrixの列の並び順（例えば、サンプルA, サンプルB, サンプルCの順）と、sample_typeの行の並び順（例えば、サンプルAの情報, サンプルCの情報, サンプルBの情報）が異なっていると、DESeq2はサンプルBのカウントデータをサンプルCの情報と関連付けてしまうなど、誤った解析を行ってしまいます。これでは正しい結果が得られません。
そこで、このコードブロックでは以下のことを行っています。
colnames(counts_matrix):
これはcounts_matrixの全ての列名（サンプルID）を順番に取り出し、文字列のベクトル（配列のようなもの）として返します。
例: counts_matrixが TCGA-T1, TCGA-T2, TCGA-N1 の順で列を持っている場合、c("TCGA-T1", "TCGA-T2", "TCGA-N1") というベクトルが返ります。
sample_type$case_id:
これはsample_typeというデータフレーム（表のようなもの）の中のcase_idという列の値を全て順番に取り出し、ベクトルとして返します。
例: sample_typeが以下のような内容だったとします。
Apply to DEG.R
          case_id sample_type
        1 TCGA-T1       tumor
        2 TCGA-N1      normal
        3 TCGA-T2       tumor
この場合、sample_type$case_idは c("TCGA-T1", "TCGA-N1", "TCGA-T2") というベクトルになります。
colnames(counts_matrix) == sample_type$case_id:
上記1と2で取得した2つのベクトルを要素ごとに比較します。
例:
colnames(counts_matrix): c("TCGA-T1", "TCGA-T2", "TCGA-N1")
sample_type$case_id: c("TCGA-T1", "TCGA-N1", "TCGA-T2")
比較結果:
"TCGA-T1" == "TCGA-T1" -> TRUE
"TCGA-T2" == "TCGA-N1" -> FALSE
"TCGA-N1" == "TCGA-T2" -> FALSE
結果として、c(TRUE, FALSE, FALSE) という論理値のベクトルが返ります。
all( ... ):
括弧の中の論理値ベクトルが全てTRUEである場合にのみTRUEを返します。一つでもFALSEが含まれていればFALSEを返します。
例: all(c(TRUE, FALSE, FALSE)) は FALSE になります。
これは、「counts_matrixの列名とsample_typeのcase_idが、その順番も含めて完全に一致しているか」を確認しています。
!all( ... ):
! は論理値を反転させる演算子です（NOT）。
!FALSE は TRUE に、 !TRUE は FALSE になります。
つまり、!all(colnames(counts_matrix) == sample_type$case_id) は、「counts_matrixの列名とsample_typeのcase_idが、順番も含めて完全に一致していない場合」に TRUE となります。
if (!all( ... )) { ... }:
このif文は、条件（!all(...)の部分）がTRUEの場合、つまりサンプルIDのリストが順番も含めて一致していない場合に、{ }の中の処理を実行します。
match(colnames(counts_matrix), sample_type$case_id):
ここが並べ替えの核となる部分です。
match(A, B) 関数は、ベクトルAの各要素が、ベクトルBの中で最初に出現する位置（インデックス番号）を返します。
A には colnames(counts_matrix)（基準としたい正しい順序のサンプルIDリスト）を指定します。
B には sample_type$case_id（現在のsample_typeのサンプルIDリスト）を指定します。
例:
A = colnames(counts_matrix): c("TCGA-T1", "TCGA-T2", "TCGA-N1")
B = sample_type$case_id: c("TCGA-T1", "TCGA-N1", "TCGA-T2") (現在のsample_typeの順)
match(A, B) の動作:
Aの1番目の要素 "TCGA-T1" は、B の中で1番目に見つかります。 -> 1
Aの2番目の要素 "TCGA-T2" は、B の中で3番目に見つかります。 -> 3
Aの3番目の要素 "TCGA-N1" は、B の中で2番目に見つかります。 -> 2
結果として、c(1, 3, 2) というインデックスのベクトルが返ります。このベクトルは、「現在のsample_typeの1行目、次に3行目、最後に2行目、という順番で取り出せば、counts_matrixの列順と一致する」という意味になります。
sample_type <- sample_type[match(...), ]:
Rでは、データフレームを行のインデックスで指定することで、行を並べ替えたり、特定の行だけを抽出したりできます。 データフレーム[行インデックス, 列インデックス] のように使います。列インデックスを省略（カンマの後に何も書かない）すると、全ての列が対象になります。
ここでは、sample_type の行を、ステップ7で得られたインデックスベクトル c(1, 3, 2) の順番に並べ替えています。
元のsample_type:
Apply to DEG.R
          case_id sample_type
        1 TCGA-T1       tumor
        2 TCGA-N1      normal  <-- 元の2行目
        3 TCGA-T2       tumor  <-- 元の3行目
sample_type[c(1, 3, 2), ] の結果、新しいsample_typeは以下のようになります:
Apply to DEG.R
          case_id sample_type
        1 TCGA-T1       tumor   (元の1行目)
        3 TCGA-T2       tumor   (元の3行目)
        2 TCGA-N1      normal  (元の2行目)
この結果、並べ替え後の sample_type$case_id は c("TCGA-T1", "TCGA-T2", "TCGA-N1") となり、colnames(counts_matrix) と完全に一致するようになります。
この処理によって、counts_matrixの各列が、sample_typeの対応する行のサンプル情報と正しく紐づけられることが保証され、DESeq2の解析が正しく行えるようになります。
もし、counts_matrixの列名に存在しないサンプルIDがsample_typeに含まれていたり、その逆があったりする場合には、このmatch関数はNA（欠損値）を返すことがあります。DESeqDataSetFromMatrixはそのような不整合があるとエラーになるため、事前にデータのクリーニング（IDの完全な一致を確認・修正する作業）も重要です。このコードはあくまで順序の不一致を修正するものです。

```

{{% /details %}}

- さすがに回答が長すぎるなと思ったらさらにAIに注文をつけましょう。

```
もう少し簡潔かつ初学者にわかりやすく解説して。
```

{{% details "AIの回答" %}}
```
もっと簡単なイメージで言うと…
生徒名簿A (counts_matrixの列名): 「山田くん、佐藤さん、田中くん」の順番で名前が並んでいます。これが「正しい順番」です。
生徒名簿B (sample_typeのcase_id): 「山田くん、田中くん、佐藤さん」の順番で名前が並んでいるとします。
if文のチェック: 「あれ？名簿Aと名簿Bで、2番目以降の名前の順番が違うぞ！」と気づきます。
match関数:
名簿Aの「山田くん」は、名簿Bの何番目？ -> 1番目
名簿Aの「佐藤さん」は、名簿Bの何番目？ -> 3番目
名簿Aの「田中くん」は、名簿Bの何番目？ -> 2番目
この結果、「1, 3, 2」という「取り出す順番」が分かります。
並べ替え: 名簿Bを「1番目の行、3番目の行、2番目の行」という順番に並べ替えます。
これで名簿Bも「山田くん、佐藤さん、田中くん」の順番になり、名簿Aと完全に一致します。
このようにして、2つのデータのサンプルが正しく対応するように並びを整えているのです。これができていないと、解析結果が全く意味のないものになってしまうため、とても大切なステップです。
```
{{% /details %}}

- 説明がシンプルすぎてもう少し詳しく教えてほしいならまた一言聞いてみます。

```
説明がどのコードに対応しているか具体的に教えて。
```

{{% details "AIの回答" %}}
```
まとめると:
if (!all(colnames(counts_matrix) == sample_type$case_id))
この行全体で、「counts_matrixのサンプルの並び順」と「sample_typeのサンプルの並び順」が異なるかどうかを確認しています。
具体的には、
colnames(counts_matrix) で counts_matrix のサンプル名リストを取得。
sample_type$case_id で sample_type のサンプル名リストを取得。
== で2つのリストを要素ごとに比較。
all() で全ての比較結果が TRUE かどうか（＝完全に一致しているか）を確認。
! でその結果を反転（＝一致していなければ TRUE になり、if文の中が実行される）。
sample_type <- sample_type[match(colnames(counts_matrix), sample_type$case_id), ]
この行全体で、sample_typeデータフレームの行の並び順を、counts_matrixの列の並び順（サンプル順）に揃える操作をしています。
match(A, B)関数が、「Aの各要素がBのどこにあるか」という「正しい順番の指示書（インデックス番号のリスト）」を作成します。
sample_type[指示書, ] とすることで、sample_typeの行をその指示書通りに並べ替えています。
この2行のコードによって、データ間のサンプルの対応関係が正しく保たれるようになります。
```
{{% /details %}}

## 3.疑問が解消できた場所から実行してみる。

- ただ文章を読むよりも実際にコードを動かして出力結果を確認してみる方が理解しやすいこともよくあります。
- そこで、疑問が解消できた場所から`Ctrl+Enter`(windows)や`command+Enter`(mac)で実行してみましょう。
- もしもエラーが起きたらエラーメッセージを張り付けて再度AIに聞いてみましょう。


## 4. 日常のタスクに使ってみる。

- 研究に携わっていると、RNA-seq解析以外にもいろいろなデータを分析したり図を作ったりすることがあります。
- 特に作図やそのためのデータ整理は論文を書くときに絶対に必要な作業です。
- 普段からAIを使いながらtidyverseによってデータを整理しggplot2を使って作図をすることで、最終的に論文を書くときにも時間を節約できます。
- 例えば以下のような作業を代わりにRを使ってやってみましょう。

例：
- Excelを使って表にまとめる
- Excelを使って棒グラフや散布図を作る
- qPCRの計算
- MTTなどの時系列データの計算