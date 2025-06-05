---
title: "データの前処理"
weight: 5
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---
# データの前処理

## 0. 関連チュートリアル

- [RNA-seqとは？]({{% ref "/docs/RNA-seq/what-is-RNA-seq/index.md" %}})
- [Cursorの使い方]({{% ref "/docs/get-started/Cursor/index.md" %}})
- [R basic grammar]({{% ref "/docs/R/R_basic_grammar/index.md" %}})
- [ggplot2]({{% ref "/docs/R/tidyverse/ggplot2/index.md" %}})
- [Tidyverse]({{% ref "/docs/R/tidyverse/index.md" %}})

## 1. データの準備(Raw Countsデータ & TPMデータ)

- 今回はすでに**tibble**形式にまとめたデータを用意してあるのでダウンロードしてください。
- ダウンロードしたファイルはプロジェクト内の`data`フォルダに保存してください。

Google Drive
- [TCGA_GTEx_colon_counts_tibble.csv](https://drive.google.com/uc?id=1qfSMqOc2pcrhflhq0BgTvAe_4I1t3ChG&export=download)
- [TCGA_GTEx_colon_metadata_tibble.csv](https://drive.google.com/uc?id=1499GCgwk5ep9IojJI-1XURvV7X0bW7Nw&export=download)



{{% hint info %}}

上記のデータは以下のサイトからダウンロードしたものです。

- [GDC Data Portal](https://portal.gdc.cancer.gov/): 大腸癌のRNA-seqデータ。
- [GTEx Portal](https://www.gtexportal.org/home/): 正常大腸のRNA-seqデータ。

※ public データのダウンロードとクリーニングが実は最難関かもしれませんので、最後に紹介します。

{{% /hint %}}




## 2. データの読み込みとlog変換


{{< highlight sh >}}
library(tidyverse)
library(here)
{{< /highlight >}}

- まずはデータを読み込んでどのようなデータなのか調べてみます。
- RNA-seqのデータはとても大きいので、全体像をただ見て把握するのは困難です。
- そこで、データを要約してくれる関数や図を使って全体像を把握します。

{{< code-tabs >}}
--- code
{{< highlight R >}}
# データの読み込み
counts <- read_csv(here("data", "TCGA_GTEx_colon_counts_tibble.csv"))
metadata <- read_csv(here("data", "TCGA_GTEx_colon_metadata_tibble.csv"))

# データの中身を確認
head(counts)
head(metadata)

# table()で各列にどのようなデータが入っているか確認
table(counts$gene_type)
table(metadata$source)
{{< /highlight >}}

--- output
{{< highlight R>}}
> counts <- read_csv(here("data", "TCGA_GTEx_colon_counts_tibble.csv"))
Rows: 61569 Columns: 103
── Column specification ────────────────────────────────────────────────────────
Delimiter: ","
chr   (3): gene_id, gene_name, gene_type
dbl (100): TCGA-D5-5540, TCGA-EI-6509, TCGA-A6-6137, TCGA-QG-A5Z2, TCGA-AA-3...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.


> metadata <- read_csv(here("data", "TCGA_GTEx_colon_metadata_tibble.csv"))
Rows: 100 Columns: 210
── Column specification ────────────────────────────────────────────────────────
Delimiter: ","
chr (210): case_id, source, project.project_id, cases.consent_type, cases.da...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.


> head(counts)
# A tibble: 6 × 103
  gene_id       gene_name gene_type `TCGA-D5-5540` `TCGA-EI-6509` `TCGA-A6-6137`
  <chr>         <chr>     <chr>              <dbl>          <dbl>          <dbl>
1 ENSG00000000… TSPAN6    protein_…          14293           4986          13086
2 ENSG00000000… TNMD      protein_…             27              2            106
3 ENSG00000000… DPM1      protein_…           6729           2769           2197
4 ENSG00000000… SCYL3     protein_…            707            565            770
5 ENSG00000000… C1orf112  protein_…           1388            543            334
6 ENSG00000000… FGR       protein_…             78            176            250
# ℹ 97 more variables: `TCGA-QG-A5Z2` <dbl>, `TCGA-AA-3489` <dbl>,
#   `HCM-CSHL-0160-C18` <dbl>, `TCGA-AD-A5EK` <dbl>, `TCGA-AA-3867` <dbl>,
#   `TCGA-AA-3975` <dbl>, `05CO007` <dbl>, `TCGA-CK-5914` <dbl>,
#   `15CO002` <dbl>, `TCGA-AZ-6598` <dbl>, `TCGA-AZ-6607` <dbl>,
#   `05CO037` <dbl>, `TCGA-D5-6539` <dbl>, `TCGA-AA-3662` <dbl>,
#   `TCGA-CK-5913` <dbl>, `TCGA-AY-6197` <dbl>, `HCM-STAN-1111-C19` <dbl>,
#   `05CO048` <dbl>, `TCGA-G4-6293` <dbl>, `HCM-CSHL-0060-C18` <dbl>, …


> head(metadata)
# A tibble: 6 × 210
  case_id     source project.project_id cases.consent_type cases.days_to_consent
  <chr>       <chr>  <chr>              <chr>              <chr>                
1 TCGA-D5-55… TCGA_… TCGA-COAD          Informed Consent   20                   
2 TCGA-EI-65… TCGA_… TCGA-READ          Informed Consent   0                    
3 TCGA-A6-61… TCGA_… TCGA-COAD          Informed Consent   0                    
4 TCGA-QG-A5… TCGA_… TCGA-COAD          Informed Consent   49                   
5 TCGA-AA-34… TCGA_… TCGA-COAD          Informed Consent   31                   
6 HCM-CSHL-0… TCGA_… HCMI-CMDC          '--                '--                  
# ℹ 205 more variables: cases.days_to_lost_to_followup <chr>,
#   cases.disease_type <chr>, cases.index_date <chr>,
#   cases.lost_to_followup <chr>, cases.primary_site <chr>,
#   demographic.age_at_index <chr>, demographic.age_is_obfuscated <chr>,
#   demographic.cause_of_death <chr>, demographic.cause_of_death_source <chr>,
#   demographic.country_of_birth <chr>,
#   demographic.country_of_residence_at_enrollment <chr>, …


> table(counts$gene_type)
                         IG_C_gene                    IG_C_pseudogene 
                                14                                  9 
                         IG_D_gene                          IG_J_gene 
                                37                                 18 
                   IG_J_pseudogene                      IG_pseudogene 
                                 3                                  1 
                         IG_V_gene                    IG_V_pseudogene 
                               145                                187 
                            lncRNA                              miRNA 
                             16901                               1881 
                          misc_RNA                            Mt_rRNA 
                              2212                                  2 
                           Mt_tRNA             polymorphic_pseudogene 
                                22                                 48 
              processed_pseudogene                     protein_coding 
                             10167                              19962 
                        pseudogene                           ribozyme 
                                18                                  8 
                              rRNA                    rRNA_pseudogene 
                                47                                497 
                            scaRNA                              scRNA 
                                49                                  1 
                            snoRNA                              snRNA 
                               943                               1901 
                              sRNA                                TEC 
                                 5                               1057 
                         TR_C_gene                          TR_D_gene 
                                 6                                  4 
                         TR_J_gene                    TR_J_pseudogene 
                                79                                  4 
                         TR_V_gene                    TR_V_pseudogene 
                               106                                 33 
  transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
                               500                                138 
transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
                               939                                  2 
 translated_unprocessed_pseudogene                 unitary_pseudogene 
                                 1                                 98 
            unprocessed_pseudogene                          vault_RNA 
                              2614                                  1 


> table(metadata$source)
   GTEx_sigmoid GTEx_transverse        TCGA_CRC 
             25              25              50 
{{</ highlight >}}
{{</ code-tabs >}}

{{% hint info %}}
### here()関数

- here()関数はファイルのパスを指定するための関数です。
- この関数を使うと、相対パスを絶対パスに変換してくれます。
- 絶対パスを自分で書くのはめんどくさいので相対パスをよく使いますが、相対パスだとうまく動かなくなる場合があります。
- パスを手書きするならそれぞれ以下のようになります。
  - 相対パス：`data/TCGA_GTEx_colon_counts_tibble.csv`
  - 絶対パス：`/home/user/proj/{プロジェクト名}/data/TCGA_GTEx_colon_counts_tibble.csv"`
- `here("data", "TCGA_GTEx_colon_counts_tibble.csv")`とだけ書いて実行し、本当に絶対パスが出力されるか確認してみましょう。

{{% /hint %}}

- gene_typeにはprotein_coding(mRNA)以外にも様々なtypeが含まれています。
- sourceはそのサンプルがどのデータベースから取得したものであるかが書かれています。
- GTExからはS状結腸と横行結腸からのサンプルが25サンプル、TCGAからは大腸癌のサンプルが50サンプル取得してあります。
- 今回はmRNAの発現量を扱うのでまずgene_typeがprotein_codingのデータを抽出します。

{{< code-tabs >}}
--- code
{{< highlight R >}}
# gene_typeがprotein_codingのデータを抽出
counts_mrna <- counts %>% 
    filter(gene_type == "protein_coding")

# データの中身を確認
table(counts_mrna$gene_type)
{{< /highlight >}}

--- output
{{< highlight R >}}
protein_coding 
         19962 
{{< /highlight >}}
{{</ code-tabs>}}

- gene_typeがprotein_codingのみになりました。
- 次に行名をgene_nameにして読みやすくします。
- さらにmatrixに変換して今後の解析に使いやすくします。

{{< code-tabs>}}
--- code
{{< highlight R >}}
# 行名をgene_nameに変更
counts_mrna_matrix <- counts_mrna %>% 
    column_to_rownames("gene_name") %>% 
    select(!c(gene_id, gene_type)) %>% 
    as.matrix()
{{< /highlight >}}

--- output
{{< highlight R >}}
Warning: non-unique values when setting 'row.names': 'ACTL10', 'AKAP17A',
'ASMT', 'ASMTL', 'CD99', 'CRLF2', 'CSF2RA', 'DHRSX', 'GTPBP6', 'IL3RA', 'IL9R',
'MATR3', 'P2RY8', 'PDE11A', 'PLCXD1', 'POLR2J3', 'PPP2R3B', 'SHOX', 'SLC25A6',
'SMIM40', 'TMSB15B', 'VAMP7', 'WASH6P', 'ZBED1'

Error in `.rowNamesDF<-`(x, value = value): 重複した 'row.names' は許されません
{{< /highlight >}}
{{</ code-tabs>}}

- ひとつのgene_nameに対して複数のgene_idがあるせいで変換できませんでした。
- 複数のアイソフォームがあるなどの理由でgene_idが複数ある事があります。
- この場合はカウントの合計をとってひとつのgene_nameにまとめます。


{{% details title="同じgene_nameをもつ行のカウント値を合計して統合し、matrixに変換" %}}
{{< highlight R >}}
# 同じ遺伝子名の値を合計して統合し、matrixに変換
counts_mrna_matrix <- counts_mrna %>% 
    group_by(gene_name) %>%
    summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
    column_to_rownames("gene_name") %>% 
    as.matrix()

{{< /highlight >}}
{{% /details %}}

{{% details title="ヒストグラムを使って発現値の分布を確認" %}}
{{< highlight R >}}
# ヒストグラムを使って発現値の分布を確認
hist(
    counts_mrna_matrix,
    main = "Distribution of expression values",
    xlab = "Expression values",
    ylab = "Frequency"
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-7-1.png" >}}

- カウントの値の頻度を見ると、ほとんどの値が非常に低い値だがごく一部非常に大きな値もあることがわかります。
- 一般にRNA-seqの発現値はこのような対数正規分布になります。
- このように値のスケールがあまりにも違いすぎるデータは理解しにくく扱いづらいので、ログ変換を行います。

{{% details title="ログ変換を行って発現値の分布を確認" %}}
{{< highlight R >}}
# ログ変換
counts_log_matrix <- log2(counts_mrna_matrix + 1)

# ログ変換後のヒストグラム
hist(
    counts_log_matrix,
    main = "Distribution of log-transformed expression values",
    xlab = "log-transformed expression values",
    ylab = "Frequency"
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-8-1.png" >}}
- このようにログ変換によって、桁違いだった値が小さくなり、ほとんどの値を0から15くらいで表せるようになりました。
- 値が0に近い遺伝子を除けば10付近を頂点とした山になっていることが分かります。
- 左右対称な山の分布は正規分布と呼び、log変換によって正規分布に近づいたといえます。


## 3. データの正規化

- サンプル間で比較ができるように**メディアン比正規化**を行います。
- [RNA-seqとは？]({{% ref "/docs/RNA-seq/what-is-RNA-seq/index.md#メディアン比正規化median-ratio-normalization" %}})で詳しく解説しています。
- 正規化はraw countsデータに対して行います。


{{< highlight R >}}
library(DESeq2)
{{< /highlight >}}

- DEGを行う代表的なパッケージであるDESeq2にはメディアン比正規化の関数が用意されています。
- 発現値が0の遺伝子が含まれるとメディアン比正規化はできませんが、これは関数内部で自動的に除外されます。
- 0でなかったとしてもあまりに低すぎると正確な結果が出ないことがあるので、平均発現値が10カウント未満の遺伝子を除外することが多いです。


{{< details title="平均発現値が10以上の遺伝子のみを保持" >}}

{{< code-tabs >}}
--- code
{{< highlight R >}}
# 平均発現値が10以上の遺伝子のみを保持
keep_genes <- rowMeans(counts_mrna_matrix) >= 10
filtered_counts_mrna_matrix <- counts_mrna_matrix[keep_genes, ]

print(paste0("除外された遺伝子数: ", nrow(counts_mrna_matrix) - nrow(filtered_counts_mrna_matrix)))
print(paste0("除外された遺伝子の割合: ", round((nrow(counts_mrna_matrix) - nrow(filtered_counts_mrna_matrix)) / nrow(counts_mrna_matrix) * 100, 2), "%"))
{{< /highlight >}}

--- output
{{< highlight R >}}
> print(paste0("除外された遺伝子数: ", nrow(counts_mrna_matrix) - nrow(filtered_counts_mrna_matrix)))
[1] "除外された遺伝子数: 3481"

> print(paste0("除外された遺伝子数: ", nrow(counts_mrna_matrix) - nrow(filtered_counts_mrna_matrix)))
[1] "除外された遺伝子の割合: 17.46%"
{{< /highlight >}}
{{</ code-tabs >}}

{{< /details >}}

{{% details title="DEseq2の関数を使ってsize factorを計算" %}}
{{< highlight R >}}
# size factorを計算
size_factors <- estimateSizeFactorsForMatrix(filtered_counts_mrna_matrix)

# 結果を確認
head(size_factors)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
     TCGA-D5-5540      TCGA-EI-6509      TCGA-A6-6137      TCGA-QG-A5Z2 
        0.6611094         1.0191326         1.0514363         1.4867364 
     TCGA-AA-3489 HCM-CSHL-0160-C18 
        1.3669717         0.6783093 
{{< /highlight >}}

{{% details title="サンプル毎に計算したsize factorsで割って正規化" %}}
{{< highlight R >}}
# サンプル毎に計算したsize factorsで割って正規化
normalized_counts <- sweep(filtered_counts_mrna_matrix, 2, size_factors, "/")

# データの中身を確認
print(normalized_counts[1:5, 1:5])
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
       TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
A1BG       1.512609     3.924906      0.00000     3.363071     3.657720
A1CF    1188.910555  2696.410738   1891.69802   434.508762   480.624449
A2M     3542.529924  6484.926700  13386.45027  2968.246389 28157.130975
A2ML1      1.512609     8.831040     16.16836     4.708299     5.852353
A4GALT    72.605225   322.823556    258.69375   122.415781   771.778986
{{< /highlight >}}


{{% details title="正規化したデータをlog変換してヒストグラムを作成" %}}
{{< highlight R >}}
# log変換
log_normalized_counts <- log2(normalized_counts + 1)

# 正規化かつlog変換後のヒストグラム
hist(
    log_normalized_counts,
    main = "Distribution of normalized expression values",
    xlab = "log-normalized expression values",
    ylab = "Frequency"
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-10-1.png" >}}

- 最終的に非常に正規分布に近い分布になりました。
- 正規分布は様々な統計解析手法で仮定される分布のため非常に扱いやすくなっています。

