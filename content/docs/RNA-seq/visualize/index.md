---
title: "可視化"
weight: 5
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---
# 可視化

## 4. heatmapで全体像をつかむ

- 無事に正規化できたのでheatmapを作成してみましょう。
- 今回は少しデータを追加できるように便利なheatmapパッケージを使用します。


{{< highlight R >}}
library(ComplexHeatmap)
{{< /highlight >}}


- すべての遺伝子を使うと計算に時間がかかります。
- 今回は大腸癌の分子サブタイプ分類に用いられる遺伝子セットを使います
- 本来は一つの遺伝子セットを用いてクラスタリングしたあとに、さらに別の遺伝子セットを用いて・・・と段階的に行う方がよいですが、今回は例としてまとめて行います。

{{< code-tabs >}}
--- code
{{< highlight R >}}
# 大腸癌の分子サブタイプ分類に用いられる遺伝子セット

## 大腸癌サブクラスタリング用の主要遺伝子セット
crc_driver_genes <- c("APC", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS", "BRAF", "TP53", "CTNNB1", "PTEN", "MSH2", "MLH1", "MSH6", "PMS2", "POLE", "POLD1")

## CMS分類関連遺伝子
cms_signature_genes <- c("CDX2", "EPHB2", "HNF4A", "KRT20", "VIL1", "CEACAM5", "EPCAM", "LGR5", "ASCL2", "OLFM4", "SOX9", "BMI1", "MSI1", "DCLK1", "CD44")

## 免疫関連遺伝子セット
immune_genes <- c("CD3E", "CD8A", "GZMA", "GZMB", "PRF1", "IFNG", "TBX21", "FOXP3", "IL10", "TGFB1", "PD1", "PDL1", "CTLA4", "LAG3")

## 間質関連遺伝子
stromal_genes <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2", "PDGFRB", "FAP", "THY1", "SPARC", "TIMP1", "VCAN")

# 遺伝子セットを結合
genes_set <- c(crc_driver_genes, cms_signature_genes, immune_genes, stromal_genes)

# log_normalized_countsに実際に存在する遺伝子のみをフィルタリング
available_genes <- intersect(genes_set, rownames(log_normalized_counts))

# 除外された遺伝子を確認
missing_genes <- setdiff(genes_set, rownames(log_normalized_counts))

print(paste0("指定した遺伝子数: ", length(genes_set)))
print(paste0("利用可能な遺伝子数: ", length(available_genes)))
print(paste0("除外された遺伝子数: ", length(missing_genes)))
print(missing_genes)
{{< /highlight >}}

--- output
{{< highlight R >}}
> print(paste0("指定した遺伝子数: ", length(genes_set)))
[1] "指定した遺伝子数: 58"

> print(paste0("利用可能な遺伝子数: ", length(available_genes)))
[1] "利用可能な遺伝子数: 56"

> print(paste0("除外された遺伝子数: ", length(missing_genes)))
[1] "除外された遺伝子数: 2"

> print(missing_genes)
[1] "PD1"  "PDL1"
{{< /highlight >}}
{{</ code-tabs >}}

- これまでのフィルタリング（10カウント未満の遺伝子を除外）により、指定した遺伝子セットの遺伝子がなくなってしまったことがわかります。
- もしもPD1やPDL1の研究をしているのであれば、これらの遺伝子を除外しないようにフィルタリング方法を変えたり特定の遺伝子だけ意図的に残したりする必要があります。


{{% details title="利用可能な遺伝子のみでheatmap用のデータを作成" %}}
{{< highlight R >}}
# 利用可能な遺伝子のみでheatmap用のデータを作成
heatmap_data <- log_normalized_counts[available_genes, ]
{{< /highlight >}}
{{% /details %}}

- 行ごとにz-scoreスケーリングを行います。
- データの平均を0、標準偏差を1にすることで、遺伝子ごとの発現パターンの相対的な変化を比較しやすくします。

{{% details title="行ごとにz-scoreスケーリングを行う" %}}
{{< highlight R >}}
# 行ごとにz-scoreスケーリング
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# スケーリングされたデータの確認（最初の数行、数列）
print(heatmap_data_scaled[1:5, 1:5])
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
       TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
APC     -0.94538403   -1.4797431 -1.104450118  -0.51245302  -0.30314542
KRAS     1.26869736    0.5109699  0.345097246   0.11027905  -0.28421717
PIK3CA  -1.68336116   -0.2764212 -0.004262841  -0.77866149  -0.04524831
FBXW7   -1.53689130   -2.5857930 -0.610560350   0.06755438   0.30084770
SMAD4   -0.01473564   -3.3935131 -0.672149570   0.35610149  -0.81051601
{{< /highlight >}}

{{% details title="ヒートマップの作成(列方向のクラスタリングのみ行う、列名は表示せず行名は表示する。)" %}}
{{< highlight R >}}
# ヒートマップの作成
ComplexHeatmap::Heatmap(
    heatmap_data_scaled,
    name = "exp",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),  
    row_names_max_width = unit(6, "cm"),
    row_names_side = "right",
    width = unit(8, "cm"),
    height = unit(14, "cm")
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-14-1.png" >}}

- これだとどれが腫瘍でどれが正常大腸なのかやどの遺伝子がどの遺伝子セットに属しているかが分かりません。
- **annotation bar**を追加することで情報を追加できます。

{{% details title="col annotation barにするための情報をmetadataから抜き出す(stageはdiagnoses.ajcc_pathologic_stage)" %}}
{{< highlight R >}}
# metadataからcol annotation barにしたい情報を抜き出す
col_annotation_data <- metadata %>% 
    select(case_id, source, diagnoses.ajcc_pathologic_stage) %>% 
    setNames(c("case_id", "source", "stage")) %>% 
    mutate(stage = if_else(is.na(stage), "Normal", stage)) %>% 
    column_to_rownames("case_id")

# データを確認
head(col_annotation_data)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
                    source      stage
TCGA-D5-5540      TCGA_CRC  Stage IIA
TCGA-EI-6509      TCGA_CRC        '--
TCGA-A6-6137      TCGA_CRC Stage IIIB
TCGA-QG-A5Z2      TCGA_CRC    Stage I
TCGA-AA-3489      TCGA_CRC   Stage II
HCM-CSHL-0160-C18 TCGA_CRC  Stage IIA
{{< /highlight >}}


{{< code-tabs >}}
--- code
{{< highlight R >}}
table(col_annotation_data$source)
table(col_annotation_data$stage)
{{< /highlight >}}

--- output
{{< highlight R >}}
> table(col_annotation_data$source)
   GTEx_sigmoid GTEx_transverse        TCGA_CRC 
             25              25              50 

> table(col_annotation_data$stage)
       '--     Normal    Stage I   Stage II  Stage IIA  Stage IIB  Stage III 
         9         50          9          4         13          1          2 
Stage IIIA Stage IIIB Stage IIIC   Stage IV  Stage IVA  Stage IVB 
         1          5          3          1          1          1 
{{< /highlight >}}
{{</ code-tabs>}}

{{% details title="RColorBrewerの色パレットを使ってcol annotation barを作成" %}}
{{< highlight R >}}
library(RColorBrewer)

# sourceの色設定
source_colors <- brewer.pal(
    n = length(unique(col_annotation_data$source)), 
    "Set1"
)
names(source_colors) <- unique(col_annotation_data$source)

# stageの色設定
stage_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(col_annotation_data$stage)))
names(stage_colors) <- unique(col_annotation_data$stage)

# col annotation barを作成
col_annotation <- HeatmapAnnotation(
    df = col_annotation_data,
    col = list(
        source = source_colors,
        stage = stage_colors
    )
)

# heatmapを作成
Heatmap(
    heatmap_data_scaled,
    name = "exp",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),  
    row_names_max_width = unit(6, "cm"), 
    row_names_side = "right",         
    width = unit(8, "cm"),
    height = unit(14, "cm"),
    top_annotation = col_annotation
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-15-1.png" >}}

- これでサンプルの情報を追加することができました。
- 続いて遺伝子がどの遺伝子セットに属しているかを追加します。


{{% details title="row annotation bar用のデータをgene setから作成して色を設定しannotation barを作成" %}}
{{< highlight R >}}
# gene setからrow annotation barを作成
row_annotation_data <- data.frame(
    set_name = c(
        rep("CRC driver", length(crc_driver_genes)),
        rep("CMS signature", length(cms_signature_genes)),
        rep("Immune", length(immune_genes)),
        rep("Stromal", length(stromal_genes))
    )
)
rownames(row_annotation_data) <- genes_set

# available_genesに一致する遺伝子のみを抽出
row_annotation_data <- row_annotation_data[available_genes, , drop = FALSE]

# gene setの色を設定
gene_set_colors <- brewer.pal(
    n = length(unique(row_annotation_data$set_name)), 
    "Set2"
)
names(gene_set_colors) <- unique(row_annotation_data$set_name)

# row annotation barを作成
row_annotation <- rowAnnotation(
    df = row_annotation_data,
    col = list(
        set_name = gene_set_colors
    )
)

# heatmapを作成
Heatmap(
    heatmap_data_scaled,
    name = "exp",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),  
    row_names_max_width = unit(6, "cm"), 
    row_names_side = "right",         
    width = unit(8, "cm"),
    height = unit(14, "cm"),
    top_annotation = col_annotation,
    left_annotation = row_annotation
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-16-1.png" >}}


- Stageの項目が多すぎて少々見づらくなっています。
- Stage I, II, III, IVというようにまとめてみましょう。


{{< code-tabs >}}
--- code
{{< highlight R >}}
# 現在のStageの項目を確認
table(col_annotation_data$stage)
{{< /highlight >}}

--- output
{{< highlight R >}}
       '--     Normal    Stage I   Stage II  Stage IIA  Stage IIB  Stage III 
         9         50          9          4         13          1          2 
Stage IIIA Stage IIIB Stage IIIC   Stage IV  Stage IVA  Stage IVB 
         1          5          3          1          1          1 
{{< /highlight >}}
{{</ code-tabs >}}

{{% details title="ステージ情報を整理" %}}
{{< highlight R >}}
# ステージ情報を整理
col_annotation_data_updated <- col_annotation_data %>% 
    mutate(stage_simplified = case_when(
        stage == "'--" ~ "Unknown",
        stage == "Normal" ~ "Normal",
        str_detect(stage, "^Stage IV") ~ "Stage_IV",
        str_detect(stage, "^Stage III") ~ "Stage_III",
        str_detect(stage, "^Stage II") ~ "Stage_II",
        str_detect(stage, "^Stage I") ~ "Stage_I",
        TRUE ~ "Other"
    ))

# 整理後のステージ分布を確認
table(col_annotation_data_updated$stage_simplified)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
   Normal   Stage_I  Stage_II Stage_III  Stage_IV   Unknown 
       50         9        18        11         3         9 
{{< /highlight >}}

{{< code-tabs >}}
--- code
{{< highlight R >}}
table(col_annotation_data$stage, col_annotation_data_updated$stage_simplified)
{{< /highlight >}}

--- output
{{< highlight R >}}
            
             Normal Stage_I Stage_II Stage_III Stage_IV Unknown
  '--             0       0        0         0        0       9
  Normal         50       0        0         0        0       0
  Stage I         0       9        0         0        0       0
  Stage II        0       0        4         0        0       0
  Stage IIA       0       0       13         0        0       0
  Stage IIB       0       0        1         0        0       0
  Stage III       0       0        0         2        0       0
  Stage IIIA      0       0        0         1        0       0
  Stage IIIB      0       0        0         5        0       0
  Stage IIIC      0       0        0         3        0       0
  Stage IV        0       0        0         0        1       0
  Stage IVA       0       0        0         0        1       0
  Stage IVB       0       0        0         0        1       0
{{< /highlight >}}
{{</ code-tabs >}}


{{% details title="ステージ情報の整理後、色設定を更新" %}}
{{< highlight R >}}
# ステージ情報の整理後、色設定を更新
stage_colors_updated <- brewer.pal(
    n = min(length(unique(col_annotation_data_updated$stage_simplified)), 11), 
    "Dark2"
)
names(stage_colors_updated) <- unique(col_annotation_data_updated$stage_simplified)

# stage列は使わないのでcol_annotation_data_updatedから削除
col_annotation_data_updated <- col_annotation_data_updated %>% 
    select(!stage)

# 色設定を確認
print(stage_colors_updated)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
 Stage_II   Unknown Stage_III   Stage_I  Stage_IV    Normal 
"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" 
{{< /highlight >}}



{{% details title="更新されたcol annotation barを作成しheatmapを作成" %}}
{{< highlight R >}}
# 更新されたcol annotation barを作成
col_annotation_updated <- HeatmapAnnotation(
    df = col_annotation_data_updated,
    col = list(
        source = source_colors,
        stage_simplified = stage_colors_updated
    )
)

# 整理されたステージ情報を使ったheatmapを作成
Heatmap(
    heatmap_data_scaled,
    name = "exp",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),  
    row_names_max_width = unit(6, "cm"), 
    row_names_side = "right",         
    width = unit(8, "cm"),
    height = unit(14, "cm"),
    top_annotation = col_annotation_updated,
    left_annotation = row_annotation
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-17-1.png" >}}



- このようにheatmapは臨床情報と遺伝子発現量の両方を可視化することができます。
- このような可視化はデータの全体像を把握するのに非常に便利です。






## 5. PCAでデータのばらつきをつかむ

- 次にデータのばらつきをつかむためにPCAを行います。
- Principal Component Analysis (PCA)はデータの次元を削減するための手法です。
- RNA-seqデータは1サンプルにつき20000遺伝子＝20000次元のデータとなります。
- このような高次元データは可視化が困難です。
- PCAはデータのばらつきを主成分軸に分解し、それぞれの軸に対してどの程度の情報を持っているかを可視化します。
- 20000次元のデータを2次元に圧縮して人間が理解できるような形にするということです。


{{% details title="PCAを計算してリザルトの構造を表示" %}}
{{< highlight R >}}
# PCAを行う（行をサンプルにする必要があるので転置する）
pca_result <- prcomp(t(log_normalized_counts))

# resultの構造を確認
str(pca_result)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
List of 5
 $ sdev    : num [1:100] 136 60 44.5 35.8 28.1 ...
 $ rotation: num [1:16457, 1:100] 0.01699 -0.01958 0.00588 0.00362 0.00951 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:16457] "A1BG" "A1CF" "A2M" "A2ML1" ...
  .. ..$ : chr [1:100] "PC1" "PC2" "PC3" "PC4" ...
 $ center  : Named num [1:16457] 4.14 7.74 14.06 3.48 9.45 ...
  ..- attr(*, "names")= chr [1:16457] "A1BG" "A1CF" "A2M" "A2ML1" ...
 $ scale   : logi FALSE
 $ x       : num [1:100, 1:100] -156.4 -120.7 -141.7 -149 -77.9 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:100] "TCGA-D5-5540" "TCGA-EI-6509" "TCGA-A6-6137" "TCGA-QG-A5Z2" ...
  .. ..$ : chr [1:100] "PC1" "PC2" "PC3" "PC4" ...
 - attr(*, "class")= chr "prcomp"
{{< /highlight >}}

- `pca_result$sdev`には各主成分の標準偏差が格納されています。これを2乗すると分散になります。
- 各主成分の分散を、すべての主成分の分散の合計で割ることで、各主成分がデータ全体のばらつきのうちどれくらいの割合を説明しているか（寄与率）がわかります。
- 計算して棒グラフで表示してみましょう。

{{% details title="主成分の寄与率を棒グラフで可視化(barplot関数)" %}}
{{< highlight R >}}
# 主成分の寄与率を棒グラフで可視化
barplot(
    pca_result$sdev^2 / sum(pca_result$sdev^2),
    main = "PCA contribution rates",
    xlab = "Principal Components",
    ylab = "Contribution Rate",
    names.arg = paste0("PC", 1:length(pca_result$sdev))
)
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-18-1.png" >}}


- PC1（最初の主成分）が最も多くの情報を持っており、PCが進むにつれて寄与率が小さくなっていることがわかります。
- 次に、最も情報量の多いPC1とPC2を使って、サンプルがどのように分布しているかを2次元の散布図で見てみましょう。
- `pca_result$x` には各サンプルが各主成分に対して持つスコアが格納されています。これを利用します。


{{% details title="PCAの結果とメタデータ(source列)を結合したデータフレームを作成(サンプルはsample_id列にする)" %}}
{{< highlight R >}}
# PCAの結果とメタデータを結合したデータフレームを作成
plot_data_pca_gg <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  source = metadata$source,
  sample_id = rownames(pca_result$x)
)

# データの中身を確認
head(plot_data_pca_gg)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
                         PC1        PC2   source         sample_id
TCGA-D5-5540      -156.38728 -82.446074 TCGA_CRC      TCGA-D5-5540
TCGA-EI-6509      -120.68225 -49.843552 TCGA_CRC      TCGA-EI-6509
TCGA-A6-6137      -141.74386   9.706455 TCGA_CRC      TCGA-A6-6137
TCGA-QG-A5Z2      -149.04311   7.815018 TCGA_CRC      TCGA-QG-A5Z2
TCGA-AA-3489       -77.93359 -33.290015 TCGA_CRC      TCGA-AA-3489
HCM-CSHL-0160-C18  -56.86226   1.675592 TCGA_CRC HCM-CSHL-0160-C18
{{< /highlight >}}

{{% details title="ggplot2で散布図をプロット" %}}
{{< highlight R >}}
# ggplot2でプロット
ggplot(plot_data_pca_gg, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    y = paste0("PC2 (", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    color = "Sample Source",
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-19-1.png" >}}


- このプロットでは、横軸にPC1、縦軸にPC2を取り、各サンプルを点で表示しています。
- 次は色分けをしてもう少しわかりやすくしてみましょう。


{{% details title="ggplot2で散布図をプロット(色分け)" %}}
{{< highlight R >}}
# ggplot2でプロット(色分け)
ggplot(plot_data_pca_gg, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = source), size = 3, alpha = 0.7) +
  scale_color_manual(values = source_colors) + 
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    y = paste0("PC2 (", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    color = "Sample Source",
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-20-1.png" >}}

- このように、大腸がんと正常大腸粘膜はかなり分離していることが分かります。
- 一方で横行結腸とS状結腸はかなり近いところにいて一部混ざっていますが、ある程度別の集団であることを示しています。
- このようなグループ分けを自動的に行う方法（クラスタリング）も試してみましょう。
- `k-means`という手法によってクラスタリングを行ってみます。


{{% details title="k-meansクラスタリングを行い、クラスタリング結果を信頼楕円としてプロット" %}}
{{< highlight R >}}
# k-meansクラスタリング (PC1とPC2を使用)
set.seed(123)
kmeans_obj <- kmeans(plot_data_pca_gg[, c("PC1", "PC2")], centers = 3, nstart = 25)
plot_data_pca_gg$cluster <- factor(kmeans_obj$cluster)

# クラスタの中心点をデータフレームとして準備
cluster_centers_gg <- as.data.frame(kmeans_obj$centers)
colnames(cluster_centers_gg) <- c("PC1", "PC2")

# クラスタを中心点から信頼楕円を書くためにクラスタIDを追加
cluster_centers_gg$cluster <- factor(1:nrow(cluster_centers_gg))

# ggplot2でプロット
ggplot(plot_data_pca_gg, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = source, shape = cluster), size = 3, alpha = 0.7) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.15, type = "t", level = 0.95, show.legend = FALSE, color = "darkgrey") + 
  scale_color_manual(values = source_colors) + 
  scale_fill_brewer(palette = "Pastel2") + 
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15)) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    y = paste0("PC2 (", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    color = "Sample Source",
    shape = "K-means Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-21-1.png" >}}

- 大腸がんと正常大腸粘膜ははっきりクラスタリングすることができました。
- 横行結腸とS状結腸の正常粘膜については完全に分けることはできませんが、集団の大部分が所属するクラスターに分けることはできました。
- 今回はすべての遺伝子を使いましたが、多くの場合ほとんどの遺伝子は大きく変動していないのでがんと正常大腸を分けることくらいしかできませんでした。
- 大腸がんをさらに細かく分けるためには、ある程度重要な遺伝子を絞り込んでからPCAを行う必要があります。
- heatmapの時に使ったgene setの遺伝子を使ってみましょう。

{{< code-tabs >}}
--- code
{{< highlight R >}}
# heatmapの時に使ったgene setのデータを確認
print(row_annotation_data)
{{< /highlight >}}

--- output
{{< highlight R >}}
             set_name
APC        CRC driver
KRAS       CRC driver
PIK3CA     CRC driver
FBXW7      CRC driver
SMAD4      CRC driver
TCF7L2     CRC driver
NRAS       CRC driver
BRAF       CRC driver
TP53       CRC driver
CTNNB1     CRC driver
PTEN       CRC driver
MSH2       CRC driver
MLH1       CRC driver
MSH6       CRC driver
PMS2       CRC driver
POLE       CRC driver
POLD1      CRC driver
CDX2    CMS signature
EPHB2   CMS signature
HNF4A   CMS signature
KRT20   CMS signature
VIL1    CMS signature
CEACAM5 CMS signature
EPCAM   CMS signature
LGR5    CMS signature
ASCL2   CMS signature
OLFM4   CMS signature
SOX9    CMS signature
BMI1    CMS signature
MSI1    CMS signature
DCLK1   CMS signature
CD44    CMS signature
CD3E           Immune
CD8A           Immune
GZMA           Immune
GZMB           Immune
PRF1           Immune
IFNG           Immune
TBX21          Immune
FOXP3          Immune
IL10           Immune
TGFB1          Immune
CTLA4          Immune
LAG3           Immune
COL1A1        Stromal
COL1A2        Stromal
COL3A1        Stromal
FN1           Stromal
VIM           Stromal
ACTA2         Stromal
PDGFRB        Stromal
FAP           Stromal
THY1          Stromal
SPARC         Stromal
TIMP1         Stromal
VCAN          Stromal
{{< /highlight >}}
{{</ code-tabs >}}


{{% details title="tibbleに変換" %}}
{{< highlight R >}}
# tibbleに変換
gene_set_tibble <- row_annotation_data |>
    rownames_to_column(var = "gene_name") |>
    as_tibble()

# データの中身を確認
print(gene_set_tibble)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
# A tibble: 56 × 2
   gene_name set_name  
   <chr>     <chr>     
 1 APC       CRC driver
 2 KRAS      CRC driver
 3 PIK3CA    CRC driver
 4 FBXW7     CRC driver
 5 SMAD4     CRC driver
 6 TCF7L2    CRC driver
 7 NRAS      CRC driver
 8 BRAF      CRC driver
 9 TP53      CRC driver
10 CTNNB1    CRC driver
# ℹ 46 more rows
{{< /highlight >}}



- これらの遺伝子だけを使って、大腸がんだけをPCAプロットしてみましょう。


{{% details title="大腸がんのサンプルの行だけを抽出" %}}
{{< highlight R >}}
# metadataから大腸がんのサンプルの行だけを抽出
tumor_samples_metadata <- metadata |>
    filter(source == "TCGA_CRC")

head(tumor_samples_metadata)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
```
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
{{< /highlight >}}



{{% details title="log_normalized_countsから大腸がんのサンプルの列だけを抽出" %}}
{{< highlight R >}}
# log_normalized_countsから大腸がんのサンプルの列だけを抽出
tumor_log_normalized_counts <- log_normalized_counts[, tumor_samples_metadata$case_id]

tumor_log_normalized_counts[1:5, 1:5]
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
       TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
A1BG       1.329186     2.300096     0.000000     2.125344     2.219624
A1CF      10.216637    11.397360    10.886229     8.766558     8.911765
A2M       11.790972    12.663097    13.708594    11.535881    14.781264
A2ML1      1.329186     3.297344     4.101680     2.513061     2.776599
A4GALT     6.201736     8.339064     8.020667     6.947383     9.593912
{{< /highlight >}}

{{< highlight R >}}
dim(tumor_log_normalized_counts)
{{< /highlight >}}

{{< highlight R >}}
[1] 16457    50
{{< /highlight >}}


{{% details title="使用するgene setの遺伝子行だけを抽出" %}}
{{< highlight R >}}
# 使用するgene setの遺伝子行だけを抽出
tumor_gene_set_log_normalized_counts <- tumor_log_normalized_counts[gene_set_tibble$gene_name, ]

tumor_gene_set_log_normalized_counts[1:5, 1:5]
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
       TCGA-D5-5540 TCGA-EI-6509 TCGA-A6-6137 TCGA-QG-A5Z2 TCGA-AA-3489
APC        9.940582     9.574465     9.831598    10.237207     10.38061
KRAS      11.824451    11.461010    11.381450    11.268820     11.07960
PIK3CA     8.961879     9.976906    10.173252     9.614568     10.14368
FBXW7      9.275326     8.763883     9.727003    10.057651     10.17140
SMAD4     11.526786     8.476560    10.933299    11.861563     10.80839
{{< /highlight >}}

{{< highlight R >}}
dim(tumor_gene_set_log_normalized_counts)
{{< /highlight >}}


{{< highlight R >}}
[1] 56 50
{{< /highlight >}}


- うまくデータを抽出できたらPCAを行ってplotしてみましょう。
- 今回はsourceは一つしかないのでStageで色分けしてみましょう。


{{% details title="PCAを行って結果をメタデータと結合したデータフレームを作成" %}}
{{< highlight R >}}
# PCAを行う（行をサンプルにする必要があるので転置する）
pca_result <- prcomp(t(tumor_gene_set_log_normalized_counts))

# PCAの結果とメタデータを結合したデータフレームを作成
plot_data_pca_gg <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stage = tumor_samples_metadata$diagnoses.ajcc_pathologic_stage,
  sample_id = rownames(pca_result$x)
)

head(plot_data_pca_gg)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
                        PC1       PC2      stage         sample_id
TCGA-D5-5540       1.041493 -4.475991  Stage IIA      TCGA-D5-5540
TCGA-EI-6509       3.471889  7.542652        '--      TCGA-EI-6509
TCGA-A6-6137       5.181010  1.625038 Stage IIIB      TCGA-A6-6137
TCGA-QG-A5Z2       6.057406 -3.254991    Stage I      TCGA-QG-A5Z2
TCGA-AA-3489      -7.073449  3.741021   Stage II      TCGA-AA-3489
HCM-CSHL-0160-C18 -3.227502  1.648932  Stage IIA HCM-CSHL-0160-C18
{{< /highlight >}}



{{% details title="ステージ情報を整理" %}}
{{< highlight R >}}
# ステージ情報を整理
plot_data_pca_gg <- plot_data_pca_gg %>% 
    mutate(stage_simplified = case_when(
        stage == "'--" ~ "Unknown",
        stage == "Normal" ~ "Normal",
        str_detect(stage, "^Stage IV") ~ "Stage_IV",
        str_detect(stage, "^Stage III") ~ "Stage_III",
        str_detect(stage, "^Stage II") ~ "Stage_II",
        str_detect(stage, "^Stage I") ~ "Stage_I",
        TRUE ~ "Other"
    ))

head(plot_data_pca_gg)
{{< /highlight >}}
{{% /details %}}

{{< highlight R >}}
                        PC1       PC2      stage         sample_id
TCGA-D5-5540       1.041493 -4.475991  Stage IIA      TCGA-D5-5540
TCGA-EI-6509       3.471889  7.542652        '--      TCGA-EI-6509
TCGA-A6-6137       5.181010  1.625038 Stage IIIB      TCGA-A6-6137
TCGA-QG-A5Z2       6.057406 -3.254991    Stage I      TCGA-QG-A5Z2
TCGA-AA-3489      -7.073449  3.741021   Stage II      TCGA-AA-3489
HCM-CSHL-0160-C18 -3.227502  1.648932  Stage IIA HCM-CSHL-0160-C18
                  stage_simplified
TCGA-D5-5540              Stage_II
TCGA-EI-6509               Unknown
TCGA-A6-6137             Stage_III
TCGA-QG-A5Z2               Stage_I
TCGA-AA-3489              Stage_II
HCM-CSHL-0160-C18         Stage_II
{{< /highlight >}}

{{< highlight R >}}
table(plot_data_pca_gg$stage, plot_data_pca_gg$stage_simplified)
{{< /highlight >}}




{{< highlight R >}}
             Stage_I Stage_II Stage_III Stage_IV Unknown
  '--              0        0         0        0       9
  Stage I          9        0         0        0       0
  Stage II         0        4         0        0       0
  Stage IIA        0       13         0        0       0
  Stage IIB        0        1         0        0       0
  Stage III        0        0         2        0       0
  Stage IIIA       0        0         1        0       0
  Stage IIIB       0        0         5        0       0
  Stage IIIC       0        0         3        0       0
  Stage IV         0        0         0        1       0
  Stage IVA        0        0         0        1       0
  Stage IVB        0        0         0        1       0
{{< /highlight >}}

- k-meansクラスタリングしてplotします。


{{% details title="k-meansクラスタリングしてplot" %}}
{{< highlight R >}}
# k-meansクラスタリング (PC1とPC2を使用)
set.seed(123)
kmeans_obj <- kmeans(plot_data_pca_gg[, c("PC1", "PC2")], centers = 3, nstart = 25)
plot_data_pca_gg$cluster <- factor(kmeans_obj$cluster)

# クラスタの中心点をデータフレームとして準備
cluster_centers_gg <- as.data.frame(kmeans_obj$centers)
colnames(cluster_centers_gg) <- c("PC1", "PC2")

# クラスタを中心点から信頼楕円を書くためにクラスタIDを追加
cluster_centers_gg$cluster <- factor(1:nrow(cluster_centers_gg))

# ggplot2でプロット
ggplot(plot_data_pca_gg, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = stage_simplified, shape = cluster), size = 3, alpha = 0.7) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.15, type = "t", level = 0.95, show.legend = FALSE, color = "darkgrey") + 
  scale_color_manual(values = stage_colors_updated) + 
  scale_fill_brewer(palette = "Pastel2") + 
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15)) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    y = paste0("PC2 (", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    color = "Stage",
    shape = "K-means Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-25-1.png" >}}


- k-meansクラスタリングは中心点の数を調整することができます。
- 4つにしてみましょう。


{{% details title="k-meansクラスタリング (PC1とPC2を使用)、中心点を4つにしてplot" %}}
{{< highlight R >}}
# k-meansクラスタリング (PC1とPC2を使用)
set.seed(123)
kmeans_obj <- kmeans(plot_data_pca_gg[, c("PC1", "PC2")], centers = 4, nstart = 25)
plot_data_pca_gg$cluster <- factor(kmeans_obj$cluster)

# クラスタの中心点をデータフレームとして準備
cluster_centers_gg <- as.data.frame(kmeans_obj$centers)
colnames(cluster_centers_gg) <- c("PC1", "PC2")

# クラスタを中心点から信頼楕円を書くためにクラスタIDを追加
cluster_centers_gg$cluster <- factor(1:nrow(cluster_centers_gg))

# ggplot2でプロット
ggplot(plot_data_pca_gg, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = stage_simplified, shape = cluster), size = 3, alpha = 0.7) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.15, type = "t", level = 0.95, show.legend = FALSE, color = "darkgrey") + 
  scale_color_manual(values = stage_colors_updated) + 
  scale_fill_brewer(palette = "Pastel2") + 
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18)) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    y = paste0("PC2 (", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1), "%)"),
    color = "Stage",
    shape = "K-means Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
{{< /highlight >}}
{{% /details %}}

{{< figure src="unnamed-chunk-26-1.png" >}}


- このように、heatmapとは違った形でサンプルのばらつきを可視化することができます。


## 6. 追加演習

- ここまでの演習では、データの前処理や可視化を行いました。
- 今回はメディアン比正規化したデータについて可視化を行いましたが、他の正規化だと結果はどのように変わるでしょうか？
- [TPM正規化したデータ](https://drive.google.com/uc?id=1PXcZ2psP1dBosYCpCSbiX_D180UpIkSv&export=download)も用意したので、同様の可視化を行って結果を比べてみましょう。
- また、例えばDEseq2パッケージにはPCAに最適化された正規化方法が実装されています(rlog, vstなど)。
- 余裕があったらこれらの正規化方法を試してみてください。
