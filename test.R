# 必要なパッケージを読み込む (未インストールの場合は install.packages() でインストール)
library(readr)
library(pheatmap)

# CSVファイル名を指定 (実際のパスに修正してください)
csv_file_path <- "your_rnaseq_counts.csv" 

# カウントデータを読み込む (1列目:遺伝子ID, それ以降:カウント値を想定)
count_data_df <- read_csv(csv_file_path)

# ヒートマップ用にデータを整形 (1列目を行名に、他を行列に)
gene_names <- count_data_df[[1]]
count_matrix <- as.matrix(count_data_df[, -1])
rownames(count_matrix) <- gene_names

# ヒートマップを生成 (行ごとにスケール、行名は非表示)
pheatmap(
  count_matrix, 
  scale = "row",
  show_rownames = FALSE
  # filename = "heatmap.png" # ファイル保存する場合
)
