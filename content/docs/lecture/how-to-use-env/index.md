---
title: "1. 解析環境を使ってみよう"
description: "HU Bioinfo Workshopで使用する解析環境とCursorの使い方を紹介します。"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# 解析環境を使ってみよう

## 0. 関連チュートリアル

- [HU Bioinfo Launcherの使い方]({{% ref "/docs/tutorials/HU-bioinfo-Launcher/index.md" %}})
- [Cursorの使い方]({{% ref "/docs/tutorials/Cursor/index.md" %}})
- [Linuxコマンドの使い方]({{% ref "/docs/tutorials/Linux-command/index.md" %}})

## 1. 解析環境のセットアップ
[HU Bioinfo Workshop 開発環境セットアップガイド]({{% ref "/docs/start-up/" %}})を参照してください。

## 2. 仮想環境について

HU Bioinfo Workshop で使用する解析環境を使う前に、仮想環境について説明します。

「仮想環境」は解析作業を行う上で重要な概念です。

仮想環境とは大まかにいうと「パソコンの中を仮想的に分割する」という技術です。  
ただしその分け方にはいくつかのレベルがあります。

1. **仮想マシン**: パソコンの中に別のパソコンを作ってその中でOSを起動する
- 全く別のパソコンを作るようなものなのでファイルを保存する場所も根本から違う
- 「`C:\Users\black\Documents`」と「`\\wsl.localhost\Ubuntu\home\shuu5`」のように一番最初の分岐から異なる。
{{< hint info >}}完全に別の家、2世帯住宅。{{< /hint >}}
2. **コンテナ**: 一つのOSの中で一部のアプリケーションやファイルを隔離する  
- マシンとOSは一緒なのでファイルを保存する大本は一緒。
- コンテナの中で利用するためにはコンテナの中にファイルを入れる作業（＝**マウント** ）が必要
{{< hint info >}}一つの家の中で部屋を分ける。部屋の中に使う資料を入れる。{{< /hint >}}
3. **言語仮想環境**: 特定のプログラミング言語に必要なファイルを隔離する。
- 特定のディレクトリ（フォルダ）に必要なものをまとめる。
- 日常的に行うフォルダによるファイル整理と全く同じ。  
{{< hint info >}}一つの部屋の中でしまう棚や引き出しを分ける。{{< /hint >}}

セットアップガイドで作成した環境を図示したものを以下に示します。

{{< figure 
    src="wsl-devcontainer-figure.drawio.svg"
    alt="仮想マシンとコンテナの関係"
    caption="仮想マシンとコンテナの関係"
>}}

{{% hint info %}}
Windows Subsystem for Linux 2 (WSL2), もしくはDocker Desktopが**仮想マシン(Virtual Machine)** を作成し、その中でUbuntuのOSを起動しています。 

さらにDocker Engine（Macの場合はDocker Desktopの中に含まれているDocker Engine）が**コンテナ(Docker Container)** を作成し、その中でUbuntuのOSを起動しています。
{{% /hint %}}

{{% hint info %}}
どちらもさらにvenvやrenvといった**言語仮想環境**を作成し、その中でPythonやRなどのプログラミング言語を起動しています。
{{% /hint %}}

仮想環境を利用することで、必要なツールやライブラリのバージョンが他のプロジェクトと衝突するのを防ぎ、安定した解析環境を維持できます。




## 3. HU Bioinfo Workshop Launcherの主な機能

HU Bioinfo Workshop Launcher（具体的には `Start bioinfo-launcher` コマンド）は、実行する際の状況によって、主に以下の2つの異なる動作をします。

#### 初回実行時：解析環境の初期セットアップと起動

初めて `Start bioinfo-launcher` コマンドを実行する場合、解析環境コンテナをゼロからセットアップします。

初期設定の詳細な手順については、[HU Bioinfo Workshop 開発環境セットアップガイド]({{% ref "/docs/start-up/" %}})の「6. HU-Bioinfo Workshop Launcher 拡張機能のセットアップと実行」を参照してください。

#### 2回目以降の実行時：既存の解析環境の起動

一度初期セットアップが完了していれば、2回目以降の `Start bioinfo-launcher` コマンド実行時には、既にセットアップ済みの環境を再利用します。

1.  **既存環境での起動:** 初回設定時に作成された `container` ディレクトリに自動的に移動し、その中の設定ファイルを使ってDev Containersで解析環境コンテナを開きます。前回の作業状態を維持したまま、すぐに作業を再開できます。

このように、Launcher は初回はセットアップを、2回目以降は既存環境での作業開始を自動化してくれる便利なツールです。

{{% hint info %}}
とにかくコマンドパレット(`Ctrl+Shift+P`, `Cmd+Shift+P` (Mac))で  
「**bioinfo-launcher: Start bioinfo-launcher**」を実行すると解析環境コンテナを開けます！
{{% /hint %}}

{{% hint info %}}
一度解析環境コンテナを開くと、Cursorを新しく立ち上げたり新しいウィンドウで開いたり（`Ctrl+Shift+N`, `Cmd+Shift+N` (Mac)）した時に、**Recent projects** が表示され、  

`project名 [開発コンテナー: bioinfo-launcher]`

のように表示されます。これをクリックするだけでも解析環境コンテナに入れます。
{{< /hint >}}



## 4. projectの考え方・作り方

### 4.1. Projectとは

{{% hint warning %}}
ここから先は**解析環境コンテナに入った状態** を前提として説明します。
{{% /hint %}}

これからいろいろな解析を行うことになりますが、その解析を行うためのファイルをまとめて管理するための方法を説明します。

解析にはいろいろな種類があり、それぞれの解析に必要なファイルやデータ、プログラムが異なります。
すべてひとまとめにしてしまうと、どれがどの解析に必要なファイルだったかが分からなくなったり、古いプログラムと新しいプログラムが混ざってエラーを起こしてしまうことがあります。

そこで、解析ごとにファイルを**ディレクトリdirectory（＝フォルダfolder）** にまとめて管理します。  
{{% hint info %}}
このようなディレクトリを**Project** と呼びます。
{{% /hint %}}

なにか一つ解析を始めようと思ったら、解析環境コンテナを立ち上げて入り、その中でProjectを作成します。

![devcontainer-projects-figure.drawio.svg](devcontainer-projects-figure.drawio.svg)

上図は二つのProjectを作成した例です。

RやPythonという言語には解析に使える**パッケージ**というものがあり、それぞれのパッケージにはそれぞれの解析に必要なファイルやデータ、作成済みのプログラムが用意されています。

解析に合わせて必要なパッケージを選択し、そのパッケージの組み合わせを記録するために**言語仮想環境** を作成します。

{{% hint info %}}
今回の例ではrenvやvenvが言語仮想環境にあたります。
{{% /hint %}}

Projectは一つの言語で作成できることが多いですが、高度な解析だと複数の言語を使うこともあります。

### 4.2. Projectの作成

今回作った環境において、Projectの作成は`prem`というコマンドによって自動化されています。

**細かいことは考えずにこのコマンドを実行しましょう！**

- `lecture1`という名前のProjectを作成したい場合
```bash
prem lecture1
```

このコマンドを実行すると`venv`(Pythonの言語仮想環境)と`renv`(Rの言語仮想環境)が自動的に作成されます。

プロジェクトディレクトリを作成できたら、Cursorの機能でこのディレクトリを開きます。

```bash
cursor lecture1
```

このコマンドを実行すると、新しいウィンドウでcursorが開き、目的のディレクトリが表示されます。

{{% hint info %}}
- コマンドパレット(`Ctrl+Shift+P`, `Cmd+Shift+P` (Mac))で`File: Open Folder`
- 左上の「ファイル」→「フォルダを開く」を選択（`Ctrl+M Ctrl+O` (Windows), `Cmd+M Cmd+O` (Mac)）

このような方法で開くこともできます。
{{% /hint %}}



## 5. R環境のデモ

次に、この `lecture1` プロジェクト内でRのパッケージをインストールし、簡単な解析を試してみましょう。

まず試しにR言語を使用するためにR script(拡張子が.Rのファイル)を作成します。

```bash
touch code/lecture1.R
```

`code`というディレクトリの中に`lecture1.R`というR言語のプログラムを記載するためのファイルが作成されます。

{{% hint info %}}
- 左のExplorer (エクスプローラー) ビューで、Projectディレクトリ（例: lecture1）を右クリックし、「New File...」を選択してファイル名を指定する。

CursorのGUIを使ってより直感的な操作も可能です。
{{% /hint %}}



### 5.1. パッケージのインストール

Rのパッケージは、特定の機能（データ読み込み、統計解析、グラフ作成など）を提供してくれる拡張機能のようなものです。`renv` を使っているプロジェクトでは、`renv::install()` 関数を使ってパッケージをインストールするのが一般的です。

ここでは例として、CSVファイルを簡単に読み込むための`readr`パッケージをインストールしてみましょう。

`lecture1.R` ファイルに以下のコードを追加してください。

{{< highlight R "linenos=inline, linenostart=1" >}}
# readr パッケージをインストール
renv::install("readr")
{{< /highlight >}}

コードを追加したら、追加した部分にカーソルを合わせて `Ctrl+Enter` (Windows), `Cmd+Enter` (Mac) を押して実行します。

インストールが完了すると、コンソールにメッセージが表示されます。`renv` はプロジェクトごとにインストールしたパッケージを管理してくれるため、他のプロジェクトに影響を与えません。

### 5.2. インストールしたパッケージを使ってみる

`readr` パッケージを使って、簡単なデータ読み込みのデモを行います。まず、試しに読み込むためのCSVファイルを作成しましょう。

Projectディレクトリ（`lecture1`）内の、`data`というディレクトリの中に`sample_data.csv` という名前で新しいファイルを作成し、以下の内容をコピーアンドペーストしてください。

{{< highlight csv "linenos=inline, linenostart=1" >}}
gene_name, sample1, sample2, sample3
gene1, 10, 20, 30
gene2, 15, 25, 35
gene3, 20, 30, 40
gene4, 25, 35, 45
gene5, 30, 40, 50
{{< /highlight >}}

次に、`lecture1.R` ファイルに以下のコードを追加します。

{{< highlight R "linenos=inline, linenostart=3" >}}
# readr パッケージを読み込みます
library(readr)

# sample_data.csv ファイルを読み込みます
data <- read_csv("sample_data.csv")
{{< /highlight >}}

こちらもコードを追加したら、追加した部分にカーソルを合わせて実行してみてください。

{{% hint info %}}
- `library(readr)` は、インストールした `readr` パッケージを現在のRセッションで使用可能にするためのコマンドです。
- パッケージを使う前には必ず`library()`で使いたいパッケージを読み込む必要があります。
{{% /hint %}}

次に読み込んだデータを表示してみます。

{{< highlight R "linenos=inline, linenostart=1" >}}
print(data)
#   gene_name sample1 sample2 sample3
#   <chr>     <dbl>   <dbl>   <dbl>
# 1 gene1        10      20      30
# 2 gene2        15      25      35
# 3 gene3        20      30      40
# 4 gene4        25      35      45
# 5 gene5        30      40      50
{{< /highlight >}}

コンソールに `sample_data.csv` の内容が表示されれば成功です。これは `readr` パッケージの `read_csv` 関数を使ってファイルが正しく読み込めたことを示しています。

{{% hint info %}}
- このような`data`を**tibble**と言います。
- **列名がついた列** と **行名がついていない行** からなるデータ形式です。
{{% /hint %}}

解析は基本的に行と列の形で整理されたデータを使います。

データの中から一部の行や列、値を取り出してみましょう。

{{< highlight R "linenos=inline, linenostart=1" >}}
# 最初の3行を表示
head(data, 3)
#   gene_name sample1 sample2 sample3
#   <chr>     <dbl>   <dbl>   <dbl>
# 1 gene1        10      20      30
# 2 gene2        15      25      35
# 3 gene3        20      30      40

# gene_name列を取り出す
print(data$gene_name)
# [1] gene1 gene2 gene3 gene4 gene5

# sample1列を取り出す
print(data$sample1)
# [1] 10 15 20 25 30

# 列番号で列を取り出す
print(data[, 1])
# [1] gene1 gene2 gene3 gene4 gene5

# 行番号で行を取り出す
print(data[1, ])
# [1] gene1 10    20    30

# 特定の値を取り出す
print(data[1, 2])
# [1] 10

# 特定の範囲の行と列を取り出す
print(data[1:2, 2:4])
#   sample1 sample2 sample3
#   <dbl>   <dbl>   <dbl>
# 1    10      20      30
# 2    15      25      35
{{< /highlight >}}

次は簡単な図を作成してみましょう。

{{< highlight R "linenos=inline, linenostart=5" >}}
# barplot 関数を使って棒グラフを作成
barplot(
    data$sample1, 
    names.arg = data$gene_name, 
    main = "Sample 1"
    )
{{< /highlight >}}

このコードを実行すると、`sample1` 列のデータを棒グラフで表示します。
