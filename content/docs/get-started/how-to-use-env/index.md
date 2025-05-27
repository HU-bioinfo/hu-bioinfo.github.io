---
title: "解析環境を使ってみよう"
description: "HU Bioinfo Workshopで使用する解析環境とCursorの使い方を紹介します。"
weight: 3
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# 解析環境を使ってみよう

## 0. 関連チュートリアル

- [HU Bioinfo Launcherの使い方]({{% ref "/docs/get-started/HU-bioinfo-Launcher/index.md" %}})
- [Cursorの使い方]({{% ref "/docs/get-started/Cursor/index.md" %}})
- [Linuxコマンドの使い方]({{% ref "/docs/get-started/linux-command/index.md" %}})

## 1. 解析環境のセットアップ
[HU Bioinfo Workshop 開発環境セットアップガイド]({{% ref "/docs/get-started/start-up/" %}})を参照してください。

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

初期設定の詳細な手順については、[HU Bioinfo Workshop 開発環境セットアップガイド]({{% ref "/docs/get-started/start-up/" %}})の「6. HU-Bioinfo Workshop Launcher 拡張機能のセットアップと実行」を参照してください。

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

- `playground`という名前のProjectを作成したい場合
```bash
prem playground
```

このコマンドを実行すると`venv`(Pythonの言語仮想環境)と`renv`(Rの言語仮想環境)が自動的に作成されます。

プロジェクトディレクトリを作成できたら、Cursorの機能でこのディレクトリを開きます。

```bash
cursor playground
```

このコマンドを実行すると、新しいウィンドウでcursorが開き、目的のディレクトリが表示されます。

{{% hint info %}}
- コマンドパレット(`Ctrl+Shift+P`, `Cmd+Shift+P` (Mac))で`File: Open Folder`
- 左上の「ファイル」→「フォルダを開く」を選択（`Ctrl+M Ctrl+O` (Windows), `Cmd+M Cmd+O` (Mac)）

このような方法で開くこともできます。
{{% /hint %}}
