---
title: "Cursor"
weight: 2
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# Cursorの使い方

[Cursor](https://www.cursor.com/)は、VS Codeというプログラミング用エディタを基に、AIによるプログラミングをサポートするように作られたエディタです。  

セットアップについては[HU Bioinfo Workshop 開発環境セットアップガイド]({{% ref "/docs/start-up/" %}})を参照してください。

基本的な使い方は分かったので、CursorのAI機能だけ知りたいという方は[CursorのAI機能](#cursorのai機能)を参照してください。

## コマンドパレット

- Cursorは**Graphical User Interface (GUI)** なのでクリックしたりドラッグしたりで操作できます。
- 一方でマウスを使わずにコマンドパレットを使ってキーボードだけで操作することもできます。
- コマンドパレットは`Ctrl+Shift+P`(win) or `Cmd+Shift+P`(mac)で開きます。

{{% hint info %}}
- このような操作方法のインターフェースを**Command Line Interface (CLI)** といいます。
- コマンドを覚えていればGUIよりも素早く操作できます。
{{% /hint %}}

## ユーザーインターフェース

{{< figure src="cursor-ui.drawio.svg" alt="Cursorのユーザーインターフェース">}}

### **1. サイドバー（左）とアクティビティバー**

{{< figure src="cursor-side.drawio.svg" alt="Cursorのサイドバー">}}

- サイドバーとその上部にあるアクティビティバー
- 画像では**エクスプローラー**を選んでいるので、今開いているディレクトリ（＝**カレントディレクトリ**）のディレクトリやファイル一覧が表示されている。
- アクティビティバーのボタンを押すと機能を切り替えられるがとりあえず**エクスプローラー**と**拡張機能**だけ知っていれば大丈夫。
- カレントディレクトリを切り替える場合はコマンドパレット(`Ctrl+Shift+P`(win) or `Cmd+Shift+P`(mac))を開いて`File: Open Folder`と入力して選択する。
- 左上の`ファイル`ボタンから操作することも可能です。

### **2. エディタエリア**

{{< figure src="cursor-editor.drawio.svg" alt="Cursorのエディタエリア">}}

- エディタエリアが実際にコードを書く場所です。
- 上部のタブで開いているファイルを切り替えることができます。
- エクスプローラーのコードファイルを`1回クリック`でプレビュー表示、`2回クリック`で編集モードに切り替わります。
- コードを編集した後は`Ctrl+S`(win) or `Cmd+S`(mac)で保存します。

### **3. ターミナルエリア**

{{< figure src="cursor-terminal.drawio.svg" alt="Cursorのターミナルエリア">}}

- ターミナルエリアはコマンドを直接打ち込んで実行する場所です。
- `ctrl+@`(win) or `cmd+@`(mac)でbashというLinux Commandの**シェル**を開くことができます。
- 言語ごとにシェルは異なります。複数開いて切り替えることも可能です。
- ターミナルにカーソルを合わせて`ctrl+K`(win) or `cmd+K`(mac)を押すと**AIによるコマンドの提案**を受けることができます。
    - 例えば" すべてのファイルの一覧を表示するコマンドは？"と聞くと`ls -la`というコマンドを提案してくれます。


### **4. AIチャットエリア**

{{< figure src="cursor-chat.drawio.svg" alt="CursorのAIチャットエリア">}}

- 右側のサイドバーは**AIチャット**です。AIに質問したり指示を出したりできます。
- `ctrl+L`(win) or `cmd+L`(mac)でチャットエリアを開くことができます。
- 詳しい使い方は[CursorのAI機能](#cursorのai機能)を参照してください。


## CursorのAI機能

### コードを書きながらリアルタイムでコードを補完(**Cursor Tab**)

{{< figure src="cursor-tab.drawio.svg" alt="CursorのCursor Tab">}}

{{< figure src="cursor-tab-example.png" alt="Cursor Tabの例">}}

- コードを書きながらリアルタイムでコードを補完する機能です。
- コードを少し書くと灰色の文字でコードの提案が出てきます。
- 提案を受け入れるには`Tab`キー、拒否するには`Esc`キーを押します。
- `ctrl+→`(win) or `cmd+→`(mac)で**単語単位で部分的に**提案を受け入れます。
- 次に修正すべき部分を自動で認識し、`Tab`で受け入れるとその場所にカーソルが移動します。


### AI chatでコードを書く

