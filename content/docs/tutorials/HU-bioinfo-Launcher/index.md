---
title: "HU Bioinfo Launcher"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# HU Bioinfo Launcher

## HU bioinfo-launcher とは？

Cursorの拡張機能で、データ解析に使うRやPythonの環境を簡単に用意してくれます。  
Docker（仮想的な開発部屋）を使って、誰でも同じ開発環境をすぐに使えるようにします。

### 機能

*   **統一環境:** RとPythonを同じ開発部屋（Dockerコンテナ）で使えます。
*   **簡単設定:** VS Codeの画面でポチポチするだけで環境ができます。
*   **キャッシュ共有:** よく使う部品（ライブラリなど）を共有して、効率的に作業できます。

### 使うには？

*   LinuxかMacで使えます。WindowsはWSL（WindowsでLinuxを使う仕組み）経由で使えます。
*   Dockerが必要です（自動インストール機能もあります）。
*   識別子：`hu-bioinfo-workshop.bioinfo-launcher` 拡張機能で検索

### 使い方の流れ

1.  Cursorのコマンドパレット（Ctrl+Shift+P または Cmd+Shift+P）を開きます。
2.  `bioinfo-launcher: Start bioinfo-launcher` を選んで実行。
3.  作業場所とGitHubの合鍵（PAT）を設定します。
4.  環境ができたら、新しいCursorウィンドウが開き、そこで作業開始！

{{% details title="GitHubの合鍵（PAT）の準備" %}}

1.  **GitHubアカウントの作成/サインイン**
    *   持っていない場合は [GitHub公式サイト](https://github.com/) で作成します。
    *   作成済みならサインインします。

2.  **PAT設定ページへ移動**
    *   右上のプロフィールアイコン → `Settings` → 左メニューの `Developer settings` → `Personal access tokens` → `Tokens (classic)` と進みます。

3.  **新しいPATを生成 (classic)**
    *   `Generate new token` ボタン → `Generate new token (classic)` を選択。
    *   (パスワード確認があれば入力)

4.  **PATの詳細設定**
    *   **Note (メモ):** 分かりやすい名前 (例: `Bioinfo Launcher WSL Token`) を入力。
    *   **Expiration (有効期限):** `30 days` や `90 days` などを選択 (無期限は非推奨)。
    *   **Select scopes (権限スコープ):** 以下にチェックを入れます。
        *   `repo` (プライベートリポジトリへのアクセス)
        *   `read:packages` (GitHub Packagesからの読み取り)
        *   `workflow` (GitHub Actionsワークフロー操作)

5.  **PATを生成して安全に保管**
    *   `Generate token` ボタンをクリック。
    *   表示されたPAT (`ghp_`で始まる文字列) は**一度しか表示されません。必ずコピーし、安全な場所に保管してください。**

このPATは後ほど使用します。パスワード同様、大切に扱ってください。

{{< hint warning >}}
**注意:** PATを紛失した、もしくは漏洩した場合は必ず **そのPATを削除** して新しく作り直してください。
{{< /hint >}}
{{% /details %}}

### ファイルの場所

設定した作業場所のフォルダに、`cache/`, `container/` といったフォルダができます。

{{% hint warning %}}
### 別の場所にコンテナを作り直すには

1.  Cursorのコマンドパレット（Ctrl+Shift+P または Cmd+Shift+P）を開きます。
2.  `bioinfo-launcher: Reset bioinfo-launcher Config` を選んで実行。
3.  新しい作業場所とGitHubの合鍵（PAT）を設定します。
4.  環境ができたら、新しいCursorウィンドウが開き、そこで作業開始！
{{% /hint %}}
