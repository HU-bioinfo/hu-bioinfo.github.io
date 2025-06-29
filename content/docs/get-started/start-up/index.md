---
title: "セットアップ"
description: "HU Bioinfoワークショップ用 解析環境セットアップ (Windows/MacOS)"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# HU Bioinfo Workshop 解析環境セットアップガイド

{{% tabs "windows-version" %}}
{{% tab "Windows 11" %}}

## Windows 11 でのセットアップ

{{< details title = "1. WSL2 (Windows Subsystem for Linux 2) の導入" >}}

Windows上でLinux環境を動かすための「WSL2」をセットアップします。これにより、Linuxベースのツールが使えるようになります。

1.  **PowerShellを管理者として起動**
    スタートメニューで "PowerShell" を検索し、「管理者として実行」を選びます。

2.  **WSLのインストール**
    ```powershell
    wsl --install
    ```
{{% hint info %}}
自動的に**Ubuntu**というLinuxディストリビューションがインストールされます。

- Linuxとは、オープンソースのOSです。  
特にその中核部分（＝**カーネル**）のことを言います。
- Linuxカーネルに様々なソフトウェアをまとめてパッケージにしたものを  
**Linuxディストリビューション**と言います。
- どのようなソフトウェアとまとめるかによって様々なディストリビューションが存在します。
- UbuntuはそのようなLinuxディストリビューションの一つです。

{{% /hint %}}

3.  **Ubuntuの初期設定**
    インストール完了後、PowerShellで `wsl` と入力してUbuntuを起動します。
    初回起動時、ユーザー名とパスワードの設定を求められます。画面の指示に従い設定し、**忘れないようにメモしておきましょう。**

4.  **Ubuntuのパッケージリストを更新・アップデート**
    開発用のツールのインストールがスムーズに進むように、パッケージリストを更新・アップデートします。

    ```bash
    sudo apt update
    sudo apt upgrade -y
    ```

{{< hint warning >}}
**Ubuntu ユーザー名とパスワード:** この後すぐに使用します。忘れた場合は以下の方法で再設定が可能です。

{{< details title = "パスワード再設定手順" >}}

1.  **WindowsのコマンドプロンプトまたはPowerShellを管理者として開きます。**
2.  **以下のコマンドでWSLをrootユーザーで起動します。** (`YourUbuntuDistro` の部分は、`wsl -l -v` で確認できる実際のディストリビューション名 (例: `Ubuntu-24.04`) に置き換えてください。)
    ```powershell
    wsl -d YourUbuntuDistro -u root
    ```
3.  **rootユーザーとしてUbuntuが起動したら、以下のコマンドでパスワードをリセットしたいユーザーのパスワードを設定します。** (`yourusername` の部分は、ご自身のUbuntuのユーザー名に置き換えてください。ユーザー名が分からない場合は、`ls /home` コマンドで確認できます。パスワード入力時は画面に文字が表示されませんが、そのまま入力してEnterを押してください。)
    ```bash
    passwd yourusername
    ```
4.  **新しいパスワードを2回入力し、設定が完了したら `exit` コマンドでrootシェルを終了します。**
5.  **Windowsのコマンドプロンプト/PowerShellで、以下のコマンドを実行し、デフォルトユーザーを元に戻します。** (`YourUbuntuDistro` と `yourusername` は適切に置き換えてください。)
    ```powershell
    YourUbuntuDistro config --default-user yourusername
    ```
    *(`Ubuntu` コマンドで起動している場合は `ubuntu config --default-user yourusername` となります。)*

これでパスワードが再設定されました。WSLを通常起動して、新しいパスワードでログインできることを確認してください。
{{% /details %}}
{{< /hint >}}

これでWSL2とUbuntuの基本的なセットアップは完了です。
{{< /details >}}

{{< details title = "2. VSCode のインストール" >}}

マイクロソフト製のコードエディタ「Visual Studio Code (VSCode)」をインストールします。

1.  **公式サイトからダウンロード**
    [VSCode公式サイト](https://code.visualstudio.com/) にアクセスし、Windows版インストーラーをダウンロードします。

2.  **インストール実行**
    ダウンロードした `.exe` ファイルをダブルクリックして実行し、画面の指示に従ってインストールを進めます。
    「このアプリがデバイスに変更を加えることを許可しますか？」と表示されたら「はい」を選択してください。

インストールが完了すると、VSCodeが使えるようになります。
{{< /details >}}

{{< details title = "3. VSCodeの日本語化" >}}

VSCodeの表示を日本語にします。

1.  **VSCodeを起動**

2.  **プライマリサイドバーを開く**
    `Ctrl + B`を押します。

3.  **拡張機能ビューを開く**
    左側のアクティビティバーにある四角いアイコン（または `Ctrl+Shift+X`）で開きます。

{{< figure src="cursor-side.drawio.svg" alt="VSCodeのサイドバー">}}

3.  **日本語言語パックを検索・インストール**
    検索ボックスに `Japanese Language Pack` と入力し、`Japanese Language Pack for Visual Studio Code` (Microsoft提供) を見つけて `Install` ボタンをクリックします。

4.  **VSCodeを再起動**
    インストール後、右下に再起動を促す通知が出たら「再起動」をクリックします。出ない場合は手動でVSCodeを再起動してください。

これでVSCodeが日本語表示になります。
{{< /details >}}

{{< details title = "4. 便利な拡張機能のインストール" >}}

開発効率を上げるため、VSCodeに以下の拡張機能をインストールします。WSL連携やコンテナ開発に役立ちます。
（日本語化の時と同じように、拡張機能ビューからIDで検索してインストールします）

#### 4.1 Remote - WSL (`ms-vscode-remote.remote-wsl`)

Windows Subsystem for Linux (WSL) と連携するための拡張機能です。VSCodeからWSL内のプロジェクトを直接扱えるようになります。

1.  拡張機能ビューで `ms-vscode-remote.remote-wsl` を検索。
2.  `Remote - WSL` (Microsoft提供) をインストール。

#### 4.2 Dev Containers (`ms-vscode-remote.remote-containers`)

Dockerコンテナを開発環境として使うための拡張機能です。プロジェクトごとに独立した環境を構築できます。

{{% hint info %}}
**Dev Containers拡張機能について:**
この拡張機能を使用して初めてコンテナ環境を作成する際に、Dockerが必要な場合は自動的にインストールプロセスが開始されます。
{{% /hint %}}

1.  拡張機能ビューで `ms-vscode-remote.remote-containers` を検索。
2.  `Dev Containers` (Microsoft提供) をインストール。

これらの拡張機能で、より高度な開発が可能になります。
{{< /details >}}

{{< details title = "5. VSCodeからWSL (Ubuntu) へ接続" >}}

インストールした `Remote - WSL` 拡張機能を使って、VSCodeからWSL上のUbuntu環境に接続します。

1.  **VSCodeを起動**

2.  **コマンドパレットを開く**
    `Ctrl+Shift+P` を押します。

3.  **WSL接続コマンドを実行**
    コマンドパレットに `WSL` と入力し、候補から `WSL: WSL に接続` (または `"WSL: Connect to WSL"`) を選択して実行します。

4.  **接続するディストリビューションを選択**
    利用可能なWSLディストリビューションのリストから `Ubuntu` を選択します。

5.  **WSL (Ubuntu) 環境へ接続完了**
    新しいVSCodeウィンドウが開き、WSL上のUbuntu環境に接続されます。
    ウィンドウ左下に `WSL: Ubuntu` のように表示されていれば成功です。
    （初回接続時は、必要なコンポーネントのインストールに少し時間がかかることがあります。）

6.  **接続確認 (任意)**
    VSCode内で新しいターミナルを開きます (`Ctrl+@` またはメニューから `ターミナル` > `新しいターミナル`)。
    以下のようなプロンプトが表示されれば、Ubuntuに接続できています。
    ```bash
    your_username@your_hostname:~$
    ```
    (`your_username` はUbuntuのユーザー名、`your_hostname` はPC名です。)

    試しに以下のコマンドを実行してみましょう。
    ```bash
    pwd       # 現在のディレクトリ (例: /home/your_username)
    ls -la    # ファイルやフォルダの一覧
    uname -a  # Linuxカーネル情報
    ```
    コマンドが正常に実行されれば、接続は完璧です。

これで、VSCodeからWSL上のUbuntu環境で開発作業を開始できます。
{{< /details >}}

{{< details title = "6. HU bioinfo launcher 拡張機能のセットアップと実行" >}}

これは、本ワークショップ用のR/Python解析環境を簡単に構築できるVSCode拡張機能です。Dockerコンテナ技術を使用します。

{{% hint warning %}}
**重要:** この拡張機能のインストールと実行は、**VSCodeがWSL (Ubuntu) に接続された状態**で行います。
{{% /hint %}}

#### 6.1 拡張機能のインストール (WSL接続環境で)

1.  **WSL (Ubuntu) 接続の確認**
    VSCodeウィンドウ左下のステータスバーが `WSL: Ubuntu` (または類似の表示) になっていることを確認します。

2.  **拡張機能ビューを開く** (`Ctrl+Shift+X`)

3.  **`hu-bioinfo-workshop.bioinfo-launcher` を検索・インストール**

#### 6.2 GitHub Personal Access Token (PAT) の準備

GitHub上のリソースにアクセスするために、Personal Access Token (PAT) が必要です。

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

#### 6.3 解析作業用ディレクトリの作成とLauncherの実行

1.  **WSL (Ubuntu) ターミナルを開く** (VSCode内で `Ctrl+@`)
    プロンプトが `your_username@your_hostname:~$` であることを確認。

2.  **作業用親ディレクトリを作成**
    ターミナルで以下のコマンドを実行し、ホームディレクトリに `BioinfoSpace` ディレクトリを作成します。
    ```bash
    mkdir ~/BioinfoSpace
    ```

3.  **`HU bioinfo launcher` を実行**
    コマンドパレット (`Ctrl+Shift+P`) を開き、`bioinfo-launcher` と入力。
    候補から `bioinfo-launcher: Start bioinfo-launcher` を選択して実行します。

{{% hint info %}}
実行するとdocker imageのダウンロードが始まりますが**初回は時間がかかります**。
Dev Containers拡張機能により、必要に応じてDockerが自動的にインストールされます。
{{% /hint %}}

4.  **作業環境ディレクトリの設定**
    親ディレクトリを選択するように求められるので、先ほど作成した `BioinfoSpace` ディレクトリを選択します。

5.  **GitHub Personal Access Token (PAT) の設定**
    次にPATの入力を求められるので、準備しておいたPATを貼り付けます。

6.  **環境構築の開始**
    設定後、`HU bioinfo launcher` がDockerを使って解析環境の構築を開始します。

7.  **コンテナ内での作業開始**
    構築が完了すると、多くの場合、新しいVSCodeウィンドウが自動で開きます。これがDockerコンテナ内の開発環境です。
    ウィンドウ左下のステータスバーがコンテナ環境を示していれば成功です。

これで、`HU bioinfo launcher` のセットアップと実行は完了です。統一された解析環境で作業を始めましょう！
{{< /details >}}

{{< details title = "7. GitHub Copilot の学生向け無料セットアップ" >}}

GitHub Copilotは、AIペアプログラミングツールで、コードの提案やオートコンプリート機能を提供します。学生は認証を行うことで、通常有料のGitHub Copilotを無料で利用できます。

#### 7.1 学生認証の手続き

1. **GitHub Educationにアクセス**
   [GitHub Education](https://education.github.com/students) にアクセスします。

2. **学生申請を開始**
   「Join GitHub Education」または「Get student benefits」ボタンをクリックします。

3. **学生情報の入力**
   - 学校名を検索・選択
   - 学術メールアドレス（.ac.jpドメインなど）を入力
   - 在学証明書類のアップロード（学生証の写真など）

4. **申請の承認を待つ**
   通常、数日以内に承認されます。承認されるとメールで通知が届きます。

#### 7.2 GitHub Copilotの有効化

1. **GitHub設定ページへアクセス**
   [GitHub Settings](https://github.com/settings/copilot) にアクセスします。

2. **Copilotを有効化**
   学生認証が完了していれば、「Enable GitHub Copilot」ボタンが表示されます。クリックして有効化します。

#### 7.3 VSCodeでGitHub Copilotを使用

1. **GitHub Copilot拡張機能のインストール**
   VSCodeの拡張機能ビューで `GitHub.copilot` を検索し、「GitHub Copilot」をインストールします。

2. **GitHubアカウントでサインイン**
   拡張機能インストール後、GitHubアカウントでのサインインを求められます。指示に従ってサインインします。

3. **Copilotの使用開始**
   サインインが完了すると、コード編集時に自動的にAIによる提案が表示されるようになります。
   - `Tab`キー: 提案を受け入れる
   - `Esc`キー: 提案を拒否する
   - `Alt + ]`: 次の提案を表示
   - `Alt + [`: 前の提案を表示

{{< hint info >}}
**GitHub Copilotの主な機能:**
- コードの自動補完
- 関数やクラスの実装提案
- コメントからのコード生成
- テストコードの生成
- ドキュメントの作成支援
{{< /hint >}}

{{< hint warning >}}
**注意事項:**
- 学生認証は定期的に更新が必要です（通常は年1回）
- 卒業後は有料プランへの移行が必要になります
- 提案されたコードは必ず確認し、理解してから使用してください
{{< /hint >}}

これで、GitHub Copilotを使った効率的なコーディングが可能になります！
{{< /details >}}

{{% /tab %}}

{{% tab "MacOS" %}}

## MacOS でのセットアップ

MacOSはUnixベースのOSであり、WSLのようなLinux仮想環境は通常不要です。

{{< details title = "1. VSCode のインストール" >}}

マイクロソフト製のコードエディタ「Visual Studio Code (VSCode)」をインストールします。

1.  **公式サイトからダウンロード**
    [VSCode公式サイト](https://code.visualstudio.com/) にアクセスし、macOS版インストーラー (通常は `.zip` ファイル) をダウンロードします。

2.  **インストール実行**
    ダウンロードした `.zip` ファイルを解凍し、VSCodeアプリアイコンを `Applications` フォルダにドラッグ＆ドロップします。
    初回起動時にセキュリティ確認が表示された場合は、「開く」を選択してください。

インストールが完了すると、VSCodeが使えるようになります。
{{< /details >}}

{{< details title = "2. VSCodeの日本語化" >}}

VSCodeの表示を日本語にします。

1.  **VSCodeを起動**

2.  **プライマリサイドバーを開く**
    `Cmd + B`を押します。

3.  **拡張機能ビューを開く**
    左側のアクティビティバーにある四角いアイコン（または `Cmd+Shift+X`）で開きます。

{{< figure src="cursor-side.drawio.svg" alt="VSCodeのサイドバー">}}

3.  **日本語言語パックを検索・インストール**
    検索ボックスに `Japanese Language Pack` と入力し、`Japanese Language Pack for Visual Studio Code` (Microsoft提供) を見つけて `Install` ボタンをクリックします。

4.  **VSCodeを再起動**
    インストール後、右下に再起動を促す通知が出たら「再起動」をクリックします。出ない場合は手動でVSCodeを再起動してください。

これでVSCodeが日本語表示になります。
{{< /details >}}

{{< details title = "3. 便利な拡張機能のインストール" >}}

開発効率を上げるため、VSCodeに以下の拡張機能をインストールします。コンテナ開発に役立ちます。
（日本語化の時と同じように、拡張機能ビューからIDで検索してインストールします）

#### 3.1 Dev Containers (`ms-vscode-remote.remote-containers`)

Dockerコンテナを開発環境として使うための拡張機能です。プロジェクトごとに独立した環境を構築できます。

{{% hint info %}}
**Dev Containers拡張機能について:**
この拡張機能を使用して初めてコンテナ環境を作成する際に、Dockerが必要な場合は自動的にインストールプロセスが開始されます。
{{% /hint %}}

1.  拡張機能ビューで `ms-vscode-remote.remote-containers` を検索。
2.  `Dev Containers` (Microsoft提供) をインストール。

この拡張機能で、より高度な開発が可能になります。
{{< /details >}}

{{< details title = "4. HU bioinfo launcher 拡張機能のセットアップと実行" >}}

これは、本ワークショップ用のR/Python解析環境を簡単に構築できるVSCode拡張機能です。Dockerコンテナ技術を使用します。

#### 4.1 拡張機能のインストール

1.  **拡張機能ビューを開く** (VSCode内で `Cmd+Shift+X`)

2.  **`hu-bioinfo-workshop.bioinfo-launcher` を検索・インストール**

#### 4.2 GitHub Personal Access Token (PAT) の準備

GitHub上のリソースにアクセスするために、Personal Access Token (PAT) が必要です。

1.  **GitHubアカウントの作成/サインイン**
    *   持っていない場合は [GitHub公式サイト](https://github.com/) で作成します。
    *   作成済みならサインインします。

2.  **PAT設定ページへ移動**
    *   右上のプロフィールアイコン → `Settings` → 左メニューの `Developer settings` → `Personal access tokens` → `Tokens (classic)` と進みます。

3.  **新しいPATを生成 (classic)**
    *   `Generate new token` ボタン → `Generate new token (classic)` を選択。
    *   (パスワード確認があれば入力)

4.  **PATの詳細設定**
    *   **Note (メモ):** 分かりやすい名前 (例: `Bioinfo Launcher Mac Token`) を入力。
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

#### 4.3 解析作業用ディレクトリの作成とLauncherの実行

1.  **ターミナルを開く** (VSCode内で `Ctrl+@` またはMacのターミナルアプリ)
    プロンプトが `your_username@your_hostname:~$` (または類似の形式) であることを確認。

2.  **作業用親ディレクトリを作成**
    ターミナルで以下のコマンドを実行し、ホームディレクトリに `BioinfoSpace` ディレクトリを作成します。
    ```bash
    mkdir ~/BioinfoSpace
    ```

3.  **`HU bioinfo launcher` を実行**
    コマンドパレット (VSCode内で `Cmd+Shift+P`) を開き、`bioinfo-launcher` と入力。
    候補から `bioinfo-launcher: Start bioinfo-launcher` を選択して実行します。

{{% hint info %}}
実行するとdocker imageのダウンロードが始まりますが**初回は時間がかかります**。
Dev Containers拡張機能により、必要に応じてDockerが自動的にインストールされます。
{{% /hint %}}

4.  **作業環境ディレクトリの設定**
    親ディレクトリを選択するように求められるので、先ほど作成した `BioinfoSpace` ディレクトリを選択します。

5.  **GitHub Personal Access Token (PAT) の設定**
    次にPATの入力を求められるので、準備しておいたPATを貼り付けます。

6.  **環境構築の開始**
    設定後、`HU bioinfo launcher` がDockerを使って解析環境の構築を開始します。

7.  **コンテナ内での作業開始**
    構築が完了すると、多くの場合、新しいVSCodeウィンドウが自動で開きます。これがDockerコンテナ内の開発環境です。
    ウィンドウ左下のステータスバーがコンテナ環境を示していれば成功です。

これで、`HU bioinfo launcher` のセットアップと実行は完了です。統一された解析環境で作業を始めましょう！
{{< /details >}}

{{< details title = "5. GitHub Copilot の学生向け無料セットアップ" >}}

GitHub Copilotは、AIペアプログラミングツールで、コードの提案やオートコンプリート機能を提供します。学生は認証を行うことで、通常有料のGitHub Copilotを無料で利用できます。

#### 5.1 学生認証の手続き

1. **GitHub Educationにアクセス**
   [GitHub Education](https://education.github.com/students) にアクセスします。

2. **学生申請を開始**
   「Join GitHub Education」または「Get student benefits」ボタンをクリックします。

3. **学生情報の入力**
   - 学校名を検索・選択
   - 学術メールアドレス（.ac.jpドメインなど）を入力
   - 在学証明書類のアップロード（学生証の写真など）

4. **申請の承認を待つ**
   通常、数日以内に承認されます。承認されるとメールで通知が届きます。

#### 5.2 GitHub Copilotの有効化

1. **GitHub設定ページへアクセス**
   [GitHub Settings](https://github.com/settings/copilot) にアクセスします。

2. **Copilotを有効化**
   学生認証が完了していれば、「Enable GitHub Copilot」ボタンが表示されます。クリックして有効化します。

#### 5.3 VSCodeでGitHub Copilotを使用

1. **GitHub Copilot拡張機能のインストール**
   VSCodeの拡張機能ビューで `GitHub.copilot` を検索し、「GitHub Copilot」をインストールします。

2. **GitHubアカウントでサインイン**
   拡張機能インストール後、GitHubアカウントでのサインインを求められます。指示に従ってサインインします。

3. **Copilotの使用開始**
   サインインが完了すると、コード編集時に自動的にAIによる提案が表示されるようになります。
   - `Cmd + ]`: 次の提案を表示
   - `Cmd + [`: 前の提案を表示
   - `Tab`キー: 提案を受け入れる
   - `Esc`キー: 提案を拒否する

{{< hint info >}}
**GitHub Copilotの主な機能:**
- コードの自動補完
- 関数やクラスの実装提案
- コメントからのコード生成
- テストコードの生成
- ドキュメントの作成支援
{{< /hint >}}

{{< hint warning >}}
**注意事項:**
- 学生認証は定期的に更新が必要です（通常は年1回）
- 卒業後は有料プランへの移行が必要になります
- 提案されたコードは必ず確認し、理解してから使用してください
{{< /hint >}}

これで、GitHub Copilotを使った効率的なコーディングが可能になります！
{{< /details >}}

{{% /tab %}}
{{% /tabs %}}

