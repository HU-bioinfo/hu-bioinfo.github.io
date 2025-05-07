---
title: "Start up"
description: "HU Bioinfoワークショップ用 開発環境セットアップ (Windows編)"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# HU Bioinfo Workshop 開発環境セットアップガイド

{{% tabs "windows-version" %}}
{{% tab "Windows 11" %}}

## Windows 11 でのセットアップ

{{< details title = "1. WSL2 (Windows Subsystem for Linux 2) の導入" >}}

Windows上でLinux環境を動かすための「WSL2」をセットアップします。これにより、Linuxベースのツールが使えるようになります。

1.  **PowerShellを管理者として起動**
    スタートメニューで "PowerShell" を検索し、「管理者として実行」を選びます。

2.  **WSL機能の有効化**
    PowerShellに以下のコマンドを入力し、Enterキーを押します。
    ```powershell
    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
    ```

3.  **仮想マシンプラットフォーム機能の有効化**
    続けて、以下のコマンドも入力します。
    ```powershell
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
    ```

4.  **PCの再起動**
    上記2つのコマンド実行後、**必ずPCを再起動してください。**

5.  **Linuxカーネルの更新**
    再起動後、[WSL2 Linuxカーネル更新プログラムパッケージ(x64マシン用)](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi) をダウンロードし、インストールします。
    (ARM64マシンの場合は[こちら](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_arm64.msi)を使用)

6.  **WSL2をデフォルトバージョンに設定**
    再度PowerShellを管理者として起動し、以下を入力します。
    ```powershell
    wsl --set-default-version 2
    ```

7.  **Ubuntu 24.04 LTSのインストール**
    PowerShellで以下のコマンドを実行し、Ubuntuをインストールします。
    ```powershell
    wsl --install -d Ubuntu-24.04
    ```
    ネットワーク環境によっては少し時間がかかります。

8.  **Ubuntuの初期設定**
    インストール完了後、PowerShellで `wsl` または `wsl -d Ubuntu-24.04` と入力してUbuntuを起動します。
    初回起動時、ユーザー名とパスワードの設定を求められます。画面の指示に従い設定し、**忘れないようにメモしておきましょう。**

これでWSL2とUbuntuの基本的なセットアップは完了です。
{{< /details >}}

{{< details title = "2. Cursorエディタのインストール" >}}

AIコーディング支援機能付きのエディタ「Cursor」をインストールします。

1.  **公式サイトからダウンロード**
    [Cursor公式サイト](https://cursor.sh/) にアクセスし、Windows版インストーラーをダウンロードします。

2.  **インストール実行**
    ダウンロードした `.exe` ファイルをダブルクリックして実行し、画面の指示に従ってインストールを進めます。
    「このアプリがデバイスに変更を加えることを許可しますか？」と表示されたら「はい」を選択してください。

インストールが完了すると、Cursorが使えるようになります。
{{< /details >}}

{{< details title = "3. Cursorの日本語化" >}}

Cursorの表示を日本語にします。

1.  **Cursorを起動**

2.  **拡張機能ビューを開く**
    左側のアクティビティバーにある四角いアイコン（または `Ctrl+Shift+X`）で開きます。

3.  **日本語言語パックを検索・インストール**
    検索ボックスに `Japanese Language Pack` と入力し、`Japanese Language Pack for Visual Studio Code` (Microsoft提供) を見つけて `Install` ボタンをクリックします。

4.  **Cursorを再起動**
    インストール後、右下に再起動を促す通知が出たら「再起動」をクリックします。出ない場合は手動でCursorを再起動してください。

これでCursorが日本語表示になります。
{{< /details >}}

{{< details title = "4. 便利な拡張機能のインストール" >}}

開発効率を上げるため、Cursorに以下の拡張機能をインストールします。WSL連携やコンテナ開発に役立ちます。
（日本語化の時と同じように、拡張機能ビューからIDで検索してインストールします）

#### 4.1 Remote - WSL (`ms-vscode-remote.remote-wsl`)

Windows Subsystem for Linux (WSL) と連携するための拡張機能です。CursorからWSL内のプロジェクトを直接扱えるようになります。

1.  拡張機能ビューで `ms-vscode-remote.remote-wsl` を検索。
2.  `Remote - WSL` (Microsoft提供) をインストール。

#### 4.2 Remote - Containers (`ms-vscode-remote.remote-containers`)

Dockerコンテナを開発環境として使うための拡張機能です。プロジェクトごとに独立した環境を構築できます。

1.  拡張機能ビューで `ms-vscode-remote.remote-containers` を検索。
2.  `Remote - Containers` (Microsoft提供) をインストール。

これらの拡張機能で、より高度な開発が可能になります。
{{< /details >}}

{{< details title = "5. CursorからWSL (Ubuntu) へ接続" >}}

インストールした `Remote - WSL` 拡張機能を使って、CursorからWSL上のUbuntu環境に接続します。

1.  **Cursorを起動**

2.  **コマンドパレットを開く**
    `Ctrl+Shift+P` を押します。

3.  **WSL接続コマンドを実行**
    コマンドパレットに `WSL` と入力し、候補から `WSL: WSL に接続` (または `"WSL: Connect to WSL"`) を選択して実行します。

4.  **接続するディストリビューションを選択**
    利用可能なWSLディストリビューションのリストから `Ubuntu-24.04` を選択します。

5.  **WSL (Ubuntu) 環境へ接続完了**
    新しいCursorウィンドウが開き、WSL上のUbuntu環境に接続されます。
    ウィンドウ左下に `WSL: Ubuntu-24.04` のように表示されていれば成功です。
    （初回接続時は、必要なコンポーネントのインストールに少し時間がかかることがあります。）

6.  **接続確認 (任意)**
    Cursor内で新しいターミナルを開きます (`Ctrl+@` またはメニューから `ターミナル` > `新しいターミナル`)。
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

これで、CursorからWSL上のUbuntu環境で開発作業を開始できます。
{{< /details >}}

{{< details title = "6. HU-Bioinfo Workshop Launcher 拡張機能のセットアップと実行" >}}

`hu-bioinfo-workshop.bioinfo-launcher` は、本ワークショップ用のR/Python解析環境を簡単に構築できるCursor拡張機能です。Dockerコンテナ技術を使用します。

**重要:** この拡張機能のインストールと実行は、**CursorがWSL (Ubuntu) に接続された状態**で行います。

#### 6.1 拡張機能のインストール (WSL接続環境で)

1.  **WSL (Ubuntu) 接続の確認**
    Cursorウィンドウ左下のステータスバーが `WSL: Ubuntu-24.04` (または類似の表示) になっていることを確認します。

2.  **拡張機能ビューを開く** (`Ctrl+Shift+X`)

3.  **`hu-bioinfo-workshop.bioinfo-launcher` を検索・インストール**
    検索ボックスに `hu-bioinfo-workshop.bioinfo-launcher` と入力し、`HU-Bioinfo Workshop Launcher` をインストールします。

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

#### 6.3 解析作業用ディレクトリの作成とLauncherの実行

1.  **WSL (Ubuntu) ターミナルを開く** (Cursor内で `Ctrl+@`)
    プロンプトが `your_username@your_hostname:~$` であることを確認。

2.  **作業用親ディレクトリを作成**
    ターミナルで以下のコマンドを実行し、ホームディレクトリに `BioinfoSpace` ディレクトリを作成して移動します。
    ```bash
    mkdir ~/BioinfoSpace
    cd ~/BioinfoSpace
    pwd # /home/your_username/BioinfoSpace と表示されることを確認
    ```
    この `BioinfoSpace` にプロジェクトファイルが格納されます。

3.  **`bioinfo-launcher` を実行**
    コマンドパレット (`Ctrl+Shift+P`) を開き、`bioinfo-launcher` と入力。
    候補から `bioinfo-launcher: Start bioinfo-launcher` を選択して実行します。

4.  **作業環境ディレクトリの設定**
    プロンプトが表示されたら、先ほど作成したディレクトリのフルパス (例: `/home/your_username/BioinfoSpace`) を入力します。
    (ターミナルの `pwd` コマンドの出力が参考になります。)

5.  **GitHub Personal Access Token (PAT) の設定**
    次にPATの入力を求められるので、準備しておいたPATを貼り付けます。

6.  **環境構築の開始**
    設定後、`bioinfo-launcher` がDockerを使って解析環境の構築を開始します。
    *   **Dockerについて:** この拡張機能はDockerが必要です。UbuntuにDocker Engineが未インストールの場合は、拡張機能がインストールを試みます。うまくいかない場合は、[Docker Engine for Ubuntu 公式サイト](https://docs.docker.com/engine/install/ubuntu/) を参考に手動でインストールしてください。
    *   環境構築には数分～数十分かかることがあります（特に初回）。

7.  **コンテナ内での作業開始**
    構築が完了すると、多くの場合、新しいCursorウィンドウが自動で開きます。これがDockerコンテナ内の開発環境です。
    ウィンドウ左下のステータスバーがコンテナ環境を示していれば成功です。
    `BioinfoSpace` ディレクトリ内に `.env` ファイルや `cache`, `container` などのディレクトリが自動生成されているはずです。

これで、`hu-bioinfo-workshop.bioinfo-launcher` のセットアップと実行は完了です。統一された解析環境で作業を始めましょう！
{{< /details >}}
{{% /tab %}}

{{% tab "Windows 10" %}}








## Windows 10 でのセットアップ

{{< details title = "1. WSL2 (Windows Subsystem for Linux 2) の導入" >}}

Windows上でLinux環境を動かすための「WSL2」をセットアップします。これにより、Linuxベースのツールが使えるようになります。
**注意:** WSL2の利用には、Windows 10 バージョン 2004 (ビルド 19041) 以降が必要です。お使いのWindowsのバージョンを確認してください。

1.  **PowerShellを管理者として起動**
    スタートメニューで "PowerShell" を検索し、「管理者として実行」を選びます。

2.  **WSL機能の有効化**
    PowerShellに以下のコマンドを入力し、Enterキーを押します。
    ```powershell
    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
    ```

3.  **仮想マシンプラットフォーム機能の有効化**
    続けて、以下のコマンドも入力します。
    ```powershell
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
    ```

4.  **PCの再起動**
    上記2つのコマンド実行後、**必ずPCを再起動してください。**

5.  **Linuxカーネルの更新**
    再起動後、[WSL2 Linuxカーネル更新プログラムパッケージ(x64マシン用)](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi) をダウンロードし、インストールします。
    (ARM64マシンの場合は[こちら](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_arm64.msi)を使用)

6.  **WSL2をデフォルトバージョンに設定**
    再度PowerShellを管理者として起動し、以下を入力します。
    ```powershell
    wsl --set-default-version 2
    ```
    このコマンドがエラーになる場合は、先にUbuntuをインストールしてから再度試すか、WSLのバージョンが古い可能性があります。

7.  **Ubuntu 24.04 LTSのインストール**
    PowerShellで以下のコマンドを実行し、Ubuntuをインストールします。
    ```powershell
    wsl --install -d Ubuntu-24.04
    ```
    ネットワーク環境によっては少し時間がかかります。
    **注:** この `wsl --install` コマンドが利用できない古いバージョンのWindows 10 (バージョン 2004より前、または必要な更新プログラムが適用されていない場合) では、Microsoft Storeを開き、「Ubuntu 24.04 LTS」を検索してインストールしてください。その後、PowerShellで `wsl --set-default-version 2` (またはディストリビューション名を指定して `wsl --set-version Ubuntu-24.04 2`) を実行してWSL2を使用するように設定してください。

8.  **Ubuntuの初期設定**
    インストール完了後、PowerShellで `wsl` または `wsl -d Ubuntu-24.04` と入力してUbuntuを起動します。
    初回起動時、ユーザー名とパスワードの設定を求められます。画面の指示に従い設定し、**忘れないようにメモしておきましょう。**

これでWSL2とUbuntuの基本的なセットアップは完了です。
{{< /details >}}

{{< details title = "2. Cursorエディタのインストール" >}}
AIコーディング支援機能付きのエディタ「Cursor」をインストールします。

1.  **公式サイトからダウンロード**
    [Cursor公式サイト](https://cursor.sh/) にアクセスし、Windows版インストーラーをダウンロードします。

2.  **インストール実行**
    ダウンロードした `.exe` ファイルをダブルクリックして実行し、画面の指示に従ってインストールを進めます。
    「このアプリがデバイスに変更を加えることを許可しますか？」と表示されたら「はい」を選択してください。

インストールが完了すると、Cursorが使えるようになります。
{{< /details >}}

{{< details title = "3. Cursorの日本語化" >}}
Cursorの表示を日本語にします。

1.  **Cursorを起動**

2.  **拡張機能ビューを開く**
    左側のアクティビティバーにある四角いアイコン（または `Ctrl+Shift+X`）で開きます。

3.  **日本語言語パックを検索・インストール**
    検索ボックスに `Japanese Language Pack` と入力し、`Japanese Language Pack for Visual Studio Code` (Microsoft提供) を見つけて `Install` ボタンをクリックします。

4.  **Cursorを再起動**
    インストール後、右下に再起動を促す通知が出たら「再起動」をクリックします。出ない場合は手動でCursorを再起動してください。

これでCursorが日本語表示になります。
{{< /details >}}

{{< details title = "4. 便利な拡張機能のインストール" >}}
開発効率を上げるため、Cursorに以下の拡張機能をインストールします。WSL連携やコンテナ開発に役立ちます。
（日本語化の時と同じように、拡張機能ビューからIDで検索してインストールします）

#### 4.1 Remote - WSL (`ms-vscode-remote.remote-wsl`)

Windows Subsystem for Linux (WSL) と連携するための拡張機能です。CursorからWSL内のプロジェクトを直接扱えるようになります。

1.  拡張機能ビューで `ms-vscode-remote.remote-wsl` を検索。
2.  `Remote - WSL` (Microsoft提供) をインストール。

#### 4.2 Remote - Containers (`ms-vscode-remote.remote-containers`)

Dockerコンテナを開発環境として使うための拡張機能です。プロジェクトごとに独立した環境を構築できます。

1.  拡張機能ビューで `ms-vscode-remote.remote-containers` を検索。
2.  `Remote - Containers` (Microsoft提供) をインストール。

これらの拡張機能で、より高度な開発が可能になります。
{{< /details >}}

{{< details title = "5. CursorからWSL (Ubuntu) へ接続" >}}
インストールした `Remote - WSL` 拡張機能を使って、CursorからWSL上のUbuntu環境に接続します。

1.  **Cursorを起動**

2.  **コマンドパレットを開く**
    `Ctrl+Shift+P` を押します。

3.  **WSL接続コマンドを実行**
    コマンドパレットに `WSL` と入力し、候補から `WSL: WSL に接続` (または `"WSL: Connect to WSL"`) を選択して実行します。

4.  **接続するディストリビューションを選択**
    利用可能なWSLディストリビューションのリストから `Ubuntu-24.04` を選択します。

5.  **WSL (Ubuntu) 環境へ接続完了**
    新しいCursorウィンドウが開き、WSL上のUbuntu環境に接続されます。
    ウィンドウ左下に `WSL: Ubuntu-24.04` のように表示されていれば成功です。
    （初回接続時は、必要なコンポーネントのインストールに少し時間がかかることがあります。）

6.  **接続確認 (任意)**
    Cursor内で新しいターミナルを開きます (`Ctrl+@` またはメニューから `ターミナル` > `新しいターミナル`)。
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

これで、CursorからWSL上のUbuntu環境で開発作業を開始できます。
{{< /details >}}

{{< details title = "6. HU-Bioinfo Workshop Launcher 拡張機能のセットアップと実行" >}}
`hu-bioinfo-workshop.bioinfo-launcher` は、本ワークショップ用のR/Python解析環境を簡単に構築できるCursor拡張機能です。Dockerコンテナ技術を使用します。

**重要:** この拡張機能のインストールと実行は、**CursorがWSL (Ubuntu) に接続された状態**で行います。

#### 6.1 拡張機能のインストール (WSL接続環境で)

1.  **WSL (Ubuntu) 接続の確認**
    Cursorウィンドウ左下のステータスバーが `WSL: Ubuntu-24.04` (または類似の表示) になっていることを確認します。

2.  **拡張機能ビューを開く** (`Ctrl+Shift+X`)

3.  **`hu-bioinfo-workshop.bioinfo-launcher` を検索・インストール**
    検索ボックスに `hu-bioinfo-workshop.bioinfo-launcher` と入力し、`HU-Bioinfo Workshop Launcher` をインストールします。

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

#### 6.3 解析作業用ディレクトリの作成とLauncherの実行

1.  **WSL (Ubuntu) ターミナルを開く** (Cursor内で `Ctrl+@`)
    プロンプトが `your_username@your_hostname:~$` であることを確認。

2.  **作業用親ディレクトリを作成**
    ターミナルで以下のコマンドを実行し、ホームディレクトリに `BioinfoSpace` ディレクトリを作成して移動します。
    ```bash
    mkdir ~/BioinfoSpace
    cd ~/BioinfoSpace
    pwd # /home/your_username/BioinfoSpace と表示されることを確認
    ```
    この `BioinfoSpace` にプロジェクトファイルが格納されます。

3.  **`bioinfo-launcher` を実行**
    コマンドパレット (`Ctrl+Shift+P`) を開き、`bioinfo-launcher` と入力。
    候補から `bioinfo-launcher: Start bioinfo-launcher` を選択して実行します。

4.  **作業環境ディレクトリの設定**
    プロンプトが表示されたら、先ほど作成したディレクトリのフルパス (例: `/home/your_username/BioinfoSpace`) を入力します。
    (ターミナルの `pwd` コマンドの出力が参考になります。)

5.  **GitHub Personal Access Token (PAT) の設定**
    次にPATの入力を求められるので、準備しておいたPATを貼り付けます。

6.  **環境構築の開始**
    設定後、`bioinfo-launcher` がDockerを使って解析環境の構築を開始します。
    *   **Dockerについて:** この拡張機能はDockerが必要です。UbuntuにDocker Engineが未インストールの場合は、拡張機能がインストールを試みます。うまくいかない場合は、[Docker Engine for Ubuntu 公式サイト](https://docs.docker.com/engine/install/ubuntu/) を参考に手動でインストールしてください。
    *   環境構築には数分～数十分かかることがあります（特に初回）。

7.  **コンテナ内での作業開始**
    構築が完了すると、多くの場合、新しいCursorウィンドウが自動で開きます。これがDockerコンテナ内の開発環境です。
    ウィンドウ左下のステータスバーがコンテナ環境を示していれば成功です。
    `BioinfoSpace` ディレクトリ内に `.env` ファイルや `cache`, `container` などのディレクトリが自動生成されているはずです。

これで、`hu-bioinfo-workshop.bioinfo-launcher` のセットアップと実行は完了です。統一された解析環境で作業を始めましょう！
{{< /details >}}
{{% /tab %}}

{{% tab "MacOS" %}}

## MacOS でのセットアップ

MacOSはUnixベースのOSであり、WSLのようなLinux仮想環境は通常不要です。

{{< details title = "1. Cursorエディタのインストール" >}}

AIコーディング支援機能付きのエディタ「Cursor」をインストールします。

1.  **公式サイトからダウンロード**
    [Cursor公式サイト](https://cursor.sh/) にアクセスし、macOS版インストーラー (通常は `.dmg` ファイル) をダウンロードします。

2.  **インストール実行**
    ダウンロードした `.dmg` ファイルを開き、Cursorアプリアイコンを `Applications` フォルダにドラッグ＆ドロップします。
    初回起動時にセキュリティ確認が表示された場合は、「開く」を選択してください。

インストールが完了すると、Cursorが使えるようになります。
{{< /details >}}

{{< details title = "2. Cursorの日本語化" >}}

Cursorの表示を日本語にします。

1.  **Cursorを起動**

2.  **拡張機能ビューを開く**
    左側のアクティビティバーにある四角いアイコン（または `Cmd+Shift+X`）で開きます。

3.  **日本語言語パックを検索・インストール**
    検索ボックスに `Japanese Language Pack` と入力し、`Japanese Language Pack for Visual Studio Code` (Microsoft提供) を見つけて `Install` ボタンをクリックします。

4.  **Cursorを再起動**
    インストール後、右下に再起動を促す通知が出たら「再起動」をクリックします。出ない場合は手動でCursorを再起動してください。

これでCursorが日本語表示になります。
{{< /details >}}

{{< details title = "3. 便利な拡張機能のインストール" >}}

開発効率を上げるため、Cursorに以下の拡張機能をインストールします。コンテナ開発に役立ちます。
（日本語化の時と同じように、拡張機能ビューからIDで検索してインストールします）

#### 3.1 Remote - Containers (`ms-vscode-remote.remote-containers`)

Dockerコンテナを開発環境として使うための拡張機能です。プロジェクトごとに独立した環境を構築できます。

1.  拡張機能ビューで `ms-vscode-remote.remote-containers` を検索。
2.  `Remote - Containers` (Microsoft提供) をインストール。

この拡張機能で、より高度な開発が可能になります。
{{< /details >}}

{{< details title = "4. Docker Desktop のインストール" >}}

コンテナ技術を利用するために、Docker Desktop for Mac をインストールします。

1.  **公式サイトからダウンロード**
    [Docker Desktop for Mac 公式サイト](https://docs.docker.com/desktop/install/mac-install/) にアクセスし、お使いのMacのチップ (Intel または Apple Silicon) に対応したインストーラーをダウンロードします。

2.  **インストール実行**
    ダウンロードした `.dmg` ファイルを開き、Dockerアプリアイコンを `Applications` フォルダにドラッグ＆ドロップします。

3.  **Docker Desktop の起動と設定**
    `Applications` フォルダからDockerを起動します。初回起動時には設定や利用規約への同意が求められる場合があります。
    画面右上のメニューバーにDockerアイコンが表示されれば起動成功です。

これでDockerが利用可能になります。
{{< /details >}}

{{< details title = "5. HU-Bioinfo Workshop Launcher 拡張機能のセットアップと実行" >}}

`hu-bioinfo-workshop.bioinfo-launcher` は、本ワークショップ用のR/Python解析環境を簡単に構築できるCursor拡張機能です。Dockerコンテナ技術を使用します。

**重要:** この拡張機能のインストールと実行は、**Docker Desktopが起動している状態**で行います。

#### 5.1 拡張機能のインストール

1.  **Docker Desktopの起動確認**
    メニューバーにDockerアイコンが表示され、起動していることを確認します。

2.  **拡張機能ビューを開く** (Cursor内で `Cmd+Shift+X`)

3.  **`hu-bioinfo-workshop.bioinfo-launcher` を検索・インストール**
    検索ボックスに `hu-bioinfo-workshop.bioinfo-launcher` と入力し、`HU-Bioinfo Workshop Launcher` をインストールします。

#### 5.2 GitHub Personal Access Token (PAT) の準備

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

#### 5.3 解析作業用ディレクトリの作成とLauncherの実行

1.  **ターミナルを開く** (Cursor内で `Ctrl+@` またはMacのターミナルアプリ)
    プロンプトが `your_username@your_hostname:~$` (または類似の形式) であることを確認。

2.  **作業用親ディレクトリを作成**
    ターミナルで以下のコマンドを実行し、ホームディレクトリに `BioinfoSpace` ディレクトリを作成して移動します。
    ```bash
    mkdir ~/BioinfoSpace
    cd ~/BioinfoSpace
    pwd # /Users/your_username/BioinfoSpace (または /home/your_username/BioinfoSpace) と表示されることを確認
    ```
    この `BioinfoSpace` にプロジェクトファイルが格納されます。

3.  **`bioinfo-launcher` を実行**
    コマンドパレット (Cursor内で `Cmd+Shift+P`) を開き、`bioinfo-launcher` と入力。
    候補から `bioinfo-launcher: Start bioinfo-launcher` を選択して実行します。

4.  **作業環境ディレクトリの設定**
    プロンプトが表示されたら、先ほど作成したディレクトリのフルパス (例: `/Users/your_username/BioinfoSpace`) を入力します。
    (ターミナルの `pwd` コマンドの出力が参考になります。)

5.  **GitHub Personal Access Token (PAT) の設定**
    次にPATの入力を求められるので、準備しておいたPATを貼り付けます。

6.  **環境構築の開始**
    設定後、`bioinfo-launcher` がDockerを使って解析環境の構築を開始します。
    *   環境構築には数分～数十分かかることがあります（特に初回）。

7.  **コンテナ内での作業開始**
    構築が完了すると、多くの場合、新しいCursorウィンドウが自動で開きます。これがDockerコンテナ内の開発環境です。
    ウィンドウ左下のステータスバーがコンテナ環境を示していれば成功です。
    `BioinfoSpace` ディレクトリ内に `.env` ファイルや `cache`, `container` などのディレクトリが自動生成されているはずです。

これで、`hu-bioinfo-workshop.bioinfo-launcher` のセットアップと実行は完了です。統一された解析環境で作業を始めましょう！
{{< /details >}}
{{% /tab %}}
{{% /tabs %}}

