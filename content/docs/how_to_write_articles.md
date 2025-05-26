---
title: "How to Write Articles"
weight: 100
draft: true
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# 記事の投稿方法

## 構成
サイトは[Hugo](https://gohugo.io/)の[Hugo-bookテーマ](https://themes.gohugo.io/themes/hugo-book/)用いて作成しました．
ディレクトリ構成は以下のとおりです．
基本的にいじるのは`content`以下のディレクトリだけです．


```
.
├── archetypes
├── config
    └── _default
        ├── hugo.toml
        └── params.toml
├── content
    ├── docs
    │   ├── tutorials
    │   │   ├── example
    │   │   │   ├── bioinfo_tutorial.webp
    │   │   │   └── index.md
    │   │   ├── test
    │   │   │   └── index.md
    │   │   ├── _index.md
    │   │   └── bioinfo_tutorial.webp
    │   ├── about.md
    │   └── members.md
    └── _index.md
├── data
├── i18n
├── layouts
├── public
├── resources
├── static
├── go.mod
└── go.sum

```

## 記事の作成
まず，[Githubリポジトリ](https://github.com/HU-bioinfo/hu-bioinfo.github.io)の`main`ブランチをローカルにcloneしてください．
自分で新しい記事を作成，修正するときは，`main`ではない別のブランチ(自分の名前とか)を作成し，そこで作業するようにしてください．


今後は基本的に`content/docs/tutorials`以下に，一記事ずつディレクトリを作成していきます．

```
├── content
    ├── docs
    │   ├── tutorials
    │   │   ├── new_article
    │   │   │   ├── images
    │   │   │   └── index.md
```
作成したい記事のディレクトリを作成し，index.mdを作り，編集します．

index.mdにはヘッダーとして以下のものをつけてください．

```markdown=
---
title: "Example"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---
```
このときtitleにはディレクトリ名と同じ名前を基本的にはつけてください．
weightはサイトの左側メニューに表示される時の順番を示しています．上に載せたいものほど小さい値をつけてください．
挿入する画像は，作成したディレクトリ以下に配置してください．

markdownの編集方法については，[Hugo-bookのサイト](https://hugo-book-demo.netlify.app/)のShortcodesや，私が作成した[記事](https://github.com/HU-bioinfo/hu-bioinfo.github.io/blob/main/content/docs/example/index.md)と[サイト](https://hu-bioinfo.github.io/docs/example/)を見比べてみてください．


## 投稿方法
記事の作成が完了したら，`git add && git commit`して，リポジトリにそのまま自分のブランチをpushしてください．

pushした自分のブランチを`main`ブランチに向けてプルリクエストを作成してください．

無事，プルリクが受理され，`main`ブランチにマージされたら，Github Actionsが走り，hugoによりサイトがビルドされます．
[北大バイオインフォのサイト](https://hu-bioinfo.github.io/)に変更がなされているか確認してみてください．


## リンク

* [Githubリポジトリ](https://github.com/HU-bioinfo/hu-bioinfo.github.io)
* [Github pages](https://hu-bioinfo.github.io/)
* [Hugo-book](https://hugo-book-demo.netlify.app/)

