---
title: "Tidyverse"
weight: 4
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# Tidyverse


## Tidyverseってなんぞや
Tidyverseは，データ分析と可視化を効率的に行うために設計されたR言語のパッケージ群です．
このパッケージの思想上の特徴として**Tidy Data**の原則のもと，一貫性のある操作体系により，データ処理の生産性を保つことを特徴としています．

とはいえ含まれる内容や思想も膨大かつ難解なので，この章では，よくする使い方をメインに書いていきます．
より詳しく知りたい方は，それぞれ詳細のページを見てください．(現在製作中です)

{{% details title="Tidyverseのパッケージ群" open=false %}}

| パッケージ名 | 役割                           | 主な機能                                                   | 
| ------------ | ------------------------------ | ---------------------------------------------------------- | 
| ggplot2      | データ可視化                   | グラフィックスの文法に基づいて、グラフを宣言的に作成できる | 
| dplyr        | データ操作                     | フィルタリング，要約，並べ替え，結合といったデータ操作     | 
| tidyr        | データ整理                     | **Tidy data**作成                                          | 
| readr        | データ読み込み                 | 矩形データ（csv、tsv、fwfなど）の高速読み込み              | 
| purrr        | 関数型プログラミング           | ベクトルやリストの操作（map関数など）                      | 
| tibble       | データフレームの現代的な再構築 | モダンなデータフレーム構造                                 | 
| stringr      | 文字列操作                     | パターンマッチングや文字列の処理                           | 
| forcats      | 因子型操作                     | 因子データのレベルの順序や値の変更など                     | 
| lubridate    | 日付・時刻データ操作           | 日付や時間の解析と処理                                     | 

{{% /details %}}


## Tidy Data

基本的にデータはExcelのように表形式となったデータを扱います．この時，ただ漫然と表を作成するのではなく，一定のルールの下に作成，加工していくことで効率や安全性が増します．具体的なルールとは以下の三つのことです.

1. 一つの変数は一つの列に（Each variable must have its own column.）
2. 一つの観測は一つの行に（Each observation must have its own row.）
3. 一つの値は一つのセルに（Each value must have its own cell.）

上記3つを満たした表形式のデータを **Tidy Data(整列データ)** といいます．

### 一つの変数は一つの列に（Each variable must have its own column.）
この原則における変数とは，データの特徴を表すものです．人を対象にしたデータだとすると名前，身長，体重，年齢などが該当します．考え方として列では，単位が揃うようなデータにするようにします．例えば，一つの列でcmやkgが同時に登場しないようにします．

{{% tabs "1st_principle" %}}
{{% tab "良い例" %}}
| 名前  | 年齢 | 身長 (cm) | 体重 (kg) |
|------|----|----------|----------|
| 太郎  | 25  | 170      | 65       |
| 花子  | 22  | 160      | 50       |
{{% /tab %}}
{{% tab "悪い例" %}}
| 名前  | 属性   | 値   |
|------|------|------|
| 太郎  | 年齢  | 25   |
| 太郎  | 身長  | 170  |
| 太郎  | 体重  | 65   |
| 花子  | 年齢  | 22   |
| 花子  | 身長  | 160  |
| 花子  | 体重  | 50   |
{{% /tab %}}
{{% /tabs %}}

### 一つの観測は一つの行に（Each observation must have its own row.）
観測とは、一つのデータのまとまりのことを指します．人を対象にしたデータだとすると名前がAlice,身長160cm，体重50kg，年齢20歳という一人のデータのことを指します．これらをまとめて一つのデータとして，行単位で追加していきます．

{{% tabs "2st_principle" %}}
{{% tab "良い例" %}}
| 名前  | 年齢 | 身長 (cm) | 体重 (kg) |
|------|----|----------|----------|
| 太郎  | 25  | 170      | 65       |
| 花子  | 22  | 160      | 50       |
{{% /tab %}}
{{% tab "悪い例" %}}
| 名前  | データ       |
|------|------------|
| 太郎  | 25歳        |
| 太郎  | 170cm      |
| 太郎  | 65kg       |
| 花子  | 22歳        |
| 花子  | 160cm      |
| 花子  | 50kg       |
{{% /tab %}}
{{% /tabs %}}

### 一つの値は一つのセルに（Each value must have its own cell.）
各セルに入る値は，一つだけにするという原則です．名前がAliceというデータについて，身長が2つや3つもあるのは不自然ですよね．
ただし，Tidyverseで提供されているtibbleというデータフレーム(Excelの表の強化版みたいなもの)は一つのセルの中に,配列を入れることができたり，さらにtibbleを格納することもできるので，必ずしもこの原則が達成されているとは限らないので注意が必要です．

{{% tabs "3st_principle" %}}
{{% tab "良い例" %}}
| 名前  | 好きな食べ物 |
|------|----------|
| 太郎  | 寿司     |
| 太郎  | ラーメン  |
| 花子  | カレー   |
| 花子  | パスタ   |
{{% /tab %}}
{{% tab "悪い例" %}}
| 名前  | 好きな食べ物        |
|------|----------------|
| 太郎  | 寿司, ラーメン    |
| 花子  | カレー, パスタ    |
{{% /tab %}}
{{% /tabs %}}

以上の内容を意識すると，統一的な記述方法でデータを解析，整形することができるので，tidy dataを意識してデータ解析すると良いでしょう．

{{% details title="余談" open=false %}}

Tidy Dataが解析しやすいデータだとすれば，それ以外のデータは解析しずらい，めんどくさいデータです．例え人にとって見やすかったり，扱いやすかったりしてもTidyでないならコンピュータには扱いづらいです．そのため，データを提供する側と解析する側では，往々にしてこの部分ですれ違いがあったりするので，何かデータを扱うときはTidyであるかどうかを気にしてみるといいかもしれません．

{{% /details %}}


## パイプライン演算子 `|>`
Tidyverseを使用していく上でパイプライン演算子`|>` を活用することで，処理を簡潔かつ解釈しやすく記述することができます．


### 関数
まず，パイプライン演算子を使う前に，プログラミングにおける関数について説明します．プログラミングにおける関数は，ある処理のまとまりを集めて，いつでも呼び出せるようにしたものです．これは数学における，ある集合からある集合への対応を示す関数とは違い，サブルーチンといった方がより適切です．そのため，処理をまとめただけの関数もあれば，ある入力に対してある出力を返す関数があったりします．

ここで入力に対応するものを引数，出力に対応するものを返り値(戻り値)といいます．

**純粋関数**という概念があります．純粋関数とは所謂，数学における関数と等価なもので，ある入力(引数)に対して，いつも同じ出力(返り値)を返す関数です．この純粋関数を組み合わせることでコーディングをしていくのがTidyverseとパイプラインを組み合わせた方法です．

具体的には純粋関数は，同じ引数に対して，特定の返り値を返すので，ある純粋関数の返り値をそのまま，別の純粋関数の引数に入力することで，処理を記述していきます．これまでのプログラミングといえばやりたいことを逐一丁寧に記述していく手続き型プログラミングといったものを行ってきました．しかし，それとは考え方，アプローチが違ってくるので戸惑う点があるかもしれません．とはいえ，基本的にプログラミングとは何をするかというと，与えられた入力(データ)を加工していく処理，作業に他なりません．その大前提を意識しておけば，どちらもやることがたいして変わりません．以下で具体例を見てみましょう．


{{% details title="余談" open=false %}}

純粋関数とは，ある引数に対して特定の返り値を返す—つまりは引数が同じ時，常に同じ返り値を返すことと，副作用が発生しない関数のことを指します．副作用とは引数以外の入力，返り値以外の出力を伴うことです．具体的には関数の外の変数を参照，変更することであったり，I/O画面やストレージと相互作用をすることを指します．

例えばprint関数であったり(文字を外部に出力するため), random関数(同じ引数を入力しても出力がランダム)といった関数が非純粋関数となります．

しかし，Rにおいてそのような関数は大多数を占めますし，本当の意味で純粋関数は少ないのですが，今回は簡単のために上記ような解説となりました．

{{% /details %}}


### パイプライン演算子を使ってみよう

初めに手続きプログラミング的にコードを書いてみます．以下のコードは，データフレームdfを宣言し，value列の値が20より大きいものを選別し，groupごとに集計し，groupごとの平均値を算出したものを出力するコードです．
```R
library(dplyr)

df <- data.frame(
    group = c("A", "B", "A", "B", "A"),
    value = c(10, 20, 30, 40, 50)
)

temp1 <- filter(df, value > 20)  # 20より大きい値をフィルタ
temp2 <- group_by(temp1, group)  # "group" でグループ化
result <- summarise(temp2, mean_value = mean(value))  # "group"ごとの平均値を計算

print(result)
```

このようなコードでは，変数temp1, temp2のように演算結果を一時的に保存するためだけの変数を宣言したりと無駄が多いです．無駄が多くなると，論理を追いづらくなったり，コードを手直しするにも大変です．

次の例では一時的な変数を作成しないようにコードを書き直してみます．
```R
result <- summarise(
    group_by(filter(df, value > 20), group),
    mean_value = mean(value)
)

print(result)
```

この書き方では，余計な変数はありませんが，関数の中に関数を書き込んだ形になってしまい，非常に可読性が悪いです．私たちは右から左，外から内に文字を読むのに，処理の流れは左から右に，内から外の順番になってしまっています．
ここでパイプライン演算子の出番です．Rのパイプライン演算子`|>`は**ある関数の返り値をそのまま，次の関数の第一引数**にしてしまいます．
```R
library(dplyr)

df |>
    filter(value > 20) |>
    group_by(group) |>
    summarise(mean_value = mean(value)) |>
    print()
```
こうすることによって，処理の流れと記述の流れが同じになりました．純粋関数を用いると同じ引数に対し，同じ返り値が帰ってくるので，このような記述が可能になります．

{{% details title="Advanced" open=false %}}

やや発展的な内容になりますが，パイプライン演算子にはRで定義されたnative pipelineと呼ばれる`|>`とTidyverseのmagrittrパッケージによって提供されている`%>%`の2種類があります．

これらの使い分けとしては，返り値を第一引数ではない場所に用いたい時に，やり方が異なります．`%>%`では第二引数などに値を渡す時`.`を使うことで可能になります(place holderといいます)．`|>`では`_`を使用しますが，`x |> f(100, y=_)`のように名前付き引数の時だけ使用できます．
```R
# place holderは`.`。引数名を指定しなくても良い。
c("apple", "pineapple", "banana") %>% grepl("apple", .)
#> TRUE  TRUE FALSE

# place holderは`_`。引数名を指定しなければいけない。
c("apple", "pineapple", "banana") |> grepl("apple", x = _)
#> TRUE  TRUE FALSE

# 引数名を指定しないとエラーになる。
c("apple", "pineapple", "banana") |> grepl("apple", _)
#> Error: pipe placeholder can only be used as a named argument
```

また複数箇所にplace holderを使用したい場合,`%>%`しか選択肢はありません．
とはいえ，基本的には`|>`の方が高速に動作し，何よりパッケージではなく，R言語そのもので定義されているので，`|>`を使用するのを個人的には推奨します．

{{% /details %}}

## Tidyverseで遊んでみよう
お待たせしました．ここからは実際にcsvデータを加工して，Tidyverseの凄さを思い知りましょう．Tidyverseには多くのパッケージが含まれますが，ここではよく使う，dplyrを中心に使っていきます．

ここで加工していくデータは,dplyr内に練習用として定義されているstarwarsのデータを用います．スターウォーズのキャラクターについての情報が詰まれたデータフレームです．
glimpse関数を使用することで，データフレームの概観を知ることができます．(そのままstarwarsと入力するだけでもいいです)

```R
starwars |> glimpse()

#> Rows: 87
#> Columns: 14
#> $ name       <chr> "Luke Skywalker", "C-3PO", "R2-D2", "Darth Vader", "L…
#> $ height     <int> 172, 167, 96, 202, 150, 178, 165, 97, 183, 182, 188, …
#> $ mass       <dbl> 77.0, 75.0, 32.0, 136.0, 49.0, 120.0, 75.0, 32.0, 84.…
#> $ hair_color <chr> "blond", NA, NA, "none", "brown", "brown, grey", "bro…
#> $ skin_color <chr> "fair", "gold", "white, blue", "white", "light", "lig…
#> $ eye_color  <chr> "blue", "yellow", "red", "yellow", "brown", "blue", "…
#> $ birth_year <dbl> 19.0, 112.0, 33.0, 41.9, 19.0, 52.0, 47.0, NA, 24.0, …
#> $ sex        <chr> "male", "none", "none", "male", "female", "male", "fe…
#> $ gender     <chr> "masculine", "masculine", "masculine", "masculine", "…
#> $ homeworld  <chr> "Tatooine", "Tatooine", "Naboo", "Tatooine", "Alderaa…
#> $ species    <chr> "Human", "Droid", "Droid", "Human", "Human", "Human",…
#> $ films      <list> <"A New Hope", "The Empire Strikes Back", "Return of…
#> $ vehicles   <list> <"Snowspeeder", "Imperial Speeder Bike">, <>, <>, <>…
#> $ starships  <list> <"X-wing", "Imperial shuttle">, <>, <>, "TIE Advance…
```

glimpse()を適用すると,行と列がひっくり返ってしまうので注意が必要です．しかし，これでstarwarsが87×14のtibbleで，それぞれの列がどのようなデータを保持しているがわかりました．例えばname列はcharacter(文字列)を保持している列ですね．

{{< button relref="/docs/tutorials/tidyverse/#問題1" >}}問題1へ{{< /button >}}

### 列の操作
まず初めに列に対する操作です．パイプ演算子`|>`を使って`select`関数を適用することで，列を抽出できます．`rename`関数では列名を新しく付け替えることができます．

```R
starwars |>
    select(name, height, homeworld) #name, height, homeworld列だけを抽出
#> # A tibble: 87 × 3
#>    name               height homeworld
#>    <chr>               <int> <chr>    
#>  1 Luke Skywalker        172 Tatooine 
#>  2 C-3PO                 167 Tatooine 
#>  3 R2-D2                  96 Naboo    
#>  4 Darth Vader           202 Tatooine 

starwars |>
    select(!starships) # starships以外の列を抽出
#> # A tibble: 87 × 13
#>    name     height  mass hair_color skin_color eye_color birth_year sex   gender
#>    <chr>     <int> <dbl> <chr>      <chr>      <chr>          <dbl> <chr> <chr> 
#>  1 Luke Sk…    172    77 blond      fair       blue            19   male  mascu…
#>  2 C-3PO       167    75 NA         gold       yellow         112   none  mascu…
#>  3 R2-D2        96    32 NA         white, bl… red             33   none  mascu…
#>  4 Darth V…    202   136 none       white      yellow          41.9 male  mascu…

starwars |>
    select(ends_with("color")) #末尾がcolorで終わる列(hair_color, skin_color, eye_color)を抽出
#> # A tibble: 87 × 3
#>    hair_color    skin_color  eye_color
#>    <chr>         <chr>       <chr>    
#>  1 blond         fair        blue     
#>  2 NA            gold        yellow   
#>  3 NA            white, blue red      
#>  4 none          white       yellow   

starwars |>
    select(where(is.character)) #値が文字列(chr)である列を抽出
#> # A tibble: 87 × 8
#>    name           hair_color skin_color eye_color sex   gender homeworld species
#>    <chr>          <chr>      <chr>      <chr>     <chr> <chr>  <chr>     <chr>  
#>  1 Luke Skywalker blond      fair       blue      male  mascu… Tatooine  Human  
#>  2 C-3PO          NA         gold       yellow    none  mascu… Tatooine  Droid  
#>  3 R2-D2          NA         white, bl… red       none  mascu… Naboo     Droid  
#>  4 Darth Vader    none       white      yellow    male  mascu… Tatooine  Human 

starwars |>
    select(homeworld) |> # homeworld列を抽出
    rename(home_planet = homeworld) # 列名homeworldをhome_planetに変更
#> # A tibble: 87 × 1
#>    home_planet
#>    <chr>      
#>  1 Tatooine   
#>  2 Tatooine   
#>  3 Naboo      
#>  4 Tatooine   
#>  5 Alderaan   
```
{{< button relref="/docs/tutorials/tidyverse/#問題2" >}}問題2へ{{< /button >}}
{{< button relref="/docs/tutorials/tidyverse/#問題4" >}}問題4へ{{< /button >}}

### 行の操作
`slice_*`系の関数は行の操作に使います．`slice`関数は指定した行番号の列を，`slice_head`,`slice_tail`はそれぞれ最初のn行と最後のn行を，`slice_min`, `slice_max`は指定して列の値が小さい順ないし大きい順にn行抽出します．

`filter`関数は何かの条件と合致する(あるいは合致しない)行のみ抽出する関数です． `arrange`関数はソート関数です．

```R
starwars |>
    select(name, homeworld) |> # name, homeworld列だけ抽出
    filter(homeworld == "Tatooine") |> # homeworldがTatooineの行のみを抽出
    slice_head(n = 4) # そのうち先頭4行だけ抽出
#> # A tibble: 4 × 2
#>   name               homeworld
#>   <chr>              <chr>    
#> 1 Luke Skywalker     Tatooine 
#> 2 C-3PO              Tatooine 
#> 3 Darth Vader        Tatooine 
#> 4 Owen Lars          Tatooine 

starwars |>
    select(name, height, homeworld) |> # name,height,homeworld列だけ抽出
    filter(homeworld != "Tatooine") |> # homeworldがTatooineでない行を抽出
    slice_max(height, n = 3) # そのうちheightが大きい順に3行抽出
#> # A tibble: 3 × 3
#>   name        height homeworld
#>   <chr>        <int> <chr>    
#> 1 Yarael Poof    264 Quermia  
#> 2 Tarfful        234 Kashyyyk 
#> 3 Lama Su        229 Kamino   

starwars |>
    select(name, hair_color, height, mass) |> # name, height, sex, homeworld列を抽出
    filter(!is.na(mass), height >= 100, hair_color == "blond") |> # massが欠損値がなくて，尚且つheightが100以上かつ，hair_colorがblondの行を抽出
    arrange(mass) # massが昇順になるようにソート
    #> # A tibble: 2 × 4
    #>   name             hair_color height  mass
#>   <chr>            <chr>       <int> <dbl>
#> 1 Luke Skywalker   blond         172    77
#> 2 Anakin Skywalker blond         188    84
```

{{< button relref="/docs/tutorials/tidyverse/#問題3" >}}問題3へ{{< /button >}}

### データの要約
`summarise`関数を使うと変数の平均値や標準偏差などの記述統計量(要約統計量)を計算できます. `group_by`関数と組み合わせることで値ごとの記述統計量を出すことができます．

```R
starwars |>
    group_by(homeworld) |> # homeworldごとに集計
    summarise(
        height_mean = mean(height, na.rm = TRUE), #欠損値を除外してheightの平均を計算
        mass_mean = mean(mass, na.rm = TRUE), #欠損値を除外してmassの平均を計算
    )
#> # A tibble: 49 × 3
#>    homeworld      height_mean mass_mean
#>    <chr>                <dbl>     <dbl>
#>  1 Alderaan              176.      64  
#>  2 Aleen Minor            79       15  
#>  3 Bespin                175       79  
#>  4 Bestine IV            180      110  
```
{{< button relref="/docs/tutorials/tidyverse/#問題7" >}}問題7へ{{< /button >}}

### データの拡張
`mutate`関数はtibble内の変数を用いて計算を行い，その結果を新しい列として追加する関数です．

```R
starwars |>
    select(name, height) |>
    mutate(height_M = height / 100) # heightを100でわりm換算したものをheight_in_mとした
#> # A tibble: 87 × 3
#>    name               height       height_M
#>    <chr>               <int>       <dbl>
#>  1 Luke Skywalker        172        1.72
#>  2 C-3PO                 167        1.67
#>  3 R2-D2                  96        0.96
#>  4 Darth Vader           202        2.02

starwars |>
    select(name, height, mass) |>
    mutate(height_M = height / 100, BMI = mass / height_M / height_M) |> # heightを100でわりm換算したものをheight_in_mとした
    mutate(is_obesity = if_else(BMI >= 25, true="obesity", false="not_obesity")) # if_else関数を用いると,条件を指定して値を指定できる
#> # A tibble: 87 × 6
#>    name               height  mass height_M   BMI is_obesity 
#>    <chr>               <int> <dbl>    <dbl> <dbl> <chr>      
#>  1 Luke Skywalker        172    77     1.72  26.0 obesity    
#>  2 C-3PO                 167    75     1.67  26.9 obesity    
#>  3 R2-D2                  96    32     0.96  34.7 obesity    
#>  4 Darth Vader           202   136     2.02  33.3 obesity    
#>  5 Leia Organa           150    49     1.5   21.8 not_obesity
```
{{< button relref="/docs/tutorials/tidyverse/#問題5" >}}問題5へ{{< /button >}}

## 練習問題
以下は練習問題です．これまで学習してきたこと＋αな内容ですので，適宜，chatGPTに聞いてみたり調べたりしながら解いてみてください．

### 問題1

`starwars` データセットの構造を確認せよ。
{{% hint info %}}
💡 ヒント: glimpse(), head(), dim() などを使う

{{< button relref="/docs/tutorials/tidyverse/#tidyverseで遊んでみよう" >}}参考{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}


```R
starwars |>
    glimpse()

starwars |>
    head()

starwars |>
    dim()
```

{{% /details %}}
### 問題2

`name`, `height`, `mass`, `homeworld` の4列のみを選択し，キャラクターを `mass` の降順で並べ替えよ．

{{% hint info %}}
💡 ヒント: `arrange`関数と，`desc`関数について調べてみよう

{{< button relref="/docs/tutorials/tidyverse/#列の操作" >}}列の操作{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}

```R
starwars |>
    select(name, heightm, mass, homeworld) |>
    arrange(desc(mass))
```
{{% /details %}}

### 問題3

`starwars`データセットのうち，`gender` が "male" かつ `mass` (体重) が 80kg 以上のキャラクターを抽出せよ

{{% hint info %}}
💡 ヒント: `filter`関数について調べてみよう

{{< button relref="/docs/tutorials/tidyverse/#行の操作" >}}行の操作{{< /button >}}

{{% /hint %}}

{{% details title="Answer" open=false %}}

```R
starwars |>
    filter(gender == "male", mass >= 80)
```

{{% /details %}}


### 問題4
`name`と列名にアンダーバー `_`を含む列のみを選択せよ

{{% hint info %}}
💡 ヒント: `select`, `contains`関数の使い方を調べてみよう

{{< button relref="/docs/tutorials/tidyverse/#列の操作" >}}列の操作{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}
```R
starwars |>
    select(name, contains("_"))
```
{{% /details %}}

### 問題5

`birth_year` が 100 より大きい場合は "Old"、それ以外は "Young" とする `age_category` という新しい列を作成せよ
{{% hint info %}}
💡 ヒント: `if_else`関数の使い方を調べてみよう

{{< button relref="/docs/tutorials/tidyverse/#データの拡張" >}}データの拡張{{< /button >}}
{{< button href="https://dplyr.tidyverse.org/reference/if_else.html" >}}if_else関数{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}
```R
starwars |>
    mutate(age_category = if_else(birth_year > 100, "Old", "Young"))
```
{{% /details %}}

### 問題6

以下のルールに基づいて `weight_category` という新しい列を作成せよ。

* `mass` が 100 kg 以上なら "Heavy"
* `mass` が 50 以上 100 未満なら "Medium"
* `mass` が 50 未満なら "Light"
* `mass` が NA の場合は "Unknown"

{{% hint info %}}
💡 ヒント: `case_when`関数の使い方を調べてみよう

{{< button href="https://dplyr.tidyverse.org/reference/case_when.html" >}}case_when関数{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}
```R
starwars |>
    mutate(weight_category = case_when(
        is.na(mass) ~ "Unknown",
        mass >= 100 ~ "Heavy",
        mass >= 50 ~ "Medium",
        TRUE ~ "Lignt"
    ))
```
{{% /details %}}

### 問題7
species ごとに、
* キャラクター数 (n())
* height の平均 (mean(height, na.rm = TRUE))
* mass の平均 (mean(mass, na.rm = TRUE))
を求めよ。

ただし、キャラクター数が 3 人未満の種族は除外せよ。

{{% hint info %}}
💡 ヒント: `group_by`, `summarise`, `filter`関数の使い方を調べてみよう

{{< button relref="/docs/tutorials/tidyverse/#データの要約" >}}データの要約{{< /button >}}
{{% /hint %}}

{{% details title="Answer" open=false %}}
```R
starwars |>
    group_by(species) |>
    summarise(
        N = n(),
        height_mean = mean(height, na.rm = TRUE),
        mass_mean = mean(mass, na.rm = TRUE)
    ) |>
    filter(N > 3)
```
{{% /details %}}


### 問題8
NA を含む数値型の列をすべて選択し、それらの NA を 0 に置き換えよ。

{{% hint info %}}
💡 ヒント: `mutate`, `across`, `where`関数の使い方を調べてみよう. Rにおける**無名関数**について調べてみよう
{{% /hint %}}

{{% details title="Answer" open=false %}}
```R
starwars |>
    mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))
```
{{% /details %}}
