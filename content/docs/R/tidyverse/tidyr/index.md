---
title: "tidyr"
weight: 1
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# tidyr

tidyrは字の如く，与えられたデータフレーム(tibble)をtidy dataに整形するのに抜群の効果を発揮します．

## wide dataとlong data

実際にtidyではないデータを用いて，tidyrの力を実感してみましょう．
以下の表は，ある生徒のテスト成績表です．

| name    | math | science | 
| ------- | ---- | ------- | 
| Alice   | 80   | 90      | 
| Bob     | 70   | 85      | 
| Charlie | 30   | 40      | 

この表は一見[ tidy data ]({{% ref "/docs/R/tidyverse/#tidy-data" %}})に見えますが，原則に合致していない部分があります．

{{% details title="tidy dataの3原則" open=false %}}
1. 一つの変数は一つの列に（Each variable must have its own column.）
2. 一つの観測は一つの行に（Each observation must have its own row.）
3. 一つの値は一つのセルに（Each value must have its own cell.）
{{% /details %}}

この表では，2と3の原則は保たれていますが，1番目の **一つの変数は一つの列に（Each variable must have its own column.）** の原則が厳密には満たされていません．
`math`の列を例に見てみると，この列は数学の**得点**の列です．単位をつけるとするならば点数の列であり，数学そのものの属性を表している列ではありません．仮に`math_score`のように列名をつけるとすると，単位と列名が合致している分，辻褄は取れています．しかし`数学の得点`と言う列は，教科と得点と言う二つの要素が混ざった列になってしまっています．

そこで以下の表のように，教科と得点と言う混ざった2要素を，`subject`と`score`と言う2つの列にそれぞれ分割します．こうすることで，元の表の内容を保ったまま，tidy dataに整形することができました．

| name    | subject | score |
| ------- | ------- | ----- |
| Alice   | math    | 80    |
| Alice   | science | 90    |
| Bob     | math    | 70    |
| Bob     | science | 85    |
| Charlie | math    | 30    |
| Charlie | science | 40    |


このように,二つに分割される前の横長の表データのことを`wide data`，分割後の縦長のデータを`long data (narrow data)`といいます．それぞれがtidy dataの原則と合致するわけではないのですが，long dataの方がtidyであることが多い印象です．また，人間の目からしてみるとwide dataの方が直感的なことが多いです．そのため，wide dataの方が人間にとっては見やすいけどlong data形式にしないと解析がしづらいといった問題が度々生じ得ます．この2つの形式を，互いに変換しあって，上手く乗りこなしていく必要があります．

### pivot_longer

それでは実際にwide dataをlong dataに書き直してみましょう．それには`tidyr`の`pivot_longer`関数を使います．
```R
library(tidyverse)

wide_data <- tibble(
    name = c("Alice", "Bob", "Charlie"),
    math = c(80, 70, 30),
    science = c(90, 85, 40),
)

long_data <- wide_data |>
                pivot_longer(
                    cols = c(math, science),
                    names_to = "subject",
                    values_to = "score"
                )

print(long_data)
```
{{< youtube wUVrsxITP4g >}}

このコードの動きは上の動画を見ていただくとわかりやすいと思います．cols引数で指定された元の表の列名が，新しい表では`subject`列に収まっています．また,元の表の値は新しく指定された`score`列に収まっています．このように，cols引数で，分割したい元の列名を指定し，それぞれ収納先の列名を新しくつけてやる必要があります．

{{% details title="colsの指定方法" open=false %}}
上記の例の他にもう一つ，colsの指定方法があります．それは，分割したくない列名に`!`をつけてやることです．`!`はの[論理型]({{% ref "/docs/R/R_basic_grammar/#logicalbool-論理型" %}})説明にもあったように否定を表す記号です．こうすることによって以下の例では`name`列以外を分割するように指定できます．

long dataに変換する際は，多くの列名を指定することが多く，その際，分割しない列名を指定した方が記述がスッキリと分かりやすくなることが多いです．
```R
df |>
    pivot_longer(
        cols = !name,
        names_to = "subject",
        values_to = "score"
    )
```
{{% /details %}}

### pivot_wider

pivot_widerはpivot_longerの逆の操作をします．つまりはlong dataからwide dataへの変換を行うことができます．

```R
library(tidyverse)

wide_data <- long_data |>
                pivot_wider(
                    names_from = subject,
                    values_from = score
                )

print(wide_data)
```
{{< youtube 4rzjJvyReFI >}}

`names_from`で指定した列から新たな列名を作り，`values_from`で指定した列から新たな値を取り出します．

このようにpivot_longerとpivot_widerは，それぞれ逆の操作となっています．

