---
title: "2. RNA-seqとは?"
weight: 2
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---
# RNA-seqとは？

## 1. 遺伝子の情報が形になるまで

私たちの体の設計図であるDNAには、様々なタンパク質を作るための情報が書かれています。しかし、DNAの情報が直接タンパク質になるわけではありません。

**セントラルドグマ**という生命の基本的な原則があります。これは、遺伝子情報が **DNA → RNA → タンパク質** という一方向の流れで伝えられる、という考え方です。

{{< mermaid >}}
graph LR
    A[DNA] -- Replication --> A
    A -- Transcription --> B[RNA]
    B -- Translation --> C[Protein]
{{< /mermaid >}}

1.  **転写 (Transcription):** DNAの必要な部分の情報が **mRNA (メッセンジャーRNA)** という分子にコピーされます。
2.  **翻訳 (Translation):** mRNAの情報をもとに、タンパク質が合成されます。

このタンパク質が、私たちの体の機能や見た目（表現型）を作り出す実働部隊です。つまり、**どの遺伝子からどれくらいの量のmRNAが作られるか**（遺伝子発現）が、細胞や個体の状態を決定する上で非常に重要なのです。

{{< hint info >}}
**補足:**
実際には、RNAにはmRNA以外にも様々な種類があり、それぞれ重要な役割を担っています。また、セントラルドグマには例外（逆転写など）も存在しますが、まずはこの基本的な流れを理解することが重要です。
{{< /hint >}}

## 2. RNA-seqとは？

**RNA-seq (RNAシーケンシング)** は、細胞の中に存在するRNA、特に **mRNA** を網羅的に読み取り、その種類と量を調べるための強力な技術です。

これにより、特定の時点や条件下で、**どの遺伝子が活発に働いているか（発現しているか）** を大規模に知ることができます。




## 3. RNA-seq実験の基本的な流れ

bulk RNA-seqは、一般的に以下のステップで進められます。

### **1. RNA抽出 (RNA Extraction):** 
- 解析したい細胞や組織からRNAを抽出します。DNAやタンパク質などの他の分子を取り除き、純粋なRNAを得ることが重要です。

{{< figure src="RNA-seq1-RNA_extraction.drawio.svg">}}

### **2. ライブラリ調製 (Library Preparation):** 
- 抽出したRNAをシーケンサーで読み取れる形に加工します。
    *   mRNAの選択(PolyA-tailing)またはrRNAの除去
    *   RNAを断片化
    *   cDNA（相補的DNA）への逆転写
    *   アダプター配列の付加（シーケンサーが認識するための目印）
    *   PCR増幅

{{< figure src="RNA-seq2-Library_preparation.drawio.svg">}}

### **3. シーケンシング (Sequencing):** 
- 調製したライブラリを次世代シーケンサー（NGS）で読み取り、大量の短い塩基配列データ（リード）を得ます。

{{< figure src="RNA-seq3-Sequencing.drawio.svg">}}

### **4. マッピング(Mapping):** 
- 得られたリードをリファレンスゲノムやトランスクリプトームにマッピングする。

{{< figure src="RNA-seq4-Mapping.drawio.svg">}}

### **5. 定量(Quantification):** 
- マッピングしたリードをもとに、遺伝子ごとにリード数をカウントすることで、遺伝子の発現量を定量する。

| 遺伝子   | Sample A | Sample B | Sample C |
| :------- | -----------: | -----------: | -----------: |
| 遺伝子1  |          150 |          200 |          180 |
| 遺伝子2  |           10 |            5 |           12 |
| 遺伝子3  |          800 |          750 |          820 |
| ...      |          ... |          ... |          ... |





## 4. bulk RNA-seq と single-cell RNA-seq

RNA-seqには、大きく分けて2つのアプローチがあります。

*   **bulk RNA-seq:** 複数の細胞の集団からRNAをまとめて抽出し、平均的な発現パターンを解析します。組織全体の傾向を見るのに適しています。
*   **single-cell RNA-seq (scRNA-seq):** 細胞を1つずつ分離し、それぞれの細胞のRNAを解析します。細胞ごとの違いや、不均一な細胞集団の内部構造を詳細に調べることができます。

もう少し詳しくデータの性質を見てみましょう。

**bulk RNA-seq** では、測定結果は組織や細胞集団全体の **平均的な遺伝子発現** を表します。通常、1人の患者さんや1つの実験条件あたり、**1つのサンプル**（RNA抽出物）からデータが得られます。

{{< mermaid >}}
graph LR
    A("組織/細胞集団 (Bulk)") --> B("まとめてRNA抽出") --> C("平均的な発現データ<br>(1サンプル)")
{{< /mermaid >}}

**Bulk RNA-seq のデータの例:**

| 遺伝子   | 症例A (サンプルA) | 症例B (サンプルB) | 症例C (サンプルC) |
| :------- | -----------: | -----------: | -----------: |
| 遺伝子1  |          150 |          200 |          180 |
| 遺伝子2  |           10 |            5 |           12 |
| 遺伝子3  |          800 |          750 |          820 |
| ...      |          ... |          ... |          ... |

（行: 遺伝子, 列: 各症例から得られたサンプル）

**scRNA-seq** では、**個々の細胞ごと** の遺伝子発現プロファイルが得られます。つまり、1人の患者さんや1つのサンプルから、**数百〜数万個の細胞（データ点）** を得ることができ、細胞間のばらつきや、異なる細胞タイプの存在を直接的に解析できます。

{{< mermaid >}}
graph LR
    D("組織/細胞集団 (sc)") --> E("1細胞ずつ分離") --> G("1細胞ごとに<br>RNA抽出/解析") --> H("1細胞ごとの発現データ<br>(多数の細胞)")
{{< /mermaid >}}

**Single-cell RNA-seq のデータの例:**

| 遺伝子   | 細胞A_1 | ... | 細胞A_1000 | 細胞B_1 | ... | 細胞B_1500 | 細胞C_1 | ... | 細胞C_2000 |
| :------- | ------: | :-: | ---------: | ------: | :-: | ---------: | ------: | :-: | ---------: |
| 遺伝子1  |       5 | ... |          8 |       7 | ... |         10 |       9 | ... |         11 |
| 遺伝子2  |       0 | ... |          1 |       0 | ... |          0 |       2 | ... |          0 |
| 遺伝子3  |      25 | ... |         22 |      30 | ... |         28 |      15 | ... |         18 |
| ...      |     ... | ... |        ... |     ... | ... |        ... |     ... | ... |        ... |

（行: 遺伝子, 列: 各症例の個々の細胞。例えば、症例Aから1000細胞、症例Bから1500細胞、症例Cから2000細胞など）

このデータの粒度の違いは、後に行うデータ解析のアプローチにも大きく影響します。



## 5. 正規化

### **転写産物長とライブラリサイズによるバイアス**

RNA-seqで得られるデータ（=**Raw Countsデータ**）は、そのメソッドの特性上、**転写産物長によるバイアス** を受けます。

つまり、転写産物が長ければその本数が少なかったとしても、断片の数が多くなるのでリードの数も多くなります。

それを図解したものが以下です。

{{< figure src="RNA-seq5-long_transcripts.drawio.svg">}}
引用: https://mbernste.github.io/posts/rna_seq_basics/

また、そもそもRNAの量が違えば、RNAの断片やその後のcDNAの量も違ってきます（=**ライブラリサイズの違い**）。

これらのバイアスを補正するために、**正規化** という手法が用いられます。
このようなバイアスを補正するために、**正規化** という手法が用いられます。

「すべての転写産物長の長さと各サンプル間のライブラリサイズが同じだったとしたらカウントはどうなるのか」を考えるのが正規化です。

### **TPM(Transcripts Per Million)**

TPMは、**転写産物長およびライブラリサイズによるバイアス** の正規化方法です。

{{% hint info %}}
- 以前はTPMではなくFPKM(Fragments Per Kilobase Million)やRPKM(Reads Per Kilobase Million)という正規化方法もありました。
- しかし、これらは発現量を正しく表せないことがあり現在はTPMが最も一般的に使用されています。
{{% /hint %}}

TPMは簡単に言うと以下の二つの計算を行っています。

-   **遺伝子の長さで割る** 
    -   まず、各遺伝子で得られたリード数を、その遺伝子の長さで割ります（長さあたりのリード数を計算）。これで、長さによる偏りが取り除かれます。
    -   短い遺伝子でも長い遺伝子でも、mRNAの数が同じなら、この「長さあたりのリード数」は同じような値になります。
-   **全体の合計値で割る（全体のバランスをとる）** 
    -   次に、全ての遺伝子について計算した「長さあたりのリード数」の合計値で各遺伝子の値を割り、100万倍します。
    -   サンプル全体での遺伝子発現の総量（長さ補正後）がサンプルによって違っても、比較しやすいようにバランスをとります。
    -   全体の長さ補正リード数の合計を100万として、その中で各遺伝子がどれくらいの割合を占めるか、という形に揃えているイメージです。


{{% details title="TPMの計算方法" %}}
TPMは、以下のステップで計算されます。

{{< katex display >}}
\begin{aligned}
i &\text{ : 特定の遺伝子} \\
G &\text{ : 全ての遺伝子の数} \\
c_i &\text{ : 遺伝子 } i \text{ にマッピングされたリード数} \\
l_i &\text{ : 遺伝子 } i \text{ の長さ（例えば、キロベース単位）} \\
\text{TPM}_i &\text{ : 遺伝子 } i \text{ のTPM値}
\end{aligned}
{{< /katex >}}

1.  **長さによる補正:**

{{< katex display >}}
\huge \frac{c_i}{l_i}
{{< /katex >}}

-   各遺伝子 \( i \) について、リード数 \( c_i \) を遺伝子長 \( l_i \) で割ります。
-   これは、遺伝子の長さあたりのリード数（RPK: Reads Per Kilobase）を計算していることになり、長さによるバイアスを取り除きます。

2.  **全体の合計値の計算:** 

{{< katex display >}}
\huge \sum_{j=1}^G \frac{c_j}{l_j}
{{< /katex >}}

-   全ての遺伝子についての長さ補正済みリード数（ステップ1で計算）を合計します。
-   これは、サンプル全体での、長さによる補正済みのリード数の合計を表します。

3.  **相対的な割合の計算:** 

{{< katex display >}}
\huge \frac{c_i / l_i}{\sum_{j=1}^G c_j / l_j}
{{< /katex >}}

-   ステップ1で得られた各遺伝子 \( i \) の値（長さあたりのリード数）を、ステップ2で計算した全体の合計値で割ります。
-   これで、長さによるバイアスを取り除いた上での、全転写産物の中でのその遺伝子の相対的な割合が得られます。

4.  **100万倍にスケーリング:** 

{{< katex display >}}
\huge \text{TPM}_i = 10^6 \times \frac{c_i / l_i}{\sum_{j=1}^G c_j / l_j}
{{< /katex >}}

-   ステップ3で得られた相対的な割合に \( 10^6 \) （100万）を掛けます。
-   これにより、全体を100万としたときの、その遺伝子の発現量の値（TPM）が得られます。

{{% /details %}}


### **TPMの利点**

TPMを使うことで、以下の利点が得られます。

*   遺伝子の長さによる影響が取り除かれるため、より正確な発現量比較が可能になります。
*   **同じサンプル内**で、異なる遺伝子の発現量（相対的な割合）を比較するのに適しています。
*   **異なるサンプル間**で、同じ遺伝子（例えば遺伝子X）の発現量（相対的な割合）を比較するのに適しています。



### **TPMの欠点**

TPMは異なるサンプル間で同じ遺伝子の相対的な発現量（サンプル全体に占める割合）を比較するのに適していますが、あくまで「相対値」であるため、**絶対的な発現量**（例: 細胞あたりのmRNA数）をサンプル間で比較することはできません。

例えば、ある遺伝子（図では青の遺伝子）の細胞当たりのmRNAの絶対数はサンプル2の方が多かったとしても、TPMを計算するとサンプル1の方がサンプル2より高かったという状況が起こり得ます。これは、サンプル間で全体のRNA量が大きく異なる場合に特に起こりやすい現象です。

{{< figure src="RNA-seq6-TPM_vs_absolute_expression.drawio.svg">}}
引用: [https://mbernste.github.io/posts/rna_seq_basics/](https://mbernste.github.io/posts/rna_seq_basics/)

RNA-seqデータでサンプル間の絶対的な発現量を比較するには、TPMのような相対値の正規化だけでは不十分であり、以下のような別の手法が必要になることがあります。

*   **スパイクインRNA:** 既知の量のRNAをサンプルに添加し、それを基準に絶対量を推定する。
*   **ハウスキーピング遺伝子:** 発現量が安定しているとされる遺伝子を基準にする。
*   **メディアン比正規化:** 多くの遺伝子の発現は変動しないという仮定に基づき、統計的にサンプル間の違いを補正する。

RNA-seqの解析結果を正しく解釈するには、TPMが相対値であるという性質を理解しておくことが重要です。


### **メディアン比正規化(Median Ratio Normalization)**

RNA-seq解析を行う目的の一つは、**二つの条件の間で変動する遺伝子を特定する**ことです。

{{% hint info %}}
-   このような解析を**遺伝子差次発現解析(DEG: Differential Expression Gene)** といいます。
{{% /hint %}}

しかし、TPMは絶対的な発現量をサンプル間で比較することはできません([TPMの欠点](#tpmの欠点))。

そこで、**メディアン比正規化(Median Ratio Normalization)** という手法が用いられます。

メディアン比正規化は、**多くの遺伝子の発現は変動しない**という仮定に基づき、統計的にサンプル間の違いを補正する手法です。

概念を図解すると以下のようになります。

{{< figure src="RNA-seq7-Median_ratio_sample.drawio.svg">}}
引用：[https://mbernste.github.io/posts/median_ratio_norm/](https://mbernste.github.io/posts/median_ratio_norm/)

-   まずすべてのサンプルを使って「Reference Sample」を作成します。

{{< figure src="RNA-seq8-Median_ratio_sample2.drawio.svg">}}
引用：[https://mbernste.github.io/posts/median_ratio_norm/](https://mbernste.github.io/posts/median_ratio_norm/)

-   次に、各サンプルの各遺伝子ごとにReferenceとの倍率変化が最も小さい遺伝子を選択します。

{{< figure src="RNA-seq9-Median_ratio_sample3.drawio.svg">}}
引用：[https://mbernste.github.io/posts/median_ratio_norm/](https://mbernste.github.io/posts/median_ratio_norm/)

-   その遺伝子をベースラインとして他の遺伝子を補正します。

このようにすると、各遺伝子のサンプル間での差はmRNAの絶対量を強く反映したものになります。

{{% details title="メディアン比正規化の計算方法" %}}
メディアン比正規化は以下のステップで計算されます。

{{< katex display >}}
\begin{aligned}
n &\text{ : サンプルの総数} \\
g &\text{ : 遺伝子の総数} \\
c_{i,j} &\text{ : サンプル } i \text{ における遺伝子 } j \text{ からのリードのカウント} \\
m_j &\text{ : 遺伝子 } j \text{ のベースライン発現値} \\
r_{i,j} &\text{ : サンプル } i \text{ における遺伝子 } j \text{ の参照サンプルとの比率} \\
s_i &\text{ : サンプル } i \text{ のサイズファクター}
\end{aligned}
{{< /katex >}}

#### 1.  **参照サンプルの発現値の計算**

{{< katex display >}}
\huge m_j := \left( \prod_{i=1}^n c_{i,j} \right)^{\frac{1}{n}}
{{< /katex >}}

-   各遺伝子のすべてのサンプルにわたるカウントの**幾何平均**を計算し、各遺伝子のベースライン発現を表す「参照サンプル」の発現値とする。
-   算術平均ではなく幾何平均が使用されるのは、外れ値に対してよりロバストであるためです。

#### 2.  **各サンプルの各遺伝子における参照サンプルとの比率の計算**

{{< katex display >}}
\huge r_{i,j} := \frac{c_{i,j}}{m_j}
{{< /katex >}}

-   各サンプルのどの遺伝子が参照サンプルの発現と一致する必要があるかを特定します。
-   各サンプルの各遺伝子について、サンプルにおける遺伝子のカウントと遺伝子のベースライン発現値の比率を計算します。
-   これは、この遺伝子の参照サンプルの発現との間の偏差（より具体的には倍率変化）を意味します。
-   仮に理想的なサンプルについて計算された比率を小さい方から並べて図示すると以下のようになります。

{{< figure src="RNA-seq10-Median_ratio_sample4.drawio.svg">}}
引用：[https://mbernste.github.io/posts/median_ratio_norm/](https://mbernste.github.io/posts/median_ratio_norm/)

-   発現減少（青）や発現増加（赤）の遺伝子も少数ありますが、**ほとんどの遺伝子（灰色）は変化しない**はずです。

#### 3.  **サンプルごとのサイズファクター（中央比率）の計算:**

{{< katex display >}}
\huge s_i := \text{median}(r_{i,1}, r_{i,2}, \ldots, r_{i,g})
{{< /katex >}}

-   2.で明らかになった比率のリストの**中央値**が、そのサンプル全体の技術的なバイアス（ライブラリサイズなど）を示す「**サイズファクター**」として使用できると仮定します。
-   このランキングの中央でベースラインからの逸脱が見られる場合、それは**ライブラリサイズによって引き起こされる**と仮定されます。
-   したがって、中央比率を、リストの中央にある比率がベースラインに近づくように、サンプルのカウントを再スケーリングするために使用できる「サイズファクター」として扱うことができます。

{{< figure src="RNA-seq11-Median_ratio_sample5.drawio.svg">}}
引用：[https://mbernste.github.io/posts/median_ratio_norm/](https://mbernste.github.io/posts/median_ratio_norm/)

{{% hint warning %}}
-   中央比率を計算する際に、すべてのサンプルで発現がゼロでない遺伝子のみを使用します。
-   検出可能であるほど十分に高い発現を持つ遺伝子のみを使用したいからです。
{{% /hint %}}

4.  **正規化されたカウントの計算:**

{{< katex display >}}
\huge \tilde{c}_{i,j} := \frac{c_{i,j}}{s_i}
{{< /katex >}}

-   最後に、計算されたサイズファクター を使用して、各サンプル内のすべてのカウントを再スケーリング（割る）することで正規化します。

{{% /details %}}




#### 参考文献
-   Figure created in the Mind the Graph Platform, available at www.mindthegraph.com
-   Veeranagouda, Y., & Didier, M. (2017). RNA fragmentation and sequencing (RF-Seq): A cost-effective approach for transcriptome analysis. Current Protocols in Molecular Biology. doi: 10.1002/cpmb.xx
-   Bernstein, M. (2021, January 7). RNA-seq: the basics. Retrieved May 11, 2025, from https://mbernste.github.io/posts/rna_seq_basics/
-   Bernstein, M. (2024, March 02). Assessing the utility of data visualizations based on dimensionality reduction. Retrieved May 11, 2025, from https://mbernste.github.io/posts/dim_reduc/