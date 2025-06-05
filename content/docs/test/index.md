---
title: "test"
weight: 99
draft: true
# bookFlatSection: false
# bookToc: true
# bookHidden: false
# bookCollapseSection: false
# bookComments: false
# bookSearchExclude: false
---

# コードとその出力例

{{< code-tabs >}}
--- test
a
--- コード
```javascript
function fibonacci(n) {
  if (n <= 1) return n;
  return fibonacci(n - 1) + fibonacci(n - 2);
}

console.log(fibonacci(10));
```
--- 出力
```
55
```
{{< figure src="test.png" >}}
--- wow1
## メー
{{< /code-tabs >}}

## Python例

{{< code-tabs >}}
Python コード
```python
def factorial(n):
    if n == 0 or n == 1:
        return 1
    return n * factorial(n - 1)

print(factorial(5))
```
---
実行結果
```
120
```
---
解説
階乗を計算する再帰関数です。
- `factorial(5)` = 5 × 4 × 3 × 2 × 1 = 120
- 時間計算量: O(n)
- 空間計算量: O(n)（再帰スタックのため）
{{< /code-tabs >}}