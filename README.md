##概要
内容
 次元格子上の効率的な SVP 解法である Lagrange 基底簡約アルゴリ
ズムの一般次元への拡張である LLL 基底簡約アルゴリズムを紹介

内容
 次元格子上の効率的な SVP 解法である Lagrange 基底簡約アルゴリ
ズムの一般次元への拡張である LLL 基底簡約アルゴリズムを紹介

##数学的定義
定義：LLL 簡約基底
n 次元格子 L の基底 {b1, . . . , bn} の GSO ベクトルを b∗1, . . . , b∗n，GSO 
係数を µi,j (1 ≤ j < i ≤ n) とする．簡約パラメータ 1/4 < δ < 1 に関する次の 2 つの条件を
満たすとき，その基底は δ に関して LLL 簡約されている（Lenstra-Lenstra-Lov´asz (LLL)-reduced）という：
(i) 基底 {b1, . . . , bn} はサイズ簡約されている
(ii) 任意の 2 ≤ k ≤ n に対して，∥b∗k∥2 ≥ (δ − µ2k,k−1)∥b∗k−1∥2.2 を満たす．

（Lov´asz 条件）ただし，各 1 ≤ ℓ ≤ n に対して



## Dependency
sagemath9.1


## Usage
ソースに記載


## Authors
私.

## References
暗号のための格子の理論
Daniele Micciancio (著), Shafi Goldwasser (著)
