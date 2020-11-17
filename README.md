

## Contents
次元格子上の効率的な SVP 解法である Lagrange 基底簡約アルゴリ
ズムの一般次元への拡張である LLL 基底簡約アルゴリズムを紹介

## Definition
定義：LLL 簡約基底

n 次元格子 L の基底 {b1, . . . , bn} の GSO ベクトルを 
b∗1, . . . , b∗n，
グラムシュミット係数を 
µi,j (1 ≤ j < i ≤ n) とする．
簡約パラメータ 1/4 < δ < 1 に関する次の 2 つの条件を
満たすとき，その基底は δ に関して LLL 簡約されている
（Lenstra-Lenstra-Lov´asz (LLL)-reduced）という：

(i) 基底 {b1, . . . , bn} はサイズ簡約されている.

(ii) 任意の 2 ≤ k ≤ n に対して，
((b∗k).norm())^2 ≥ (δ − µ^2_(k,k−1) ) * ((b∗_k−1).norm())^2 を満たす．

## theorem
定理： LLL 簡約基底の性質
n 次元格子 L の基底 {b1, . . . , bn} が簡約パラメータ 1
4 < δ < 1 に関して
LLL 簡約されているとする．このとき，
α = 4 / (4δ − 1)
とおくと，次が成り立つ：
(1) 任意の 1 ≤ j ≤ i ≤ n に対して，∥b∗j∥^2 ≤ α^(i−j)∥b∗i∥^2 が成り立つ

(2) 不等式 ∥b1∥ ≤ α^((n−1)/4) * vol(L)1n が成り立つ

(3) 任意の 1 ≤ i ≤ n に対して，∥bi∥ ≤ α^((n−1)/2) *λi(L) が成り立つ
(4) 不等式 Yni=1∥bi∥ ≤ α^(n(n−1)/4) *vol(L) が成り立つ

## Dependency
 SageMath version 9.1, Release Date: 2020-05-20                     │
│ Using Python 3.7.3.


## Usage
ソースに記載
簡約したい行列を変えることで簡約可能.
Deep_LLL_methodを使うとSageMathでも
40次元の行列は
簡約化可能です。

## Authors
私.

## References
暗号のための格子の理論
Daniele Micciancio (著), Shafi Goldwasser (著)
