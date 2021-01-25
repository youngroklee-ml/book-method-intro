# Discrete choice model {#discrete-choice-model}


```r
library(extraDistr)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(gt)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(tidyr)
library(rlang)
library(purrr)
```

```
## 
## Attaching package: 'purrr'
```

```
## The following objects are masked from 'package:rlang':
## 
##     %@%, as_function, flatten, flatten_chr, flatten_dbl, flatten_int,
##     flatten_lgl, flatten_raw, invoke, list_along, modify, prepend,
##     splice
```

```
## The following object is masked from 'package:extraDistr':
## 
##     rdunif
```

```r
library(mlogit)
```

```
## Loading required package: dfidx
```

```
## 
## Attaching package: 'dfidx'
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

우리의 삶에는 수많은 선택의 순간들이 있다. 크게는 집을 구매할 때 어떤 집을 선택할지부터, 소소하게는 짬뽕과 짜장면 중 어떤 음식을 먹을지 선택해야한다. 물론, 가진 자원이 많다면 여러 채의 집을 구매할 수도 있고, 짬뽕과 짜장면 둘 다 주문에서 먹을 수도 있다. 하지만, 제한된 자원에서는 선택지 중 한 가지를 선택하면 나머지는 포기해야 하는 상황에 부딪힌다. 이 때, 사람들이 어떤 선택을 하는지를 모델링하는 방법이 있다. 본 장에서는 **discrete choice model**이라 불리는 방법에 대해 알아보자.


## Data

아래와 같이 데이터를 정의해보자.

### Choice

- $j = 1, \ldots, J$: 선택할 수 있는 $J$개의 서로 다른 대안(alternative). 예를 들어, 짬뽕($j = 1$)과 짜장면($j = 2$).
- $i = 1, \ldots, N$: 선택하여야 하는 상황(choice situation). 예를 들어, 100명의 손님이 있다면, 각각의 손님($i = 1, \ldots, 100$)은 짬뽕과 짜장면 중 한 가지 음식을 선택하여야 한다.
- $y_{ij}$: $i$번째 선택 상황에서 $j$번째 대안이 선택되었다면 1, 아니라면 0.

이 때, $y_{ij}$는 모든 $i$에 대해 다음과 같은 제약을 만족한다.

\[
\sum_{j = 1}^{J} y_{ij} = 1\\
y_{ij} \in \{0, 1\}
\]

### Covariate

선택 변수 $y_{ij}$의 값은 선택 상황에 대한 설명 변수들과 상관관계가 있을 수 있다. 예를 들어, 손님의 나이가 어릴수록 짜장면을 선택할 확률이 높다거나, 날씨가 비가 오면 짬뽕을 선택할 확률이 높다거나, 음식 가격이 비쌀수록 해당 음식을 선택할 확률이 줄어들 수 있다. 이러한 변수의 벡터를 아래 $\mathbf{x}_{ij}$와 $\mathbf{z}_{i}$ 두 가지로 나누어 정의하자.

- $\mathbf{x}_{ij}$: $i$번째 선택 상황에서 $j$번째 대안에 대한 설명변수 관측치 (음식 가격 등)
- $\mathbf{z}_{i}$: $i$번째 선택 상황에 대한 설명변수 관측치 (손님 나이, 날씨 등)

위에서 $\mathbf{x}_{ij}$는 이후에 다시 두 가지 벡터로 분리된다. 


### Availability

$i$번째 선택 상황에서 $j$번째 대안이 선택 가능하지 않는 경우들이 있다. 예를 들어, 재료수급 문제로 짬뽕이 품절이었거나, 손님이 가진 현금이 짬뽕을 주문하기 부족한 경우, 짬뽕은 선택 가능한 대안에서 제외되게 된다. 이를 표현하기 위해 변수 $a_{ij}$를 아래와 같이 정의하자.

- $a_{ij}$: $i$번째 선택 상황에서 $j$번째 대안이 선택 가능했다면 1, 아니라면 0.

이 때, $a_{ij}$ 와 $y_{ij}$간에는 아래와 같은 관계가 성립한다.

\[
a_{ij} \geq y_{ij}
\]

따라서,

\[
\sum_{j = 1}^{J} a_{ij} \geq \sum_{j = 1}^{J} y_{ij} = 1
\]

즉, 적어도 하나의 대안은 선택 가능하여야 한다.

$i$번째 선택 상황에서 선택 가능한 대안의 집합을 $A_i$라 하자.

\[
A_i = \left\{j : a_{ij} = 1 \right\}
\]

이 때, $y_{ij}$값은 아래와 같은 조건을 만족한다. 

\[
\sum_{j \in A_i} y_{ij} = 1\\
y_{ij} \in \{0, 1\}, \; \forall j \in A_i\\
y_{ij} = 0, \; \forall j \notin A_i\\
\]


### "No choice" alternative

앞서 availability 측면에서 고려했던, 짬뽕이 선택 가능하지 않고 짜장면만 선택 가능한 경우, $\sum_{j = 1}^{J} y_{ij} = 1$에 따라 무조건 짜장면을 선택하게 된다. 하지만, 짬뽕이 선택 가능하지 않았을 때, 짜장면을 먹는 대신 다른 음식점을 찾아가는 경우를 예상해볼 수 있다. 이 때, 손님은 선택 가능한 대안 중 어떠한 선택도 하지 않은 것이 되어 $\sum_{j = 1}^{J} y_{ij} = 1$을 만족하지 못한다.

이 때, 해당 선택 상황을 분석에서 제외하는 대신, 모든 선택 상황에서 $J$개의 대안 중 어떠한 대안도 선택하지 않는 $J + 1$째의 대안을 고려할 수 있다. 이에 대한 구체적은 예는 이후 다시 다루도록 하자.


### 데이터 예


```r
tribble(
  ~손님ID, ~손님나이, ~음식, ~가격, ~선택여부,
  1, 40, "짬뽕", 8000, 1,
  1, 40, "짜장면", 7000, 0,
  2, 10, "짬뽕", 5000, 0,
  2, 10, "짜장면", 5000, 1,
  3, 70, "짜장면", 6500, 1
) %>%
  knitr::kable()
```



| 손님ID| 손님나이|음식   | 가격| 선택여부|
|------:|--------:|:------|----:|--------:|
|      1|       40|짬뽕   | 8000|        1|
|      1|       40|짜장면 | 7000|        0|
|      2|       10|짬뽕   | 5000|        0|
|      2|       10|짜장면 | 5000|        1|
|      3|       70|짜장면 | 6500|        1|



## Utility

$i$번째 선택 상황에서 $j$번째 대안에 대한 효용성의 값을 $u_{ij} \in \mathbb{R}$이라 하자. 이 때, 합리적인 선택은 효용성이 가장 큰 대안을 선택하는 것이다. 본 장에서 다룰 기본적인 discrete choice model들은 이러한 합리적인 선택을 가정한다.

\[
y_{ij} = \begin{cases}
1 & j = \arg\max_{j \in A_i} u_{ij}\\
0 & \text{otherwise}
\end{cases}
\]

이 때, $u_{ij}$는 관측되지 않는 값이다. 이를 두 개의 부분으로 아래와 같이 나눈다.

\[
u_{ij} = v_{ij} + \varepsilon_{ij}
\]

이 때, $v_{ij}$는 관측된 설명변수 $\mathbf{x}_{ij}$와 $\mathbf{z}_i$에 대한 함수이며, $\varepsilon_{ij}$는 설명변수로는 설명되지 않는 부분이다.

\[
v_{ij} = f(\mathbf{x}_{ij}, \mathbf{z}_i)
\]

이 때, 함수 $f()$의 형태나 $\varepsilon_{ij}$에 대한 가정에 따라 다양한 형태의 discrete choice model이 존재한다. 이후 본 장에서는 가장 간단한 두 가지 모형들(multinomial logit model, nested logit model)을 살펴보기로 하자.


### Observable utility (Representative utility)

$v_{ij}$는 observable utility 혹은 representative utility라 한다. $v_{ij}$를 설명변수에 대한 affine 함수로 아래와 같이 정의하자.

\[
v_{ij} = \alpha_j + \boldsymbol{\beta}^{\top} \mathbf{x}_{ij} + \boldsymbol{\gamma}_j^{\top} \mathbf{z}_i + \boldsymbol{\delta}_j \mathbf{x}_{ij}
\]

이 때, $\boldsymbol{\beta}$는 모든 대안 $j$에 동일하게 적용되는 계수벡터이며, $\boldsymbol{\gamma}_j$와 $\boldsymbol{\delta}_j$는 대안 $j$에 따라 다른 값을 지니는 계수벡터이다. 변수 벡터 $\mathbf{x}_{ij}$내의 각 변수는 $\boldsymbol{\beta}$와 $\boldsymbol{\delta}_j$ 중 한 가지의 회귀계수에만 해당한다고 가정하고, 설명변수 벡터 $\mathbf{w}_{ij}$를 추가로 정의하자.

\[
v_{ij} = \alpha_j + \boldsymbol{\beta}^{\top} \mathbf{x}_{ij} + \boldsymbol{\gamma}_j^{\top} \mathbf{z}_i + \boldsymbol{\delta}_j \mathbf{w}_{ij}
\]

여기에서, 각각의 변수에 대한 예는 아래와 같이 생각해볼 수 있다.

- $\mathbf{x}_{ij}$: (예: 음식가격) 음식 가격은 각 메뉴별로 다르나, 손님이 느끼는 효용성에 가격이 미치는 영향은 단위 가격당 동일한 경우, 음식 가격에 적용되는 계수 $\beta$의 값은 대안 $j$와 상관없이 동일할 수 있다.
- $\mathbf{z}_{i}$: (예: 나이) 나이가 어릴수록 짜장면을 선호하고 짬뽕을 덜 선호하는 경향이 커지는 경우, 나이에 적용되는 계수 $\gamma$의 값은 짬뽕($j = 1$)에 대한 계수값보다 짜장면($j = 2$)에 대한 계수값이 크다고 볼 수 있다 ($\gamma_1 < \gamma_2$).
- $\mathbf{w}_{ij}$: (예: 음식가격) 만약 음식에 느끼는 효용성에 가격이 미치는 영향이 메뉴별로 다를 경우, 음식 가격은 $\mathbf{x}_{ij}$가 아닌 $\mathbf{w}_{ij}$에 해당하는 변수라 할 수 있다. 짜장면 가격이 500원 증가할 때 느끼는 효용성의 감소가 짬뽕 가격이 500원 증가할 때 느끼는 효용성의 감소보다 크다면, 즉 효용성 가격 민감도가 짜장면의 경우가 더 크다면, $\delta_2 < \delta_1 < 0$라 볼 수 있다 ($|\delta_2| > |\delta_1|$).

위의 세 가지 변수 종류 중 본 장에서는 두 가지 변수 $\mathbf{x}_{ij}$와 $\mathbf{z}_{i}$를 중심으로 살펴보자. 아래와 같이 보다 단순한 형태의 $v_{ij}$를 고려하도록 하자.

\[
v_{ij} = \alpha_j + \boldsymbol{\beta}^{\top} \mathbf{x}_{ij} + \boldsymbol{\gamma}_j^{\top} \mathbf{z}_i
\]

이 때, 주어진 데이터로부터 회귀계수들의 추정치 $\hat{\alpha}_j, \hat{\boldsymbol{\beta}}, \hat{\boldsymbol{\gamma}}_j$을 구하는 것이 discrete choice model 추정이다.


### Unobservable utility

$\varepsilon_{ij}$는 효용성 중 covariate으로 설명되지 않는 부분이다. 같은 나이의 두 손님이 동일한 가게에서 동일한 시각에 각자 식사를 하더라도, 한 명은 짬뽕을 선택하고 다른 한 명은 짜장면을 선택할 수 있다. 이는 설명변수(나이, 가격 등)로 설명되지 않는 각 손님의 개인적인 취향에 기인한 것으로 생각할 수 있다. 

임의의 $i$번째 선택 상황에서, $j \in A_i$에 대해 $v_{ij}$의 참값은 알고 있고 $y_{ij}$값은 관측되지 않았다 할 때, $y_{ij}$에 대한 기대값은 아래와 같이 표현할 수 있다.

\begin{eqnarray*}
E[y_{ij}] &=& \prod_{k \in A_i \\ j} P(u_{ij} > u_{ik})\\
&=& \prod_{k \in A_i \\ j} P(v_{ij} + \varepsilon_{ij} > v_{ik} + \varepsilon_{ik})\\
&=& \prod_{k \in A_i \\ j} P(\varepsilon_{ij} - \varepsilon_{ik} > v_{ik} - v_{ij})
\end{eqnarray*}

여기서, $y_{ij}$의 기대값을 $p_{ij}$라 하면 ($p_{ij} = E[y_{ij}]$), $p_{ij}$는 아래의 식들을 만족한다.

\[
\sum_{j \in A_i} p_{ij} = 1\\
p_{ij} \in [0, 1]
\]

$y_{ij}$ ($j = 1, \ldots, J$)는 확률 $p_{ij}$ ($j = 1, \ldots, J$)로 정의된 multinomial distribution으로부터 얻어진 관측치라 하자.

\[
(y_{i1}, \ldots, y_{iJ}) \sim multinom(p_{i1}, \ldots, p_{iJ})
\]

위에서, $\varepsilon_{ij}$의 분포에 대한 가정에 따라 $p_{ij}$의 값이 다르게 추정된다. 따라서 discrete choice model 추정에서는 observable utility에 대한 model specification 뿐만 아니라 unobservable utility에 대한 분포 가정 또한 중요하다. 다음 장에서 살펴볼 multinomial logit model과 nested logit model은 $\varepsilon_{ij}$의 분포에 대해 서로 다른 가정을 지닌다.



## Multinomial logit model

Discrete choice model을 얘기할 때 가장 기본으로 다루는 모형이다. 

### Utility

$$
u_{ij} = v_{ij} + \varepsilon_{ij}
$$

에서, unobservable utility $\varepsilon_{ij}$가 standard Gumbel distribution으로부터 얻어진다고 가정한다.

$$
\varepsilon_{ij} \overset{i.i.d.}{\sim} Gumbel(0, 1)
$$

이 때, standard Gumbel distribution은 아래와 같다.

$$
F(\varepsilon_{ij}) = \exp\left(-\exp\left(- \varepsilon_{ij}\right) \right)
$$


### 데이터 모델 (예)

네 가지 메뉴가 있는 식당을 고려해보자. 각 선택 상황에 따라 가격과 메뉴 주문 가능 여부는 매번 달라진다고 가정하자.

- 대안
  - $j = 1$: 짬뽕
  - $j = 2$: 짜장면
  - $j = 3$: 삼선짬뽕
  - $j = 4$: 삼선짜장
- 설명변수
  - $x_{ij}$: $i$번째 손님이 $j$번째 음식에 대해 지불해야 할 가격
$$
x_{i1} \sim U(5000, 7000)\\
x_{i2} \sim U(5000, x_{i1})\\
x_{i3} \sim x_{i1} + U(2000, 4000)\\
x_{i4} \sim x_{i2} + (x_{i3} - x_{i1}) + U(-500, 0)
$$
  - $z_{i}$: $i$번째 손님의 나이
$$
z_i \sim U(5, 95)
$$
- Observable utility
\begin{eqnarray*}
v_{i1} &=& 0 - 0.001 * x_{i1} + 0.1 * z_{i}\\
v_{i2} &=& 0 - 0.001 * x_{i2} - 0.2 * z_{i}\\
v_{i3} &=& 3 - 0.001 * x_{i3} + 0.1 * z_{i}\\
v_{i4} &=& 3 - 0.001 * x_{i4} - 0.2 * z_{i}
\end{eqnarray*}
- Unobservable utility
  - $\varepsilon_{ij} \overset{i.i.d.}{\sim} Gumbel(0, 1)$






### Independent from irrelevant alternatives (IIA)

서로 다른 두 개의 대안 $j$와 $k$간의 상대적인 선택확률($p_{ij} / p_{ik}$)은 또 다른 대안 $l \notin \{j, k\}$을 선택할 확률 $p_{il}$에 영향을 받지 않는다.




## Nested logit model

