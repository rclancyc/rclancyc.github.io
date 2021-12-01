---
title: "Robust Least Squares"
header:
  overlay_image: /assets/images/swift_current_peak.jpg
  caption: "Swift Current Peak, Galcier NP"
permalink: /projects/robustLS/
---

This section is based on the work published in [Signal Processing](https://www.sciencedirect.com/science/article/abs/pii/S0165168420302541)[^1]. A preprint can be found [here](https://arxiv.org/pdf/2003.12004). 

Due to the limitations of classical methods, we approach the problem through the lens of robust optimization. The primary goal is to recover unknown parameters from a noisy observation and an uncertain linear operator. In particular, our mathematical model for robust least squares is

$$
\min_{\mathbf x} \, \left\{ \max_{ \boldsymbol \Delta  \in \,  \mathcal U }\; \| (\mathbf A+ \boldsymbol \Delta ) \mathbf x - \mathbf y \|^2 \right\}
$$

which we refer to as our robust optimization (RO) problem. The Euclidean norm is denoted by $$\| \cdot \|$$ and $$\mathcal{U}$$ is the uncertainty set from which perturbations in $$\mathbf A$$ are drawn. The above RO formulation is motivated by two situations. In both cases, we let $$\bar{\mathbf A}$$ and $$\bar{\mathbf x}$$ represent the _true and unknown_ data matrix and parameter vector, respectively. We model $$\mathbf y = \bar{ \mathbf A} \bar{ \mathbf x} + \boldsymbol \eta$$, with $$\boldsymbol  \eta$$ i.i.d. Gaussian, and we only have knowledge of $$\mathbf A = \bar{ \mathbf A }- \boldsymbol \Delta$$. We don't know $$\boldsymbol \Delta$$ explicitly but can make inferences based on the problem. In the first situation, we consider a data matrix subject to quantization or round-off error. Suppose the observed matrix $$\mathbf A$$ has elements rounded to the hundredth place. Our uncertainty set can be written as $$\mathcal U= \{ \boldsymbol \Delta  \in \mathbb R^{m \times n } : \| \boldsymbol \Delta \|_{\infty} \le \delta \}$$ with $$\delta = 0.005$$. If $$\mathbf A_{i,j} = 0.540$$, then we know the true $$\bar {\mathbf A}_{i,j} \in (0.535,\, 0.545]$$, hence $$\boldsymbol \Delta _{i,j} \in (-0.005, 0.005]$$. The norm $$\| \cdot \|_{\infty}$$ takes the maximum absolute value of any element in the matrix. 

The second problem considers a data matrix with uncertainty proportional to the magnitude of each entry, i.e. 
$$\mathcal U = \{ \boldsymbol \Delta  \in \mathbb R^{m \times n} : \boldsymbol \Delta_{i,j} \in (-p |\mathbf A_{i,j}|, p | \mathbf A_{i,j}| ] \}$$. Here, $$p$$ denotes a proportionality constant. Data subject to $$\pm 1 \%$$ uncertainty would have $$p = 0.01$$. The two cases cover the effects of finite-precision in fixed and floating point representations, respectively. In both problems, the uncertainty sets are specified element-wise allowing us to decouple along rows.

In this project, we consider box constraints over entries of the data matrix. Our main contribution is the formulation of a robust objective to handle quantization and floating point error in least squares problems and the presentation of methods to solve it. Although the proposed method requires selection of an uncertainty parameter, it is chosen in a natural and theoretically appropriate way based on the observed extent of quantization. This is in contrast to ridge regression and MLE based methods that require involved parameter tuning or *a priori* knowledge of the probability distributions from which model uncertainty is drawn. We anticipate our method to be most effective under moderate to heavy quantization where fidelity loss is greater than $$0.1\%$$.

To solve the problem, it must be tractable. In our paper, we show that the following are equivalent

$$
\min_{\mathbf x} \bigg\{ \  \underset{f(\mathbf x)}{\underbrace{ \,  \max_{ |\boldsymbol \Delta | \le \mathbf D}\; \| (\mathbf A+ \boldsymbol \Delta ) \mathbf x - \mathbf y \|^2 }} \ \bigg\} \\
\ \qquad \qquad   \Updownarrow \\
\min_{\mathbf x} \left\{ \big\| \mathbf A \mathbf x-\mathbf y \big\|^2 + 2 \langle \, | \mathbf A \mathbf x- \mathbf y|, \mathbf D|\mathbf x| \, \rangle + \big\| \, \mathbf D| \mathbf x| \, \big\|^2 \right \}
$$

where $$|\cdot |
$$ is the absolute value acting component-wise on vectors/matrices and $$\mathbf D$$ is a matrix filled with positive entries constraining the elements of $$\boldsymbol  \Delta$$.

By letting $$\mathbf D = \delta \mathbf 1_m^{} \mathbf 1_n^T$$, we recover the fixed-point error uncertainty set and $$f( \mathbf x)$$ simplifies to

$$
f( \mathbf x) = \bigg\{ \|\mathbf A \mathbf x- \mathbf y \|^2 + 2 \delta \|\mathbf x\|_1 \| \mathbf A \mathbf x- \mathbf y \|_1 + m \delta^2 \| \mathbf x \|_1^2 \bigg \}.
$$

A challenge to overcome is that $$f(\mathbf x)$$ is not differentiable and finding elements of the subgradient can be challenging. Another contribution is that we provide an element of the subgradient that is easy to compute for all $$\mathbf x \in \mathbb R^n$$. We do so by recognizing that since $$\| (\mathbf A+ \boldsymbol \Delta ) \mathbf x - \mathbf y \|^2$$ is convex in $$\mathbf x$$ for fixed $$\boldsymbol  \Delta$$ and the maximum of convex functions is also convex, then $$f(\mathbf x)$$ is convex as well. The conditions of [Danskin's theorem](https://en.wikipedia.org/wiki/Danskin's_theorem) are met in our problem which gives an element of the subgradient,

$$
f'(\mathbf x) := 2 (\mathbf A+\boldsymbol \Delta _{\mathbf x})^T\big[(\mathbf A+\boldsymbol \Delta _{\mathbf x}) \mathbf x - \mathbf y \big] \  \in \ \partial f(\mathbf x),
$$

where $$\boldsymbol \Delta _{\mathbf x} = \mathbf D \odot \text{sign}[ \mathbf x \, (\mathbf A \mathbf x - \mathbf y)^T]$$ is the optimal $$\boldsymbol \Delta $$ for a given $$\mathbf x$$. The symbol $$\odot$$ denotes the Hadamard or element-wise product. Equipped with a subgradient, we can efficiently solve the problem.

[^1]: Becker, Stephen, and Richard J. Clancy. "Robust least squares for quantized data matrices." Signal Processing 176 (2020): 107711.