---
title: "Uncertain Regression"
usemathjax: true
header:
  overlay_image: /assets/images/vatican_staircase.jpg
  caption: "The Vatican"
permalink: /projects/regression/
toc: true
---

# Background

When the regressors used to fit a statistical model are uncertain, ordinary least squares (the classic regression workhorse) is no longer theoretically appropriate. Unfortunately, accounting for design matrix uncertainty with linear models is often intractable. In what follows we consider two approachs, one based on robust optimization and the other on an "approximate" maximum likelihood framework.


## Generative model

An exceptionally useful technique in engineering and the sciences is to model a response variable as an affine function of related input data. Indeed, simple linear regression serves as an introductory example of model fitting for high school students across the country. Assuming knowledge of $$\mathbf{A} \in \mathbb{R}^{m \times n}$$ and $$\mathbf{y} \in \mathbb{R}^m$$ and using the generative model

$$
\mathbf{y} = \mathbf{A} \mathbf{x} + \boldsymbol{\eta},
$$

with $$\mathbf{x} \in \mathbb{R}^n$$ and $$\boldsymbol{\eta} \in \mathbb{R}^m$$, the goal of regression is to infer the model parameters $$\mathbf{x}$$ that best explain the observations $$\mathbf{y}$$. We focus on the over-determined case where $$m>n$$. In so using this model to fit data, we can easily understand the relationship between regressors and response variables in a simple and intuitive way.

## Ordinary least squares
This problem has been studied extensively with the most common method being ordinary least squares (OLS) which can be written as

$$
\mathbf{\hat x}_{\text{OLS}} = \text{argmin}_{\mathbf{x}} \| \mathbf{A x - y}\|^2.
$$

For OLS, we minimize sum of residuals squared, that is, we make the mismatch between the observed data and modeled data as small as possible. The solution can be found in easily found and is given by the closed form solution $$\mathbf{\hat x}_{\text{OLS}} = (\mathbf{A}^T \mathbf{A})^\dagger \mathbf{A}^T \mathbf{y}$$ where $$\dagger$$ denotes the Moore-Penrose pseudo-inverse. Although the problem is intuitively appealing and its solution is simple, OLS suffers from several draw backs. 

First, the problem is often poorly conditioned which a numerical analysts way of saying the solution is sensitive to noise in the data. The normal equation solution provided above exacerbates problems with conditioning. Second, the standard formulation assumes that the data or design matrix $$\mathbf A$$ is known precisely which is often a poor assumption. Typical causes of operator uncertainty are sampling error, measurement error, human error, modeling error, or rounding error.

## Uncertain design matrix
 A classical method for addressing this uncertainty is by using total least squares (TLS) which solves the problems

$$
\min_{\mathbf{x}, \boldsymbol{\eta, \Delta}}  \quad \big\| [\boldsymbol{ \Delta, \ \ \eta}] \big\|^2_F \quad \\
\text{subject to} \quad  (\mathbf{A} + \boldsymbol{\Delta}) \mathbf{x = y} + \boldsymbol{\eta}, \\
$$


and has the solution $$\mathbf{\hat x}_{\text{TLS}} = (\mathbf{A}^T \mathbf{A} - \sigma_{n+1}^2 \mathbf{I})^\dagger \mathbf{A}^T \mathbf{y}$$ where $$\sigma_{n+1}$$ is the smallest singular value of the matrix $$[\mathbf{A, y}]$$. $$F$$ denotes the Frobenius norm. Unfortunately, the TLS solution is even worse condition that OLS. There is also an implicit assumption of Gaussian uncertainty in TLS as well which might be undesirable. In what follows, we consider two approaches for handling uncertainty in the design matrix. In particular, we focus on a robust least squares formulation and another method that forms and optimizes an approximate maximum likelihood function.


# Robust Least Squares
This subsection is based on the work published in [Signal Processing](https://www.sciencedirect.com/science/article/abs/pii/S0165168420302541)[^1]. A preprint can be found [here](https://arxiv.org/pdf/2003.12004). 


## Robust formulation
Due to the limitations of classical methods, we approach the problem through the lens of robust optimization. The primary goal is to recover unknown parameters from a noisy observation and an uncertain linear operator. In particular, our mathematical model for robust least squares is

$$
\min_{\mathbf x} \, \left\{ \max_{ \boldsymbol \Delta  \in \,  \mathcal U }\; \| (\mathbf A+ \boldsymbol \Delta ) \mathbf x - \mathbf y \|^2 \right\}
$$

which we refer to as our robust optimization (RO) problem. The Euclidean norm is denoted by $$\| \cdot \|$$ and $$\mathcal{U}$$ is the uncertainty set from which perturbations in $$\mathbf A$$ are drawn. The above RO formulation is motivated by two situations. In both cases, we let $$\bar{\mathbf A}$$ and $$\bar{\mathbf x}$$ represent the _true and unknown_ data matrix and parameter vector, respectively. We model $$\mathbf y = \bar{ \mathbf A} \bar{ \mathbf x} + \boldsymbol \eta$$, with $$\boldsymbol  \eta$$ i.i.d. Gaussian, and we only have knowledge of $$\mathbf A = \bar{ \mathbf A }- \boldsymbol \Delta$$. We don't know $$\boldsymbol \Delta$$ explicitly but can make inferences based on the problem. In the first situation, we consider a data matrix subject to quantization or round-off error. Suppose the observed matrix $$\mathbf A$$ has elements rounded to the hundredth place. Our uncertainty set can be written as $$\mathcal U= \{ \boldsymbol \Delta  \in \mathbb R^{m \times n } : \| \boldsymbol \Delta \|_{\infty} \le \delta \}$$ with $$\delta = 0.005$$. If $$\mathbf A_{i,j} = 0.540$$, then we know the true $$\bar {\mathbf A}_{i,j} \in (0.535,\, 0.545]$$, hence $$\boldsymbol \Delta _{i,j} \in (-0.005, 0.005]$$. The norm $$\| \cdot \|_{\infty}$$ takes the maximum absolute value of any element in the matrix. 

The second problem considers a data matrix with uncertainty proportional to the magnitude of each entry, i.e. 
$$\mathcal U = \{ \boldsymbol \Delta  \in \mathbb R^{m \times n} : \boldsymbol \Delta_{i,j} \in (-p |\mathbf A_{i,j}|, p | \mathbf A_{i,j}| ] \}$$. Here, $$p$$ denotes a proportionality constant. Data subject to $$\pm 1 \%$$ uncertainty would have $$p = 0.01$$. The two cases cover the effects of finite-precision in fixed and floating point representations, respectively. In both problems, the uncertainty sets are specified element-wise allowing us to decouple along rows.

In this project, we consider box constraints over entries of the data matrix. Our main contribution is the formulation of a robust objective to handle quantization and floating point error in least squares problems and the presentation of methods to solve it. Although the proposed method requires selection of an uncertainty parameter, it is chosen in a natural and theoretically appropriate way based on the observed extent of quantization. This is in contrast to ridge regression and MLE based methods that require involved parameter tuning or *a priori* knowledge of the probability distributions from which model uncertainty is drawn. We anticipate our method to be most effective under moderate to heavy quantization where fidelity loss is greater than $$0.1\%$$.

## Solving the robust optimization problem
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







# Approximate Maximum Likelihood 
The work is this subsection is based on our preprint that can be found [here](https://arxiv.org/pdf/2104.03307)[^2].

## Regression as maximum likelihood estimation
A reasonable question to ask is ``why should we consider maximum likelihood estimation for regression?''. It turns out that when we consider additive noise alone, i.e., $$\boldsymbol \eta = \mathbf{y - Ax}$$ is a random variable, the OLS solution coincides with the maximum likelihood estimator (MLE) for $$\boldsymbol{\eta} \sim \mathcal{N} (\boldsymbol{0}, \sigma^2 \mathbf{I})$$. Similarly the solutions to $$\min_{\mathbf{x}} \| \mathbf{A x - y} \|_{\infty}$$ (minimax regression) and $$\min_{\mathbf{x}} \| \mathbf{A x - y} \|_1$$ (least deviation regression) coincide with the MLEs for uniform and Laplacian noise, respectively. This relationship between regression and MLE for different noise models suggests using MLE to infer model parameters in more complicated settings such as when the design matrix uncertain.

## MLE formulation
We consider the general case where the design matrix $$\mathbf{A}$$ is a random variable (not necessarily Gaussian) as well. To amplify this fact, we change the design matrix to be a random $$\mathbf G$$ instead of a fixed $$\mathbf A$$. MLE regression problems rely on knowledge of the vector $$\mathbf y$$'s joint probability density function (PDF). Note that each component of $$\mathbf y$$ in the generative model is the sum of scaled random variables, i.e., $$y_i = \mathbf g_i^T \mathbf x + \eta_i = \sum_{i=1}^n G_{ij}x_j + \eta_i$$ where $$\mathbf g_i^T$$ is the $$i^\text{th}$$ row of $$\mathbf G$$ and subscripts denote the component of the corresponding vector. Despite the innocuous form, sums of random variables are difficult to work with: individual PDFs must be convolved to obtain a PDF for their sum. For general noise models, forming an exact MLE is intractable and necessitates alternate techniques for solving the MLE problem.

To avoid the difficulties associated with convolving many densities, we use properties of moment generating functions (which we denote by $$M_A(t)$$ for a random variable $$A$$), then draw on ideas from complex analysis by using the saddle point approximation to estimate a density for $$\mathbf y$$.

The saddle point method or [method of steepest descent](https://en.wikipedia.org/wiki/Method_of_steepest_descent) is a generalization of Laplace's method and was first used in statistics by Daniels in his seminal paper[^3] to estimate the PDF of sample means. An accessible overview of the method for the purposes here was written by [Goutis and Casellas](https://www.tandfonline.com/doi/abs/10.1080/00031305.1999.10474463?casa_token=WdBajCdgWu8AAAAA%3ATfkJRURS_IiM34aU321ekY6n3JTeiB80j9FpZxZWkqC2GpyI_fYfghhyDu1qhF_BfwgvTKLgvJZJ&)[^4]. Assuming the elements of $$\mathbf G$$ and $$\boldsymbol \eta$$ are independent and that their distributions are known as $$\mathcal  P_{\mathbf G}$$ and $$\mathcal  P_{\boldsymbol \eta}$$, respectively, we can write the approximate density as

$$
p(\mathbf y \, | \, \mathbf x, \mathcal  P_{\mathbf G},\mathcal  P_{\boldsymbol \eta}) \approx \prod_{i=1}^m  \left( \frac{\exp\left\{ \sum_{i=1}^m K_{\mathbf g_i^T \mathbf x+ \eta_i}(t_i) - t_i y_i\right\}}{\sqrt{2\pi \, K_{\mathbf g_i^T \mathbf x+ \eta_i}''(t_i)}}  \right) \ ,
$$

where $$t_i$$ is the solution $$q_i(\mathbf x, t_i) = K_{\mathbf g_i^T \mathbf x + \eta_i}'(t_i) - y_i = 0$$ and $$K_{\mathbf g_i^T \mathbf x + \eta_i}(t_i) = \ln M_{\mathbf g_i^T \mathbf x + \eta_i}(t_i)$$ is the **cumulant generating function** (CGF) of $$y_i$$. By taking the logarithm of the PDF and dropping unnecessary constants we quickly arrive at our element-wise approximate log-likelihood function

$$
\ell(\mathbf x) = \sum_{i=1}^m  \left[ \rule{0cm}{.75cm} \right.  K_{\eta_i}(t_i) \ +\  \left(\sum_{j=1}^n K_{G_{ij}}( t_i x_j) \right) \qquad \\ \qquad \qquad
\ -\  t_i y_i - \frac{1}{2} \ln\left( K''_{\eta_i} (t_i) + \sum_{j=1}^n K''_{G_{ij}}(t_i x_j) \right)  \left. \rule{0cm}{.75cm} \right].
$$

We note that the terms of the log-likelihood function are simple univariate CGFs and can be retrieve from most mathematical statistics texts. For notational simplicity, it is helpful to consider the matrix/vector version given by

$$
\ell(\mathbf x) = \boldsymbol 1^T \left( K_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \frac 1 2 \ln \left( K''_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t)   \right) \right) - \mathbf t^T \mathbf y,
$$

where $$\mathbf 1$$ is a vector of ones and $$\mathbf t$$ is the solution to $$\mathbf q(\mathbf x, \mathbf t) = K'_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t) - \mathbf y = \mathbf 0$$ which can be solved with a bracketed root-finding algorithm. We use the notation that $$K^{(i)}_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t) = [K^{(i)}_{\mathbf g_1^T \mathbf x + \eta_1}(t_1), \, \dots, K^{(i)}_{\mathbf g_m^T \mathbf x + \eta_m}(t_m)]^T$$. 

## MLE as an optimization problem
We can cast the approximate MLE problem as

$$
\text{argmax}_{\mathbf x, \mathbf t} \quad \ell(\mathbf x, \mathbf t)
$$

$$
\text{subject to} \quad K'_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \mathbf y= \boldsymbol 0.
$$

Interestingly, there appears to be little effort made towards solving this problem for more that a few parameters in the literature. A contribution of this project is that we find an exact gradient for $$\ell(\mathbf x)$$ that can be used in ``off-the-shelf'' first-order algorithms. Specifically, we have

$$
\nabla_{\mathbf x} \ell = \frac{\partial \ell}{\partial \mathbf x}  - \left( \frac{\partial \ell}{\partial \mathbf t} \right) \left( \frac{\partial \boldsymbol q}{\partial \mathbf t} \right)^{-1}  \left( \frac{\partial \boldsymbol q}{\partial \mathbf x} \right).
$$

with the different terms and factors given by

$$
\frac{\partial \ell}{\partial \mathbf x} \ \ = \ \  \boldsymbol 1^T \left(   \frac{\partial K_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t)}{\partial \mathbf x} - \frac 1 2 \left\{ \text{diag}\left( K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) \right)  \right\}^{-1} \frac{\partial K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t)}{\partial \mathbf x}    \right)_{_{}},
$$

$$
\frac{\partial \ell}{\partial \mathbf t} \ \ = \ \  \left( K'_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \frac 1 2 \left( K'''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) \oslash K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) \right) - \mathbf y \right)^T _{_{}},
$$

$$
\frac{\partial \boldsymbol q}{\partial \mathbf t_{_{}}}  \ \ = \ \   \text{diag}\left( K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) \right)_{_{}} ,
$$

$$
\frac{\partial \boldsymbol q}{\partial \mathbf x} \ \ = \ \   \frac{\partial}{\partial \mathbf x}  K'_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t).
$$

The symbol $$\oslash$$ denotes element-wise matrix division. Although these derivatives can typically be calculated by hand, it is advantageous to use automatic differentiation software if it is readily available.

Equipped with an approximate MLE function and a method for computing the gradient, we are capable of using a variety of first order methods to solve the regression problems subject to uncertain data. 

[^1]: Becker, Stephen, and Richard J. Clancy. "Robust least squares for quantized data matrices." Signal Processing 176 (2020): 107711.

[^2]: Clancy, Richard J., and Stephen Becker. "Approximate maximum likelihood estimators for linear regression with design matrix uncertainty." arXiv preprint arXiv:2104.03307 (2021).

[^3]: Daniels, Henry E. "Saddlepoint approximations in statistics." The Annals of Mathematical Statistics (1954): 631-650.
    
[^4]: Goutis, Constantino, and George Casella. "Explaining the saddlepoint approximation." The American Statistician 53.3 (1999): 216-224.
    

