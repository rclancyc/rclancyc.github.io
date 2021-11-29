---
title: "Approximate Maximum Likelihood Estimation"
usemathjax: true
header:
  overlay_image: /assets/images/mammoth.jpg
  caption: "Mammoth Hot Springs, Yellowstone NP"
permalink: /projects/amle/
---
# Generative model 1

An exceptionally useful technique in engineering and the sciences is to model a response variable as an affine function of related input data. Indeed, simple linear regression serves as an introductory example of model fitting for high school students across the country. Assuming knowledge of $\mathbf{A}$ and $\mathbf{y}$ and using the generative model

$$
\mathbf{y} = \mathbf{A} \mathbf{x} + \boldsymbol{\eta}, \qquad \text{with} \qquad \mathbf{A} \in \mathbb{R}^{m \times n}, \ \mathbf{x} \in \mathbb{R}^n, \ \mathbf{y}, \  \boldsymbol{\eta} \in \mathbf{R}^m,
$$

the goal of regression is to infer the model parameters $\mathbf{x}$ that best explain the observations $\mathbf{y}$. We focus on the over-determined case where $m>n$. In so using this model to fit data, we can easily understand the relationship between regressors and response variables in a simple and intuitive way.

# Regression as maximum likelihood estimation

This problem has been studied extensively with the most common method being ordinary least squares (OLS) which can be written as

$$
\mathbf{\hat x}_{\text{OLS}} = \argmin_{\mathbf{x}} \| \mathbf{A x - y}\|^2,
$$

and is the maximum likelihood estimator (MLE) for $\mathbf{x}$ when $\boldsymbol{\eta} = \mathbf{y-A x} \sim \mathcal{N} (\boldsymbol{0}, \sigma^2 \mathbf{I})$. An enticing feature of OLS is that it has a closed form solution given by $\mathbf{\hat x}_{\text{OLS}} = (\mathbf{A}^T \mathbf{A})^\dag \mathbf{A}^T \mathbf{y}$ where $\dag$ denotes the Moore-Penrose psuedo-inverse. Unfortunately, the method is known to suffer from poor conditioning. There are other formulations such as $\min_{\mathbf{x}} \| \mathbf{A x - y} \|_{\infty}$ (minimax regression) or $\min_{\mathbf{x}} \| \mathbf{A x - y} \|_1$ (least deviation regression) that coincide with the MLEs for uniform and Laplacian noise, respectively. This relationship between regression and MLE for different noise models suggests using MLE to infer model parameters in more complicated settings.

# Uncertain design matrix

All the above problems make the strong and often unrealistic assumption that the design matrix or operator $\mathbf{A}$ is known with certainty. Typical causes of operator uncertainty are sampling error, measurement error, human error, modeling error, or rounding error. A classical method for addressing this uncertainty is by using total least squares (TLS) which solves the problems

$$
\min_{\mathbf{x}, \boldsymbol{\eta, \Delta}}  \quad \big\| [\boldsymbol{ \Delta, \ \ \eta}] \big\|^2_F \quad \text{subject to} \quad  (\mathbf{A} + \boldsymbol{\Delta}) \mathbf{x = y} + \boldsymbol{\eta}
$$

and has the solution $\mathbf{\hat x}_{\text{TLS}} = (\mathbf{A}^T \mathbf{A} - \sigma_{n+1}^2 \mathbf{I})^\dag \mathbf{A}^T \mathbf{y}$ where $\sigma_{n+1}$ is the smallest singular value of the matrix $[\mathbf{A, y}]$. Unfortunately, the TLS solution is even worse condition that OLS. There is also an implicit assumption of Gaussian uncertainty in TLS as well which might be undesirable. 

We consider the general case where the design matrix $\mathbf{A}$ is a random variable (not necessarily Gaussian) as well. To amplify this fact, we change the design matrix to be a random $\mathbf G$ instead of a fixed $\mathbf A$. MLE regression problems rely on knowledge of the vector $\mathbf y$'s joint probability density function (PDF). Note that each component of $\mathbf y$ in the generative model is the sum of scaled random variables, i.e., $y_i = \mathbf g_i^T \mathbf x + \eta_i = \sum_{i=1}^n G_{ij}x_j + \eta_i$ where $\mathbf g_i^T$ is the $i^\text{th}$ row of $\mathbf G$ and subscripts denote the component of the corresponding vector. Despite the innocuous form, sums of random variables are difficult to work with: individual PDFs must be convolved to obtain a PDF for their sum. For general noise models, forming an exact MLE is intractable and necessitates alternate techniques for solving the MLE problem.

# Approximate maximum likelihood estimation

To avoid the difficulties associated with convolving many densities, we use properties of moment generating functions (which we denote by $M_A(t)$ for a random variable $A$), then draw on ideas from complex analysis by using the saddle point approximation to estimate a density for $\mathbf y$.

The saddle point method or [method of steepest descent]](https://en.wikipedia.org/wiki/Method_of_steepest_descent) is a generalization of Laplace's method and was first used in statistics by Daniels in his seminal paper[^1] to estimate the PDF of sample means. An accessible overview of the method for the purposes here was written by [Goutis and Casellas](https://www.tandfonline.com/doi/abs/10.1080/00031305.1999.10474463?casa_token=WdBajCdgWu8AAAAA%3ATfkJRURS_IiM34aU321ekY6n3JTeiB80j9FpZxZWkqC2GpyI_fYfghhyDu1qhF_BfwgvTKLgvJZJ&)[^2]. Assuming the elements of $\mathbf G$ and $\boldsymbol \eta$ are independent and that their distributions are known as $\mathcal  P_{\mathbf G}$ and $\mathcal  P_{\boldsymbol \eta}$, respectively, we can write the approximate density as

$$
p(\mathbf y \, | \, \mathbf x, \mathcal  P_{\mathbf G},\mathcal  P_{\boldsymbol \eta}) \approx \prod_{i=1}^m \left[ \left(2\pi \, K_{\mathbf g_i^T \mathbf x+ \eta_i}''(t_i)\right)^{-1/2}\right]  \ \exp\left\{ \sum_{i=1}^m K_{\mathbf g_i^T \mathbf x+ \eta_i}(t_i) - t_i y_i\right\},
$$

where $t_i$ is the solution $q_i(\mathbf x, t_i) = K_{\mathbf g_i^T \mathbf x + \eta_i}'(t_i) - y_i = 0$ and $K_{\mathbf g_i^T \mathbf x + \eta_i}(t_i) = \ln M_{\mathbf g_i^T \mathbf x + \eta_i}(t_i)$ is the **cumulant generating function** (CGF) of $y_i$. By taking the logarithm of the PDF and dropping unnecessary constants we quickly arrive at our element-wise approximate log-likelihood function

$$
\ell(\mathbf x) = \sum_{i=1}^m  \left[ \rule{0cm}{.75cm} \right.  K_{\eta_i}(t_i) \ +\  \left(\sum_{j=1}^n K_{G_{ij}}( t_i x_j) \right)
\ -\  t_i y_i - \frac{1}{2} \ln\left( K''_{\eta_i} (t_i) + \sum_{j=1}^n K''_{G_{ij}}(t_i x_j) \right)  \left. \rule{0cm}{.75cm} \right].
$$

We note that the terms of the log-likelihood function are simple univariate CGFs and can be retrieve from most mathematical statistics texts. For notational simplicity, it is helpful to consider the matrix/vector version given by

$$
\ell(\mathbf x) = \boldsymbol 1^T \left( K_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \frac 1 2 \ln \left( K''_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t)   \right) \right) - \mathbf t^T \mathbf y,
$$

where $\mathbf 1$ is a vector of ones and $\mathbf t$ is the solution to $\mathbf q(\mathbf x, \mathbf t) = K'_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t) - \mathbf y = \mathbf 0$ which can be solved with a bracketed root-finding algorithm. We use the notation that $K^{(i)}_{\mathbf G \mathbf x + \boldsymbol \eta}(\mathbf t) = [K^{(i)}_{\mathbf g_1^T \mathbf x + \eta_1}(t_1), \, \dots, K^{(i)}_{\mathbf g_m^T \mathbf x + \eta_m}(t_m)]^T$. We can cast the approximate MLE problem as

$$
\argmax_{\mathbf x, \mathbf t} \quad \ell(\mathbf x, \mathbf t)
$$

$$
\text{subject to} \quad K'_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \mathbf y= \boldsymbol 0.
$$

Interestingly, there appears to be little effort made towards solving this problem for more that a few parameters in the literature. A contribution of this project is that we find an exact gradient for $\ell(\mathbf x)$ that can be used in ``off-the-shelf'' first-order algorithms. Specifically, we have

$$
\nabla_{\mathbf x} \ell = \frac{\partial \ell}{\partial \mathbf x}  - \left( \frac{\partial \ell}{\partial \mathbf t} \right) \left( \frac{\partial \boldsymbol q}{\partial \mathbf t} \right)^{-1}  \left( \frac{\partial \boldsymbol q}{\partial \mathbf x} \right).
$$

with the different terms and factors given by

$$
\frac{\partial \ell}{\partial \mathbf x} \ \ = \ \  \boldsymbol 1^T \left(   \frac{\partial}{\partial \mathbf x}K_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) - \frac 1 2 \left\{ \text{diag}\left( K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t) \right)  \right\}^{-1} \frac{\partial}{\partial \mathbf x}K''_{\mathbf G \mathbf x+ \boldsymbol \eta}(\mathbf t)    \right)_{_{}},
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

The symbol $\oslash$ denotes element-wise matrix division. Although these derivatives can typically be calculated by hand, it is advantageous to use automatic differentiation software if it is readily available.

Equipped with an approximate MLE function and a method for computing the gradient, we are capable of using a variety of first order methods to solve the regression problems subject to uncertain data. More details on this project can be found in our [preprint](https://arxiv.org/pdf/2104.03307)[^3].

[^1]: Daniels, Henry E. "Saddlepoint approximations in statistics." The Annals of Mathematical Statistics (1954): 631-650.
    
[^2]: Goutis, Constantino, and George Casella. "Explaining the saddlepoint approximation." The American Statistician 53.3 (1999): 216-224.
    
[^3]: Clancy, Richard J., and Stephen Becker. "Approximate maximum likelihood estimators for linear regression with design matrix uncertainty." arXiv preprint arXiv:2104.03307 (2021).
