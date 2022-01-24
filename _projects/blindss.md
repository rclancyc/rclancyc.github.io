---
title: "Blind Source Separation"
header:
  overlay_image: /assets/images/chicago_skyline2.jpg
  caption: "Chicago, IL"
permalink: /projects/blindss/
toc: true
classes: wide
---

Blind source separation (BSS) focuses on extracting meaningful components from a signal mixture which is often represented as time course data from a sensor array. The sensors are often microphones, radar receivers, or antenna arrays. My primary motivativation for studying BSS is to localize brain activity in neuroimaging, but the field is rich and worthy of study in its own right. Although machine learning often focuses on deep neural and convolutional networks, BSS methods are popular for unsupervised learning tasks and enjoy broad applicability. On this page, we introduce the basic premise of the problem through example then discuss two methods for separating the constituent signals.

# Background: The Cocktail Party Problem

Consider the classic "cocktail party problem" with several groups in a room all having separate conversations. An observer at the bar hears the din of conversation and can skillfully eavesdrop on a particular conversation. This is a simple task for humans, but difficult to formulate mathematically.

Assume there are $$ N $$ separate conversations and $$M $$ microphones placed throughout the room. Each microphone records a combination of conversations. Groups that are closer to a microphone or speak loudly are heavily represented in that recording or mixture. 

Each conversation can be represented by time series data or a stochastic process given by audio signals from the corresponding group which we write as $$s_1(t), \, \dots, \, s_N(t)$$. The subscript indicates the group and $$ t $$ represents the time of a reading. Similarly, each microphone picks up a time series and is written as $$y_1(t), \, \dots, \, y_M(t)$$. Since each $$y_i(t)$$ for $$i \in \{1, \dots, M\}$$ is a linear combination of $$s_1(t), \, \dots, \, s_N(t)$$, we let $$\mathbf a_i$$ represent the corresponding weights such that $$y_i(t) = \sum_{j=1}^N (\mathbf a_i)_j s_j(t)$$. We can write $$\mathbf y(t) = \mathbf A \mathbf s(t)$$ as the vector of microphone readings as a function of the source signal vector at time $$ t$$. For notational simplicity we drop the $$t $$ when the meaning being clear from context and write the collection of all time points in matrix form as $$\mathbf{Y  = A S}.$$ The rows of $$\mathbf{Y, A}$$ and $$\mathbf S$$ correspond to the time course for each mircophone, the weights of each source for the corresponding recording, and the time courses for each each source, respectively. The colums of $$\mathbf Y$$ and $$\mathbf S$$ correspond to different times. Explicitly,

$$ \mathbf{Y = AS} \\
\Updownarrow
\\
 \begin{bmatrix} y_1(t_1) & \dots & y_1(t_T) \\  \vdots & \ddots & \vdots \\ y_M(t_1) & \dots & y_M(t_T)  \end{bmatrix} =  \begin{bmatrix} (\mathbf a_1)_1 & \dots & (\mathbf a_1)_N \\  \vdots & \ddots & \vdots \\ (\mathbf a_M)_1 & \dots & (\mathbf a_M)_N  \end{bmatrix}  \begin{bmatrix} s_1(t_1) & \dots & s_1(t_T) \\  \vdots & \ddots & \vdots \\ s_N(t_1) & \dots & s_N(t_T)  \end{bmatrix}.$$

The goal of blind source separation is to decompose the signal mixture matrix $$\mathbf Y$$ into a mixing matrix $$\mathbf A$$ and a source matrix $$\mathbf S$$.

# Independent Component Analysis

## The Central Limit Theorem
To understand independent component analysis or ICA, it is helpful to introduce the central limit theorem (CLT). In its simplest form, the CLT states that sums of identically distributed indpendent (i.i.d.) random variables converge to a normal distribution. More precisely, supposing that $$X_1, \dots, X_N$$ is a random sample drawn from a distributuon of mean $$\mu$$ and bounded variance $$\sigma^2$$, then 

$$
  Z = \lim_{N \rightarrow \infty} \sqrt{N} \left(\frac{\bar X_N - \mu}{\sigma}\right) \sim \mathcal{N}(0,1),
$$

where $$\bar X_N$$ is the sample mean. The Lyapunov CLT relaxes the condition that the random variables must be identically distributed but still demands independence. There are other results with different statements and looser restrictions, but for our purposes here, it suffices to say that the distribution for a sum of random variables is more Gaussian than the individual random variables (assuming the random variables are not all Gaussian). By assuming that sources are indpendent, we can exploit the CLT to find an "unmixing matrix" that yields sources signals (or components in the vernacular) that are as indpendent as possible. 

## Measuring Independence

To use ideas from probability theory, we treat source signals and their mixtures as random processes with associated, but unknown, distributions. The source processes are assumed to be independent of each other (e.g. $$s_1(t)$$ does not depend on $$s_2(t)$$) and non-Gaussian. Recall that a mixture is a weighted sum of source signals, i.e., $$y_i(t) = \sum_{j=1}^N (\mathbf a_i)_j s_j(t)$$ or $$\mathbf{Y = AS}$$. Because all sensor readings are sums of the same source signals, then the readings from different microphones will be statistically _dependent_. However, based on assumption, we know that the sources themselves are _independent_. 

Our ultimate goal is to find the source signals that gave rise to the mixtures. We can do this by finding weight vectors $$\mathbf w_1, \dots, \mathbf w_M$$ such that $$z_p(t) = \mathbf w_p^T \, \mathbf y(t)$$ and $$z_q(t) = \mathbf w_q^T \, \mathbf y(t)$$ are as close to being statistically independent as possible for $$p \ne q$$. By stacking the different $$z_p(t)$$ and weight vectors, we write $$\mathbf z(t) = \mathbf W \mathbf y(t)$$ where $$\mathbf W = [\mathbf w_1, \dots, \mathbf w_M]^T$$. We call $$\mathbf W$$ the unmixing matrix. Concatenating the vectors for each time point (such that $$\mathbf Z = [\mathbf z(t_1), \dots, \mathbf z(t_T)]$$) gives the matrix equation $$\mathbf Z = \mathbf {W Y}$$. Since we want $$\mathbf W$$ to unmix the combined signals, the hope is that $$\mathbf Z = \mathbf {W Y} \approx \mathbf S$$. It is easy to see that if $$\mathbf A$$ is invertible, then $$\mathbf A^{-1}$$ is the unmixing matrix since $$\mathbf{Z = W Y = A^{-1} Y = A^{-1} (AS) = S}.$$

So we want to find a matrix, $$ \mathbf W$$, that makes the components of $$\mathbf z(t)$$ independent, but how exactly do we measure independence? We just stated that sums of source signals will be dependent. From the previous section we know that sums of random variables are more normal than their constituent parts. Consquently, we can use Gaussianity as proxy for independence, that is, the more normal, the more likely it is to be a mixture rather than a pure signal. So we use Gaussianity to measure independence but are left with a new question: how do we measure Gaussianity?

## Measuring Gaussianity
To determine how normal or Gaussian a distribution is, we use [cumulants](https://en.wikipedia.org/wiki/Cumulant) which are closely related to moments. The first, second, third, and fourth cumulants correspond to mean, variance, skew, and kurtosis. A full treatment of the cumulants is beyond the scope of this page, but it is important to note that normal distributions are uniquely defined by their mean and variance with all higher order cumulants equaling zero (using appropriate convention for kurtosis). If the sample skew or kurtosis is very large for a random variable, then it is unlikely to be normal. 

Since skew is zero for all distributions symmetric about their mean, kurtosis is typically used to measure Gaussianity. The further the sample kurtosis is from zero, the less Gaussian it is. _A good unmixing matrix $$\mathbf W$$ is one that drives the sample kurtosis for the rows of $$\mathbf{Z=WY}$$ far from zero_.


## Separating Source Signals Using Kurtosis
Source separation can be done in series (one signal at a time) or in parallel (all at once). For clarity of exposition, we consider the former. Accordingly, we focus on finding a single row, $$\mathbf w_i^T$$, of the matrix $$ \mathbf W$$ that gives $$z_i = \mathbf w_i^T \mathbf y$$ the largest possible magnitude value for kurtosis possible. The logic is that if kurtosis is large, the random variable is very non-Gaussian implying that a pure source signal has been recovered rather than a mixture. Since we don't know the true distribution of $$\mathbf y$$ and therefore $$\mathbf z$$, we cannot calculate its true kurtosis directly.  We do on the other hand have time course data which can be used to compute sample kurtosis. Letting $$ k_S(z_i) $$ be the sample kurtosis of $$ z $$, we can find an unmixing vector $$\mathbf w_i$$ by solving the following optimization problem

$$
  \underset{\| \mathbf w_i \| = 1}{ \text{argmax} } \quad  \big| k_S(\mathbf {w}_i^T \mathbf y) \big|.
$$

We must constrain the unmixing vector for their to exist a solution, otherwise $$k_S(\mathbf w_i^T \mathbf{y})$$ is unbounded. So our goal is to solve for the unit vector $$\mathbf w_i$$ that maximizes the sample kurtosis of $$\mathbf w_i^T \mathbf y$$.

Although there might not be a closed for solution, we can solve the problem numerically using gradient based methods. In particualr, since $$k_S$$ has a functional form, we can compute the gradient (or subgradient when $$k_S(\mathbf {w}^T \mathbf{y}) = 0$$) with respect to $$ \mathbf w$$ and use (stochastic) gradient descent. Once a maximizing $$\mathbf w_i$$ is recovered, we can deflate the problem by projecting out the component $$\mathbf w_i^T \mathbf Y$$ from the signal mixture matrix $$\mathbf Y$$. We then find  the next unmixing vector that maximizes the magnitude of sample kurtosis for the delated matrix $$\mathbf Y_d$$. 

## Variations 
Using kurtosis as a proxy for independence is relatively simple and intuitively appealing, but is far from the only approach. Other constrast functions for measuring independence can be derived. Popular alternatives are based on information theoretic (negentropy, mutual information) and distributional (maximum likelihood estimation) frameworks. The derivations of different objectives to use in the ICA optimization problem can be found in the excellent text, _Independent Component Analysis_[^1], and draws on many nice branches of mathematics and statistics. In fact, the use of probability density function approximation through Edgeworth expansions in ICA inspired our work on [approximate maximum likelihood estimation](/projects/regression/#approximate-maximum-likelihood). A well know algorithm called FastICA[^2] has an easy to use implementations in Sci-Kit Learn that quickly factorizes a mixture matrix and is a nice place to start for those interested in blind source separation. 


# Beamformers
ICA focuses on the statistical properties of a source to help separate its constituent parts. Beamformers take a markedly different approach and use spatial filtering to seperate sources. A patial filter is similar to a bandpass filter in the frequency domain; instead of isolating frequencies, spatial filters attenuate/amplify signals from certain locations. More details to come soon!

[^1]: Hyvärinen, Aapo; Karhunen, Juha; Oja, Erkki (2001). Independent component analysis (1st ed.). New York: John Wiley & Sons. ISBN 978-0-471-22131-9.
[^2]: Hyvärinen, Aapo, and Erkki Oja. "A fast fixed-point algorithm for independent component analysis." Neural computation 9.7 (1997): 1483-1492.