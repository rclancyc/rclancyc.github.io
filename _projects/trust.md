---
title: "Trust Region Methods"
header:
  overlay_image: /assets/images/positano1.jpg
  caption: "Positano, Italy"
permalink: /projects/trust/
---

# Background
Optimization is concerned with the minimization of an objective or cost function, $$f: \mathbb R^n \rightarrow \mathbb R$$, and is cast mathematically as 

$$
\min_{x \in \mathcal C} f(x) 
$$

where $$\mathcal C$$ is a generic constraint set. For many problems, the objective and/or constaints are too cumbersome for a solution to be computed analytically. In such cases, we turn to numerical optimization methods.

In numerical optimization, an objective or cost function is typically minimized iteratively. This is done by stepping through the function domain by choosing a descent direction at each point. Gradient descent, as the name suggests, does so by moving in the direction opposite the gradient at each iterate. The step-size is specified by the ``learning rate'' which is usually set by the user or some other heuristic. Although the negative gradient is guaranteed to be a descent direction when it is exists, it is possible to step too aggressively and realize a function increase. Barring strong assumptions like Lipschitz continuity, it is difficult to choose the step-size in a principled way. 

There are two broad approaches to address this: line-search and trust regions. Line-search methods choose a descent direction then evaluate the objective at different points along the ray until a sufficient reduction in the cost is observed. The methods is intuitively simple and fairly easy to implement. Unfortunately, in some applications such as PDE constrained optimization, a single function evaluation can be prohibitively expensive. Evaluating the objective numerous times for a single iteration is undesirable.

Trust region methods circumvent costly function evaluations by maintaining a set over which an easily computed model serves a good approximation. At the current iterate $$x_k$$, a model objective function is formed $$m_k(\mathbf p)$$ where $$\mathbf p$$ is the displacement vector from the current iterate. A popular choice of a model is the second order Taylor expansion, i.e., $$f(\mathbf x) \approx m_k(\mathbf p) = f(\mathbf x_k) + \mathbf p^T \nabla f(\mathbf x_k) + (1/2) \mathbf p^T \nabla^2 f(\mathbf x_k) \mathbf p$$. This model is optimized constrained by the trust region. This is known as the trust region subproblem. A simple model and combined with a simple constraint results in an easy to solve subproblem. If a sufficient decrease is realized, the next iterate is chosen to be $$\mathbf{x_{k+1} = x_k + p_k}$$.

<p align="center">
  <img title="https://optimization.mccormick.northwestern.edu/index.php/File:Trust-Region_Method_Overview.png" alt="https://optimization.mccormick.northwestern.edu/index.php/File:Trust-Region_Method.png" src="/assets/images/Trust-Region_Method_Overview.png">
</p>
<p align = "center">
  <em> <small> Figure 1: Trust region progress over several iterations. Courtesty of Ye Wenhe </small> </em>
</p>

This is performed iteratively until a minimizer is found. The algorithm is given by: 

__Trust region algorithm__
- Input $$\mathbf x_0, \ \eta_{good}, \ \eta_{great}, \Delta_0, \gamma_{inc}, \gamma_{dec}$$
- For $$k = 0, 1, 2, ..$$, loop until converged
    - Solve the trust region subproblem $$\mathbf p_k = \text{argmin}_{\mathbf p} \ \ m_k(\mathbf p) \ \   \text{such that }\| \mathbf p \| \le \Delta_k$$
    - Calculate $$\rho_k = \frac{f(\mathbf x_k) - f(\mathbf{x_k + p_k})}{m_k(\mathbf 0) - m_k(\mathbf{p_k})}$$
      - If $$\rho_k > \eta_{good}$$
        - Accept step: $$\mathbf {x_{k+1} = x_{k} + p_k}$$
        - If $$\rho_k > \eta_{great}$$
            - Expand trust region: $$\Delta_{k+1} = \gamma_{inc} \Delta_k$$
        - Else 
            - Use last iterate: $$\mathbf{x_{k+1} = x_k}$$
            - Keep trust region fixed: $$\Delta_{k+1} = \Delta_k$$
      - Else
        - Shrink trust region: $$\Delta_{k+1} = \gamma_{dec} \Delta_k$$


- next $$k$$