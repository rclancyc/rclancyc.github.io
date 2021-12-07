---
title: "Trust Region Methods"
header:
  overlay_image: /assets/images/positano1.jpg
  caption: "Positano, Italy"
permalink: /projects/trust/
---
Optimization is concerned with the minimization of an objective or cost function, $$f: \mathbb R^n \rightarrow \mathbb R$$, and is cast mathematically as 

$$
\min_{x \in \mathcal C} f(x) 
$$

where $$\mathcal C$$ is a generic constraint set. For many problems, the objective and/or constaints are too cumbersome for a solution to be computed analytically. In such cases, we turn to numerical optimization methods.

A standard approach in numerical optimization involves computing a direction to step based on a model that approximates the true objective then stepping some distance in that direction. For those familiar with gradient descent, the direction is chosen in the opposite direction of the gradient and the step-size is chosen by the learning rate. Chosing an appropriate learning rate is often difficult and very frustrating when convergence is proves to be painfully slow. More complex methods aim to ``automate'' this parameter tuning. There are two methods frequently used in practice; line search and trust region methods.

Line search methods seek a step-size by evaluating the objective multiple times along the search direction. The distance to step is determined when a minimum function decrease condition is satisfied. This method is simple and frequently used but potentially requires many function evaluations which can prohibitive if the objective is costly to evaluate. Trust region methods, on the other hand, form a simplified model of the objective (Taylor Series expansions or interpolation) then step to the minimizer of the surrogate. If the minimizer is beyond the "trust region", we step the boundary of the trust region and choose where to expand or keep fixed the trust region radius if the surrogate model is accurate, or shrink it if it performs poorly. 

