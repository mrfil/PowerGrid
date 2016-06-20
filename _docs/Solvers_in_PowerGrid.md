---
layout: project
title: Solvers in PowerGrid
category: docs
order_page: 2
---

## Solvers
{: .content-subhead }


#### solve_pwls_pcg

Implements a penalized weighted least square solver using preconditioned conjugate gradient. A penalty object i.e Robject can be used with this solver.

#### solve_grad_desc

Implements a gradient descent solver without a line search. This is not a practical solver but a pedagogic example of how a solver can be implemented in PowerGrid. We strongly encourage you to use this as a boilerplate template if you wish to implement more complex solvers other than the preconditioned conjugate gradient solver above.
