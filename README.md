Peter's code for solving variable coefficient scattering problems in 2D. The exact form of these scattering problems are below. This code uses the fact that the Lippman-Schwinger equation for each problem looks like an identity plus the sum of translationally invariant operators:

$$ V_0(x) \mu(x) + \sum_{i=1}^N V_i(x) \int_{\Omega} K_i (x-y) \mu(y)  d y  = f(x) $$ 

where $\Omega \subset \mathbb{R}^2$.

There is currently support for:

- Helmholtz

$$ \Delta u + k^2(1+V(x)) u = 0 $$
