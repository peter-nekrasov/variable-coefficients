Peter's code for solving variable coefficient scattering problems in 2D. The exact form of these scattering problems are below. This code uses the fact that the Lippman-Schwinger equation for each problem looks like an identity term plus the sum of diagonal times translationally invariant operators:

$$ \mu(x) + \sum_{i=1}^N V_i(x) \int_{\Omega} K_i (x-y) \mu(y) \ d y  = f(x) $$ 

where $\Omega \subset \mathbb{R}^2$.

There is currently support for:

- Helmholtz

$$ \Delta u + k^2(1+V(x)) u = 0 $$

- Flexural-Gravity

$$  \Delta( \alpha \Delta) \partial_z \phi + (1 - \nu) \left( 2 \partial_{xy} \alpha \ \partial_{xy} - \partial_x^2 \alpha \ \partial^2_y  - \partial_y^2 \alpha \ \partial_x^2  \right) \partial_z \phi - \beta \partial_z \phi + \gamma \phi = 0  $$
