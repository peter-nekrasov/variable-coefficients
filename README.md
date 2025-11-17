Peter's repository for solving variable coefficient scattering problems in 2D. The exact form of these scattering problems is below. This code uses the fact that the (adjoint) Lippman-Schwinger equation for each problem looks like an identity term plus the sum of diagonal times translationally invariant operators:

$$ \mu(x) + \sum_{i=1}^N V_i(x) \int_{\Omega} K_i (x-y) \mu(y) \ d y  = f(x) $$ 

where $\Omega \subset \mathbb{R}^2$. The code stores both the coefficients $V_i$ and the kernels $K_i$ as rank 3 tensors and sums them up in a systematic way. The quadrature used here is Zeta-corrected trapezoidal rule (see Wu & Martinson, 2021).

There is currently support for:

- Helmholtz:

$$ \Delta u + k^2(1+V(x)) u = 0 $$

- Simplified Flexural (sf):

$$ \Delta^2 u - k^4(1+V(x)) u = 0 $$

- Simplified Flexural-Gravity (sfg):

$$ \frac12 ( \alpha \Delta^2 - \beta (1+V) ) u  + \gamma \int_{\Omega} \frac{1}{4\pi | x - y \ |} u(y) dy = 0 $$

- Flexural-Gravity (flex-gravity):

$$  \frac12 \left( \Delta( \alpha \Delta) + (1-\nu) \left( 2 \partial_{x_1 x_2} \alpha \ \partial_{x_1 x_2} - \partial_{x_1}^2 \alpha \ \partial^2_{x_2}  - \partial_{x_2}^2 \alpha \ \partial_{x_1}^2  \right)  - \beta \right) u + \gamma \int_{\Omega} \frac{1}{4\pi | x - y \ |} u(y) dy = 0  $$



These scattering problems are solved using either FFT + GMRES or using recursive skeletonization with proxy surfaces. The former tends to be faster, but the latter may be better for geometries with resonances where GMRES tends to suffer. 
