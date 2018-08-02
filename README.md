# Omega Equation

## About

The Omega equation is a time independent partial differential which governs the vertical
velocity of water. Normally, equations governing fluids are time independent, but a particular
approximation, known as the quasigeostrophic approximation, allows one to reduce the complicated
Navier-Stokes equations into this simpler form.

This approximation has been hisotrically used in the open ocean, but recently, Libe Washburn and Chris Gotschalk
discovered the presence of phytoplankton well below the depth which sunlight allows them to live. Using data
that they gathered in the Santa Barbara channel, one can solve the Omega equation for the vertical velocity
of the water column as a possible explination for the surprising depth of the phytoplankton.

## Running the code  
Simply clone the repo and run the Matlab scipt omega_eqn_zero_neumann.m.

## Technical Information
The *.mat files are the necessary data in order to run the code. The omega equation is a non-constant coefficient
second order elliptic pde, so we solve it iteratively by GMRES. The domain is rectangular prism, with constant step
size in each axis direction.

The two matlab files correspond to different boundary conditions:

omega_eqn_zero_dirichlet_fft_z.m corresponds to zero Dirichlet boundary conditions on all faces of the box. This is 
not physically resonable, and is included for testing purposes, and for comparing against the academic literature.

omega_eqn_zero_neumann.m is a mixed type Neumann-Dirichlet boundary condition.  On the sides of the box we use a natural
boundary condition. This amounts to choosing a Neumann boundary condition with flux coming from the Q vector. On the
bottom of the box we use a zero Neumann boundary condition, while on the top we use a zero Dirichlet boundary conditions


## References
Omega Equation: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96JC01144

GMRES: https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
