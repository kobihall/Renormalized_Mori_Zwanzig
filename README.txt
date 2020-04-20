This folder contains Matlab and Mathematica files used to construct and simulate reduced order models of partial differential equations on periodic domains. Each subfolder here contains a README that explains the contents of the subdirectories.

%%%%%%%%%%%%%%%%
%SUBDIRECTORIES%
%%%%%%%%%%%%%%%%

Burgers: Matlab files for constructing and simulating reduced order models of Burgers’ equation u_t + uu_x = 0.

Euler: Matlab files for constructing and simulating reduced order models of the three-dimensional Euler’s equations u_t + u dot grad u = -grad p, grad dot u = 0.

KdV: Matlab files for constructing and simulating reduced order models of the Korteweg-de Vries equation u_t + epsilon * u_{xxx} + uu_x = 0.

NLSE: Matlab files for constructing and simulating reduced order models of the Nonlinear Schrödinger Equation i * u_t+u_{xx}-kappa * |u|^2 * u=0

Mathematica: Mathematica files for deriving the functional form of reduced order models for the equations contained in this directory. This includes both the BCH approximation and complete memory approximation.