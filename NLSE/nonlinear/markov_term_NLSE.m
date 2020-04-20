function [nonlin0,u_full] = markov_term_NLSE(u,M,N,kappa)
%
%Computes the Markov term of NLSE for a given state vector and size of full 
%model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u        =  a state vector (positive modes only)
%
%  M        =  size of the "full" model upon which we will be basing 
%              calculations
%
%  N        =  size of reduced model
%
%  kappa    =  coefficient on the nonlinear term
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin0  =  the Markov term
%
%  u_full   =  the full state vector (positive and negative modes)

u_full = zeros(2*M+1,1);
u_full(1:N+1) = u(1:N+1);
u_full(2*M-N+2:2*M+1) = u(N+2:2*N+1);

%compute first convolution
nonlin0 = -convolution_sum_NLSE(u_full,u_full,u_full,kappa);