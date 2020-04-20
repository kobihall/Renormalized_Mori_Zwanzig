function [nonlin1,uuu_tilde] = tmodel_term_NLSE(u_full,nonlin0,kappa,F_modes)
%
%Computes the t-model term of KdV for a given state vector
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin0  =  the result of a convolution of u with itself (markov
%              term)
%
%  kappa    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin1  =  the t-model term
%
%  uu_star  =  the unresolved part of the markov convolution

%eliminate resolved modes of first convolution
uuu_tilde = -nonlin0;
uuu_tilde(F_modes)=0;

%compute t model term
nonlin1 = 2*convolution_sum_NLSE(u_full,u_full,uuu_tilde,kappa)+convolution_sum_NLSE(u_full,uuu_tilde,u_full,kappa);