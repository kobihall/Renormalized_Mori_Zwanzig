function [nonlin2, CA, CD] = t2model_term_NLSE(u_full,nonlin0,uuu_tilde,kappa,F_modes,G_modes,k)
%
%Computes the complete t^2-model term of NLSE for a given state vector
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  nonlin0  =  the result of a convolution of u_full with itself (markov
%              term)
%
%  uuu_tilde  =  the unresolved part of the markov convolution
%
%  kappa    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%  G_modes  =  vector of which modes in u_full correspond to unresolved
%              modes
%
%  k        =  vector of wavenumbers corresponding to entries of u_full
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin2  =  t^2-term
%
%  CA

uuu = -nonlin0;
uuu(G_modes) = 0;

A = 2*convolution_sum_NLSE(uuu_tilde,u_full,uuu_tilde,kappa);

B = 4*convolution_sum_NLSE(u_full,uuu_tilde,uuu_tilde,kappa);

CA = 1i*k.^2.*uuu_tilde;
CB = convolution_sum_NLSE(u_full,u_full,-2i*k.^2.*u_full - 2*uuu + 4*uuu_tilde,kappa);
CB(F_modes) = 0;
CC = convolution_sum_NLSE(u_full,-2i*k.^2.*u_full - 2*uuu + 2*uuu_tilde,u_full,kappa);
CC(F_modes) = 0;
CD = convolution_sum_NLSE(-1i*k.^2.*u_full-uuu,u_full,u_full,kappa);
CD(F_modes) = 0;
C = convolution_sum_NLSE(u_full,u_full, 2*CA+CB+CC+2*CD ,kappa);

DB = convolution_sum_NLSE(u_full,u_full,-1i*k.^2.*u_full - uuu + 2*uuu_tilde,kappa);
DB(F_modes) = 0;
DC = convolution_sum_NLSE(u_full,-1i*k.^2.*u_full - uuu + uuu_tilde,u_full,kappa);
DC(F_modes) = 0;
D = convolution_sum_NLSE(u_full, CA+DB+DC+CD ,u_full,kappa);

nonlin2 = A + B + C + D;