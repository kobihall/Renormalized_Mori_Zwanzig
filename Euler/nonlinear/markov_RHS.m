function du_dt = markov_RHS(u_full,a,b,k,a_tilde,N)
%
% Computes the RHS for every mode in the markov model for 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%        a  =  indices of positive resolved modes 1:M
%
%        b  =  indices of negative resolved modes -M:-1
%
%        k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%  a_tilde  =  indices of unresolved modes
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  derivative of each mode

% the full model is a simple convolution Ck(u,u)
du_dt = markov_term(u_full,a,b,k,a_tilde);

du_dt = u_squishify(du_dt,N);