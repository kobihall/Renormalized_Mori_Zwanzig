function nonlin=renormalized_complete_1st_order_NLSE(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the complete 
%t-model of the NLSE equation based upon a "full" model with M positive 
%modes (M>N)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  current state vector
%
%  t  =  current time
%
%  simulation_params: structure containing details of simulation
%
%      kappa    =  coefficient on the nonlinear term
%
%      F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%      G_modes  =  vector of which modes in u_full correspond to unresolved
%                  modes
%
%      k        =  vector of wavenumbers corresponding to entries of u_full
%
%      coeffs   =  2x1 vector of coefficients for terms in memory expansion
%
%      N        =  resolution of ROM
%
%      M        =  resolution of "full" model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin  =  nonlinear part of RHS according to this model



%gather parameters needed for simulation
F_modes = simulation_params.F_modes;
G_modes = simulation_params.G_modes;
k = simulation_params.k;
coeffs = simulation_params.coeffs;
N = simulation_params.N;
M = simulation_params.M;
kappa =simulation_params.kappa;

%compute Markov term
    [nonlin0,u_full] = markov_term_NLSE(u,M,N,kappa);

%compute t-model term
    [nonlin1,~] = tmodel_term_NLSE(u_full,nonlin0,kappa,F_modes);

%compute nonlinear part of right hand side of t-model for NLSE
if simulation_params.time_dependence == 1
    
    nonlin = nonlin0([1:N+1,2*M-N+1:2*M]) + t*coeffs(1)*nonlin1([1:N+1,2*M-N+1:2*M]);
    
elseif simulation_params.time_dependence == 0
    
    nonlin = nonlin0([1:N+1,2*M-N+1:2*M]) + coeffs(1)*nonlin1([1:N+1,2*M-N+1:2*M]);%1:2*N+1
    
end