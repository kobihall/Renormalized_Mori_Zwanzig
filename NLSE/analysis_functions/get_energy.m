function energy = get_energy(kappa,u,N)
%
%Calculates the energy in the first N modes at each time step for the 
%solution matrix u
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  N x length(t_list) array of Fourier modes at all timesteps
%
%  N  =  maximal mode to include when computing energy
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  energy  =  the energy in modes 1:N at all times
k = [0:N,-N:-1].';

energy = sum(1/2 * k.^2 .* (abs(u(1:2*N+1,:))).^2 + 1/2 * kappa * (abs(u(1:2*N+1,:))).^4).';