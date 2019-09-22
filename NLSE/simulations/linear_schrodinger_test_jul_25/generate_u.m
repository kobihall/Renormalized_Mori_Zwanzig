%when kappa=0 the NLSE reduces to the linear schrodinger equation in free
%space (i.e. no potential energy). Exact solutions for this case are known
%and computationally trivial. This is a test to check the validity of the
%underlying code and PDE solver.

clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis_functions

N = 32;

simulation_params.kappa = 0;  %coefficient on nonlinear term in NLSE
simulation_params.dt = 1e-4;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 10;   %how often to save state vector
simulation_params.blowup = 1;     %1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = N;          %number of positive modes to simulate
simulation_params.order = 1;        %use the t-model
simulation_params.initial_condition = @(x) 1j*1.8*exp(-5*(x-pi).^2);
simulation_params.initialization = @(x) complete_init_NLSE(x);  %full simulation

[t_list,u_ROM_freq] = PDE_solve(simulation_params);
save t_list t_list
save u_ROM_freq u_ROM_freq

[x,u_ROM_real] = make_real_space(u_ROM_freq,N);
u_ROM_energy = get_energy(u_ROM_freq,N);
save x x
save u_ROM_real u_ROM_real
save u_ROM_energy u_ROM_energy

x = linspace(0,2*pi*(2*N-1)/(2*N),2*N+1);
u_init = fft_norm(simulation_params.initial_condition(x).');
u_exact_freq = zeros(2*N+1,length(t_list));
k = [0:N,-N:-1];
for i = 1:2*N+1
    u_exact_freq(i,:) = u_init(i)*exp(-1j*(k(i))^2*t_list(:));
end
save u_exact_freq u_exact_freq

[~,u_exact_real] = make_real_space(u_exact_freq,N);
u_exact_energy = get_energy(u_exact_freq,N);
save u_exact_real u_exact_real
save u_exact_energy u_exact_energy