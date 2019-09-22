%a simple test of just the markov and t terms

clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis_functions

simulation_params.kappa = -1;  %coefficient on nonlinear term in NLSE
simulation_params.time_dependence = 1;
simulation_params.dt = 1e-4;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 10;   %how often to save state vector
simulation_params.blowup = 1;     %1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 128;          %number of positive modes to simulate
simulation_params.order = 1;        %use the t-model
simulation_params.initial_condition = @(x) 1j*1.8*exp(-(x-pi).^2);
simulation_params.initialization = @(x) complete_init_NLSE(x);  %full simulation

[t_list,u_list] = PDE_solve(simulation_params);
save t_list t_list
save u_list u_list

[x,u_real] = make_real_space(u_list,simulation_params.N);
E = get_energy(u_list,simulation_params.N);
save x x
save u_real u_real
save u_energy E