%a simple test of just the markov and t terms

clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis_functions

dt = 1e-3;
endtime = 2;
howoften = 10;

simulation_params.kappa = -1;  %coefficient on nonlinear term in NLSE
simulation_params.time_dependence = 1;
simulation_params.dt = dt;      %timestep
simulation_params.endtime = endtime;   %end of simulation
simulation_params.howoften = howoften;   %how often to save state vector
simulation_params.blowup = 1;     %1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 32;          %number of positive modes to simulate
simulation_params.order = 2;        %use the t^2-model
simulation_params.initial_condition = @(x) 1j*1.8*exp(-(x-pi).^2);
simulation_params.initialization = @(x) complete_init_NLSE(x);  %full simulation

simulation_params.coeffs = zeros(2,1);
simulation_params.coeffs(1) = 1;

numcoeff = 50;
t_list = 0:dt*howoften:endtime;
E_list = zeros(numcoeff,length(t_list));

for i = 1:numcoeff
    simulation_params.coeffs(2) = -exp(-i);%-0.792320392542639*kappa^-3.783222616952010*N^-5.825426679579797;
    
    [t_list,u_list] = PDE_solve(simulation_params);
    [x,u_real] = make_real_space(u_list,simulation_params.N);
    E = get_energy(u_list,simulation_params.N);
    for j = 1:length(E)
        E_list(i,j) = E(j);
    end
end