clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis_functions

N = 1024;
dt = 1e-3;
endtime = 2;
howoften = 10;

simulation_params.kappa = 1;  %coefficient on nonlinear term in NLSE
simulation_params.dt = dt;      %timestep
simulation_params.endtime = endtime;   %end of simulation
simulation_params.howoften = howoften;   %how often to save state vector
simulation_params.blowup = 1;     %1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = N;          %number of positive modes to simulate
simulation_params.order = 1;        %use the exact ODE
simulation_params.initial_condition = @(x) x;%1j*1.8*exp(-5*(x-pi).^2);
simulation_params.initialization = @(x) complete_init_NLSE(x);  %full simulation

[t_list,ue_list] = PDE_solve(simulation_params);
save ue_list ue_list

[x,ue_real] = make_real_space(ue_list,simulation_params.N);
E = get_energy(ue_list,simulation_params.N);
% t_list = 0:dt*howoften:endtime;
% 
% x = linspace(0,2*pi*(2*N-1)/(2*N),2*N+1);
% u_exact = zeros(2*N+1,length(t_list));
% for i = 1:2*N+1
%     u_exact(i,:) = u_init(i)*exp(-1j*(k(i))^2*t_list(:));
% end
% 
% save u_exact u_exact