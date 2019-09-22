function simulation_params = complete_init_NLSE(simulation_params)
%
%simulation_params = complete_init_NLSE(simulation_params)
%
%Takes simulation_params structures and initializes them for a complete
%memory approximation with known effective coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type is saved
%
%           N         =  number of positive resolved modes
%
%           kappa   =  degree of nonlinearity
%
%           dt        =  timestep
%
%           coeffs    =  4x1 vector of coefficients for memory terms
%                        (computed from scaling laws if not provided)
%
%           order     =  number of terms to use in memory (2 or 4, 4 by
%                        default)
%
%    time_dependence  =  if 0, algebraically decaying renormalization
%                        coefficients (default). If 1, constant
%                        renormalization coefficients
%
%  initial_condition  =  initial condition as a function of x
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type is saved
%
%       M        =  size of the full model
%
%       F_modes  =  a cell array of indices for resolved modes in the full
%                   model
%
%       G_modes  =  a cell array of indices for unresolved modes in the full
%                   model
%
%       k        =  a cell array of wave numbers in full model
%
%       u        =  state vector of positive modes of full model
%
%       A        =  matrix of linear part of KdV model
%
%       B        =  function handle of nonlinear part of RHS in full model
%
%       As       =  stiff portion of A to be handled implicitly
%
%       s        =  indices of stiff portion of A
%
%       ns       =  indices of nonstiff portions of A
%
%       cutoff   =  cutoff for stiff portions of A

%create shorthand for ROM system size
N = simulation_params.N;

%general simulation parameters
kappa = simulation_params.kappa;
dt = simulation_params.dt;


%define the ordinates in real space
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N+1);

%define the initial condition
u = fft_norm(simulation_params.initial_condition(x).');

%if not already specified, generates the size of the full model
if ~isfield(simulation_params,'M')
    simulation_params.M = 3*N;
end

%initialize cells indicating index information, and populate them
simulation_params.F_modes = [1:N,2*N:4*N+3,5*N+3:6*N+1];
simulation_params.G_modes = N+1:5*N+2;
simulation_params.k = [0:3*N,-3*N:-1].';

%compute cutoff for different regimes of behavior
cutoff = ceil((2.8/(dt*kappa^2))^(1/3));

if cutoff>N
    ns  = 1:2*N+1;
    s = [];
else
    ns = [1:cutoff,2*N-1-cutoff:2*N+1];
    s = cutoff+1:2*N-2-cutoff;
end

%define the linear and nonlinear portions of the right hand side
A=diag(-1j*[0:N,-N:-1].^2);

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%save data into simulation_params
simulation_params.u = u;

%default order = 4
if ~isfield(simulation_params,'order')
    simulation_params.order = 4;
end

%default time dependence is none
if ~isfield(simulation_params,'time_dependence')
    simulation_params.time_dependence = 0;
end

%if coefficients are not provided, compute them from function
if ~isfield(simulation_params,'coeffs')
    
    if simulation_params.order == 4
        simulation_params.coeffs = zeros(4,1);
        
        simulation_params.coeffs(2) = -1.200557496049101*kappa^-3.701319784646586*N^-5.738411475366152;
        simulation_params.coeffs(4) = -0.331765896671333*kappa^-7.405248509731031*N^-11.471487980232354;
        
        
    elseif simulation_params.order == 2
        simulation_params.coeffs = zeros(2,1);
        simulation_params.coeffs(2) = -0.792320392542639*kappa^-3.783222616952010*N^-5.825426679579797;
        simulation_params.coeffs(1) = 1;
        
    elseif simulation_params.order == 1
        simulation_params.coeffs = zeros(1,1);
        simulation_params.coeffs(1) = 1;
    end
    
end

%save relevant parameters into structure
simulation_params.A = A;

if simulation_params.order == 4
    simulation_params.B = @(x,t) renormalized_complete_4th_order_NLSE(x,t,simulation_params);
elseif simulation_params.order == 2
    simulation_params.B = @(x,t) renormalized_complete_2nd_order_NLSE(x,t,simulation_params);
elseif simulation_params.order == 1
    simulation_params.B = @(x,t) renormalized_complete_1st_order_NLSE(x,t,simulation_params);
end

simulation_params.As = As;
simulation_params.s = s;
simulation_params.ns = ns;
simulation_params.cutoff = cutoff;