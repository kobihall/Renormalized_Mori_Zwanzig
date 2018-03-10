function resolve_data(N)
%
% A function to compute the hypothetical flow of energy out of the full
% model of size N
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%         N  =  resolution of full model  
%
%
%%%%%%%%%
%OUTPUS:%
%%%%%%%%%
%
%  tmodel_size_listN.dat  =  saved vector of the time derivative of the
%                            t-model of the full model at associated times

% load relevant folders into the path
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

if exist(sprintf('tmodel_size_list%i.mat',N),'file') == 2
    
    % if there is already data for this resolution, load it and continue the
    % simulation from the previous end time up to the proposed end time
    load(sprintf('tmodel_size_list%i.mat',N))
    completed = length(tmodel_size_list);
    
else
    
    % otherwise, start at time 0
    completed = 0;
    tmodel_size_list = [];
    
end

% load data
load(sprintf('u%i.mat',N))
load(sprintf('t%i.mat',N))

if length(t) > completed
    
    % gather the new data
    u_new = u(:,:,:,:,:,completed+1:end);
    t_new = t(completed+1:end);
    
    % use the resolve_array function to resolve each term
    tmodel_size_list_new = resolve_array(u_new,t_new);
    
    % append the results and save it to the directory
    tmodel_size_list = [tmodel_size_list;tmodel_size_list_new];
    save(sprintf('tmodel_size_list%i.mat',N),'tmodel_size_list')
    
end