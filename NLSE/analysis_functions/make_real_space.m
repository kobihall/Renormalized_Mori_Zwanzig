function [x,u_real] = make_real_space(u,N)
%
%Takes output from a simulation and creates an array of the real space
%solution for plotting and error purposes.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  length(t_list) x resolution array of positive Fourier modes of solutions at
%        every time step
%
%  N  =  number of modes to use in creating real space solution
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  x       =  x coordinates of real space solution
%
%  u_real  =  length(x) x length(t_list) array of real space solution at
%             every time step

%find how many time steps there are
[M,num_times] = size(u);

%augment array with zeros if not using full results
if 2*N+1~=M
    u(N+1:end,:) = 0;
end

u_real = zeros(size(u));

%take inverse fourier transform
for i = 1:num_times
    u_real(:,i) = ifft_norm(u(:,i));
end

%construct x-coordinate system
x = linspace(0,2*pi*(2*M-1)/(2*M),M);