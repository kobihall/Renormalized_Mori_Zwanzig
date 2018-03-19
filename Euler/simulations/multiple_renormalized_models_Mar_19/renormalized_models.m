function renormalized_models(N,end_time)

format long

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

M = 3*N;

% uniform grid
x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% create initial condition
eval = taylor_green(X,Y,Z);
u_full = fftn_norm(eval);
u = u_squishify(u_full,N);

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

params.k = k;
params.N = N;
params.M = M;
params.func = @(x) tmodel_RHS(x);
params.coeff = scaling_law(N,1); 
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
params.print_time = 1;
params.no_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t1,u_raw1] = ode45(@(t,u) RHS(u,t,params),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array1 = zeros([size(u) length(t1)]);
for i = 1:length(t1)
    u_array1(:,:,:,:,:,i) = reshape(u_raw1(i,:),[N,N,N,3,4]);
end

save(sprintf('t1_%i',N),'t1');
save(sprintf('u_array1_%i',N),'u_array1');



params.func = @(x) t2model_RHS(x);
params.coeff = scaling_law(N,2); 

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t2,u_raw2] = ode45(@(t,u) RHS(u,t,params),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array2 = zeros([size(u) length(t2)]);
for i = 1:length(t2)
    u_array2(:,:,:,:,:,i) = reshape(u_raw2(i,:),[N,N,N,3,4]);
end

save(sprintf('t2_%i',N),'t2');
save(sprintf('u_array2_%i',N),'u_array2');




params.func = @(x) t3model_RHS(x);
params.coeff = scaling_law(N,3); 

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t3,u_raw3] = ode45(@(t,u) RHS(u,t,params),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array3 = zeros([size(u) length(t3)]);
for i = 1:length(t3)
    u_array3(:,:,:,:,:,i) = reshape(u_raw3(i,:),[N,N,N,3,4]);
end

save(sprintf('t3_%i',N),'t3');
save(sprintf('u_array3_%i',N),'u_array3');


params.func = @(x) t4model_RHS(x);
params.coeff = scaling_law(N,4); 

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t4,u_raw4] = ode45(@(t,u) RHS(u,t,params),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array4 = zeros([size(u) length(t4)]);
for i = 1:length(t4)
    u_array4(:,:,:,:,:,i) = reshape(u_raw4(i,:),[N,N,N,3,4]);
end

save(sprintf('t4_%i',N),'t4');
save(sprintf('u_array4_%i',N),'u_array4');







% plot the energy in some modes
energy1 = get_3D_energy(u_array1,N);
energy2 = get_3D_energy(u_array2,N);
energy3 = get_3D_energy(u_array3,N);
energy4 = get_3D_energy(u_array4,N);
figure(1)
hold off
plot(log(t1),log(energy1),'linewidth',2)
hold on
plot(log(t2),log(energy2),'r','linewidth',2)
plot(log(t3),log(energy3),'k','linewidth',2)
plot(log(t4),log(energy4),'c','linewidth',2)
legend(sprintf('ROM order 1, N = %i',N),sprintf('ROM order 2, N = %i',N),sprintf('ROM order 3, N = %i',N),sprintf('ROM order 4, N = %i',N),'location','southwest')
title('Energy in resolved modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,sprintf('energy%i',N),'png')




% plot the helicity
w1 = helicity(u_array1);
w2 = helicity(u_array2);
w3 = helicity(u_array3);
w4 = helicity(u_array4);
figure(2)
hold off
plot(t1,w1,'linewidth',2)
hold on
plot(t2,w2,'r','linewidth',2)
plot(t3,w3,'k','linewidth',2)
plot(t4,w4,'c','linewidth',2)
legend(sprintf('ROM order 1, N = %i',N),sprintf('ROM order 2, N = %i',N),sprintf('ROM order 3, N = %i',N),sprintf('ROM order 4, N = %i',N),'location','southwest')
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,sprintf('helicity%i',N),'png')