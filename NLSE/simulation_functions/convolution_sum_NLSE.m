function convo=convolution_sum_NLSE(u,v,w,kappa)
%
%convo = convolution_sum_NLSE(u,v,alpha)
%
%computes the convolution of three vectors in frequency space C(u,v,w)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u, v, w  = frequency space vectors to be convolved
%
%  kappa  =  degree of nonlinearity
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  convo  = the convolution sum of u, v, and w as C(u,v,w)


%compute length of input vectors for use in generating wavenumber vector
L = length(u);

%flip v to be negative
v = ifftshift(flip(fftshift(v)));

%compute the convolution sum by multiplying in real space, and returning to
%Fourier space
convo = fft_norm(ifft_norm(u).*ifft_norm(conj(v)).*ifft_norm(w));

%implement the derivative and the constant
convo = 1j*kappa*convo;