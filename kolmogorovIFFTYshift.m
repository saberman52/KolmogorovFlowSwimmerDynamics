%KOLMOGOROVIFFFT performs inverse Fourier transform of probability
%density Fourier coefficients. Upsamples the data with Fourier
%interpolation (i.e., padding the Fourier coefficients with zeros).
% Shifts Fourier coefficients so y range is from -pi to pi.
function [P,y,theta] = kolmogorovIFFTYshift(Pmn)
M = (size(Pmn,1) - 1)/2;
N = (size(Pmn,2) - 1)/2;
maxF = 350; % target number of total Fourier modes in each dimension
Mp = maxF - M; % number of zeros to add on + and - sides of y Fourier coefficients
Np = maxF - N; % number of zeros to add on + and - sides of theta Fourier coefficients


%% shift Fourier modes so y is going from -pi to pi
m = -M:M;
n = -N:N;
[~,m_] = meshgrid(n,m);
Pmn = (-1).^m_.*Pmn;

%% compute ifft
Pmn_t = zeros(2*maxF+1,2*maxF+1); % new set of Fourier modes
Pmn_t(Mp+1:(Mp+2*M+1),Np+1:(Np+2*N+1)) = Pmn; % fill in given Fourier coefficients; zeros elsewhere
P =  (2*maxF+1)*(2*maxF+1)*ifft2(ifftshift(Pmn_t),'symmetric');
% y and theta domain corresponding to the Fourier coefficients
y = 2*pi/(2*maxF+1)*(0:2*maxF) - pi;
theta = 2*pi/(2*maxF+1)*(0:2*maxF);

end