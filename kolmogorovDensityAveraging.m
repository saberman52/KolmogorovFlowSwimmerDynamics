% KOLMOGOROVDENSITYAVERAGING computes the phase space density of swimmer in
% Kolmogorov flow using the output of the averaging technique.
% Inputs:
% v0 - swim speed
% alpha - swimmer shape parameter
% psi - vector of Psi values for which averaging calculations were
% performed
% p0 - psi density according to averaging technique at each psi value
% Tp - orbit period at each Psi value
%
% Outputs:
% P - phase space density matrix
% y - grid of y values
% theta - grid of theta values
function [P,y,theta,g] = kolmogorovDensityAveraging(v0,alpha,psi,p0,Tp)
% set up grid
n = 2*350+1;
y = linspace(-pi,pi,n+1);
y = y(1:end-1);
theta = linspace(0,2*pi,n+1);
theta = theta(1:end-1);
[TH,Y] = meshgrid(theta,y);

% construct g function, orbit density factor
inIsland = abs(psi) > abs(Psi(0,0,v0,alpha));
g = 2.^(inIsland).*p0./Tp;

% assign density values
psiVals = Psi(Y,TH,v0,alpha);
% p0Mat = interp1(psi,p0,psiVals,'linear','extrap');
% Tmat = interp1(psi,Tp,psiVals,'linear','extrap');
gMat = interp1(psi,g,psiVals,'linear','extrap');

P = gMat ./ (1 - alpha*cos(2*TH));

end