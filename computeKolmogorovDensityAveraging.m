% COMPUTEKOLMOGOROVDENSITYAVERAGING calculates the steady-state density
% P(y,\theta) predicted by the averaging model in the sigma -> 0 limit for
% the Kolmogorov flow, assuming rotational diffusion only.
%
% Inputs:
% v0 - swim speed
% alpha - swimmer shape parameter
% nPsi - number of psi values to use in initial psi grid (default = 500)
%
% Outputs:
% P - the probability density matrix P(y,\theta)
% y - y values where P is evaluated
% theta - theta values where P is evaluated
% psi - psi values where averaging calculations performed
% p0 - psi probability density p_0(\Psi)
% g - orbit weighting function in phase space g(\psi) = 2.^(inIsland)*p0./T
function [P,y,theta,psi,p0,g] = computeKolmogorovDensityAveraging(v0,alpha,nPsi)
if nargin < 3
    nPsi = 500;
end

% compute p0 and T
if v0 < kolmogorovBifurcation(alpha) % select function based on bifurcation
    densityFunc = @computePsiDensityBelow2;
else
    densityFunc = @computePsiDensityAbove2;
end
[psi,p0,~,~,Tp] = densityFunc(v0,alpha,nPsi);


% compute density
[P,y,theta,g] = kolmogorovDensityAveraging(v0,alpha,psi,p0,Tp);

end