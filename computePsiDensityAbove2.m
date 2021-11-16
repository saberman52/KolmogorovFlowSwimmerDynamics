% COMPUTEPSIDENSITYABOVE calculates the Psi distribution under the adiabatic
% averaging model for a diffusing swimmer in the Kolmogorov flow, for
% parameter values (v0,alpha) above the bifurcation.
% Diffusivity is averaged directly, rather than averaging square root of
% diffusivity.
%
% Inputs:
% v0 - swim speed
% alpha - swimmer shape parameter
% nPsi - number of (positive) psi values to sample (optional, default =
% 300)
% 
% Outputs: 
% psi - range of psi values where adiabatic averaging calculation succeeded
% P_psi - normalized psi density at every point
% fPsi - average drift vector
% D - average diffusivity
% Tp - period of orbit with given Psi
function [psi,P_psi,fPsi,D,Tp] = computePsiDensityAbove2(v0,alpha,nPsi)
if nargin < 3
    nPsi = 300;
end
maxPsi = Psi(0,pi,v0,alpha); % value of psi at elliptic equilibrium point
psi = linspace(0,maxPsi,nPsi+1);
% resample the end, where f/D varies a lot 
psi = psi(1:end-1);
psiEnd = psi(end);
endThresh = 0.1;
endPoints = 100;
psi = psi(abs(psi - psi(end)) > endThresh);
psiEdge = linspace(psi(end),psiEnd,endPoints+1);
psi = [psi,psiEdge(2:end)];
nPsi = length(psi);

y0 = 0;
% convert psi to th0
th0 = thetaOrbit(y0,psi);

% storage
dPsi = zeros(1,nPsi);
DInt = zeros(1,nPsi);
Tp = zeros(1,nPsi);
Y = cell(1,nPsi);
TH = cell(1,nPsi);
T = cell(1,nPsi);

% perform averaging over orbits
for i = 1:nPsi
    [dPsi(i),DInt(i),Tp(i),Y{i},TH{i},T{i}] = kolmogorovAdiabaticDriftAbove(y0,th0(i),v0,alpha);
    if Tp(i) == 0
        disp(['Psi = ' num2str(psi(i))])
        disp('orbit never completed')
    end
end
% add on equilibrium values
% Tp(end) = 2*pi/sqrt(0.5*v0*(1-alpha));
% dPsi(end) = Tp(end) * v0*cos(pi)*(1-2*alpha+alpha*cos(2*pi))/(1-alpha*cos(2*pi)).^2;
% DInt(end) = Tp(end) *  0.5*(2*abs(v0*sin(pi)/(1 - alpha*cos(2*pi)))).^2;

% compute averaged drift and diffusion strength
psi = psi(Tp>0);
fPsi = dPsi(Tp>0)./Tp(Tp>0); % drift
D = DInt(Tp>0)./Tp(Tp > 0);
Tp = Tp(Tp>0);


% map onto full psi interval using symmetry
if psi(1) == 0 % don't double count psi0 if it remained after averaging code
    indFlip = 2;
else
    indFlip = 1;
end
psi = [-fliplr(psi(indFlip:end)) psi]; % odd
fPsi = [-fliplr(fPsi(indFlip:end)) fPsi]; % odd function
D = [fliplr(D(indFlip:end)) D]; % even function
Tp = [fliplr(Tp(indFlip:end)) Tp];

% compute unnormalized psi distribution
P_psi = exp(cumtrapz(psi,fPsi./D))./D;
% normalize
P_psi = P_psi/trapz(psi,P_psi);

    function th = thetaOrbit(y_,psi_)
        if alpha > 0
            th = acos(sqrt((1+alpha)/(2*alpha))*tanh(sqrt(2*alpha*(1+alpha))/(2*v0)*(cos(y_) - psi_)));
        elseif alpha == 0
            th = acos(1/(2*v0)*(cos(y_) - psi_));
        else
            th = acos(sqrt((1+alpha)/(2*abs(alpha)))*tan(sqrt(2*abs(alpha)*(1+alpha))/(2*v0)*(cos(y_) - psi_)));
        end
    end

end