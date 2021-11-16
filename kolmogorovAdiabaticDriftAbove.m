% KOLMOGOROVADIABATICDRIFTABOVE computes the integral of the noise-induced drift
% of a swimmer with rotational diffusion in Kolmogorov flow around a 
% periodic orbit of the deterministic system.
% Also computes the diffusivity by averaging.
% Assumes parameters are above phase space bifurcation
%
% Inputs:
% (y0,th0) - initial condition
% v0 - swim speed
% alpha - swimmer shape parameter
%
% Outputs:
% dPsi - change in Psi (per unit noise intensity sigma^2) after one period
% DInt - integrated diffusivity 1/2 |d Psi/dtheta|^2 after one period
% Tp - period of orbit (zero if the orbit was never completed in time tf)
% Y - Y(t)
% TH - TH(t)
% T - time steps
function [dPsi,DInt,Tp,Y,TH,T] = kolmogorovAdiabaticDriftAbove(y0,th0,v0,alpha)
Q0 = [y0,th0,0,0];
tf = 1400; % maximum integration time
options = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events',@eventFun);

% variables for stopping condition
psi = Psi(y0,th0,v0,alpha);
psi00 = Psi(0,0,v0,alpha);

[T,Q,Tp,~,~] =  ode45(@velocity,[0,tf],Q0,options);

Y = Q(:,1);
TH = Q(:,2);
dPsi = Q(end,3);
DInt = Q(end,4);

if isempty(Tp)
    Tp = 0; % return zero if the event was never triggered
elseif abs(T(end) - tf) < 1e-2
    Tp = 0;
else
    Tp = Tp(end); % only want time of final event, which lead to termination
end
    
    
    % phase space velocity function, including dPsi evolution
    function dQ = velocity(t,Q)
        y = Q(1);
        th = Q(2);
        dQ = zeros(4,1);
        dQ(1) = v0*sin(th);
        dQ(2) = 0.5*sin(y)*(1-alpha*cos(2*th));
        dQ(3) = v0*cos(th)*(1-2*alpha+alpha*cos(2*th))/(1-alpha*cos(2*th)).^2;
        dQ(4) = 0.5*(2*abs(v0*sin(th)/(1 - alpha*cos(2*th)))).^2;
    end


    function [value,isterminal,direction] = eventFun(t,Q)
        yy = Q(1);
        insideSep = abs(psi) > abs(psi00); % condition to be inside separatrix
        if insideSep % if in the jet region
            value = yy-y0;
        else
            value = yy - (y0+2*pi);
        end
        direction = 1;
        if insideSep && t < 1
            isterminal = 0;
        else
            isterminal = 1;
        end
    end

end