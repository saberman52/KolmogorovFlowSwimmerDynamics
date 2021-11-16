function plotKolmogorovSeparatrix(v0,alpha)
line = 'LineWidth';
lw = 2;
lineStyle = 'r:';
% set T function, s.t. Psi  = cos(y) - T(theta)
if alpha > 0
    psi_ = @psiKolmogorov1;
    T = @(th) 2*v0*atanh( sqrt(2*alpha/(1+alpha)) * cos(th) ) / sqrt(2*alpha*(1+alpha));
elseif alpha == 0
    psi_ = @psiKolmogorov0;
    T = @(th) 2*v0*cos(th);
else
    psi_ = @psiKolmogorov_;
    T = @(th) 2*v0*atan( sqrt(2*abs(alpha)/(1+alpha)) * cos(th) ) / sqrt(2*abs(alpha)*(1+alpha));
end

if v0 < kolmogorovBifurcation(alpha)
    theta = linspace(0,2*pi,1000);
    y = acos(psi_(0,0,v0,alpha) + T(theta));
    y2 = acos(psi_(pi,pi,v0,alpha) + T(theta));
    % island around central elliptic fixed point
    plot(theta,y,lineStyle,line,lw)
    hold on
    plot(theta,-y,lineStyle,line,lw)
    % other island
    plot(theta,y2,lineStyle,line,lw)
    hold on    
    plot(theta,-y2,lineStyle,line,lw)

else
    y = linspace(0,pi,500);
    y = y(1:end-1);
    psiStar = Psi(pi,pi,v0,alpha);
    if alpha > 0
        theta = acos(sqrt((1+alpha)/(2*alpha))*tanh(sqrt(2*alpha*(1+alpha))/(2*v0)*(cos(y) - psiStar)));
    elseif alpha == 0
        theta = acos(1/(2*v0)*(cos(y) - psiStar));
    else
        theta = acos(sqrt((1+alpha)/(2*abs(alpha)))*tan(sqrt(2*abs(alpha)*(1+alpha))/(2*v0)*(cos(y) - psiStar)));
    end
    % island around central fixed point
    plot(theta,y,lineStyle,line,lw)
    hold on
    plot(theta,-y,lineStyle,line,lw)
    plot(2*pi-theta,y,lineStyle,line,lw)
    plot(2*pi-theta,-y,lineStyle,line,lw)
    % other island
    plot(pi-theta,pi-y,lineStyle,line,lw)
    plot(pi-theta,y-pi,lineStyle,line,lw)
    plot(pi+theta,pi-y,lineStyle,line,lw)
    plot(pi+theta,y-pi,lineStyle,line,lw)


    
end

end