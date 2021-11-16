function vals = Psi(y,th,v0,alpha)

if alpha > 0
    vals = cos(y) - 2*v0*atanh(sqrt(2*alpha/(1+alpha))*cos(th))/sqrt(2*alpha*(1+alpha));
elseif alpha == 0
    vals = cos(y) - 2*v0*cos(th);
else
    vals = cos(y) - 2*v0*atan(sqrt(2*abs(alpha)/(1+alpha))*cos(th))/sqrt(2*abs(alpha)*(1+alpha));
end

end