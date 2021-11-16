% PSIKOLMOGOROV_ computes the conserved quantity of kolmogorov flow for
% alpha < 0.

function psi = psiKolmogorov_(y,th,v0,alpha)

psi = cos(y) - 2*v0*atan( sqrt(2*abs(alpha)/(1+alpha)) * cos(th) ) / sqrt(2*abs(alpha)*(1+alpha));
 
end