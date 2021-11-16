% PSIKOLMOGOROV1 computes the conserved quantity of kolmogorov flow for
% alpha > 0.

function psi = psiKolmogorov1(y,th,v0,alpha)

psi = cos(y) - 2*v0*atanh( sqrt(2*alpha/(1+alpha)) * cos(th) ) / sqrt(2*alpha*(1+alpha));
 
end