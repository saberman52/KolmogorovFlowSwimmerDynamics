% PSIKOLMOGOROV0 computes the conserved quantity of kolmogorov flow for
% alpha = 0.

function psi = psiKolmogorov0(y,th,v0,alpha)

psi = cos(y) - 2*v0*cos(th);
 
end