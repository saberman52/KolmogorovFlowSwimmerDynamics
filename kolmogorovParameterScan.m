clear
fs = 'FontSize';
font = 18;
int = 'Interpreter';
la = 'latex';

% run Fourier solve over all possible combinations in (v0,alpha,sigma)
% space
v0 = 0.01:0.03:0.61;
alpha = [-0.9:0.1:-0.1, 0.1:0.1:0.9];
sigma = 0.05:0.0125:0.25;
M = 50;
N = 50;
nv = length(v0);
na = length(alpha);
ns = length(sigma);
% create meshgrid of alpha x v0 x sigma
[ALPHA,V0,SIGMA] = meshgrid(alpha,v0,sigma);
Pmn = cell(nv,na,ns); % cell to store all Fourier coefficients
solveTime = zeros(nv,na,ns); % store linear solve times

for i = 1:nv
    disp(['i = ' num2str(i)])
    for j = 1:na
        disp(['j = ' num2str(j)])
        for k = 1:ns
            disp(['k = ' num2str(k)])
            [Pmn{i,j,k},~,~,solveTime(i,j,k)] = kolmogorovFourierSolve(v0(i),alpha(j),sigma(k),M,N);
        end
    end
end

save('kolmogorov/fourier_parameterscan2')
