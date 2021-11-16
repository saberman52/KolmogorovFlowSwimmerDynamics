% KOLMOGOROVFOURIERSOLVE solves the steady-state Fokker-Planck equation for
% the (y,theta) variables of a rotationally diffusing swimmer in the
% Kolmogorov flow. Solves a linear system of equations for the Fourier
% coefficients Pmn of the stationary distribution P(y,theta).
% 
% Inputs:
% v0 - swim speed
% alpha - swimmer shape parameter
% sigma - noise intensity
% M - maximum index of Fourier coefficients in y direction
% N - maximum index of Fourier coefficients in theta direction
%
% Outputs:
% Pmn - matrix of Fourier coefficients

function [Pmn,A,matTime,solveTime] = kolmogorovFourierSolve(v0,alpha,sigma,M,N)
    % number of unknowns
    Ntot = (2*M+1)*(2*N+1);
    % righthand side
    b = zeros(Ntot,1);
    b((2*M+1)*N+M+1) = 1/(4*pi^2); % (0,0) equation for P00, fixed by normalization
    % matrix of coefficients for each equation
    A = zeros(Ntot);
    
    mm = -M:M; % y fourier index (rows)
    nn = -N:N; % theta fourier index (columns)
    
    sz = [2*M+1,2*N+1]; % size of underlying matrices of unknowns and equations
    
    tic;
    % fill A matrix
    for i = 1:Ntot
        % get (m,n)
        [mInd,nInd] = ind2sub(sz,i);
        m = mm(mInd);
        n = nn(nInd);
        
        % For (m,n) = (0,0), do normalization condition
        if m == 0 && n == 0
            A(i,i) = 1;
        else
            % assign coefficients of equations, provided they aren't in
            % excluded range
            if m-1 >= -M
                % (m-1,n-2)
                if n-2 >= -N 
                    ind = sub2ind(sz,mInd-1,nInd-2);
                    A(i,ind) = alpha*n/8;
                end
                % (m-1,n)
                ind = sub2ind(sz,mInd-1,nInd);
                A(i,ind) = -n/4;
                % (m-1,n+2)
                if n+2 <= N
                    ind = sub2ind(sz,mInd-1,nInd+2);
                    A(i,ind) = alpha*n/8;
                end
            end
            % (m,n-1)
            if n-1 >= -N
                ind = sub2ind(sz,mInd,nInd-1);
                A(i,ind) = -v0*m/2;
            end
            %(m,n)
            A(i,i) = -sigma.^2*n.^2/2;
            % (m,n+1)
            if n+1 <= N
                ind = sub2ind(sz,mInd,nInd+1);
                A(i,ind) = v0*m/2;
            end
            if m+1 <= M
                % (m+1,n-2)
                if n-2 >= -N
                    ind = sub2ind(sz,mInd+1,nInd-2);
                    A(i,ind) = -alpha*n/8;
                end
                % (m+1,n)
                ind = sub2ind(sz,mInd+1,nInd);
                A(i,ind) = n/4;
                % (m+1,n+2)
                if n+2 <= N
                    ind = sub2ind(sz,mInd+1,nInd+2);
                    A(i,ind) = -alpha*n/8;
                end
            end
        end
    end
    matTime = toc;
    
    % solve
    tic;
    Pmn = A\b;
    solveTime = toc;
    % return as matrix
    Pmn = reshape(Pmn,sz);
end