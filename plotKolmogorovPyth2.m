% script to plot phase space densities for \alpha > 0
clear
sigma = [0.05, 0.15*ones(1,5)];
v0 = [0.02,0.2,0.02,0.02,0.8,0.6];
alpha = [0.9,0.9,-0.7,0.9,0.9,-0.7];
n = length(sigma);

line = 'LineWidth';
lw = 3;
font = 22;
fs = 'FontSize';
int = 'Interpreter';
la = 'latex';

letters = {'(a)','(b)','(c)','(d)','(e)','(f)'};

N = 50; % number of Fourier coefficients in theta direction
% storage
Pmn = cell(1,n);
P = cell(1,n);


iter = 1;
figure
for i = 1:n
    if sigma(i) < 0.1 && v0(i) < 0.1
        M = 100;
    else
        M = 50;
    end
    
    [Pmn{i},~,~,~] = kolmogorovFourierSolve(v0(i),alpha(i),sigma(i),M,N);
    [P{i},y,theta] = kolmogorovIFFTYshift(Pmn{i});
    
    subplot(2,3,i)
    % plot density
    imagesc(theta,y,P{i})
    axis xy
    hold on
    plotKolmogorovSeparatrix(v0(i),alpha(i))
    if i > 3
        xlabel('$\theta$',int,la)
    end
    if i == 1 || i == 4
        ylabel('$y$',int,la)
    end
    colorbar('northoutside')
    % add text on plot
    text(0.1,2.5,[letters{i} ' $\sigma = ' num2str(sigma(i)) ',\,\,v_0 = ' num2str(v0(i)) ',\,\,\alpha = ' num2str(alpha(i)) '$'],...
        fs,font,int,la,'Color',ones(1,3))
    %         colorbar
    set(gca,fs,font)
    
    drawnow;
end

save('kolmogorov/paper_data/Pyth_alpha2')