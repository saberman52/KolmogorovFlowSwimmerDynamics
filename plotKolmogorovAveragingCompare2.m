clear
load('kolmogorov/paper_data/Pyth_alpha2')
[TH,Y] = meshgrid(theta,y);

figure

idx = [1,3,5];
iter = 1;

for i = idx
    % plot true distribution multiplied by 1 - alpha*cos(2*th)
    subplot(2,3,iter)
    imagesc(theta,y,P{i}.*(1-alpha(i)*cos(2*TH)))
    axis xy
    hold on
    plotKolmogorovSeparatrix(v0(i),alpha(i))
    if iter == 1
        ylabel('$y$',int,la)
    end
    colorbar('northoutside')
    % add text on plot
    text(0.1,2.5,[letters{iter} ' $\sigma = ' num2str(sigma(i)) ',\,\,v_0 = ' num2str(v0(i)) ',\,\,\alpha = ' num2str(alpha(i)) '$'],...
        fs,font,int,la,'Color',ones(1,3))
        set(gca,fs,font)
    drawnow;

    % plot distribution according to averaging multiplied by 1 - alpha*cos(2*th)
    [Pk,yk,thetak,psi,p0,g] = computeKolmogorovDensityAveraging(v0(i),alpha(i));
    [THk,Yk] = meshgrid(thetak,yk);
    subplot(2,3,iter+3)
    imagesc(thetak,yk,Pk.*(1-alpha(i)*cos(2*THk)))
    axis xy
    hold on
    plotKolmogorovSeparatrix(v0(i),alpha(i))
    xlabel('$\theta$',int,la)
    if iter == 1
        ylabel('$y$',int,la)
    end
    colorbar('northoutside')
    % add text on plot
    text(0.1,2.5,[letters{iter+3} ' $v_0 = ' num2str(v0(i)) ',\,\,\alpha = ' num2str(alpha(i)) '$'],...
        fs,font,int,la,'Color',ones(1,3))
    set(gca,fs,font)
    iter = iter+1;
end