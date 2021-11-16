% KOLMOGOROVBIFURCATION returns the critical speed vBif at which phase
% space bifurcates for swimmer with given alpha in Kolmogorov flow
function vBif = kolmogorovBifurcation(alpha)
vBif = zeros(size(alpha));

if ~isempty(vBif(alpha < 0))
    vBif(alpha < 0) = sqrt(2*abs(alpha(alpha < 0)).*(1+alpha(alpha < 0)))./(2*atan(sqrt(2*abs(alpha(alpha < 0))./(1+alpha(alpha < 0)))));
end
if ~isempty(vBif(alpha > 0))
    vBif(alpha > 0) = sqrt(2*alpha(alpha > 0).*(1+alpha(alpha > 0)))./(2*atanh(sqrt(2*alpha(alpha > 0)./(1+alpha(alpha > 0)))));
end
if ~isempty(vBif(alpha == 0))
    vBif(alpha == 0) = 0.5;
end

end