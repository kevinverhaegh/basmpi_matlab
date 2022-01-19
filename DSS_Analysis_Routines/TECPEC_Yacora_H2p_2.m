function PECv = TECPEC(N, Ne, Tev)

%THis function determines the total emission coefficients and the photon
%emission coefficients of recombination and excitation given the deuterium
%ADAS data

load([matlab_home,'/data_files/Yacora/H2p_2/Yacora_H2p_2'],'PEC','Te','nel')
Tev = Tev(:);
Ne = Ne(:);

Tev(Tev<min(Te)) = min(Te);
Tev(Tev>max(Te)) = max(Te);
Ne(Ne<min(nel)) = min(nel);
Ne(Ne>max(nel)) = max(nel);

TEC = zeros(numel(N), numel(Ne));
PECv = 0.*TEC;

for i=1:numel(N)
    PECv(i,:) = interp2L(nel, Te, squeeze(PEC(:,:,N(i)))', Ne, Tev);
end

function alpha = interp2L(x,y,V,xN,yN)

%performes 2D surface interpolation (spline) on a logaritmic grid

F = griddedInterpolant({log10(x), log10(y)}, log10(V'), 'spline');
alpha = 10.^F(log10(xN), log10(yN));
