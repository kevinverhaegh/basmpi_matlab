function PECv = TECPEC(N, Ne, Tev)

%THis function determines the total emission coefficients and the photon
%emission coefficients of recombination and excitation given the deuterium
%ADAS data

F=load([matlab_home,'/data_files/Yacora/H2_2/Yacora_H2_2'],'PEC','Te','nel');

Ne = Ne(:);
Tev = Tev(:);

Tev(Tev<min(F.Te)) = min(F.Te);
Tev(Tev>max(F.Te)) = max(F.Te);
Ne(Ne<min(F.nel)) = min(F.nel);
Ne(Ne>max(F.nel)) = max(F.nel);

TEC = zeros(numel(N), numel(Ne));
PECv = 0.*TEC;

for i=1:numel(N)
    PECv(i,:) = interp2L(F.nel, F.Te, squeeze(F.PEC(:,:,N(i)))', Ne, Tev);
end

function alpha = interp2L(x,y,V,xN,yN)

%performes 2D surface interpolation (spline) on a logaritmic grid

%F = griddedInterpolant({log10(x), log10(y)}, log10(V'), 'linear','none');
%alpha = 10.^F(log10(xN), log10(yN));
alpha = 10.^interp2(log10(x),log10(y),log10(V),log10(xN),log10(yN),'spline');