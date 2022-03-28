function [PECi] = TECPEC(N, Ne, TeP, TiHmP, TiHpP)

%THis function determines the total emission coefficients and the photon
%emission coefficients of recombination and excitation given the deuterium
%ADAS data

%uses updated Yacora data (Yacora 1.6.0); which uses different (fundamental) data for the
%H+ + H- reaction. Furthermore, it uses 'fittedmodel' to include the
%influence of TiHmP and TiHpP instead of a n-d interpolation

load([matlab_home,'/data_files/Yacora/HmHp_N_2/Yacora_HmHp_2'],'PEC','Te','nel','fittedmodel')
%TEC = zeros(numel(N), numel(Ne));
PECi = zeros(numel(N),numel(Ne(:)));
%PECexc = 0.*TEC;
TiHmP = TiHmP(:)*11600;
TiHpP = TiHpP(:)*11600;
Ne = Ne(:);
TeP = TeP(:);

%fittedmodel has been interrogated using a TiHmP and TiHpP from 1000 to
% 57900 K; which has maximum of 1 (at 1000, 1000 K) and minimum of 0.4 (at
% 57900, 57900 K). If extrapolative behaviour is used cap the minimum and
% maximum.
fittedmodel_max = 1;
fittedmodel_min = 0.4; 

% TiHmP(TiHmP<1e3) = 1e3;
% TiHpP(TiHpP<1e3) = 1e3;
% TiHmP(TiHmP<11600) = 11600;
% TiHpP(TiHpP<11600) = 11600;
% TiHmP(TiHmP>57900) = 57900;
% TiHpP(TiHpP>57900) = 57900;
% TiHmP(TiHmP<min(TiHm)) = min(TiHm);
% TiHpP(TiHpP<min(TiHp)) = min(TiHp);

TeP(TeP<min(Te)) = min(Te);
TeP(TeP>max(Te)) = max(Te);
Ne(Ne<min(nel)) = min(nel);
Ne(Ne>max(nel)) = max(nel);

%interrogate fittedmodel
if license('test', 'curveFitter')==1
    AdjF = fittedmodel(TiHmP,TiHpP);
    AdjF(AdjF<fittedmodel_min) = fittedmodel_min;
    AdjF(AdjF>fittedmodel_max) = fittedmodel_max;
else
    A = load([matlab_home,'/data_files/Yacora/HmHp_N_2/Yacora_HmHp_2'],'AdjF','TiHmP','TiHpP');
    [xx,yy] = meshgrid(A.TiHmP, A.TiHpP);
    AdjF = griddata(xx,yy,A.AdjF,TiHmP,TiHpP);
    AdjF(TiHmP>max(A.TiHmP,[],'omitnan')) = fittedmodel_min;
    AdjF(TiHpP>max(A.TiHpP,[],'omitnan')) = fittedmodel_min;
    AdjF(TiHmP<min(A.TiHmP,[],'omitnan')) = fittedmodel_max;
    AdjF(TiHpP<min(A.TiHpP,[],'omitnan')) = fittedmodel_max;
end

%for j=1:numel(Ne)
            for i=1:numel(N)
                PECi(i,:) = AdjF.*interp2L(nel, Te, reshape(PEC(:,:,N(i)), [18,9])', Ne, TeP);
           end
      %      CoolingdQdL(j) = sum(TEC(:,j).*(h.*c)./(lambdaBalmerExcitationTn2'.*1e-9)); %Radiated power per m
        %CoolingdQdL(j) = 0;
     %disp(100*j/numel(Ne))
%end

function alpha = interp2L(x,y,V,xN,yN)

%performes 2D surface interpolation (spline) on a logaritmic grid

alpha = 10.^interpn(log10(x), log10(y), log10(V), log10(xN), log10(yN), 'spline');
