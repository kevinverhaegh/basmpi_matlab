function input = DSS_Analysis_Calc_RecIonFRec(input)

%Version 4 - final matlab version

%Function to determine FRec, recombination and ionisation rates. Both the
%main calculation and MC analysis is done

%% Disable auto parpool creation
%  ps = parallel.Settings;
%   ps.Pool.AutoCreate = false;
%  delete(gcp('nocreate')); %close off parpool

%% Analyse input, get necessary info from input structure
username=char(java.lang.System.getProperty('user.name'));

disp('Analysing input....')

floc = ['/home/kver/DSS_Analysis/', num2str(input.shot), '_', num2str(input.n1), num2str(input.n2), '.mat' ];
save(floc, 'input', '-v7.3');

n1 = input.n1; %First balmer line index
n2 = input.n2; %second Balmer line index
if n1>n2
    error('n2 should correspond to the higher Balmer line -> please reverse')
end
n1Int = input.n1Int; %First Balmer line intensity
if ~isfield(input, 'n1IntErr')
    input.n1IntErr = input.n1Int.*input.AbsErr;
end
n1IntErr = input.n1IntErr; %1-sigma error First Balmer line intensity (random)
n2Int = input.n2Int; %second Balmer line intensity
RelErr = input.RelErr;
none = input.none; %neutral fraction
noneL = input.noneL; %lower limit neutral fraction
noneH = input.noneH; %upper limit neutral fraction
Den = input.Den; %electron density
DenErr = input.DenErr; %random error electron density
DL = input.DL; %Divertor leg length
DLL = input.DLL; %lower estimate divertor leg length
DLH = input.DLH; %upper estimate divertor leg length
R = input.R; %radius 
Z = input.Z; %z-position
Time = input.Time; %time vector
shot = input.shot; %shot number

%% Get optional info from input structure

%number of random parameters
if isfield(input, 'Iter')
    Iter = input.Iter;
else
    Iter = 500;
end

%ADAS
if isfield(input, 'ADASInfo')
    ADASInfo = input.ADASInfo;
else
    ADASInfo = load([matlab_home,'/data_files/DeuteriumADAS']);
end

%Minimum density
if isfield(input, 'DenMin')
    DenMin = input.DenMin;
else
    DenMin = 5e18;
end

%Correction algorithm for Frec
if isfield(input, 'FRecForceTrend')
    DoFRecForceTrend = input.FRecForceTrend;
else
    DoFRecForceTrend = 1;
end

%Te vector used for FRec determination
if isfield(input, 'TeRestrict')
    TeRestrict = input.TeRestrict;
else
    TeRestrict = 1;
end

%Deprecated
if isfield(input, 'ExtremaSol')
    ExtremaSol = input.ExtremaSol;
else
    ExtremaSol = 0;
end

%Quantile probability calculation (Deprecated)
if isfield(input, 'Quant')
    Quant = input.Quant;
else
    Quant = linspace(0, 1, 100);
end

%Parallel execution yes/no
% if isfield(input, 'DoParFor')
%     if input.DoParFor == 1
%         parpool(28)
%     end
% end

%Correction algorithm for increase in TeE
if ~isfield(input, 'TeE_filter')
    input.TeE_filter = 1;
end

%Time at which the correction algorithm starts
if ~isfield(input, 'TeE_filter')
    if ~isfield(input, 'TeE_filterTime')
        input.TeE_filterTime = 1;
    end
end

%rejects line ratios which deviate too much from atomic data
if ~isfield(input,'FRec_reject')
    input.FRec_reject = 1.2;
end

FRec_reject = input.FRec_reject;

%select molecular emission model
if ~isfield(input,'MolEmissModel')
    input.MolEmissModel = 'YACORA';
end

MolEmissModel = input.MolEmissModel;

%Log-uniform distribution (as opposed to the default uniform) for no/ne
if isfield(input,'none_loguniform')
    none_loguniform = input.none_loguniform;
else
    none_loguniform = 1;
end

%regenration of NaNs during iteration
if isfield(input,'RegenNaNs')
    RegenNans = input.RegenNans;
else
    RegenNans = 0;
end

if ~isfield(input,'pMolG')
    input.pMolG = [-1.6876   17.2156];
end
pMolG = input.pMolG;

if ~isfield(input,'THm')
    input.THm = repmat(reshape(0.5.*rand(Iter,1), [1, 1, Iter]), [size(input.n1Int),1]);
end
THm = input.THm;

if ~isfield(input,'n1MolMCFilter')
    n1MolMCFilter = 1;
else
    n1MolMCFilter = input.n1MolMCFilter;
end

if ~isfield(input,'D2pDm_model')
    D2pDm_model = 2;
else
    D2pDm_model = input.D2pDm_model;
end

%Inclusion of uncertainties in the ADAS & YACORA coefficients
if isfield(input, 'CoeffUncertainty')
    CoeffUncertainty = input.CoeffUncertainty;
    if ~isfield(input,'Unc')
    %random vector for rate uncertainty
    input.Unc.U_R1 = abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_E1 = abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_R2 = abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_E2 =abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_DAR =abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_DAE =abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_DBR =abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_DBE =abs(1 + CoeffUncertainty.*(rand(Iter,1)-0.5));
    for i=1:4
    input.Unc.U_D2(i,:)=abs(1 + 2.*CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_D2p(i,:)=abs(1 + 2.*CoeffUncertainty.*(rand(Iter,1)-0.5));
    input.Unc.U_Dm(i,:)=abs(1 + 2.*CoeffUncertainty.*(rand(Iter,1)-0.5));
    end
    end
else
    if ~isfield(input,'Unc')
    %random vector for rate uncertainty
    input.Unc.U_R1 = ones(Iter,1);
    input.Unc.U_E1 = ones(Iter,1);
    input.Unc.U_R2 = ones(Iter,1);
    input.Unc.U_E2 =ones(Iter,1);
    input.Unc.U_DAR = ones(Iter,1);
    input.Unc.U_DAE = ones(Iter,1);
    input.Unc.U_DBR = ones(Iter,1);
    input.Unc.U_DBE =ones(Iter,1);
    for i=1:4
    input.Unc.U_D2(i,:)=ones(Iter,1);
    input.Unc.U_D2p(i,:)=ones(Iter,1);
    input.Unc.U_Dm(i,:)=ones(Iter,1);
    end
    end
end

U_E1 = input.Unc.U_E1;
U_E2 = input.Unc.U_E2;
U_R1 = input.Unc.U_R1;
U_R2 = input.Unc.U_R2;
U_DAR = input.Unc.U_DAR;
U_DAE = input.Unc.U_DAE;
U_DBR = input.Unc.U_DBR;
U_DBE = input.Unc.U_DBE;
U_D2 = input.Unc.U_D2;
U_D2p = input.Unc.U_D2p;
U_Dm = input.Unc.U_Dm;

DSS_Aux_structvars(input.Unc);

X.input = input;
input = X;
%% Initialise random values
if ~isfield(input.input,'inputMC')
disp('Intialising matrices and random values....')

%make empty cells

n1IntMC = zeros([size(n1Int),Iter])+NaN;
n2IntMC = zeros([size(n1Int),Iter])+NaN;
noneMC = zeros([size(n1Int),Iter])+NaN;
DenMC = zeros([size(n1Int), Iter]) + NaN;
DLMC = zeros([size(n1Int),Iter])+NaN;
TiTeR = zeros([size(n1Int),Iter])+NaN;

%initialise random vectors

noneRand = rand(Iter,1); %random vector none
randn_DenMC = randn(Iter, 1);
randn_TiTeR = randn(Iter, 1);
randn_DL = randn(Iter,1);

if ~isfield(input.input,'DaMea')
N = 2;
else
    N=4;
end
%sample multivariate distribution
u1 = 1;
u2 = 1;
s = nanmean(n1IntErr(:)./n1Int(:));
r = RelErr;

correl_coeff = 0.5*((u1/u2)^2 + 1 - (r*u2/s)^2)/(u1/u2); 

covariance = ((s.^2).*ones(N,N)).*(correl_coeff*(~eye(N)) + eye(N)*1);

R = mvnrnd(ones(N,1),covariance,Iter);

%adjust randn_DenMC
%% utilize random vectors to assign the MC values for the density using rejection samplign to account for minimum density
for j=1:numel(n1Int(1,:)) 
    for i=1:numel(n1Int(:,1))
        if ~isnan(Den(i,j))&& ~isnan(none(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(DL(i,j)) %parameter check, if parameters not valid return NaN

            DenMC(i,j,:) = Den(i,j) + DenErr(i,j).*randn_DenMC(:);
            NOK = 1;
            while NOK==1
               
                Slice = squeeze(DenMC(i,j,:));
                V = find(Slice<DenMin);
                if ~isempty(V)
                Vr = randn(numel(V),1);
                Vn = Den(i,j) + DenErr(i,j).*Vr; 
                randn_DenMC(V) = Vr;
                DenMC(i,j,(DenMC(i,j,:)<DenMin)) = Vn;
                else
                    NOK=0;
                end
            end

        end
    end
end
%DenMC = zeros([size(n1Int),Iter])+NaN;
%% utilize random vectors to assign all the other MC values
for j=1:numel(n1Int(1,:)) 
    for i=1:numel(n1Int(:,1))
        if ~isnan(Den(i,j))&& ~isnan(none(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(DL(i,j)) %parameter check, if parameters not valid return NaN
            n1IntMC(i,j,:) = n1Int(i,j).*R(:,1);
            n2IntMC(i,j,:) = n2Int(i,j).*R(:,2);

            if ~none_loguniform
            noneMC(i,j,:)  = (noneH(i,j) - noneL(i,j))*noneRand + noneL(i,j);
            else
            noneMC(i,j,:) = DSS_Aux_LogUniform_Sample(noneRand, noneL(i,j), noneH(i,j));     
            end
            %DenMC(i,j,:) = Den(i,j) + DenErr(i,j).*randn_DenMC(:);
            
            TiTeR(i,j,:) = (1.5-0.8).*randn_TiTeR + 0.8; %assume that the Ti is random between 0.8 and 1.5 times Te
            
            randV = randn_DL(:);
            randV(randV<0) = (DL(i,j) - DLL(i,j)).*randV(randV<0);
            randV(randV>0) = (DLH(i,j) - DL(i,j)).*randV(randV>0);
            DLMC(i,j,:) = DL(i,j) + randV;
        end
    end
end

n1IntTotMC = n1IntMC.*(1/1.0);
n2IntTotMC = n2IntMC.*(1/1.0);

else
    DSS_Aux_structvars(input.input.inputMC);
end
input.inputMC = DSS_Aux_structvars({'fieldNames','n1IntMC','DLMC','n2IntMC','noneMC','DenMC','n1IntTotMC','n2IntTotMC','TiTeR'});


%% Make empty matrices for storing calculation results

%monte carlo solutions
FRec1MC = 0.*n1IntMC + NaN;
FRec2MC = 0.*n1IntMC + NaN;
%TeEMC1 = 0.*n1IntMC + NaN;
TeEMC1 = 0.*n1IntMC + NaN;
TeRMC2 = 0.*n1IntMC + NaN;
dimen = size(FRec1MC);
%% Make interpolation objects for atomic data
n1Exc = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerExcitation(n1-2,:,:))'));
n2Exc = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerExcitation(n2-2,:,:))'));
n1Rec = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerRecombination(n1-2,:,:))'));
n2Rec = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerRecombination(n2-2,:,:))'));

%% Start FRec calculation

Te = 0.5:0.01:100; 
disp('Starting FRec calculation....')
tic
parfor i=1:dimen(1)
    FRec1MCt = zeros(dimen(2), numel(DenMC(1,1,:)));
    FRec2MCt = FRec1MCt;
    [ii] = find(nansum(n1Int,1)>0);
    
    DenMCs = reshape(DenMC(i,:,:),[numel(n1Int(1,:)) Iter]);
    noneMCs = reshape(noneMC(i,:,:), [numel(n1Int(1,:)) Iter]);
    LRs = reshape(n2IntMC(i,:,:)./n1IntMC(i,:,:), [numel(n1Int(1,:)) Iter]);
    for j=ii
        disp([i j])
        PE2 = repmat(U_E2,[1 numel(Te)]).*10.^n2Exc(log10(repmat(squeeze(DenMCs(j,:)'),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMCs(j,:)')) 1])));
        PR2 = repmat(U_R2,[1 numel(Te)]).*10.^n2Rec(log10(repmat(squeeze(DenMCs(j,:)'),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMCs(j,:)')) 1])));
        PE1 = repmat(U_E1,[1 numel(Te)]).*10.^n1Exc(log10(repmat(squeeze(DenMCs(j,:)'),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMCs(j,:)')) 1])));
        PR1 = repmat(U_R1,[1 numel(Te)]).*10.^n1Rec(log10(repmat(squeeze(DenMCs(j,:)'),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMCs(j,:)')) 1])));
        LR = (squeeze(noneMCs(j,:)') .* PE2 + PR2) ./ (squeeze(noneMCs(j,:)') .* PE1 + PR1);
        FRec1 = PR1./(squeeze(noneMCs(j,:)')*ones(size(Te)).*PE1 + PR1);
        FRec2 = PR2./(squeeze(noneMCs(j,:)')*ones(size(Te)).*PE2 + PR2);
        %clear PE2 PR2 PE1 PR1
        [Mi,Ii] = nanmin(LR,[],2);
        [Ma,Ia] = nanmax(LR,[],2);
        FRec1MCtt = zeros(numel(noneMCs(j,:)),1)+NaN;
        FRec2MCtt = zeros(numel(noneMCs(j,:)),1)+NaN;
        for k=1:numel(noneMCs(j,:))
            I = Ia(k):Ii(k);
            if ~isnan(LRs(j,k))
            FRec1MCtt(k) = interp1(LR(k,I),FRec1(k,I),LRs(j,k));
            FRec2MCtt(k) = interp1(LR(k,I),FRec2(k,I),LRs(j,k));
            end
        end
        Ii2 = (squeeze((LRs(j,:)')<Mi) & squeeze((LRs(j,:)')>Mi/FRec_reject));
        Ia2 = (squeeze((LRs(j,:)')>Ma) & squeeze((LRs(j,:)')<Ma.*FRec_reject));
        P = [find(Ii2) Ii(Ii2)];
        %FRec1MCt = zeros(numel(P(:,1)));
        %FRec2MCt = zeros(numel(P(:,1)));
        for k=1:numel(P(:,1))
        ki = P(k,1);
        FRec1MCtt(ki) = rand(size(Mi(P(k,1)))).*FRec1(P(k,1), P(k,2));
        FRec2MCtt(ki) = rand(size(Mi(P(k,1)))).*FRec2(P(k,1), P(k,2));
        end
        P = [find(Ia2) Ia(Ia2)];
        for k=1:numel(P(:,1))
            ki = P(k,1);
        FRec1MCtt(ki) = rand(size(Ma(P(k,1)))).*(1-FRec1(P(k,1),P(k,2))) + FRec1(P(k,1),P(k,2));
        FRec2MCtt(ki) = rand(size(Ma(P(k,1)))).*(1-FRec2(P(k,1),P(k,2))) + FRec2(P(k,1),P(k,2));
        end
        FRec1MCt(j,:) = FRec1MCtt;
        FRec2MCt(j,:) = FRec2MCtt;
    end
    FRec1MC(i,:,:) = FRec1MCt;
    FRec2MC(i,:,:) = FRec2MCt;
end

tp = toc;
disp(['FRec calculation finished - ', num2str(tp), ' seconds elapsed'])
%% Perform post-processing on FRec calculation to clean-up the FRec vaectors
EVal = ones(size(FRec1MC));
EVal(isnan(FRec1MC))=0;
FRecF = EVal;
%first force FRec trend - if requested (correction algorithm to ensure
%correct FRec values at high FRec.
if DoFRecForceTrend
    disp('Performing post-processing on FRec - enforcing increasing FRec trend....')
    %makie backups of calculated data
    FRec1MC_o = FRec1MC;
    FRec2MC_o = FRec2MC;
    %Te_o = Te;
    for i=1:numel(n1Int(:,1))
        for j=1:Iter
            FRec1MC(i,:,j) = FRecForceTrend(FRec1MC(i,:,j)); %apply post-process FRec
            FRec2MC(i,:,j) = FRecForceTrend(FRec2MC(i,:,j));
            %TeMC(i,:,j) = FRecForceTrend(TeMC(i,:,j));
        end
        %deprecated
    end
    disp('FRec post-processing finished')
    
    %deprecated
%  FRec1 = DSS_Aux_smooth2a(FRec1, 2,2); %smooth post-processed FRec
%  FRec1 = imadjust(FRec1, [min(FRec1(:)) max(FRec1(:))], [min(FRec1_o(:)) max(FRec1_o(:))]); %adjust post-smoothed FRec to realign with previous minima/maxima
% % 
%  FRec2 = DSS_Aux_smooth2a(FRec2, 2,2);
%  FRec2 = imadjust(FRec2, [min(FRec2(:)) max(FRec2(:))], [min(FRec2_o(:)) max(FRec2_o(:))]); %adjust post-smoothed FRec to realign with previous minima/maxima
% % 
%  for i=1:Iter %now apply smoothing of post-processed FRec on the Monte Carlo execution grid
%      MiP = min(min(min(FRec1MC(:,:,i))));
%      MaP = max(max(max(FRec1MC(:,:,i))));
%      
%      FRec1MC(:,:,i) = FRec1MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
%      
%      FRec1MC(:,:,i) =  DSS_Aux_smooth2a(squeeze(FRec1MC(:,:,i)), 2,2);   %Apply smoothing
%      FRec1MC(:,:,i) = imadjust(squeeze(FRec1MC(:,:,i)), [min(min(min(FRec1MC(:,:,i)))) max(max(max(FRec1MC(:,:,i))))], [MiP MaP]); %ajust smoothed FRec to make it align with previous minima/maxima
%      
%      FRec1MC(:,:,i) = FRec1MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
%      
%      MiP = min(min(min(FRec2MC(:,:,i))));
%      MaP = max(max(max(FRec2MC(:,:,i))));
%      
%      FRec2MC(:,:,i) = FRec2MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
%      
%      FRec2MC(:,:,i) =  DSS_Aux_smooth2a(squeeze(FRec2MC(:,:,i)), 2,2);  %Apply smoothing 
%      FRec2MC(:,:,i) = imadjust(squeeze(FRec2MC(:,:,i)), [min(min(min(FRec2MC(:,:,i)))) max(max(max(FRec2MC(:,:,i))))], [MiP MaP]); %ajust smoothed FRec to make it align with previous minima/maxima
%     
%      FRec2MC(:,:,i) = FRec2MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
%    
%  end
 
    
end
% 

%% Perform excitation/reombination temperature calculation
%make general input settings
disp('Starting recombination / ionisation rate calculations for both Balmer lines...')
tic
looper=0;

parfor i=1:numel(n1Int(:,1))
    TeEMC1t = zeros(numel(n1Int(i,:)), Iter)+NaN;
    %TeEMC2t = TeEMC1t;
    TeRMC2t = TeEMC1t;    
    [ii] = find(nansum(n1Int,1)>0);
    DenMCs = reshape(DenMC(i,:,:),[numel(n1Int(1,:)) Iter]);
    DLMCs = reshape(DLMC(i,:,:),[numel(n1Int(1,:)) Iter]);
    noneMCs = reshape(noneMC(i,:,:), [numel(n1Int(1,:)) Iter]);
    n1IntMCs = reshape(n1IntMC(i,:,:), [numel(n1Int(1,:)) Iter]);
    n2IntMCs = reshape(n2IntMC(i,:,:), [numel(n1Int(1,:)) Iter]);
    FRec1MCs = reshape(FRec1MC(i,:,:), [numel(n1Int(1,:)) Iter]);
    FRec2MCs = reshape(FRec2MC(i,:,:), [numel(n1Int(1,:)) Iter]);
    for j=ii
        if ~isnan(n1Int(i,j).*Den(i,j))
        disp([i j])
        %PE2 = 10.^n2Exc(log10(repmat(squeeze(DenMC(i,j,:)),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMC(i,j,:))) 1])));
        PR2 = repmat(U_R2,[1 numel(Te)]).*10.^n2Rec(log10(repmat(squeeze(DenMCs(j,:)),[numel(Te) 1]))', log10(repmat(Te', [size(squeeze(DenMCs(j,:))) 1]))');
        PE1 = repmat(U_E1,[1 numel(Te)]).*10.^n1Exc(log10(repmat(squeeze(DenMCs(j,:)),[numel(Te) 1]))', log10(repmat(Te', [size(squeeze(DenMCs(j,:))) 1]))');
        %PR1 = 10.^n1Rec(log10(repmat(squeeze(DenMC(i,j,:)),[1 numel(Te)])), log10(repmat(reshape(Te, [1 numel(Te)]), [size(squeeze(DenMC(i,j,:))) 1])));
            for k=1:Iter
                if ~isnan(FRec2MCs(j,k))
                TeEMC1t(j,k) = interp1(DLMCs(j,k).*noneMCs(j,k).*1e-6.*PE1(k,:).*(DenMCs(j,k).^2),Te,n1IntMCs(j,k).*(1-FRec1MCs(j,k)));
                %TeEMC2t(j,k) = interp1(DLMCs(j,k).*noneMCs(j,k).*1e-6*PE2(k,:).*(DenMCs(j,k).^2),Te,n2IntMCs(j,k).*(1-FRec2MCs(j,k)));
                TeRMC2t(j,k) = interp1(DLMCs(j,k).*1e-6.*PR2(k,:).*(DenMCs(j,k).^2),Te,n2IntMCs(j,k).*FRec2MCs(j,k));
                end
            end   
        end
    end
    TeEMC1(i,:,:) = TeEMC1t;
    %TeEMC2(i,:,:) = TeEMC2t;
    TeRMC2(i,:,:) = TeRMC2t;
end

tp = toc;
disp(['Rec/Ion calculation finished - ', num2str(tp), ' seconds elapsed'])

%% employ TeE filter
if input.input.TeE_filter
    disp('TeE data filter selected, filtering out erroneous TeE points, assuming a continuous decrease of TeE')
    
    if isfield(input.input, 'TeE_filterTime')
        tbeg = input.input.TeE_filterTime; 
        else
        tbeg=0.8; %beginning of TeE checks
    end
    %filter assumes at the end of the discharge a decrease in TeE
    TeEMC1S = TeEMC1;

    %EVal = zeros(size(TeEMC1))+1;
    
    for i=1:numel(n1Int(:,1)) %for each ROI
    %check if TeE1 is decreasing
        
    for j=1:numel(Iter)
    TeEMC1S(i,:,j) = DSS_Aux_smooth2a(squeeze(TeEMC1(i,:,j)),4); %smooth TeE to obtain rough TeE trend
    end
    
    V = diff(squeeze(TeEMC1(i,:,:)))>0; %Check where TeE increases
    V(input.input.Time < tbeg,:) = 0;
    V(2:end+1,:) = V;
    TeEMC1(i,V) = NaN;
    TeRMC2(i,V) = NaN;
    %turn results to NaN
    EVal(i,V==1) = 0;
    end 
end

input.ResultMC = DSS_Aux_structvars({'fieldNames','TeEMC1','TeRMC2','EVal','FRec1MC','FRec2MC','FRecF'});

disp('Filter done! Erroneous TeE entries overwritten with NaN')
  
%% employ molecular analysis

if isfield(input.input,'DaMea')
    
   %Get excitation/recombination emission
    BIntExc_n1 = EVal.*(1-FRec1MC).*n1IntMC;
    BIntRec_n2 = EVal.*(FRec2MC).*n2IntMC;
    
    naExc = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerExcitation(1,:,:))'));
    naRec = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerRecombination(1,:,:))'));
    
    LRaExc = repmat(reshape(U_DAE,[1 1 Iter]),[size(n1Int) 1]).*(10.^naExc(log10(DenMC.*EVal), log10(TeEMC1.*EVal)))./(repmat(reshape(U_E1,[1 1 Iter]),[size(n1Int) 1]).*(10.^n1Exc(log10(DenMC.*EVal), log10(TeEMC1.*EVal))));
    LRaRec = repmat(reshape(U_DAR,[1 1 Iter]),[size(n1Int) 1]).*(10.^naRec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)))./(repmat(reshape(U_R2,[1 1 Iter]),[size(n1Int) 1]).*(10.^n2Rec(log10(DenMC.*EVal), log10(TeRMC2.*EVal))));
    
    %extrapolate atomic to |Da
    BIntAtom_Da = EVal.*(BIntExc_n1.*LRaExc + BIntRec_n2.*LRaRec);
    if ~isfield(input.input,'inputMC')
    DaTot = input.input.DaMea; %Get total Da - correct for window
    %DaTot = input.DaMea;

    %Process uncertainties
    ind = find(~isnan(sum(n1Int)) & sum(n1Int)>0,1);
    AbsErr = n1Int(:,ind)*ones(1,Iter)./squeeze(n1IntTotMC(:,ind,:)); %Get absolute uncertainty
    %RelErr = ones(size(input.n1Int(:,ind)))*(1+randn(1,input.Iter).*DaRelErr); %relative uncertainty
    DaTotMC = zeros(size(n1IntMC))+NaN;
    for i=1:numel(input.input.DaMea(:,1))
    for j=1:numel(input.input.DaMea(1,:)) %apply uncertainties
        DaTotMC(i,j,:) = DaTot(i,j).*R(:,3);
        %DaTotMC(j,i,:) = DaTot(j,i).*squeeze(n1IntTotMC(j,i,:)./n1Int(j,i)).*(1 + RelErr.*randn_DaMea);%.*RelErr(j,:);
    end
    end
    else
        DaTotMC = input.input.inputMC.DaTotMC;
    end
    input.inputMC.DaTotMC = DaTotMC;
    %DaMin = minDa.*rand(size(DaTotMC)).*DaTotMC;

    %Get Molecular Da
%     DaMolMC = DaTotMC - BIntAtom_Da;
%     DaMolMC(DaMolMC<0) = 0;

    if isfield(input.input,'DbMea')
        nbExc = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerExcitation(2,:,:))'));
        nbRec = griddedInterpolant({log10(ADASInfo.AdasNe), log10(ADASInfo.AdasTe)}, log10(squeeze(ADASInfo.BalmerRecombination(2,:,:))'));

        LRbExc = repmat(reshape(U_DBE,[1 1 Iter]),[size(n1Int) 1]).*(10.^nbExc(log10(DenMC.*EVal), log10(TeEMC1.*EVal)))./(repmat(reshape(U_E1,[1 1 Iter]),[size(n1Int) 1]).*(10.^n1Exc(log10(DenMC.*EVal), log10(TeEMC1.*EVal))));
        LRbRec = repmat(reshape(U_DBR,[1 1 Iter]),[size(n1Int) 1]).*(10.^nbRec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)))./(repmat(reshape(U_R2,[1 1 Iter]),[size(n1Int) 1]).*(10.^n2Rec(log10(DenMC.*EVal), log10(TeRMC2.*EVal))));
        %extrapolate atomic to Db
        BIntAtom_Db = EVal.*(BIntExc_n1.*LRbExc + BIntRec_n2.*LRbRec);
        
        if ~isfield(input.input,'inputMC') 
        DbTot = input.input.DbMea; %Get total Db - correct for window
       
        %Process uncertainties - Db
        ind = find(~isnan(sum(n1Int)) & sum(n1Int)>0,1);
        AbsErr = n1Int(:,ind)*ones(1,Iter)./squeeze(n1IntTotMC(:,ind,:)); %Get absolute uncertainty
        %RelErr = (n2Int(:,ind)./n1Int(:,ind))*ones(1,Iter)./squeeze(n2IntMC(:,ind,:)./n1IntMC(:,ind,:)); %relative uncertainty
        DbTotMC = zeros(size(n1IntMC))+NaN;
        for i=1:numel(input.input.DaMea(:,1))
        for j=1:numel(input.input.DaMea(1,:)) %apply uncertainties
            DbTotMC(i,j,:) = DbTot(i,j).*R(:,4);
            %DbTotMC(j,i,:) = DbTot(j,i).*AbsErr(j,:).*(1 + RelErr.*randn_DbMea);%.*RelErr(j,:);
        end
        end
        else
            DbTotMC = input.input.inputMC.DbTotMC;
        end
        
        input.inputMC.DbTotMC = DbTotMC;
        
        %Get Molecular Da
%         DbMolMC = DbTotMC - BIntAtom_Db;
%         DbMolMC(DbMolMC<0) = 0;
    end
    if strcmp(MolEmissModel,'YACORA')
    %interrogate PECs
    PEC_D2T = TECPEC_Yacora_H2_2( [2 3 input.input.n1-1 input.input.n2-1],DenMC(:).*EVal(:), TeEMC1(:).*EVal(:));
    PEC_D2T = reshape(PEC_D2T, [4 numel(DenMC(:))]);
    PEC_D2pT = TECPEC_Yacora_H2p_2([2 3 input.input.n1-1 input.input.n2-1],DenMC(:).*EVal(:), TeEMC1(:).*EVal(:));%rand(size(TeEMC1(:))).*EVal(:));
    PEC_D2pT = reshape(PEC_D2pT, [4 numel(DenMC(:))]);
    PEC_DmT = TECPEC_Yacora_HmHp_2([2 3 input.input.n1-1 input.input.n2-1],DenMC(:).*EVal(:), TeEMC1(:).*EVal(:),THm(:).*EVal(:),TiTeR(:).*TeEMC1(:).*EVal(:));
    PEC_DmT = reshape(PEC_DmT, [4 numel(DenMC(:))]);
        
    %separate Da emission channels
    elseif strcmp(MolEmissModel,'AMJUEL')
        
        %load('/home/verhaegh/Public files/get_amjuel')
        AMJ = get_amjuel();
        Nv=[2 3 input.input.n1-1 input.input.n2-1];
        PEC_D2T = zeros(numel(Nv),numel(DenMC(:)));
        PEC_DmT = PEC_D2T;
        PEC_D2pT = PEC_D2T;
        %get AMJUEL emission rates
        indexA = 'bacde'; %array of indexes for 1:6
        Ai = [4.7e8 4.4114e7 8.4217e6 2.5311e6 9.73465e5]; %Einstein coefficients for 1:6
        for i=1:numel(Nv)
            if Nv(i)>numel(indexA)
                continue
            end
        PEC_D2pT(i,:) = Ai(Nv(i)).*eval(['amjuel_tables(''h12'',AMJ.H12_2_2_14',indexA(Nv(i)),'.table,DenMC(:).*EVal(:),TeEMC1(:).*EVal(:))'])./(DenMC(:).*EVal(:));
        PEC_DmT(i,:) = Ai(Nv(i)).*eval(['amjuel_tables(''h12'',AMJ.H12_7_2',indexA(Nv(i)),'.table,DenMC(:).*EVal(:),TeEMC1(:).*EVal(:))'])./(DenMC(:).*EVal(:));
        PEC_D2T(i,:) = Ai(Nv(i)).*eval(['amjuel_tables(''h12'',AMJ.H12_2_2_5',indexA(Nv(i)),'.table,DenMC(:).*EVal(:),TeEMC1(:).*EVal(:))'])./(DenMC(:).*EVal(:));
        %PEC_D3p(i,:) = Ai(i).*eval(['amjuel_tables(''h12'',AMJ.H12_2_2_15',indexA(i),'.table,ne(:),Te(:))'])./ne(:);
        end

    end
    
    for i=1:4
        dum = repmat(reshape(U_D2(i,:),[1 1 Iter]),[size(n1Int) 1]);
        PEC_D2T(i,:) = dum(:)'.*PEC_D2T(i,:);
        dum = repmat(reshape(U_D2p(i,:),[1 1 Iter]),[size(n1Int) 1]);
        PEC_D2pT(i,:) = dum(:)'.*PEC_D2pT(i,:);
        dum = repmat(reshape(U_Dm(i,:),[1 1 Iter]),[size(n1Int) 1]);
        PEC_DmT(i,:) = dum(:)'.*PEC_DmT(i,:);
    end
    
    PEC_D2 = PEC_D2T(1,:);
    PEC_D2 = reshape(PEC_D2, size(DenMC));
    PEC_D2p = PEC_D2pT(1,:);
    PEC_D2p = reshape(PEC_D2p, size(DenMC));
    PEC_Dm = PEC_DmT(1,:);
    PEC_Dm = reshape(PEC_Dm, size(DenMC));
    
    %Use transport to include D2 contribution in D2+, D- emission values 
    
   
%     %estimate range of possible D2 emission 
    if ~isfield(input.input,'inputMC') 
        MolG = rand(size(n1IntMC))-0.5;
    else
        MolG = input.input.inputMC.MolG;
    end
        
%     if ~isfield(input.input,'inputMC')
%        nomol = 1;
%     else
%         if ~isfield(input.input.inputMC,'n1IntTotMC')
%             nomol = 1;
%         else
%             nomol = 0;
%         end
%     end
%     if nomol == 1
%         n1IntTotMC = n1IntMC;
%         n2IntTotMC = n2IntMC;
%     else
%        n1IntTotMC = input.input.inputMC.n1IntTotMC;
%        n2IntTotMC = input.input.inputMC.n2IntTotMC;
%     end
    
    input.inputMC.n1IntTotMC = n1IntTotMC;
    input.inputMC.n2IntTotMC = n2IntTotMC;
    input.inputMC.MolG=MolG;
    
    MolG = 10.^(MolG + polyval(pMolG, log10(TeEMC1)));
    MolG(isinf(MolG)) = NaN;
        
    %subtract D2 emission from all lines
    n1NoD2 = n1IntTotMC - reshape(PEC_D2T(3,:),size(DenMC)).*(MolG).*DenMC;
    n1NoD2(n1NoD2<0) = NaN;
    n2NoD2 = n2IntTotMC - reshape(PEC_D2T(4,:),size(DenMC)).*(MolG).*DenMC;
    n2NoD2(n2NoD2<0) = NaN;
    DaNoD2 = DaTotMC - reshape(PEC_D2T(1,:),size(DenMC)).*(MolG).*DenMC;
    DaNoD2(DaNoD2<0) = NaN;

    if isfield(input.input,'DbMea')
        DbNoD2 = DbTotMC - reshape(PEC_D2T(2,:),size(DenMC)).*(MolG).*DenMC;
        DbNoD2(DbNoD2<0) = NaN;
        EVal(isnan(n1NoD2+n2NoD2+DaNoD2+DbNoD2))=0;
    else
        EVal(isnan(n1NoD2+n2NoD2+DaNoD2))=0;
    end
    
    %DbNoD2(DbNoD2<BIntAtom_Db)=NaN;
    %DaNoD2(DaNoD2<BIntAtom_Da)=NaN;
    
    %each of the Balmer lines without D2 consist of emission due to D2+,
    %D-, atomic excitation, atomic recombination; use analytic solutions
    
   % MultA = (1-FRec1MC).*(10.^naExc(log10(DenMC.*EVal), log10(TeEMC1.*EVal)))./(10.^n1Exc(log10(DenMC.*EVal), log10(TeEMC1.*EVal))) + ...
   %     FRec1MC.*(10.^naRec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)))./(10.^n1Rec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)));
   % MultB = (1-FRec1MC).*(10.^nbExc(log10(DenMC.*EVal), log10(TeEMC1.*EVal)))./(10.^n1Exc(log10(DenMC.*EVal), log10(TeEMC1.*EVal))) + ...
   %     FRec1MC.*(10.^nbRec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)))./(10.^n1Rec(log10(DenMC.*EVal), log10(TeRMC2.*EVal)));
    
   % Aa = 1 - reshape(PEC_D2pT(3,:)./PEC_D2pT(1,:),size(n1NoD2)).*MultA;
   % Ab = 1 - reshape(PEC_D2pT(3,:)./PEC_D2pT(2,:),size(n1NoD2)).*MultB;
   % Ba = 1 - reshape(PEC_DmT(3,:)./PEC_DmT(1,:),size(n1NoD2)).*MultA;
   % Bb = 1 - reshape(PEC_DmT(3,:)./PEC_DmT(2,:),size(n1NoD2)).*MultB;   
   % Ca = MultA.*n1NoD2;
   % Cb = MultB.*n1NoD2;
    
%     Ca(DaNoD2-Ca<0) = DaNoD2(DaNoD2-Ca<0);
%     Cb(DaNoD2-Ca<0) = DbNoD2(DaNoD2-Ca<0);
%      X = (DaNoD2.*(Ba - Bb.*(DaNoD2./DbNoD2)) + (DaNoD2./DbNoD2).*(Bb.*Ca - Ba.*Cb))./...
%          (Ab.*Ba - Aa.*Bb).*(DaNoD2./DbNoD2);
%      Y = (DaNoD2.*(Aa - Ab.*(DaNoD2./DbNoD2)) + (DaNoD2./DbNoD2).*(Ab.*Ca - Aa.*Cb))./...
%          (Ab.*Ba - Aa.*Bb).*(DaNoD2./DbNoD2);
    
    %disp('test')
    %DaNoD2 = Da.*(fD2+ + fDm + fDA); f
%     
%     D2Estimn1 = PEC_D2n1.*(MolG).*DLMC.*DenMC;
%     D2Estimn2 = PEC_D2n2.*(MolG).*DLMC.*DenMC;
%     
%     TR = DaMolMC - PEC_D2.*(MolG).*DLMC.*DenMC;
%     TR(TR<0) = 0;
%     fTR = TR./DaMolMC; %fTR is fraction of DaMol due to Dm, D2p
%     
    if ~isfield(input.input,'DbMea') %no dbeta measurement
        if D2pDm_model==2
            AMJUELloc1 = [matlab_home, '/data_files/amjuel_h12-e.txt'];
            AMJUELloc2 = [matlab_home, '/data_files/amjuel_h11-e.txt'];
            R1=Interrogate_AMJUEL_DoubleFit(AMJUELloc1,'2.0c',TeEMC1.*EVal,EVal.*DenMC)*0.95;%.
            R2=Interrogate_AMJUEL_SingleFit(AMJUELloc2,'7.0a', TeEMC1.*EVal)*0.7;
            RmRH2p = R2./R1;
        elseif D2pDm_model==1
            RmRH2p = 0;
        end
        
        fRmRH2pR = EVal.*(PEC_Dm.*RmRH2p)./(PEC_Dm.*RmRH2p + PEC_D2p); %fraction of emission (of that from Dm and D2p) due to Dm
    else
        %Get line ratios emission coefficients Db/Da for Hm and H2p
        LR_H2p = reshape(PEC_D2pT(2,:)./PEC_D2pT(1,:),size(DenMC));
        LR_Hm = reshape(PEC_DmT(2,:)./PEC_DmT(1,:),size(DenMC));
        
        %Real line ratio
        LR_R = (DbNoD2-BIntAtom_Db)./(DaNoD2-BIntAtom_Da);
                       
        %get effective emission fraction H2+ (1- this)=Hm fraction
        fRH2pRmR = (LR_R-LR_Hm)./(LR_H2p - LR_Hm);
        %EMF_H2p = 0.*Ne + NaN;
        %set overshoots to 0, 1
        fRH2pRmR(LR_R<LR_Hm) = 0;
        fRH2pRmR(LR_R>LR_H2p) = 1;
    
        fRH2pRmR(DbNoD2<BIntAtom_Db) = 1;
        fRH2pRmR(DaNoD2<BIntAtom_Da) = 1;
    
        %Get Hm fraction
        fRmRH2pR = 1 - fRH2pRmR;
    end  
        %Now we thus know the division between Da[H2+, H-] while accounting
        %for how changes in this contribution leads to changes in the
        %medium-n Balmer lines which influence the atomic part of Da
        %(through Ab./Aa and Bb./Ba)
        
        %Contamination due to H2p in n1Int,n2Int
        %Molecular part Da [D2+,D-]:
        Da_H2p_Hm = (DaNoD2 - BIntAtom_Da);%.*(1 + (1-Aa).*(1-fRmRH2pR) + (1-Ba).*(fRmRH2pR));
        Da_H2p_Hm(Da_H2p_Hm<0) = 0;
        
        %atomic part Da
        Da_H2 = DaTotMC - DaNoD2;
        DaMolMC = (Da_H2 + Da_H2p_Hm).*EVal;
        fTR = (1-(Da_H2./DaMolMC)).*EVal;
        
        %apply EVal filter
        fRmRH2pR = fRmRH2pR.*EVal;
        
        %nanmax(fTR(:))
        %nanmin(fTR(:))
        
        %get contamination Balmer lines
    %end

    %now that all emission channels are known; extrapolate these to n1, n2
    %Balmer lines
    
    n1MolMC = EVal.*DaMolMC.*((1-fTR).*reshape(PEC_D2T(3,:)./PEC_D2T(1,:),size(DenMC)) + ...
        fTR.*fRmRH2pR.*reshape(PEC_DmT(3,:)./PEC_DmT(1,:),size(DenMC)) + ...
        fTR.*(1 - fRmRH2pR).*reshape(PEC_D2pT(3,:)./PEC_D2pT(1,:),size(DenMC)));
    
    n2MolMC = EVal.*DaMolMC.*((1-fTR).*reshape(PEC_D2T(4,:)./PEC_D2T(1,:),size(DenMC)) + ...
        fTR.*fRmRH2pR.*reshape(PEC_DmT(4,:)./PEC_DmT(1,:),size(DenMC)) + ...
        fTR.*(1 - fRmRH2pR).*reshape(PEC_D2pT(4,:)./PEC_D2pT(1,:),size(DenMC)));
    
    if n1MolMCFilter
    n1MolMC(n1MolMC>n1IntTotMC) = NaN;
    n2MolMC(n2MolMC>n2IntTotMC) = NaN;
    end
    
    n1MolMCD2 = (n1MolMC./n1MolMC).*DaMolMC.*((1-fTR).*reshape(PEC_D2T(3,:)./PEC_D2T(1,:),size(DenMC)));
    n1MolMCDm = (n1MolMC./n1MolMC).*DaMolMC.*fTR.*fRmRH2pR.*reshape(PEC_DmT(3,:)./PEC_DmT(1,:),size(DenMC));
    n1MolMCD2p = (n1MolMC./n1MolMC).*DaMolMC.*fTR.*(1 - fRmRH2pR).*reshape(PEC_D2pT(3,:)./PEC_D2pT(1,:),size(DenMC));
    
%     [ii] = find(nansum(n1Int,1)>0);
%     for i=1:ii
%        Cn1MolMC(:,i,:) = DSS_Aux_inpaint_nans(squeeze(Cn1MolMC(:,i,:)),4); 
%        Cn2MolMC(:,i,:) = DSS_Aux_inpaint_nans(squeeze(Cn2MolMC(:,i,:)),4);
%     end
    
    input.ResultMolMC = DSS_Aux_structvars({'fieldNames','DaMolMC','n1MolMC','n2MolMC','fRmRH2pR','fTR','n1MolMCD2','n1MolMCD2p','n1MolMCDm'});
       
end



function FRecC = FRecForceTrend(FRec)
%When calculating FRec, during a density ramp experiment, FRec should be increasing, resulting in a hyperbolic tangent curve as function of time.
%However, due to the issues with converging to a unique FRec number at low/high FRec, this sometimes results in a 'S-like' behaviour.

%This function tends to rectify this behaviour

Window=0.1;

FRecC = FRec;

for i=1:numel(FRec(:,1))
    %for each wavelength
    Validity = ~isnan(FRec(i,:));
    if sum(Validity)>1
        %FRecTime = smooth(FRec(i,:),Window,'loess').*Validity';
        FRecTime = FRec(i,:);
        FRecTime(~Validity)=NaN;

        %get min/max of FrecTime
        [~,imin]=min(FRecTime);
        [~,imax] = max(FRecTime);

        indxmin = find(Validity);
        indxmin=indxmin(indxmin<imin);
        indxmax = find(Validity);
        indxmax=indxmax(indxmax>imax);

        FRecC(i,indxmin) = FRec(i,imin);
        FRecC(i,indxmax) = FRec(i,imax);
    
    end
    
end
