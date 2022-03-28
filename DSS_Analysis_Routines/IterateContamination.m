%Script for using DSS_Analysis_Calc_RecIonFRec in an iterative way to
%include plasma-molecule interaction

if isfield(settings,'starting_contamin')
    starting_contamin_n1 = settings.starting_contamin.n1;
    starting_contamin_n2 = settings.starting_contamin.n2;
else
    starting_contamin_n1 = 0.*inputBase.input.n1Int;
    starting_contamin_n2 = 0.*inputBase.input.n1Int;
end

inputBase.input.n1Int = inputBase.input.n1Int .* (1-starting_contamin_n1);
inputBase.input.n2Int = inputBase.input.n2Int .* (1-starting_contamin_n2);

%analysis without assuming plasma-mol. influences higher-n Balmer line
inputBase = DSS_Analysis_Calc_RecIonFRec(inputBase.input);

inputBase.input.n1Int = inputBase.input.n1Int ./ (1-starting_contamin_n1);
inputBase.input.n2Int = inputBase.input.n2Int ./ (1-starting_contamin_n2);
inputBase.inputMC.n1IntTotMC = inputBase.inputMC.n1IntTotMC ./ repmat((1-starting_contamin_n1),1,1,inputBase.input.Iter);
inputBase.inputMC.n2IntTotMC = inputBase.inputMC.n2IntTotMC ./ repmat((1-starting_contamin_n2),1,1,inputBase.input.Iter);

%save intermediate result
save('tempBase','inputBase','-v7.3');

%convergence criteria; ConvCritM (median of changes to contamination should
%be between + and - this
ConvCritM = 2e-3;
ConvCritB = 2e-2; % 68 conf. interval should be between plus/minus this
Epss = 2e-4; %numerical accuracy to distinguish from zero
ConvDura = 4; %how many analysis runs have to be converged according to the above criteria

%do first iteration
clear input Cn1MolMCd
input{1} = inputBase; %set up input structure
input{1}.input.inputMC = inputBase.inputMC; %keep randomization
input{1}.input.inputMC.n1IntMC = input{1}.input.inputMC.n1IntTotMC - input{1}.ResultMolMC.n1MolMC; %apply molecular estimates of higher-n Balmer line emission
input{1}.input.inputMC.n2IntMC = input{1}.input.inputMC.n2IntTotMC - input{1}.ResultMolMC.n2MolMC;
%input{1}.input.inputMC = RegenerateNaN(input{1}.input); %apply regeneration of 'kicked out' data
input{1} = DSS_Analysis_Calc_RecIonFRec(input{1}.input); %analyse
Cn1MolMCd(1,:,:,:) = (input{1}.ResultMolMC.n1MolMC-inputBase.ResultMolMC.n1MolMC)./input{1}.ResultMolMC.n1MolMC; %relative change in molecular contamination estimates for the lowest-n Balmer line used in the atomic analysis
D = Cn1MolMCd(1,:,:,:);
Conv(1,:)= DSS_Aux_wprctile(D(~isnan(D)),[0.5,0.16,0.84]); %calculate quantiles of the above changes in molecular contamination
temp = input{1};
save('temp_1','temp','-v7.3');

Ctimer = 0;

for i=2:10 %loop for remaining iterations
   %as above
   input{i} = input{i-1};
   input{i-1} = [];
   input{i}.input.inputMC = inputBase.inputMC;
   input{i}.input.inputMC.n1IntMC = input{i}.input.inputMC.n1IntTotMC - input{i}.ResultMolMC.n1MolMC;
   input{i}.input.inputMC.n2IntMC = input{i}.input.inputMC.n2IntTotMC - input{i}.ResultMolMC.n2MolMC;
   %if i<RegenNM %apply regeneration for a certain number of iterations
    %input{i}.input.inputMC = RegenerateNaN(input{i}.input);
   %end
   input{i} = DSS_Analysis_Calc_RecIonFRec(input{i}.input);
   %figure
   Cn1MolMCd(i,:,:,:) = (input{i}.ResultMolMC.n1MolMC-(input{i}.inputMC.n1IntTotMC-input{i}.inputMC.n1IntMC))./input{i}.ResultMolMC.n1MolMC; %relative change in contamination
   D = Cn1MolMCd(i,:,:,:);
   Conv(i,:)= DSS_Aux_wprctile(D(~isnan(D)),[0.5,0.16,0.84]); %quantiles based on this
   %apply convergence criteria based on quantiles
    if (abs(Conv(i,1)) < ConvCritM) && (Conv(i,2) < Epss) && (Conv(i,3) >-Epss) && (Conv(i,2) > -ConvCritB) && (Conv(i,3) < ConvCritB)
        Ctimer = Ctimer+1; %converged; add +1 to convergence timer
     else
         Ctimer = 0; %not converged; reset convergence timer
    end
     
    %save temporary outputs
    temp = input{i};
    save(['temp_',num2str(i)],'temp','-v7.3');

    %if convergence succeeded for a suffisient duration; end iteration
     if Ctimer >ConvDura
        disp('Convergence succeeded')
        i
        break
     end
    
end

%convergence ended (or iteration number exceeded)
input = input{i};

%get data from latest analysis
TeEMC1 = input.ResultMC.TeEMC1;
TeRMC2 = input.ResultMC.TeRMC2;
FRec1MC = input.ResultMC.FRec1MC;
FRec2MC = input.ResultMC.FRec2MC;
n1IntMC = input.inputMC.n1IntMC;
n2IntMC = input.inputMC.n2IntMC;
DaMolMC = input.ResultMolMC.DaMolMC;
n1MolMC = input.ResultMolMC.n1MolMC;
n2MolMC = input.ResultMolMC.n2MolMC;
fRmRH2pR = input.ResultMolMC.fRmRH2pR;
fTR = input.ResultMolMC.fTR;
%get booleans for where the data is available
TeEMC1I = int16(isnan(input.ResultMC.TeEMC1));
TeRMC2I = int16(isnan(input.ResultMC.TeRMC2));
FRec1MCI = int16(isnan(input.ResultMC.FRec1MC));
FRec2MCI = int16(isnan(input.ResultMC.FRec2MC));
n1IntMCI = int16(isnan(input.inputMC.n1IntMC));
n2IntMCI = int16(isnan(input.inputMC.n2IntMC));
DaMolMCI = int16(isnan(input.ResultMolMC.DaMolMC));
n1MolMCI = int16(isnan(input.ResultMolMC.n1MolMC));
n2MolMCI = int16(isnan(input.ResultMolMC.n2MolMC));
fRmRH2pRI = int16(isnan(input.ResultMolMC.fRmRH2pR));
fTRI = int16(isnan(input.ResultMolMC.fTR));

for i=(numel(Conv(:,1))-ConvDura+1):(numel(Conv(:,1))-1) %iterate over the last ConvDura discharges and average over these
    load(['temp_',num2str(i)]) %load data
    tmp = cat(4,TeEMC1,temp.ResultMC.TeEMC1); TeEMC1 = sum(tmp,4,"omitnan");  %concatenate all data together with all the other data; then sum the data
    tmp = cat(4,TeEMC1I,isnan(temp.ResultMC.TeEMC1)); TeEMC1I = sum(tmp,4,"omitnan");
    tmp = cat(4,TeRMC2,temp.ResultMC.TeRMC2); TeRMC2 = sum(tmp,4,'omitnan');
    tmp = cat(4,TeRMC2I,isnan(temp.ResultMC.TeRMC2)); TeRMC2I = sum(tmp,4,'omitnan');
    tmp = cat(4,FRec1MC,temp.ResultMC.FRec1MC); FRec1MC = sum(tmp,4,'omitnan');
    tmp = cat(4,FRec1MCI,isnan(temp.ResultMC.FRec1MC)); FRec1MCI = sum(tmp,4,'omitnan');
    tmp = cat(4,FRec2MC,temp.ResultMC.FRec2MC); FRec2MC = sum(tmp,4,'omitnan');
    tmp = cat(4,FRec2MCI,isnan(temp.ResultMC.FRec2MC)); FRec2MCI = sum(tmp,4,'omitnan');
    tmp = cat(4,n1IntMC,temp.inputMC.n1IntMC); n1IntMC = sum(tmp,4,'omitnan');
    tmp = cat(4,n1IntMCI,isnan(temp.inputMC.n1IntMC)); n1IntMCI = sum(tmp,4,'omitnan');
    tmp = cat(4,n2IntMC,temp.inputMC.n2IntMC); n2IntMC = sum(tmp,4,'omitnan');
    tmp = cat(4,n2IntMCI,isnan(temp.inputMC.n2IntMC)); n2IntMCI = sum(tmp,4,'omitnan');
    tmp = cat(4,DaMolMC,temp.ResultMolMC.DaMolMC); DaMolMC = sum(tmp,4,'omitnan');
    tmp = cat(4,DaMolMCI,isnan(temp.ResultMolMC.DaMolMC)); DaMolMCI = sum(tmp,4,'omitnan');
    tmp = cat(4,n1MolMC,temp.ResultMolMC.n1MolMC); n1MolMC = sum(tmp,4,'omitnan');
    tmp = cat(4,n1MolMCI,isnan(temp.ResultMolMC.n1MolMC)); n1MolMCI = sum(tmp,4,'omitnan');
    tmp = cat(4,n2MolMC,temp.ResultMolMC.n2MolMC); n2MolMC = sum(tmp,4,'omitnan');
    tmp = cat(4,n2MolMCI,isnan(temp.ResultMolMC.n2MolMC)); n2MolMCI = sum(tmp,4,'omitnan');
    tmp = cat(4,fRmRH2pR,temp.ResultMolMC.fRmRH2pR); fRmRH2pR = sum(tmp,4,'omitnan');
    tmp = cat(4,fRmRH2pRI,isnan(temp.ResultMolMC.fRmRH2pR)); fRmRH2pRI = sum(tmp,4,'omitnan');
    tmp = cat(4,fTR,temp.ResultMolMC.fTR); fTR = sum(tmp,4,'omitnan');
    tmp = cat(4,fTRI,isnan(temp.ResultMolMC.fTR)); fTRI = sum(tmp,4,'omitnan');
end

%calculate over how many converged points the sum took place
input.ResultMC.TeEMC1 = TeEMC1./(ConvDura - TeEMC1I);
input.ResultMC.TeEMC1(isinf(input.ResultMC.TeEMC1)) = NaN; %set non-physical values to NaN
input.ResultMC.TeRMC2 = TeRMC2./(ConvDura - TeRMC2I);
input.ResultMC.TeRMC2(isinf(input.ResultMC.TeRMC2)) = NaN;
input.ResultMC.FRec1MC = FRec1MC./(ConvDura - FRec1MCI);
input.ResultMC.FRec1MC(isinf(input.ResultMC.FRec1MC)) = NaN;
input.ResultMC.FRec2MC = FRec2MC./(ConvDura - FRec2MCI);
input.ResultMC.FRec2MC(isinf(input.ResultMC.FRec2MC)) = NaN;
input.inputMC.n1IntMC = n1IntMC./(ConvDura - n1IntMCI);
input.inputMC.n1IntMC(isinf(input.inputMC.n1IntMC)) = NaN;
input.inputMC.n2IntMC = n2IntMC./(ConvDura - n2IntMCI);
input.inputMC.n2IntMC(isinf(input.inputMC.n2IntMC)) = NaN;
input.ResultMolMC.DaMolMC = DaMolMC./(ConvDura - DaMolMCI);
input.ResultMolMC.DaMolMC(isinf(input.ResultMolMC.DaMolMC)) = NaN;
input.ResultMolMC.n1MolMC = n1MolMC./(ConvDura - n1MolMCI);
input.ResultMolMC.n1MolMC(isinf(input.ResultMolMC.n1MolMC)) = NaN;
input.ResultMolMC.n2MolMC = n2MolMC./(ConvDura - n2MolMCI);
input.ResultMolMC.n2MolMC(isinf(input.ResultMolMC.n2MolMC)) = NaN;
input.ResultMolMC.fRmRH2pR = fRmRH2pR./(ConvDura - fRmRH2pRI);
input.ResultMolMC.fRmRH2pR(isinf(input.ResultMolMC.fRmRH2pR)) = NaN;
input.ResultMolMC.fTR = fTR./(ConvDura - fTRI);
input.ResultMolMC.fTR(isinf(input.ResultMolMC.fTR)) = NaN;
