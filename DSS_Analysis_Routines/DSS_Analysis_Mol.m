function input = DSS_Analysis_Mol(input)

%Script for using DSS_Analysis_Calc_RecIonFRec in an iterative way to
%include plasma-molecule interaction.  Requires input structure

%set-up parpool
%parpool()
settings = input.settings;
fname = settings.fname;

%file directory for storing info
if ~exist(fname, 'dir')
   mkdir(fname)
else
    error('fname already exists!')
end

eval(['cd ',fname])

%check input structures for any references to randomization; remove those
if isfield(input,'inputMC')
    input = rmfield(input,'inputMC');
end
if isfield(input,'THm')
    input = rmfield(input,'THm');
end
if isfield(input,'Unc')
    input = rmfield(input,'Unc');
end

if ~isfield(settings,'MaxIter')
    settings.MaxIter = 500;
end

if input.Iter > settings.MaxIter
        N = ceil(input.Iter/settings.MaxIter);
        input.Iter = settings.MaxIter;
else
    N = 1;
end

BaseInput = input; %get basic input
save('BaseInput','BaseInput')
clear input

%divide the analysis into smaller analysis sets with less randomization
%points and run a loop of those. Total iteration points is N*Iter. This
%speeds up the initial analysis and enables one to get intermediate results
for ii=1:N
    clear inputBase
    clear Conv
    inputBase.input = BaseInput; %put input in the correct form
    IterateContamination %compute
    save(['TempF_',num2str(ii)],'input','-v7.3') %save
end

input = Read_TempF(fname);

save([fname,'output'],'input','-v7.3')