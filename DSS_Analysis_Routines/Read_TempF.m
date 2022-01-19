function input = Read_TempF(inputloc)

%combines input structures in the form TempF

flist = dir([inputloc,'TempF_*']);

load([inputloc,flist(1).name]);

inputB = input;
inputB.inputMC.THm = inputB.input.THm;
dum = inputB.input.Unc;
f = fieldnames(dum);
for i=1:length(f)
    inputB.inputMC.(f{i}) = dum.(f{i});
end
inputB.input = rmfield(inputB.input,'THm');
inputB.input=rmfield(inputB.input,'Unc');

inputMC_list = strcat('.inputMC.',fieldnames(inputB.inputMC));
ResultMC_list = strcat('.ResultMC.',fieldnames(inputB.ResultMC));
ResultMolMC_list = strcat('.ResultMolMC.',fieldnames(inputB.ResultMolMC));

namelists = [inputMC_list; ResultMC_list; ResultMolMC_list;];

for i=2:numel(flist)
    load([inputloc,flist(i).name]);
    input.inputMC.THm = input.input.THm;
    dum = input.input.Unc;
    f = fieldnames(dum);
    for j=1:length(f)
        input.inputMC.(f{j}) = dum.(f{j});
    end
    for j=1:numel(namelists)
    eval(['inputB',namelists{j},' = cat(',num2str(sum(size(eval(['inputB',namelists{j}]))>1)),',inputB',namelists{j},',input',namelists{j},');']);
    end
end

input = inputB;