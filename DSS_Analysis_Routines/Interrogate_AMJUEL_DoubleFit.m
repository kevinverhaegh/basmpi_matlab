function V = Interrogate_AMJUEL_DoubleFit(AMJUELloc,Reac,temp,den)
%% This function, provided a temperature (in eV) and another parameter (usually 
%density in m^-3), interrogates AMJUEL e data fits using a double
%polynomial interrogation

%Kevin Verhaegh v. 0.1 - University of York / CCFE

%AMJUELloc = '/home/verhaegh/Public files/amjuel_h12-e.txt'; %AMJUEL file location
%% Find appropriate AMJUEL reaction line

%first find appropriate reaction in AMJUEL datafile
fid = fopen(AMJUELloc); %open file
SAMJUEL = textscan(fid,'%s','delimiter','\n'); %make string of file
fclose(fid); %close file

Header= find(~cellfun(@isempty,strfind(SAMJUEL{1},[Reac, ' &']))); % find row by finding which cell (each cell is a different line) contains the string specified

if isempty(Header)
    error(['Reaction: ', Reac, ' not found in AMJUEL file:', AMJUELloc])
end

%% Find appropriate AMJUEL matrix
%below this line, there is a 8x8 matrix in SAMJUEL; this needs to be loaded
%in

fid = fopen(AMJUELloc); %open file
M = textscan(fid,'%f%f%f%f%f%f%f%f%f%[^\n\r]', 9,'delimiter','&', 'HeaderLines', Header, 'EndOfLine', '\r\n'); %Get matrix
fclose(fid); %close file

M(end) = []; %clean up last index in the cell - which contrains the uncertainty of the fit
M = cell2mat(M); % convert the matrix from a cell strucutre provided by textscan into a matrix

%% Interrogate AMJUEL matrix using a double polynomial fit

den = den / 1e6; %convert in per cubic metre
den = den/ 1e8; %convert to 1e8 per cubic metre (AMJUEL default)
V = Interp_DoubleFit(M,temp,den); %reaction rate in cm^3 / s
%V = V*1e6; %convert reaction rate to m^3 /s


function V=Interp_DoubleFit(M, T, E)
%double fit interrogration for AMJUEL-like data of the form: sum (log(T).^i
%* sum(log(E).^j * c_ij)).

V = 0.*T; %output value

for i=1:numel(M(:,1))
    Temp = 0.*T; %temporary
    for j=1:numel(M(1,:))
        Temp = Temp + M(i,j).*log(E).^(j-1);
    end
    V = V + Temp.*log(T).^(i-1);
end

V=exp(V);