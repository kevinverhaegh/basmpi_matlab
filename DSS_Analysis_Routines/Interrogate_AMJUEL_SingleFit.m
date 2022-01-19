function V = Interrogate_AMJUEL_SingleFit(AMJUELloc,Reac,temp)
%% This function, provided a temperature (in eV) and another parameter (usually 
%density in m^-3), interrogates AMJUEL e data fits using a single
%polynomial interrogation

%Kevin Verhaegh v. 0.1 - University of York / CCFE
%% Find appropriate AMJUEL reaction line

%first find appropriate reaction in AMJUEL datafile
fid = fopen(AMJUELloc); %open file
SAMJUEL = textscan(fid,'%s','delimiter','\n'); %make string of file
fclose(fid); %close file

Header= find(~cellfun(@isempty,strfind(SAMJUEL{1},[Reac, ' &']))); % find row by finding which cell (each cell is a different line) contains the string specified

if isempty(Header)
    error(['Reaction: ', Reac, ' not found in AMJUEL file:', AMJUELloc])
end

%% Find the length of the appropriate polynomial

%% Find appropriate AMJUEL matrix
%below this line, there is a 8x8 matrix in SAMJUEL; this needs to be loaded
%in

fid = fopen(AMJUELloc); %open file
M = textscan(fid,'%*s%*s%*s%*s%*f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 1,'delimiter','&', 'HeaderLines', Header-1, 'EndOfLine', '\r\n'); %Get matrix
fclose(fid); %close file

M(end) = []; %clean up last index in the cell - which contrains the uncertainty of the fit
M = cell2mat(M); % convert the matrix from a cell strucutre provided by textscan into a matrix
M(isnan(M)) = []; %remove NaNs

%% Interrogate AMJUEL matrix using a  single polynomial fit
V = polyval(fliplr(M),log(temp));
V = exp(V);
