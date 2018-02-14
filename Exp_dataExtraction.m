%% Read in optimization text from MATLAB, and then create a data set to store the 
% 1) number of OPT cycles
% 2) value of meritfunction in each cycle
% 3) curvature variables in each cycle

%%
% ********************Format of the data**********************
%                           
%                           COLUMN
%                
%                 LINE    Cycle number
%                         MF
%                         var1
%                         var2
%                         .
%                         .
%                         .
%                         varN
%
% *************************************************************

%%


% Surface number of the variables 
SurNum = {' 1',' 2',' 3',' 4',' 5',' 6',' 8',' 9','10','11'};  % an extra blank space is put for the integer smaller than 10 in order to match the number of space.  %%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE %%%%%%%%%%%%%%%%%%%%%%%%%
%Number of variables 
NumVar = size(SurNum,2); 


disp('Load text file with optimization information');
fp = uigetfile({'*.txt'}, 'Select the txt file');
text = fileread(fp);
disp('txt file loaded');

% Index of the Cycle Number
CycInd = strfind(text,'CYCLE NUMBER');
CycInd = CycInd(1:2:end-1); % 'CYCLE NUMBER' appears twice for each cycle, we only choose the first one. 

%Get the number of cycle
LastCycInd = CycInd(end);
if text(LastCycInd+15)==':'
    NumCyc = str2double(text(LastCycInd+13:LastCycInd+14))+1;
else
    NumCyc = str2double(text(LastCycInd+13:LastCycInd+15))+1;
end

EFInd = strfind(text,'ERR. F.  ');
EFInd = EFInd(1:2:end-1); 

%Create index for the variables
VarInd = zeros(NumVar,NumCyc);

for n=1:NumVar
    strtag = ['  ',SurNum{n},':'];  % two blank spaces are added for the search of numbers 
    VarInd(n,:) = strfind(text,strtag);
end


TrackData = zeros(NumVar+2,NumCyc);
TrackData(1,1:end)=0:NumCyc-1; % First line gives the number of cycle 

for n=1:NumCyc
    TrackData(2,n) = str2double(text(EFInd(n)+17:EFInd(n)+27));
end


for n=3:NumVar+2
    for m=1:NumCyc
        TrackData(n,m) = str2double(text(VarInd(n-2,m)+10:VarInd(n-2,m)+20));
    end
end



%% Process variables - convert curvature into coordinate information

% index of the curvature variables for plotting
varInd = [3,4,5]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---VARIABLE---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cRef = -0.134775624347017; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---VARIABLE---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors
stepX = zeros(1,NumCyc);
stepY = zeros(1,NumCyc);
stepZ = zeros(1,NumCyc);

for n=1:NumCyc
    c1 = TrackData(varInd(1)+2,n);
    c2 = TrackData(varInd(2)+2,n);
    stepX(n) = (c1-(c2+cRef)/2)*sqrt(2);
    stepY(n) = (c2-cRef)/sqrt(2/3);
    stepZ(n) = TrackData(2,n);
end

%figure;
hold on;
scatter3(stepX,stepY,stepZ,'filled','y');
pline = line(stepX,stepY,stepZ); % connecting all the points in the scatter point plot.
pline.Color = 'green';
