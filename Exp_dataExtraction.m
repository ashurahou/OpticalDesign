%% Function of this code
% 
% This code reads the text file created by CODEV, with optimization
% information stored. It extracts the following information from the text
% 1) number of optimization cycles
% 2) value of the merit function in each cycle
% 3) curvature variables in each cycle
%
% This code is only for the case of SPCS, where the starting point of the optimization is
% near a SPCS saddle point. 

%%


%%
% Input
%
% text file genereated by CODEV with optimization history will be read 
%
% Surface number of the variables
In_SurNum = {' 1',' 2',' 3',' 4',' 5',' 6',' 8',' 9','10','11'};  % !!!an extra blank space is put for the integer smaller than 10 in order to match the number of space.!!! 

% Index of the variable within the SurNum array, which will be used for 
% plotting the 2D plane.  
In_varInd = [4,5,6];  % Example: If the variable surface are {' 2',' 3',' 4',' 6',' 7',' 8',' 9'}, In_varInd = [2,3,4] means that Surface 3,4,6 will be used to calculated the trajectory in the 2D plane. 

% Value of the reference curvature
In_cRef = -0.132831298308587; 
%
%
%% Output
%
%
% Variable in the workspace
%
    %  TrackData; Storing the information about Cycle number, MF and
    %  variables. It follows the following format.
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
%
%  
% Plot:
% The trajectory of the optimzation will be plot w.r.t. the three selected
% variables. This plot should be combined with the plot of the basins of
% attraction, which is the SPCS saddle point of the same system. The
% basins of attraction plot can be generated with Exp_plotMFSPCS.m. With
% the last plot from Exp_plotMFSPCS.m active, this macro will plot the
% trajectory on the basin of attraction plot.
%
%
%% Dependance
% 
% 1. text file genereated by CODEV with optimization history will be read 
%
% 2. plots generated by Exp_plotMFSPCS.m are used to analyse the trajectory
% plotted by this macro.
%
%
%% Version
%
% V 1.0 FEB 2018   
% Zhe HOU
% zhe.hou2011@gmail.com
%

%% Extraction of data from text

% Surface number of the variables 
SurNum = In_SurNum;
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



%% Plot on the SPCS 2D plane 

% index of the curvature variables for plotting
varInd = In_varInd; 
cRef = In_cRef; 

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


% Plot:
% The trajectory of the optimzation will be plot w.r.t. the three selected
% variables. This plot should be combined with the plot of the basins of
% attraction, which is the SPCS saddle point of the same system. The
% basins of attraction plot can be generated with Exp_plotMFSPCS.m. With
% the last plot from Exp_plotMFSPCS.m active, this macro will plot the
% trajectory on the basin of attraction plot.

%figure;
hold on;
scatter3(stepX,stepY,stepZ,'filled','g');
pline = line(stepX,stepY,stepZ); % connecting all the points in the scatter point plot.
pline.Color = 'black';
