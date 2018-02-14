%% This code fulfill the following function
% 
% ! This code plots the Merit Function (MF) height around a saddle point, which is constructed 
% ! with the Saddle-Point Construction (SPC) method for the special case when the curvatures 
% ! of the added null-element meniscus have the same values as the curvature of the reference 
% ! surface and the glasses of the null element and the reference surface are the same. 
% ! The plot is done in a 2D cut of the parameter space in which the two equimagnitude MF lines 
% ! cross. The coordinates in this plane are given by xp and yp. 
%
%
%% Input
%  Input variables. They are modified for different systems.
%
%
  % !!!Lens file and meritfunction macro will be loaded. !!!
%
  % variable name in CODEV 
  In_varName = ['CUY S2';'CUY S3';'CUY S4'];
  
  % tag for the reference surface in the varName array. Example:  ['CUY S2';'CUY S3';'CUY S4'], reference surface is S2, then In_reftag = 1. 
  In_reftag = 1; 
  
  % resolution for the scanning, dim x dim will be the resolution of the scan plot
  In_dim = 81;                        
  
  % lower boundary for the scanning range
  In_bound_down = -0.002;                 
  % upper boundary for the scanning range
  In_bound_up = 0.002;                    
%
%% Output
%
% Six plots will be generated:
%   1) Basins of attraction after executing the optimization macro once. (simulating one optimization click)
%   2) Basins of attraction after executing two times of optimization
%   marcro. (simulating two optimization clicks).
%   3) Merit function (MF) landscape without optimization (Surf plot).
%   4) Merit function landscape without optimization (contour plot).
%   5) Line plot of meritfucntion in X and Y direction. This will show the
%   ascending and descending 
%   6) Surf plot with basins of attractions + contour plot of meritfunction
%   without optimization + indication of equal MF value line
% 
%   
%
%
%
%
% 
%% Dependance
%  1. CODEV lens file in .seq or .len format
%  2. CODEV optimization macro in .seq format
%
%% Version

% V 1.0 FEB 2018
% ZHE HOU
% zhe.hou2011@gmail.com


%%

% Defining variables
global CVS
global CVsys
global CVmeritFunction
global varName
global varIndex
global varRefValue
global CVoptionSet
global convPrec

% Define CVoptionSet
CVoptionSet = cell(8,2);
CVoptionSet{1,1} = 'CVversion';
CVoptionSet{1,2} = '105';
CVoptionSet{2,1} = 'CVpath';
CVoptionSet{2,2} = [];   % will be defined later
CVoptionSet{3,1} = 'CVsys';
CVoptionSet{3,2} = [];  % will be defined later, in the process loop
CVoptionSet{4,1} = 'CVtimeoutSet';
CVoptionSet{4,2} = 5000000;
CVoptionSet{5,1} = 'CVbuffersizeSet';
CVoptionSet{5,2} = 300000;
CVoptionSet{6,1} = 'CVmeritFunction';
CVoptionSet{6,2} = []; % will be defined later['"',get(handles.CVmf,'String'),'"','    ',num2str(0),'    ','"N"'];
CVoptionSet{7,1} = 'CVautOpt';
CVoptionSet{7,2} = []; % omitted here //['"',get(handles.CVmf,'String'),'"','    ',num2str(100),'    ','"N"'];
CVoptionSet{8,1} = 'CVautOptFull';
CVoptionSet{8,2} = []; % omitted here //['"',get(handles.CVmf,'String'),'"','    ',num2str(1000),'    ','"N"'];


%% Load lens files and ERF files

disp('Load lens seq files.');
[CVname, CVpath, ~] = uigetfile({'*.seq';'*.len'}, 'Select CodeV seq files for the systems');
disp('Lens files loaded');


disp('Load Optimization Macro seq files.');
[ERFname, ~ , ~] = uigetfile([CVpath '*.seq'], 'Select a CodeV seq file as the MF');   % Error function should be put at the same folder of all the systems
disp('Optimization Macro loaded');


if ischar(CVname)
    SysNum = 1;
else
    SysNum = size(CVname,2);  % number of systems loaded
end

CVoptionSet{2,2} = CVpath;   % load CVpath
CVoptionSet{6,2} = ['"',ERFname,'"']; % Define Error function

cd(CVpath); % set the lens file folder as the current folder

%% Load a specific system and excute CodeV macro

if ischar(CVname)
    CVoptionSet{3,2} = ['"',CVname,'"']; % define system (only one system is loaded)
    CVsys = CVname(1:end-4);  % define the system used in a single loop & remove extension
else
    CVoptionSet{3,2} = ['"',CVname{1,n},'"'];   % define the system
    CVsys = CVname{1,n};
    CVsys = CVsys(1:end-4); 
end

varName = In_varName;
varIndex = true(3,1);
convPrec = '%.15f'; %accurancy of the number

reftag = In_reftag;

CVstart; % start connection with CODEV

[varRefValue] = CVreadSys(CVS, varName);

CurRef = varRefValue(reftag);

%% Plot 2D variation

% % MODE 1, scan around a point
% p = 0.2;
% t = 11;  % use an odd number so that the original system sits at the center of the square
% tempArrayA = linspace(varRefValue(1)*(1-p),varRefValue(1),t);
% tempArrayB = linspace(varRefValue(1)*(1+p),varRefValue(1),t);
% dim = 2*t-1;
% 
% if varRefValue(1) > 0
%     x = [tempArrayA,fliplr(tempArrayB(1:end-1))];
% else
%     x = [tempArrayB,fliplr(tempArrayA(1:end-1))];
% end

% MODE 2, scan in a designated area
dim = In_dim;                       
bound_down = In_bound_down;               
bound_up = In_bound_up;                    

xp = linspace(bound_down,bound_up,dim);

yp = xp;
[XP,YP] = meshgrid(xp,yp);
colon='; ';
colonArray=colon(ones(sum(varIndex,1),1),:);
whiteSpace=' ';
whiteSpaceArray=whiteSpace(ones(sum(varIndex,1),1),:);

ERFopt = zeros(dim,dim);
ERFopt_two = zeros(dim,dim);

CVcomm(CVS,'SAV temp.len');
for m = 1:dim
    for n = 1:dim      
        CVcomm(CVS, 'RES temp');

        c1 = CurRef + XP(m,n)/sqrt(2) + YP(m,n)/sqrt(6);
        c2 = CurRef + YP(m,n)*sqrt(2/3);
        c3 = CurRef - XP(m,n)/sqrt(2) + YP(m,n)/sqrt(6);
        cur1 = num2str(c1,convPrec);
        cur2 = num2str(c2,convPrec);
        cur3 = num2str(c3,convPrec);
        
        %to make strings equal size
        maxlength = max([numel(cur1),numel(cur2),numel(cur3)]);
        if numel(cur1)<maxlength, cur1=pad(cur1,maxlength); end
        if numel(cur2)<maxlength, cur2=pad(cur2,maxlength); end
        if numel(cur3)<maxlength, cur3=pad(cur3,maxlength); end
                
        commandString = [varName(varIndex,:), whiteSpaceArray, [cur1;cur2;cur3], colonArray];  
        commandStringRearr = commandString';
        CVcommandString = [(commandStringRearr(:))' '; ' 'in ' CVmeritFunction];
        [CVresponse,index1,index2] = CVcomm(CVS,CVcommandString);
        text = '=';
        index = strfind(CVresponse,text);
        ERFopt(m,n) = str2double(CVresponse(index+2:index+2+9));
        %trigger the second optimization
        CVcommandString = ['in ' CVmeritFunction];
        [CVresponse,index1,index2] = CVcomm(CVS,CVcommandString);
        text = '=';
        index = strfind(CVresponse,text);
        ERFopt_two(m,n) = str2double(CVresponse(index+2:index+2+9));
        
    end
end
        
ERFopt(ERFopt>9e9) = 0; % put ray failed one to zero
ERFopt_two(ERFopt_two>9e9) = 0;

figure;
surf(XP,YP,ERFopt);
xlabel('X');ylabel('Y');zlabel('Optimized MF value');title('Basin of attractions (one optimization click)');colorbar;
figure;
surf(XP,YP,ERFopt_two);
xlabel('X');ylabel('Y');zlabel('Optimized MF value');title('Basin of attractions (two optimization clicks)');colorbar;


ERF = zeros(dim,dim);
for m = 1:dim
    for n = 1:dim      
        CVcomm(CVS, 'RES temp');
        
        c1 = CurRef + XP(m,n)/sqrt(2) + YP(m,n)/sqrt(6);
        c2 = CurRef + YP(m,n)*sqrt(2/3);
        c3 = CurRef - XP(m,n)/sqrt(2) + YP(m,n)/sqrt(6);
        cur1 = num2str(c1,convPrec);
        cur2 = num2str(c2,convPrec);
        cur3 = num2str(c3,convPrec);
        
        %to make strings equal size
        maxlength = max([numel(cur1),numel(cur2),numel(cur3)]);
        if numel(cur1)<maxlength, cur1=pad(cur1,maxlength); end
        if numel(cur2)<maxlength, cur2=pad(cur2,maxlength); end
        if numel(cur3)<maxlength, cur3=pad(cur3,maxlength); end
        
        commandString = [varName(varIndex,:), whiteSpaceArray, [cur1;cur2;cur3], colonArray];  
        commandStringRearr = commandString';
        CVcommandString = [(commandStringRearr(:))' '; ' 'in ' CVmeritFunction ' 0 0'];  %read the MF without optimization
        [CVresponse,index1,index2] = CVcomm(CVS,CVcommandString);
        text = '=';
        index = strfind(CVresponse,text);
        ERF(m,n) = str2double(CVresponse(index+2:index+2+9));
    end
end
        
ERF(ERF>9e9) = 0; % put ray failed one to zero

figure;
surf(XP,YP,ERF);
xlabel('X');ylabel('Y');zlabel('MF value');title('MF landscape (without optimization)');colorbar;
figure;
[C,ct] = contour(XP,YP,ERF,100);
clabel(C,ct,'Color','black');
xlabel('X');ylabel('Y');zlabel('MF value');title('MF landscape (without optimization)');colorbar;
%% Plot line variation

% % MODE 1, scan around a point
% line_p = 0.2;
% line_t = 21;
% tempArrayA = linspace(varRefValue(1)*(1-line_p),varRefValue(1),line_t);
% tempArrayB = linspace(varRefValue(1)*(1+line_p),varRefValue(1),line_t);
% line_dim = 2*line_t-1;
% 
% if varRefValue(1) > 0
%     x = [tempArrayA,fliplr(tempArrayB(1:end-1))];
% else
%     x = [tempArrayB,fliplr(tempArrayA(1:end-1))];
% end

% MODE 2, scan on a designated line
% p1 = -0.1702;
% p2 = -0.165;
% p3 = -0.12;
line_dim = dim;   
xp = linspace(bound_down,bound_up,line_dim);
yp = xp;

colon='; ';
colonArray=colon(ones(sum(varIndex,1),1),:);
whiteSpace=' ';
whiteSpaceArray=whiteSpace(ones(sum(varIndex,1),1),:);

ERFlineX = zeros(line_dim,1);


for m = 1:line_dim      
    CVcomm(CVS, 'RES temp');
    
    %yp = 0
    c1 = CurRef + xp(m)/sqrt(2);
    c2 = CurRef ;
    c3 = CurRef - xp(m)/sqrt(2);
    cur1 = num2str(c1,convPrec);
    cur2 = num2str(c2,convPrec);
    cur3 = num2str(c3,convPrec);   
    
    %to make strings equal size
    maxlength = max([numel(cur1),numel(cur2),numel(cur3)]);
    if numel(cur1)<maxlength, cur1=pad(cur1,maxlength); end
    if numel(cur2)<maxlength, cur2=pad(cur2,maxlength); end
    if numel(cur3)<maxlength, cur3=pad(cur3,maxlength); end
    
    commandString = [varName(varIndex,:), whiteSpaceArray, [cur1;cur2;cur3], colonArray];  
    commandStringRearr = commandString';
    CVcommandString = [(commandStringRearr(:))' '; ' 'in ' CVmeritFunction ' 0 0'];  %read the MF without optimization
    [CVresponse,index1,index2] = CVcomm(CVS,CVcommandString);
    text = '=';
    index = strfind(CVresponse,text);
    ERFlineX(m) = str2double(CVresponse(index+2:index+2+9));
end
        
ERFlineX(ERFlineX>9e9) = 0; % put ray failed one to zero


figure;
plot(xp,ERFlineX,'DisplayName','along the X line');
xlabel('X');ylabel('MF value');title('MF landscape (without optimization)');



%% plot along the acending and decending lines



ERFlineY = zeros(line_dim,1);
for m = 1:line_dim    
    CVcomm(CVS, 'RES temp');
    
    % xp = 0
    c1 = CurRef + yp(m)/sqrt(6);
    c2 = CurRef + yp(m)*sqrt(2/3);
    c3 = CurRef + yp(m)/sqrt(6);
    cur1 = num2str(c1,convPrec);
    cur2 = num2str(c2,convPrec);
    cur3 = num2str(c3,convPrec);   
    
    %to make strings equal size
    maxlength = max([numel(cur1),numel(cur2),numel(cur3)]);
    if numel(cur1)<maxlength, cur1=pad(cur1,maxlength); end
    if numel(cur2)<maxlength, cur2=pad(cur2,maxlength); end
    if numel(cur3)<maxlength, cur3=pad(cur3,maxlength); end
            
    commandString = [varName(varIndex,:), whiteSpaceArray, [cur1;cur2;cur3], colonArray];  
    commandStringRearr = commandString';
    CVcommandString = [(commandStringRearr(:))' '; ' 'in ' CVmeritFunction ' 0 0'];  %read the MF without optimization
    [CVresponse,index1,index2] = CVcomm(CVS,CVcommandString);
    text = '=';
    index = strfind(CVresponse,text);
    ERFlineY(m) = str2double(CVresponse(index+2:index+2+9));
end
        
ERFlineY(ERFlineY>9e9) = 0; % put ray failed one to zero

hold on;
% plot(xp,ERFlineY);
% xlabel('X');ylabel('MF value');title('MF landscape (without optimization) along the scanning line');
plot(xp,ERFlineY,'DisplayName','along the Y line');
legend('show');

%% Plot a overlay with equal-MF line, double optimized results and contour profile without optimization


figure;
h = surf(XP,YP,ERFopt_two);
xlabel('X');ylabel('Y');zlabel('Optimized MF value');title('Basin of attractions (two optimization clicks)');colorbar; %shading interp;
z_max = max(max(get(h,'Zdata')));
%shading interp  %remove the grid
hold on;
pl1 = line(xp,sqrt(3)*xp,z_max*ones(1,size(xp,2)));
pl2 = line(xp,-sqrt(3)*xp,z_max*ones(1,size(xp,2)));
pl1.Color = 'red';
pl2.Color = 'red';
pl1.LineWidth = 2;
pl2.LineWidth = 2;
[C,ct] = contour(XP,YP,ERF,20);
ct.ContourZLevel = z_max;
clabel(C,ct,'Color','white');
view([0,90]);
axis([bound_down bound_up bound_down bound_up]);
pbaspect([1 1 1]);
