%**************************************************************************
% Export grain shape, grain mean-orientation from MTEX
%  to make text files containing Euler angles and volume fraction for each
%  mineral
%  
% MTEX 5.2.7
% MATLAB R2019a
% 
% Eunyoung Kim Jan. 1, 2020
    %
    % Original script from David Mainprice 
    % *************************************************************************
    % Export data from mtex
    % *************************************************************************
    %
    %  1) script imports a EBSD data set called forsterite.
    %  2) models the EBSD orientation pixels as grains using calcGrains using
    %  3) only indexed points ebsd('indexed') and the segmentation angle
    %  4) of grain boundaries, the minimum number indexed points per grain to 
    %  5) reject poorly characterized grains
    %  6) From selected grains a new EBSD is created by ebsd = ebsd(grains).
    %  7) and second grain model using calcGrains =ebsd(grains)=>ebsd('indexed')
    %  8) Fit Ellipse to grains 'area' option
    %  9) Plot map with Ellipse
    % 
    % MTEX 5.2.beta3
    %
    % David Mainprice 12/02/2019
%**************************************************************************
% clearvars; close all; clc
%**************************************************************************
%% Setup directories
dir_figs = '/Users/ekim/Documents/MATLAB/AnisEulerSC/';
% Setup filename
Title = 'Demo_MTEX_Flexible_Data_Export_Lherzolite';

% load EBSD variable containing the data
load('Original_Forsterite.mat');

% Font Size
setMTEXpref('FontSize',25);
setMTEXpref('figSize','large');

%**************************************************************************
%% Grain modelling prototcole - using indexed points and segmentation angle
% and re-calculate grain model to cleanup grain boundaries
%**************************************************************************
% segmentation angle typically 10 to 15 degrees
seg_angle = 10;
% minimum indexed points per grain between 5, 10 or 20
min_points = 20;
% Restricted to indexed points only 
[grains,ebsd.grainId,ebsd.mis2mean] = ...
    calcGrains(ebsd('indexed'),'angle',seg_angle*degree);
% Remove small grains with less than min_points indexed points 
grains = grains(grains.grainSize > min_points);
% Re-calculate grain model to cleanup grain boundaries with less than minimum indexed points 
% used ebsd points within grains having the minium indexed number of points (e.g. 10 points)
ebsd = ebsd(grains);
[grains,ebsd.grainId,ebsd.mis2mean] = ...
    calcGrains(ebsd('indexed'),'angle',seg_angle*degree);
% Plot grains
plot(grains)
drawNow(gcm,'figSize','large')
set(gcf,'renderer','zBuffer')
% % Save figure
% pname = [dir_figs,Title,'_Whole_Map_Phase_Fo_En_Di_map'];
% saveFigure ([pname,'.png'])

%**************************************************************************
%% Plot and export data 
%**************************************************************************
% Remove staircase effect from grain boundaries due to square data grid
grains = smooth(grains,5);

mineral = {'Forsterite','Enstatite','Diopside'};


nmin = length(mineral);

for i = 1:nmin
    minName = cell2mat(mineral(i));    
    %**********************************************************************
    % Fit Ellipse to grains 'area' option
    %**********************************************************************
    % a = long axis
    % b = short axis
    % omega = angle of long axis to x-axis in radians
    % use area for the fit
    [omega,a,b] = principalComponents(grains(minName),'area');
    close all
    
    figure
    plot(grains)
    hold on
    plot(grains.boundary,'lineWidth',2)
    hold on
    plotEllipse(grains(minName).centroid,a,b,omega,'lineColor','y','linewidth',2)
    hold off
    %drawNow(gcm,'figSize','large')
    set(gcf,'renderer','zBuffer')    
%     pname = [dir_figs,Title,'_',minName,'_map_with_fitted_Ellipses'];
%     saveFigure ([pname,'.png'])

    
    %**********************************************************************
    % Rose diagram using polarhistogram (0-360)
    %**********************************************************************
    % Plot X to East
    plotx2east;
    % 4:3 aspect ratio width:height
    height=800; width=round(4/3*height);
    % [left bottom width height]
    close all
    
    figure ('position',[0 0 width height])
    % 72 bins = 5 degree bins over 360
    % recent MATLAB function - Introduced in R2016b
    % polarhistogram (specify the values in radians)
    ph = polarhistogram(omega,72,'FaceColor','m','FaceAlpha',0.5);
    
    % axis properties
    ax_ph = gca;
    ax_ph.LineWidth = 2;
    ax_ph.FontSize = 40;
    ax_ph.Title.String = 'Angles of specimen X-axis at 0\circ to grain long-axis';
    ax_ph.Title.FontSize = 40;
    ax_ph.Title.FontWeight = 'normal';    
    ax_ph.ThetaLim = [0 180];
    ax_ph.RAxis.Label.String = 'Number of grains';
    
    % text properties
    tx = text;
    tx.String = minName;
    tx.FontSize = 44;
    tx.Position = [15*degree ax_ph.RLim(2)*0.5 0];

    %drawNow(gcm,'figSize','large')
    set(gcf,'renderer','zBuffer')

%     % Save figure
%     pname = [dir_figs,Title,'_',minName,'_polarhistogram_of_grain_long_axes'];
%     saveFigure ([pname,'.png'])
  
    
    %**********************************************************************
    % Export file name 
    %**********************************************************************
    fname = [dir_figs,minName,'_Demo_export.txt'];
    % List of mean orientations for one mineral
    export(grains(minName).meanOrientation,fname,'degree','Bunge')
    % Add area to S.area as addition data to your export
    S.area = grains(minName).area; 
    % Fit ellipses (S.omega= angle of S.a to X, S.a=long-axis ,S.b=short-axis)
    [S.omega,S.a,S.b] = fitEllipse(grains(minName)); % this works
    % Convert to degrees
    S.omega = S.omega./degree;
    
    %**************************************************************************
    % Export file structure
    %**************************************************************************
    %   phi1         Phi        phi2        area       omega           a           b
    %  39.91     131.882     337.894       85356     177.047     263.951     102.935
    % 23.6968     62.9256      272.38      178238    0.244368     364.628     155.596
    % 155.228     118.256     239.688     76745.2     100.294     169.707     143.947
    %
    export(grains(minName).meanOrientation,fname,S)
        
    %**********************************************************************    
    % Export only Euler angles and volume fraction for SC calculation
    %**********************************************************************    
    delimiterIn = ' ';
    headerlinesIn = 1;
    EulerAreaS = importdata(fname,delimiterIn,headerlinesIn);

    % Calculate area fraction of each grain
    SumArea = sum(EulerAreaS.data(:,4));
    % Preallocation for Euler angles and area fraction
    EulerVf = [];
    EulerVf.data = zeros(length(EulerAreaS.data(:,4)),4);
    EulerVf.colheaders = cell(1,4);
    % Euler angles (phi1 Phi phi2)
    EulerVf.data(:,1:3) = EulerAreaS.data(:,1:3);
    % Volume (area) fraction 
    EulerVf.data(:,4)   = EulerAreaS.data(:,4)./SumArea;
    % Column headers
    EulerVf.colheaders(1,1:3) = EulerAreaS.colheaders(1,1:3); % phi1,Phi,phi2
    EulerVf.colheaders(1,4)   = {'area fraction'};
    % Print Euler angles and volume (area) fraction on text file
    dlmwrite([dir_figs,minName,'_EulerVf.txt'],EulerVf.data,'delimiter','\t')    
    
    %**********************************************************************    
    % Plot grain shapes and orientations
    %**********************************************************************      
    % Calculate aspect ratio      
    grainShape.(minName).a = S.a; % long axis
    grainShape.(minName).b = S.b; % short axis
    grainShape.(minName).AspectRatio = S.a./S.b;    
    grainShape.(minName).omega = S.omega; % degree
    grainShape.(minName).area = S.area; % area fraction

    % Plot aspect ratio and orientation
    figure ('position',[0 0 width 2/3*height])    
    % 0 <= a < 1000
    [row,~] = find(grainShape.(minName).a >= 0 & ...
                   grainShape.(minName).a < 1000);               
    apGSoAl1 = plot(grainShape.(minName).omega(row,1),...
                    grainShape.(minName).AspectRatio(row,1),...
                   'Marker','o','MarkerSize',10,...
                   'MarkerFaceColor',[0.49,0.73,0.91],... % color: Aero 
                   'MarkerEdgeColor','k',...
                   'LineStyle','none');
    hold on    
    % 1000 <= a < 3500
    [row,~] = find(grainShape.(minName).a >= 1000 & ...
                   grainShape.(minName).a < 3500);  
    apGSoAl2 = plot(grainShape.(minName).omega(row,1),...
                    grainShape.(minName).AspectRatio(row,1),...
                   'Marker','o','MarkerSize',10,...
                   'MarkerFaceColor',[0,0.19,0.56],... % AirForceBlue
                   'MarkerEdgeColor','k',...
                   'LineStyle','none');
    hold on 
    
    % axis properties
    ax_apGSoAl = gca;
    ax_apGSoAl.LineWidth = 1.5;
    ax_apGSoAl.FontSize = 40;
    ax_apGSoAl.XAxis.Label.String = '\omega (angle of long a-axis to x-axis)';
    ax_apGSoAl.YAxis.Label.String = 'Aspect Ratio (a/b)';
    ax_apGSoAl.XAxis.Label.FontSize = 46;
    ax_apGSoAl.YAxis.Label.FontSize = 46;
    ax_apGSoAl.XAxis.Limits = [0 180];
    ax_apGSoAl.YAxis.Limits = [1 6];
    
    % legend properties
    lg = legend([apGSoAl1 apGSoAl2],...
                '0 < a < 1000',...
                '1000 <= a < 3500');
    lg.Location = 'northwest';    
    
    % text properties
    tx = text;
    tx.String = minName;
    tx.FontSize = 44;
    tx.Position = [120 5.2];
                
%     % Save figure
%     pname = [dir_figs,Title,'_',minName,'_AspectRatio_omega_colored'];
%     saveFigure ([pname,'.png'])
    
end

% % Save grains
% save([dir_figs,Title,'_EBSD_grains.mat'],'grains');
% save([dir_figs,Title,'_EBSD_ebsd.mat'],'ebsd');
