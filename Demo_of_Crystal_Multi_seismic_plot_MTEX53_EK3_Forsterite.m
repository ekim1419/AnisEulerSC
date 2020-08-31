%% Demo of Multi_seismic_plot_MTEX53_EK3.m
% Seismic velocities and anisotropy - single crystal plots
% Tested with MTEX 5.0.3 and MATLAB 2016a  10/05/2018 David Mainprice
% Tested with MTEX 5.3 and MATLAB R2016b & R2019a 10/06/2020 Eunyoung Kim
% by convention they are upper hemisphere for single crystals
%    and depending in the hemisphere of your polefigure plots
%    upper or lower to match your polycrytalline data
%
% C = stiffness tensor C
%
% labels for x,y,z tensor axes
% You can change these labels to fit your application
% Label_X='X';Label_Y=;'Y',Label_Z='Z';
%
% hemisphere ='upper' or hemisphere ='lower'
%
% symm = 1;% one max min marker
% symm = 2;% for generating symmetrise max and min values on the plot
%
% David Mainprice 10/05/2018
% Eunyoung Kim 10/06/2020

%% clear all
clearvars; close all; clc;
%**************************************************************************
% Define elastic constants
%**************************************************************************
% crystal symmetry - Orthorhombic mmm
cs_tensor = crystalSymmetry('mmm', [4.8 10 6],...
    'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]);

% *Import 4th rank tensor as 6 by 6 matrix*
%
% Forsterite elastic stiffness (Cij) tensor in GPa
% Reference of elastic constants
%     Isaak, D. G., Anderson, O. L., Goto, T., & Suzuki, I. (1989), 
%     Elasticity of single-crystal forsterite measured to 1700 K, 
%     J. Geophys. Res. 94, 5895-5906.
% Reference of crystal structure
%     Tommasi, A., Vauchez, A., & Ionov, D. A. (2008),
%     Deformation, static recrystallization, and reactive melt transport 
%     in shallow subcontinental mantle xenoliths 
%     (Tok Cenozoic volcanic field, SE Siberia), 
%     Earth and Planetary Science Letters, 272(1), 65-77.
%
% Enter tensor as 6 by 6 matrix,M line by line.
M = [[  328.00   69.00   69.00    0.00    0.00    0.00];...
    [    69.00  200.00   73.00    0.00    0.00    0.00];...
    [    69.00   73.00  235.00    0.00    0.00    0.00];...
    [     0.00    0.00    0.00   66.70    0.00    0.00];...
    [     0.00    0.00    0.00    0.00   81.30    0.00];...
    [     0.00    0.00    0.00    0.00    0.00   80.90]];

% Define density (g/cm3)
rho = 3.221;

% Define tensor object in MTEX
% Cij -> Cijkl - elastic stiffness tensor
C = stiffnessTensor(M,cs_tensor,'density',rho);
% Run tensor symmetry check
if(C.checkSymmetry == false)
warning('MTEX:tensor','Tensor does not have the correct symmetry')
else
disp('MTEX:Tensor symmetry is correct')    
end
%**************************************************************************
%% plot single crystal - Enstatite
%**************************************************************************
% plot x to east
plotx2east
% plot a to east
plota2east
%
% You can change these labels to fit your application
Label_X='a[100]'; Label_Y='b[010]';Label_Z='c[001]';
% Label_X='X'; Label_Y='Y';Label_Z='Z';
%
% 'upper'; or 'lower'
hemisphere ='upper';
%
% symm = 1 one max min marker
% symm = 2 for generating symmetrise max and min values on the plot 
% 
symm = 1;% one max min marker
%
% plot font_size = set your fontsize  25,20,15  etc
font_size = 38;
% marker_size (max min)
marker_size = 30;
%
setMTEXpref('FontSize',font_size)
%
Multi_seismic_plot_MTEX53_EK3(C,Label_X,Label_Y,Label_Z,hemisphere,symm,font_size)
%
drawNow(gcm,'figSize','medium')
set(gcf,'renderer','zBuffer')

function [  ] = Multi_seismic_plot_MTEX53_EK3(C,Label_X,Label_Y,Label_Z,hemisphere,symm,font_size)
%% Modified by Ralf 13/03/2019
%
% Tested with MTEX 5.1.1 and MATLAB 2016a, David Mainprice 13/03/2019
% Tested with MTEX 5.3 and MATLAB R2016b & R2019a, Eunyoung Kim 10/06/2020
%
% Plot matrix of anisotropic seismic properties
%        1 Vp        2 AVs      3 S1 polarizations  4 Vs1
%        5 Vs2       6 dVs      7 Vp/Vs1            8 Vp/Vs2
%
% INPUT
%
% MTEX tensor defined with density
% C = stiffnessTensor(M,cs_tensor,'density',rho)
% User defined labels at MTEX X,Y, and Z positions reference frame
% on Vp plot (top right plot)
% You can define the labels to suit your application single crystal or
% polycrystalline stiffness tensor with triclinic symmetry('1')
%
% X_Label = for example 'a', 'X1' or 'X'
% Y_Label = for example 'b', 'X2' or 'Y'
% Z_Label = for example 'c', 'X3' or 'Z'
%
% hemisphere = 'upper' or 'lower'
%
% New feature is generating symmetically equivalent max and min Vp,Vs1,vs2
% etc using the harmonic method for symmetry of single crystals
%
% symm = 1 only one max and min values on the plot
% symm = 2 for generating symmetrise max and min values on the plot
%
% plot font_size = set your fontsize  25,20,15  etc
%
%**************************************************************************
% Remember to set the position of x-axis in specimen coordinates
% before calling this function
% plot2xeast,  plot2xnorth, plot2xwest or plot2xsouth
% to suit your single crystal or polycrystalline elastic stiffness tensor
%**************************************************************************
%
%
% Compute seismic velocities as functions
% using option velocity(C,[]) or C.velocity 
% BOTH generate velocities on sphere using harmonic method
% OR using the traditional velocity(C,x) where x is the propagation
% directions , many directions or a grid.
% such as XY_grid = equispacedS2Grid('upper','resolution',1*degree)
% [vp,vs1,vs2,pp,ps1,ps2] = velocity(C,XY_grid)
% Using a user defined grid use classical method not harmonic method.

% Seismic velocities as functions on the sphere with harmonic method
% is recommended for a smooth representation that actuately respects
% the symmetry of single crystals and can be plotted in 3D.
%
[vp,vs1,vs2,pp,ps1,ps2] = C.velocity('harmonic');
%
% you can sample the velocity values in any x,y,z direction by using
% vp_xyz = eval(vp,vector3d(x,y,z)) or
% vp_rho_theta = eval(vp,vector3d('rho',45*degree,'theta',90*degree))
% etc
%
%**************************************************************************
% Plotting section
%**************************************************************************

% set colour map to seismic color map : red2blueColorMap
% Red = slow to Blue = fast in Seismology
%
setMTEXpref('defaultColorMap',red2blueColorMap)

% some options
blackMarker = {'Marker','s','MarkerSize',10,'antipodal',...
  'MarkerEdgeColor','white','MarkerFaceColor','black','doNotDraw'};
whiteMarker = {'Marker','o','MarkerSize',10,'antipodal',...
  'MarkerEdgeColor','black','MarkerFaceColor','white','doNotDraw'};

% some global options for the titles
titleOpt = {'visible','on','color','k','FontSize',font_size};

% Setup multiplot
% define plot size [origin X,Y,Width,Height]
mtexFig = mtexFigure('position',[0 0 1000 1000]);
% set up spacing between subplots default is 10 pixel
%mtexFig.innerPlotSpacing = 20;

% Standard Seismic plot with 8 subplots in 3 by 3 matrix
%
% Plot matrix layout 3 by 3
%--------------------------------------------------------------------------
%        1 Vp        2 AVs      3 S1 polarizations
%        4 Vs2       5 Vp/Vs1   6 Vp/Vs2
%--------------------------------------------------------------------------
%**************************************************************************
% 1 Vp : Plot P-wave velocity (km/s)
%**************************************************************************

% Plot P-wave velocity (km/s)
plot(vp,'contourf','complete',hemisphere)
mtexColorbar
mtexTitle('Vp (km/s)',titleOpt{:})

% extrema
[maxVp, maxVpPos] = max(vp);
[minVp, minVpPos] = min(vp);

% percentage anisotropy
AVp = 200*(maxVp-minVp)./(maxVp+minVp);

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
  plot(maxVpPos.symmetrise,blackMarker{:})
  plot(minVpPos.symmetrise,whiteMarker{:})
else
  plot(maxVpPos(1),blackMarker{:})
  plot(minVpPos(1),whiteMarker{:})    
end
hold off

% subTitle
xlabel(['*Vp Anisotropy = ',num2str(AVp,'%6.1f')],titleOpt{:})
%**************************************************************************
% 2 AVs : Plot S-wave anisotropy percentage for each propagation direction
% defined as AVs = 200*(Vs1-Vs2)/(Vs1+Vs2)
%**************************************************************************

% create a new axis
nextAxis

% Plot S-wave anisotropy (percent)
AVs = 200*(vs1-vs2)./(vs1+vs2);
plot(AVs,'contourf','complete',hemisphere);
mtexColorbar
mtexTitle('S-wave anisotropy (%)',titleOpt{:})

% Max percentage anisotropy
[maxAVs,maxAVsPos] = max(AVs);
[minAVs,minAVsPos] = min(AVs);

xlabel(['*Max Vs Anisotropy = ',num2str(maxAVs,'%6.1f')],titleOpt{:})

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
  plot(maxAVsPos.symmetrise,blackMarker{:})
  plot(minAVsPos.symmetrise,whiteMarker{:})
else
  plot(maxAVsPos(1),blackMarker{:})
  plot(minAVsPos(1),whiteMarker{:})
end
hold off

% mark crystal axes - Label_X,Label_Y,Label_Y
t = text([xvector,yvector,zvector],{Label_X,Label_Y,Label_Z},...
  'backgroundcolor','w','doNotDraw','FontSize',17);
% (a*(100)) properties
t(1).Position = [1.3 0 0]; 

%**************************************************************************
% 3 Vs1 : Plot Vs1 velocities (km/s) with polarization
%**************************************************************************
% create a new axis
nextAxis

plot(vs1,'contourf','doNotDraw','complete',hemisphere);
mtexColorbar
mtexTitle('Vs1 (km/s)',titleOpt{:})

% Percentage anisotropy
[maxS1,maxS1pos] = max(vs1);
[minS1,minS1pos] = min(vs1);
AVs1=200*(maxS1-minS1)./(maxS1+minS1);

xlabel(['*Vs1 Anisotropy = ',num2str(AVs1,'%6.1f')],titleOpt{:}) 

hold on
plot(ps1,'linewidth',2,'color','black')

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
   plot(maxS1pos.symmetrise,blackMarker{:})
   plot(minS1pos.symmetrise,whiteMarker{:})
else
   plot(maxS1pos(1),blackMarker{:})
   plot(minS1pos(1),whiteMarker{:})
end
hold off
%**************************************************************************
% 4 Vs2 : Plot Vs2 velocities (km/s)
%**************************************************************************
% create a new axis
nextAxis

plot(vs2,'contourf','doNotDraw','complete',hemisphere);
mtexColorbar
mtexTitle('Vs2 (km/s)',titleOpt{:})

% Percentage anisotropy
[maxS2,maxS2pos] = max(vs2);
[minS2,minS2pos] = min(vs2);
AVs2=200*(maxS2-minS2)./(maxS2+minS2);
xlabel(['*Vs2 Anisotropy = ',num2str(AVs2,'%6.1f')],titleOpt{:})

hold on
plot(ps2,'linewidth',2,'color','black')

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
  plot(maxS2pos.symmetrise,blackMarker{:})
  plot(minS2pos.symmetrise,whiteMarker{:})
else
  plot(maxS2pos(1),blackMarker{:})
  plot(minS2pos(1),whiteMarker{:})
end

hold off
%**************************************************************************
% 5 Vp/Vs1 : Plot Vp/Vs1 ratio (no units)
%**************************************************************************
% create a new axis
nextAxis

vpvs1 = vp./vs1;
plot(vpvs1,'contourf','complete',hemisphere);
mtexColorbar
mtexTitle('Vp/Vs1',titleOpt{:})

% Percentage anisotropy
[maxVpVs1,maxVpVs1Pos] = max(vpvs1);
[minVpVs1,minVpVs1Pos] = min(vpvs1);
AVpVs1=200*(maxVpVs1-minVpVs1)/(maxVpVs1+minVpVs1);

xlabel(['*Vp/Vs1 Anisotropy = ',num2str(AVpVs1,'%6.1f')],titleOpt{:})

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
  plot(maxVpVs1Pos.symmetrise,blackMarker{:})
  plot(minVpVs1Pos.symmetrise,whiteMarker{:})
else
  plot(maxVpVs1Pos(1),blackMarker{:})
  plot(minVpVs1Pos(1),whiteMarker{:})
end
hold off
%**************************************************************************
% 6 Vp/Vs2 : Plot Vp/Vs2 ratio (no units)
%**************************************************************************
% create a new axis
nextAxis

vpvs2 = vp./vs2;
plot(vpvs2,'contourf','complete',hemisphere);
mtexColorbar
mtexTitle('Vp/Vs2',titleOpt{:})

% Percentage anisotropy
[maxVpVs2,maxVpVs2Pos] = max(vpvs2);
[minVpVs2,minVpVs2Pos] = min(vpvs2);
AVpVs2=200*(maxVpVs2-minVpVs2)/(maxVpVs2+minVpVs2);

xlabel(['*Vp/Vs2 Anisotropy = ',num2str(AVpVs2,'%6.1f')],titleOpt{:})

% mark maximum with black square and minimum with white circle
hold on
if(symm == 2)
  plot(maxVpVs2Pos.symmetrise,blackMarker{:})
  plot(minVpVs2Pos.symmetrise,whiteMarker{:})
else
  plot(maxVpVs2Pos(1),blackMarker{:})
  plot(minVpVs2Pos(1),whiteMarker{:})
end
hold off
%**************************************************************************
% add colorbars to all plots and save plot
%**************************************************************************
mtexColorbar
% mtexColorbar('Location','manual','Position',[0.4 0.5 0.5 0.3])
drawNow(gcm,'figSize','large')
set(gcf,'renderer','zBuffer')
saveFigure('Multi_seismic_plot_MTEX53_EK2.eps')
%**************************************************************************
% reset to default MTEX colormap
%**************************************************************************
setMTEXpref('defaultColorMap',WhiteJetColorMap);
end



