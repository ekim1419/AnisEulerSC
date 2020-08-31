%
% MATLAB version tested 
%     with 1) MTEX 5.3
%          2) MATLAB R2019a
% 
% by Eunyoung Kim (6/07/2020)
% E-mail: brilliant@snu.ac.kr
%
%**************************************************************************
% AnisEulerSC.m
%
% Seismic properties from CPO
% using Bunge(1982) Euler angles g = (phi1,PHI,phi21)
% wave propagation in anisotropic polycrystals
% (Plane wave and Group velocities)
%
% Self-consistent version - mechanical and scattering methods
%
%
%             LAST REVISED FORTRAN CODE : 7/06/2017
%                  subroutine Error changed
%
%
%                       DAVID MAINPRICE,
%                      Geosciences Montpellier,
%            Universite Montpellier,Place E.Bataillon,
%                34095 MONTPELLIER cedex 05,FRANCE.
%                 e-mail : David.Mainprice@gm.univ-montp2.fr
%
% The orientation of the elastic constants cartesian reference frame 
%   axes X1,X2 and X3 are to be given in geographic coordinates 
%   by phi1,PHI,phi2.
%
% N.B. For low symmetry minerals the reference frame for
%   the Euler angles and the Elastic constants may not be the same.
%
% Elastic constant right-handed reference frame
%   conventions for axes X1,X2,X3
%   cubic,tetragonal,othorhombic:A,B,C
%   hexagonal,trigonal:A1,M=(C x A1),C
%   monoclinic:A,B,C* (alkali feldspar)
%              A*,B,C (diopside,augite,hornblende) *used here*
%   triclinic :A,B,C* (plagiclase) no general convention
%
% Input
%        Single crystal elastic constants in GPa
%        Orientation data in Bunge Euler angles
% Output files
%        *-Print   : Printed results
%        *-GPa     : Rock elastic constants
%        *-SC-Cij  : Rock elastic constants for every SC iteration
%
%**************************************************************************

% Load input data file ('Title'_Inputs_ModelSC.mat)
load('AnisEulerSC_Example_Lherzolite_Inputs_ModelSC.mat');

% Preallocation 
phi1max = 0; 
PHImax  = 0; 
phi2max = 0; 
 
xi    = zeros(1,3); cxr   = zeros(1,3); x1c   = zeros(1,3);
xr    = zeros(1,3); cyr   = zeros(1,3); x3c   = zeros(1,3); 
zr    = zeros(1,3); czr   = zeros(1,3); eaxis = zeros(1,3);
xhkl1 = zeros(1,3);  
xhkl2 = zeros(1,3); 
xhkl3 = zeros(1,3); 

r      = zeros(3,3); ee    = zeros(3,3);
rbacki = zeros(3,3);  
rback  = zeros(3,3);

C     = zeros(6,6); Cg    = zeros(6,6);
Csc   = zeros(6,6); CscM  = zeros(6,6); CscR = zeros(6,6);
EC    = zeros(6,6); ES    = zeros(6,6); 
CMV  = zeros(6,6); SMR  = zeros(6,6); CVRH = zeros(6,6);
CSUM = zeros(6,6); SSUM = zeros(6,6);

Am    = zeros(6,6); 
CiAm  = zeros(6,6); sumt  = zeros(6,6);
sum1  = zeros(6,6); sumtm = zeros(6,6);
sum2  = zeros(6,6); tmatm = zeros(6,6);
 
C4 = zeros(3,3,3,3); 
 
sn1 = zeros(1,3); vg1 = zeros(1,3); pm1 = zeros(1,3);
sn2 = zeros(1,3); vg2 = zeros(1,3); pm2 = zeros(1,3);
sn3 = zeros(1,3); vg3 = zeros(1,3); pm3 = zeros(1,3);

column = cell(1,22); 
column(1,:) = {'No',...
               'C11 SC','C12 SC','C13 SC','C14 SC','C15 SC','C16 SC',...
                        'C22 SC','C23 SC','C24 SC','C25 SC','C26 SC',...
                                 'C33 SC','C34 SC','C35 SC','C36 SC',...
                                          'C44 SC','C45 SC','C46 SC',...
                                                   'C55 SC','C56 SC',...
                                                            'C66 SC'};


format320 = ['\n ',repmat('%9.4f',1,6),' \n'];
rad = pi/180;
%**************************************************************************
%                         INPUT SECTION (User defined)                    *
%**************************************************************************



% Variables for SC calculation
dir_Input  = ModelSC.dir_Input;
dir_Output = ModelSC.dir_Output;

% Specify output file name (excluding file extension)
Title   = ModelSC.Title;

% SC scheme:
%   1 = Mechanical - Willis 1977 (most stable)
%   2 = Scattering - Gubernatis & Krumhansl 1975 
ioptsc  = ModelSC.ioptsc;

% Initial value for SC scheme:
%   1 = Voigt (best when fluids present)
%   2 = Reuss
%   3 = Voigt-Reuss-Hill
isch    = ModelSC.isch;

% Euler file: 
%   1 = phi1,PHI,phi2
%   2 = no,phi1,PHI,phi2
%   3 = phi1,PHI,phi2,volume fraction
%   4 = no,phi1,PHI,phi2,volume fraction
ioptn   = ModelSC.ioptn;

% Specify origin of Euler angles:
%   1 = Universal stage (UMII)
%   2 = EBSD (Channel+)
ieuler  = ModelSC.ieuler;

% Xs of azimuth (in range 0 to 360 degrees)
% if ieuler == 2, EBSD (Channel+)
% 
%     Orientation of Euler angle reference frame
%     (Xs Ys Zs) in sample coordinates
%     Zs is upward and normal from the specimen surface
%     and Xs Ys and Zs are orthogonal then
%     only the position of Xs needs to be defined
%     by the its azimuth(0-360 degrees clockwise
%     from Top or North)
%     
%               Montpellier EBSD geometry
%     ITop/North
%     I__________________________________
%     I       Top when tilted 70 degrees I
%     I                                  I
%     I   Ys              North=0        I    _______
%     I   I                 I            I   I       I
%     I   I      West=270<--I-->East=90  I   I  Si   I
%     I   I                 I            I   IcrystalI
%     I   I____Xs         South=180      I   I_______I
%     I                                  I
%     I Thin Section in Al sample holder I
%     I    Bottom when tilted 70 degrees I
%     I__________________________________I
%     I
%     IBottom/South
%     (you can type a value different from 90 degrees
%     to correct for sample misalignment
%     80 will give anticlockwise 10 degree correction)
%     A negative value of azimuth of Xs will give access
%     to the HELP screen (e.g.-90)
%     *** normally Montpellier users should type 90 ***
%     *** or 0 when pole to foliation is to the North***
% 
%     if XsAz < 0
%         
%         * extra screen information *
%         
%         Euler angles define the relative orientation of
%         a crystal (XcYcZc) and a sample (XsYsYs) right-
%         handed reference frames.The sample XsYsYs used to
%         define the Euler angles may have any orientation.
%         Commonly the sample reference frame is defined as:
%         Xs in the North (azimuth 0),
%         Ys West (azimuth 270)
%         Zs (up) for example see
%         Bunge 1982 Texture Analysis in Material Science,
%         Butterworths fig.2.22 p28
%         or
%         Xs in the East (azimuth 90),
%         Ys North (azimuth 0)
%         Zs (up) for example see
%         Kocks,Tome and Wenk 1998 Texture and Anisotropy,
%         Cambridge University Press fig.15 p64  
%
%         You may also use an arbitary value of the azimuth
%         of Xs to rotate pole figures to the desired
%         plotted orientation
XsAz    = ModelSC.XsAz;

% Max. Number SC iterations (e.g. 40);
nsc     = ModelSC.nsc;

% Max. %error for Greens Tensor (e.g. 0.1);
errmax  = ModelSC.errmax;

% Max. %error for SC Convergence (e.g. 0.1);
conmax  = ModelSC.conmax;



% Create text files to save outputs
%
% FILES: Pfile (Print)
%        Rfile (GPa, Rock elastic constants file)
%        Cfile1 (for SC analysis of all minerals)
Pfile  = [dir_Output,deblank(Title),'-SC-Print.txt'];
Rfile  = [dir_Output,deblank(Title),'-SC-GPa.txt'];
Cfile  = [dir_Output,deblank(Title),'-SC-Cij.txt'];

fid12 = fopen(Pfile,'w+');
fid14 = fopen(Rfile,'w+');
fid17 = fopen(Cfile,'w+');

fprintf(fid12,['\n\n  Seismic Anisotropy Analysis',...
    repmat(' ',1,8),'%s\n\n'],Title);
for i = 1:length(column)
    fprintf(fid17,'%8s \t',cell2mat(column(i))); 
end
fprintf(fid17,' \n');



% Preallocation for all minerals
% rho    = zeros(nmin,1);
% M      = zeros(nmin,6,6);
% abc    = zeros(nmin,3); % axis a, b, c 
% abg    = zeros(nmin,3); % alpha, beta, gamma
% index1 = zeros(nmin,1);

% Specify the name of minerals 
mineral = ModelSC.PhaseList; % list of minerals
nmin    = ModelSC.nPhase; % total number of minerals

% Specify the volume faction of an aggregate in Range 0 to 1
vf = zeros(nmin,1);
% fprintf(' Volume Fraction Total = %6.4f',sum(modelSC.vf));

% FILE Ofile: Euler Data file
Ofile = cell(nmin,1);

isymm  = zeros(nmin,1);
ixoptm = zeros(nmin,1);
izoptm = zeros(nmin,1);
x1e    = zeros(nmin,3);
x3e    = zeros(nmin,3);
ioptel = zeros(nmin,1);
ioptoe = zeros(nmin,1);
axis   = zeros(nmin,3);
al     = zeros(nmin,1);
xl     = zeros(nmin,1);
af     = zeros(nmin,1);
xf     = zeros(nmin,1);

% Variables for SC calculation
for i = 1:nmin
    minName = cell2mat(mineral(i));
    minNo = ModelSC.Phase.(minName).No;
    
    vf(minNo)     = ModelSC.Phase.(minName).vf;
    
    Ofile(minNo)  = {ModelSC.Phase.(minName).Ofile};
    
    %**********************************************************************    
    % MINERAL information
    %**********************************************************************
    % Specify Euler space for crystal- *triclinic sample symmetry*
    %  11 Proper point groups (Schoenflies - International)
    %   Crystal symmetry                  phi1 : PHI : phi2
    %   0 = Triclinic (C1 - 1)             360 : 180 : 360
    %   1 = Monoclinic (C2 - 2)            360 : 180 : 180
    %   2 = Orthorhombic (D2 - 222)        180 : 180 : 180
    %   3 = Trigonal(D3 - 32) (2-fold//Y)  180 : 180 : 120
    %   4 = Trigonal(C3 - 3)               360 : 180 : 120
    %   5 = Tetragonal (D4 - 422)          180 : 180 :  90
    %   6 = Tetragonal (C4 - 4)            360 : 180 :  90
    %   7 = Hexagonal  (D6 - 622)          180 : 180 :  60
    %   8 = Hexagonal  (C6 - 6)            360 : 180 :  60
    %   9 = Cubic (O - 432)                180 : 180 :  90
    %  10 = Cubic (T - 23)                 360 : 180 :  90
    %
    %  N.B. 0 = No crystal symmetry operations
    isymm(minNo)  = ModelSC.Phase.(minName).isymm;    
    
    % Elastic reference directions X1,X2 and X3 for elastic constants
    %
    % Cij elastic reference directions X1, X2 and X3 for elastic constants
    % Cubic/Othorhombic/Tetragonal Xc=100 Yc=010 Zc=001
    % Hexagonal/Trigonal Xc=a(2-1.0) Yc=m(01.0) Zc=(00.1) (*elastic*)
    % Monoclinic Xc=a* Yc=010 Zc=001 (*clinopyroxenes*)
    % Triclinic Xc=100 Yc=010 Zc=c*  (*plagioclase*)
    % ixoptm or izoptm=1 for uvw e.g.[100]
    % ixoptm or izoptm=2 for hkl e.g.(100)
    %
    % Specify Cij elastic reference directions X1,X2 and X3
    %    user defined Crystal Reference
    %    Directions(X1,X2,X3)
    %    As X1,X2 and X3 must be orthogonal
    %    only X1 and X3 need be defined
    %
    % X1
    %   1 = dirn.[UVW]
    %   2 = pole normal (HKL)
    ixoptm(minNo) = ModelSC.Phase.(minName).ixoptm;
    
    % if ixoptm == 1
    %   X1 direction indices - U,V,W 
    % else
    %   X1 pole normal indices - H,K,L
    x1e(minNo,:)  = ModelSC.Phase.(minName).x1e;
    
    % X3
    %   1 = dirn.[UVW]
    %   2 = pole normal (HKL)
    izoptm(minNo) = ModelSC.Phase.(minName).izoptm;
    
    % if izoptm == 1
    %   X3 direction indices - U,V,W,
    % else
    %   X3 pole normal indices - H,K,L
    x3e(minNo,:)  = ModelSC.Phase.(minName).x3e;  
    
    % Specify grain shape definition
    % Ellipsoid semi-axes a1,a2,a3
    %   sphere            a1=a2=a3
    %   oblate ellipsoid  a1=a2>a3
    %   prolate ellipsoid a1=a2<a3
    % 1 = Sphere
    % 2 = Ellipsoid
    ioptel(minNo) = ModelSC.Phase.(minName).ioptel;
    
    % Specify orientation of ellipsoid semi-axes
    % 1 = // tensor axes X1 X2 X3 of each crystal
    %         A1 // X1   (elastic tensor axes of crystal)
    %         A2 // X2    
    %         A3 // X3
    %         X1=100 olivine, X1=a in calcite etc.
    % 2 = // user defined orientation in PF coordinates (e.g. fluids)
    ioptoe(minNo,1) = ModelSC.Phase.(minName).ioptoe;
    
    % Ellipsoid - parallel to tensor axes X1 X2 X3
    % Specify Ellipsoid Semi-axes lengths: [A1,A2,A3]
    %
    % Spherical grains for testing : result close VRH and Geometric mean
    axis(minNo,:) = ModelSC.Phase.(minName).axis;
    
    % Specify Foliation with pole given by
    % AZ (0-360) and INC (0-90)
    %   al: A1 AZ Azimuth of ellipsoid semi-axis
    %   xl: A1 INC Inclination of ellipsoid semi-axis
    %   af: A3 AZ Azimuth of ellipsoid semi-axis
    %   xf: A3 INC Inclination of ellipsoid semi-axis
    al(minNo,1)   = ModelSC.Phase.(minName).al;
    xl(minNo,1)   = ModelSC.Phase.(minName).xl;
    af(minNo,1)   = ModelSC.Phase.(minName).af;
    xf(minNo,1)   = ModelSC.Phase.(minName).xf;                 
                    
end
 

%**************************************************************************
%                                                                         *
%     MINERAL LOOP                                                        *
%                                                                         *
%**************************************************************************
% Preallocation for minerals
ngrains = zeros(nmin,1); 
remin   = zeros(nmin,3,3);
ECM     = zeros(nmin,6,6); 
Cmin    = zeros(nmin,6,6); 
Smin    = zeros(nmin,6,6);
vstore  = zeros(nmin,500000); % ** Dimension for 500000 points **
rstore  = zeros(nmin,500000,3,3);
sc      = zeros(500000,6,6);

C_Voigt = [];
C_Reuss = [];
C_Hill  = [];
    
for minn = 1:nmin
    minName = cell2mat(mineral(minn));
    minNo = ModelSC.Phase.(minName).No;
    ECM(minNo,:,:) = ModelSC.Phase.(minName).M;
        
    if ioptel == 1
    % Sphere
        axis(minNo,:) = [1,1,1];
        ioptoe(minNo,1) = 1;
    else
        ioptel = 2;
        % Ellipsoid - parallel to tensor axes X1 X2 X3
        % ioptoe(minNo,1) == 1
        if ioptoe(minNo,1) == 2
        % Ellipsoid - parallel to geographic (pole figure) specimen direction
            disp(' User defined ellipsoid orientation');
            disp(' * same for all grains of this');
            disp(' composition in sample PF coordinates *');
            disp(' A1 should be 90 degrees to A3');
            disp(' ');
            disp(' *************** Foliation normal to Z');
            disp('               North=0 AZ =0');
            disp('                     .');
            disp('                    /I\');
            disp('                     Z');
            disp('                   **A1**');
            disp('                *         *');
            disp('              *              *');
            disp('             *                 *');
            disp('            *                    *');
            disp('           *                      *  AZ=90');
            disp('East=270--*ffffffffff*Y*fffffffff A3*-X');
            disp('           *                     * West=90');
            disp('            *                   *');
            disp('             *                *');
            disp('               *             *');
            disp('                 *         *');
            disp('                   *****');
            disp('              South=180 AZ=180');
            disp(' f = trace of foliation');
            disp(' A1 should be 90 degrees to A3');
            disp(' -----------------------------------------');
            disp(' A1 = Pole to foliation  AZ =90   INC = 0');
            disp(' A3 = Lineation          AZ = 0   INC = 0');
            disp(' -----------------------------------------');
            disp(' ');
            disp(' ***************** Foliation normal to X');
            disp(' ');
            disp('                    North=0 AZ =0');
            disp('                          .');
            disp('                         /I\');
            disp('                          Z');
            disp('                        **A1**');
            disp('                     *    f      *');
            disp('                  *       f         *');
            disp('                *         f          *');
            disp('               *          f           *');
            disp('              *                       * AZ=90');
            disp('  East=270  --*          *Y*          A3*-X');
            disp('              *                       *West=90');
            disp('               *          f          *');
            disp('                *         f         *');
            disp('                  *       f       *');
            disp('                     *    f     *');
            disp('                        *****');
            disp('                       South=180 AZ=180');
            disp(' f = trace of foliation');
            disp(' A1 should be 90 degrees to A3');
            disp(' -----------------------------------------');
            disp(' A1 = Pole to foliation  AZ = 0  INC = 0');
            disp(' A3 = Lineation          AZ =90  INC = 0');
            disp(' -----------------------------------------');
            disp(' ');
            disp(' ***************** Foliation normal to Y');
            disp('');
            disp('                    North=0 AZ =0');
            disp('                          .');
            disp('                         /I\');
            disp('                          Z');
            disp('                        **A1**');
            disp('                     *    f   f   *');
            disp('                  *    f  f         *');
            disp('                *   f     f    f   f  *');
            disp('               *          f            *');
            disp('              *   f   f   f     f   f  *AZ=90');
            disp('  East=270  --*          *Y*           A3*-X ');
            disp('              *    f      f    f       *West=90');
            disp('               *      f   f      f    *');
            disp('                *   f     f   f      *');
            disp('                  *       f        *');
            disp('                     *    f      *');
            disp('                        *****');
            disp('                   South=180 AZ=180');
            disp(' f = trace of foliation');
            disp(' A1 should be 90 degrees to A3');
            disp(' -----------------------------------------');
            disp(' A1 = Pole to foliation  AZ =90  INC =  0');
            disp(' A3 = Lineation          AZ =90  INC = 90');
            disp(' -----------------------------------------');
            disp(' ');
            disp(' Foliation with pole given by');
            disp(' AZ (0-360) and INC (0-90) ');
            disp(' ');
            % X=L Z=FN  Y= Z x X
            [RE] = azr(al(minNo),xl(minNo),af(minNo),xf(minNo));

            disp('Shape ellipsoid rotation matrix');
            disp(RE);
            remin(minNo,:,:) = RE(:,:);
            disp(' ');
        end
    end

    % FILE = Ofile: Euler Data file
    fid10 = fopen(cell2mat(Ofile(minNo)),'r');

    % Bond transformation matrices TC and CT
    [CT,TC] = bond(ModelSC.Phase.(minName).abc,ModelSC.Phase.(minName).abg,-1);        
    
    % Euler symmetry 
    [isym,phi1max,PHImax,phi2max] = ...
        eulsym(phi1max,PHImax,phi2max,isymm(minNo));
    
    % Crystal Reference Directions(Xc,Yc,Zc) for Euler angle determination
    switch ieuler
        case 1
        [xr,ixoptr,zr,izoptr] = refustage(isym,xr,zr);
        case 2
        [xr,ixoptr,zr,izoptr] = refebsd(xr,zr);
    end   
    
    % Convert Euler reference to cartesian (ortho-normal)
    [cxr] = ortho(CT,TC,ixoptr,xr);
    [czr] = ortho(CT,TC,izoptr,zr);
    
    % Vector cross-product 3 x 1 = 2 right-handed
    [cyr] = xprod(czr,cxr);
    
    % Convert elastic reference to cartesian (ortho-normal)
    [x1c] = ortho(CT,TC,ixoptm(minNo),x1e(minNo,:));
    [x3c] = ortho(CT,TC,izoptm(minNo),x3e(minNo,:));

    % Calculate angles --> XH1,XK1,XL1 (elastic X1) pole in cartesian  
    %                      Euler ref Xr,Yr,Zr
    % Xref to Phkl
    [ang] = angleAB(cxr,x1c);
    xhkl1(1) = cos(ang*rad);
    
    % Yref to Phkl
    [ang] = angleAB(cyr,x1c);
    xhkl1(2) = cos(ang*rad);
    
    % Zref to Phkl
    [ang] = angleAB(czr,x1c);
    xhkl1(3) = cos(ang*rad);
    
    %  Convert to direction cosines
    xm = sqrt(xhkl1(1)^2 + xhkl1(2)^2 + xhkl1(3)^2);
    xhkl1 = xhkl1/xm;

    % Calculate angles --> XH3,XK3,XL3 (elastic X3) pole in cartesian 
    %                      Euler ref Xr,Yr,Zr
    % Xref to Phkl
    [ang] = angleAB(cxr,x3c);
    xhkl3(1) = cos(ang*rad);
    
    % Yref to Phkl
    [ang] = angleAB(cyr,x3c);
    xhkl3(2) = cos(ang*rad);
    
    % Zref to Phkl
    [ang] = angleAB(czr,x3c);
    xhkl3(3) = cos(ang*rad);
    
    % convert to direction cosines
    xm = sqrt(xhkl3(1)^2 + xhkl3(2)^2 + xhkl3(3)^2);
    xhkl3 = xhkl3/xm;
    
    [xhkl2] = xprod(xhkl3,xhkl1);
    
    % Write to screen
     fprintf([' EULER-ELASTIC Transform matrix \n',...
             repmat('  %10.4f',1,3),' \n',...
             repmat('  %10.4f',1,3),' \n',...
             repmat('  %10.4f',1,3),' \n'],xhkl1,xhkl2,xhkl3);        
         
    % Fill Transform matrix by Columns    
    for i = 1:3
        ee(i,1) = xhkl1(i);
        ee(i,2) = xhkl2(i);
        ee(i,3) = xhkl3(i);
    end
    
    % Invert Cij to Sij
    EC = squeeze(ECM(minNo,:,:));
    ES = inv(EC);

    % Read Euler angle file
    frewind(fid10);
    ngrain = 0;
    
    %**********************************************************************
    %                                                                     *
    %               EULER ANGLE LOOP                                      *
    %                                                                     *
    %**********************************************************************

    % READ euler angles
    %
    % Euler file contents into array: 
    switch ioptn
        case 1
            EulerArray = fscanf(fid10,'%f %f %f',[3 inf]);
            EulerArray = EulerArray';
            [EArow,EAcol] = size(EulerArray);
        case 2
            EulerArray = fscanf(fid10,'%d %f %f %f',[4 inf]);
            EulerArray = EulerArray';    
            [EArow,EAcol] = size(EulerArray);
        case 3
            EulerArray = fscanf(fid10,'%f %f %f %f',[4 inf]);
            EulerArray = EulerArray';
            [EArow,EAcol] = size(EulerArray);
        case 4
            EulerArray = fscanf(fid10,'%d %f %f %f %f',[5 inf]);
            EulerArray = EulerArray';    
            [EArow,EAcol] = size(EulerArray);
    end          
    
    nEArow = 0;        
    while EArow > 0   
        nEArow = nEArow + 1;
        switch ioptn
            case 1
                phi1 = EulerArray(nEArow,1);
                PHI  = EulerArray(nEArow,2);  
                phi2 = EulerArray(nEArow,3); 
                vfm  = 1; 
            case 2
                nb   = EulerArray(nEArow,1);
                phi1 = EulerArray(nEArow,2);
                PHI  = EulerArray(nEArow,3);    
                phi2 = EulerArray(nEArow,4); 
                vfm  = 1; 

        % VF = volume fractions of given mineral phase, sum to 1.0
            case 3
                phi1 = EulerArray(nEArow,1);
                PHI  = EulerArray(nEArow,2);
                phi2 = EulerArray(nEArow,3);
                vfm  = EulerArray(nEArow,4);

            case 4
                nb   = EulerArray(nEArow,1);
                phi1 = EulerArray(nEArow,2);
                PHI  = EulerArray(nEArow,3);
                phi2 = EulerArray(nEArow,4);
                vfm  = EulerArray(nEArow,5);
        end
        EArow = EArow - 1;
        
        ngrain = ngrain + 1;
        %******************************************************************
        % data with origin Xs at an azimuth not equal 0 degrees (north)
        %******************************************************************
        phi1 = phi1 - XsAz;
        if phi1 < 0
            phi1 = phi1 + 360;
        end
        
        % Bunge Euler matrix
        [g] = gmatrix(phi1,PHI,phi2);
        
        % Rij = Gki * EEkj = transpose(G) * EE
        for i = 1:3
            for j = 1:3
                r(i,j) = 0;
                for k = 1:3
                    r(i,j) = r(i,j) + g(k,i)*ee(k,j);
                end
            end
        end
        
        % Store values
        vstore(minNo,ngrain)     = vfm;
        rstore(minNo,ngrain,:,:) = r(:,:);
        
        % *****************************************************************
        % ROTATE ELASTIC CONSTANTS FROM CRYSTAL FRAME TO GEOGRAPHIC       *
        % *****************************************************************
        [CMV] = crot(r,EC);
        [SMR] = srot(r,ES);
        CSUM = CSUM + CMV*vfm;
        SSUM = SSUM + SMR*vfm;

    end

    fclose(fid10);

    % Store no. grains in mineral phase
    ngrains(minNo,1) = ngrain;

    %*****************************************************************
    %     CALCULATE AVERAGE ELASTIC CONSTANTS FOR MINERAL            *
    %*****************************************************************
    % normalise elastic stiffnesses by volume fractions
    if ioptn <= 2
        CMV = CSUM./ngrains(minNo,1);
        SMR = SSUM./ngrains(minNo,1);
    elseif ioptn >= 3
        CMV = CSUM;
        SMR = SSUM;
    end       
    
    % Voigt average - CMV    Reuss average - 1/SMR = CMR
    CMR = inv(SMR);
    
    % Voigt-Ruess-Hill average = (CMV + CMR)/2
    CVRH(:,:) = (CMV(:,:)+CMR(:,:))/2;

    % Store for aggregate total
    Cmin(minNo,:,:) = CMV(:,:);
    Smin(minNo,:,:) = SMR(:,:);

    %*****************************************************************
    %     PRINT ELASTIC CONSTANTS
    %*****************************************************************

    fprintf(fid12,['\n Phase %4i: %s',...
                   '\n  Single Crystal Elastic Stiffness Matrix (GPa)',...
                   '\n   %s\n   %s\n'],minNo,minName,...
                   cell2mat(ModelSC.Phase.(minName).reference));
    for i = 1:6
        fprintf(fid12,format320,EC(i,1:6));
    end

    fprintf(fid12,'\n  Single Crystal Elastic Compliance Matrix (1/GPa)\n');

    for i = 1:6
        fprintf(fid12,format320,ES(i,1:6));
    end

    fprintf(fid12,'\n  Density = %7.4f g/cm3 \n',ModelSC.Phase.(minName).rho);
    fprintf(fid12,'\n  Voigt average elastic stiffness (GPa)\n');    
    for i = 1:6
        fprintf(fid12,format320,CMV(i,1:6)); 
    end

    fprintf(fid12,'\n  Reuss average elastic stiffness (GPa)\n');    
    for i = 1:6
        fprintf(fid12,format320,CMR(i,1:6)); 
    end    

    fprintf(fid12,'\n  Voigt-Reuss-Hill average elastic stiffness (GPa)\n');    
    for i = 1:6
        fprintf(fid12,format320,CVRH(i,1:6)); 
    end    

    fprintf(fid12,['\n   No. of Grains = %6i \n',...
                   '   Grain shape A1:A2:A3 = ',repmat('%8.2f',1,3),'\n',...
                   '   Volume fraction = %8.4f \n\n',],...
                    ngrains(minNo),...
                    axis(minNo,1),axis(minNo,2),axis(minNo,3),...
                    vf(minNo));

end
%**************************************************************************
%                                                                         *
%     END OF MINERAL LOOP                                                 *
%                                                                         *
%**************************************************************************

%*****************************************************************
%
%     Volume Loop
% 
%*****************************************************************

% Rock density and volume averaged elastic constants
drock = 0;
CMV(:,:) = 0;
SMR(:,:) = 0;

% Calculate rock density and averages of aggregate
for minn = 1:nmin
    minName = cell2mat(mineral(minn));
    minNo = ModelSC.Phase.(minName).No;
    drock = drock + vf(minNo)*ModelSC.Phase.(minName).rho;
    % Voigt average: CMV
    CMV(:,:) = CMV(:,:) + vf(minNo)*squeeze(Cmin(minNo,:,:));
    SMR(:,:) = SMR(:,:) + vf(minNo)*squeeze(Smin(minNo,:,:));    
end

% Reuss average: 1/SMR = CMR
[CMR] = inv(SMR);

% Voigt-Ruess-Hill average: (CMV + CMR)/2
CVRH(:,:) = (CMV(:,:) + CMR(:,:))/2;


%*****************************************************************
%     PRINT ELASTIC CONSTANTS
%*****************************************************************

fprintf(fid12,'\n\n  Voigt aggregate elastic stiffness (GPa) \n');    
for i = 1:6
    fprintf(fid12,format320,CMV(i,1:6)); 
end

fprintf(fid12,'\n  VRH aggregate elastic stiffness (GPa) \n');    
for i = 1:6
    fprintf(fid12,format320,CVRH(i,1:6)); 
end    

fprintf(fid12,'\n  Reuss aggregate elastic stiffness (GPa) \n');  
for i = 1:6
    fprintf(fid12,format320,CMR(i,1:6)); 
end    

fprintf(fid12,'\n  Aggregate density = %7.4f g/cm3 \n\n',drock);

if isch == 1
    C(:,:) = CMV(:,:); % GPa
elseif isch == 2
    C(:,:) = CMR(:,:); % GPa
elseif isch == 3
    C(:,:) = CVRH(:,:); % GPa
end


%*****************************************************************
%
%     SC LOOP - Calculation in GPa
%
%*****************************************************************
fprintf(fid12,'  Xs of azimuth = %g degree \n\n',XsAz);
fprintf(fid12,'  Self-consistent scheme:');
if ioptsc == 1
    fprintf(fid12,'  Mechanical - Willis 1977 \n');
elseif ioptsc == 2
    fprintf(fid12,'  Scattering - Gubernatis & Krumhansl 1975 \n');
end
if isch == 1
    fprintf(fid12,'  Initial value = Voigt \n');
elseif isch == 2
    fprintf(fid12,'  Initial value = Reuss \n');
elseif isch == 3
    fprintf(fid12,'  Initial value = Voigt-Reuss-Hill \n');
elseif isch == 4
    fprintf(fid12,'  Initial value = Geometric mean \n');
end
fprintf(fid12,'  Max. %%error for Greens Tensor = %2.1f \n',errmax);
fprintf(fid12,'  Max. number SC iterations = %g \n',nsc);
fprintf(      '  Max. %%error for SC Convergence = %g \n\n',conmax);
fprintf(fid12,'  Max. %%error for SC Convergence = %g \n\n',conmax);
fprintf(fid12,'  Volume fractions \n');

for Min = 1:nmin
    fprintf(fid12,'  %g  VF = %4.2f    N grains = %g \n',...
                    Min,vf(Min),ngrains(Min));
end


% sum total no. grains for polyphase aggregate
ngrain = 0;
for i = 1:nmin
    ngrain = ngrain + ngrains(i);
end

% initial SC estimate
sc(1,:,:) = C(:,:);

%**********************************************************************
%     SC iterative loop
%**********************************************************************
iskip = 0;
iter  = 0;

disp('  SC iterative loop');

for  isc = 1:nsc
    % test convergence
    if isc > 1 && iskip == 0
        diag = 0;
        % upper diagonal of symmetric matrix
        for i = 1:6
            for j = i:6
                % sum diagonal terms only to avoid rounding errors
                % on off-diagonal terms in isotropic case
                if i == j
                    diag = diag + ...
                        abs(sc(isc-1,i,j)-sc(isc,i,j))/abs(sc(isc,i,j));
                end
            end
        end
        % average difference of diagonal terms as percentage
        diagav = 100*diag/6;
        if diagav < conmax
            iskip = 1;
        end
        fprintf('  iteration = %g   convergence %% = %g \n\n',iter,diagav);
    end
    % converged - skip the rest
    if iskip == 1
        continue;
    end
    % iteration counter
    iter = iter + 1;
    sum1(:,:)  = 0;
    sum2(:,:)  = 0;
    sumtm(:,:) = 0;
    % estimate for loop
    Csc(:,:) = sc(isc,:,:);
    
    %**********************************************************************
    %     Mineral loop
    %**********************************************************************
    vftotal = 0;
    for minn = 1:nmin
        % ellipsidal axes for mineral
        eaxis(1,:) = axis(minn,:);
        % elastic constants for mineral
        EC(:,:) = ECM(minn,:,:);
        % error analysis for Greens Tensor
        % mineral form in background medium Cscij
        [nb,x,w,ntcheb] = error4G(Csc,errmax,eaxis);
        fprintf('  Mineral No. %g  No. Greens tensor iterations = %g\n',...
            minn,nb);


        %******************************************************************
        %     Grain loop for mineral
        %******************************************************************
        for  n = 1:ngrains(minn)
            % grain volume in aggragate
            if ioptn <= 2
                vfg = vf(minn)/ngrains(minn);
            elseif ioptn >= 3
                vfg = vf(minn)*vstore(minn,n);
            end
            vftotal = vftotal + vfg;
            % rotation matrix R
            for i = 1:3
                for j = 1:3
                    % inclusion/grain
                    r(i,j) = rstore(minn,n,i,j);
                    % background - shape parallel to crystal elastic axes
                    if ioptoe(minn) == 1
                        rback(i,j) = rstore(minn,n,i,j);
                    end
                    % background - shape parallel specimen axes
                    if ioptoe(minn) == 2
                        rback(i,j) = remin(minn,i,j);
                    end
                end
            end
            % inverse rotation matrix for (background) matrix RBACKI
            rbacki(:,:) = rback';            
            % rotate inclusion/grain
            [Cg] = crot(r,EC);
            % rotate SC (background) matrix C*
            [CscR] = crot(rback,Csc);                         
            % solve Eshelby inclusion problem in (background) matrix for 
            %   Green's tensor Gij
            [Gtr] = green3(CscR,eaxis,x,w,ntcheb);                   
            % rotate Gij to (background) matrix axes
            [Gt] = crot(rbacki,Gtr);             
            
            if ioptsc == 1
                % interaction matrix Am(g) and C(g)*Am(g)
                [Am,CiAm] = inter(Csc,Cg,Gt);
                
                % summation
                sum1(:,:) = sum1(:,:) + vfg*CiAm(:,:);
                sum2(:,:) = sum2(:,:) + vfg*Am(:,:);                                
            else
                % scattering t matrix
                [tmatm] = scater(Csc,Cg,Gt);
                % summation
                sumtm(:,:) = sumtm(:,:) + vfg*tmatm(:,:);                                
            end
        end
    end
    if ioptsc == 1
        % invert sum2ij
        [sum2] = inv(sum2);        
        % CscMatrix(ij) = sum1ik*inv(sum2kj)
        for i = 1:6
            for j = 1:6
                CscM(i,j) = 0;
                for k = 1:6
                    CscM(i,j) = CscM(i,j) + sum1(i,k)*sum2(k,j);
                end
            end
        end   
        
        % convert from matrix to Voigt tensor
        [Csc] = Kelvin2Voigt(CscM);  
        
        % check if Voigt tensor is symmetric (Cij = Cji)
        Check_Voigt_Matrix_Symmetry_EK(Csc)
        
        % Make Voigt tensor to be symmetric (Cij = Cji)
        Csc = (Csc + Csc')/2;
              
    else
        % transform sum t-matrix to Voigt tensor
        [sumt] = Kelvin2Voigt(sumtm);
        
        % check if Voigt tensor is symmetric (Cij = Cji)
        Check_Voigt_Matrix_Symmetry_EK(sumt)
        
        % Make Voigt tensor to be symmetric (Cij = Cji)
        sumt = (sumt + sumt')/2;
                
        Csc = Csc + sumt;
    end

    % Store current estimate
    sc(isc+1,:,:) = Csc(:,:);
    
    % write to screen and file
    xn = isc;
    fprintf(fid17,['%g',repmat('\t %12.4f',1,21),'\n'],...
                    xn,Csc(1:6,1),Csc(2:6,2),Csc(3:6,3),...
                       Csc(4:6,4),Csc(5:6,5),Csc(6,6)); 
end

% CLOSE OUTPUT FILE
fclose(fid17);
% Save output Csc as a mat file
ModelSC.Csc = Csc;
save([dir_Output,deblank(Title),'-SC-Csc.mat'],'Csc');
save([dir_Output,deblank(Title),'_SC_ModelSC.mat'],'ModelSC');
%*****************************************************************
%
%     END OF SC LOOP
%
%*****************************************************************



%****************************************************************
%% Write Rock Elastic Constants to FILE 14 and FILE 12
%*****************************************************************
C = Csc; % GPa

% FILE 14
fprintf(fid14,['Self Consistent elastic stiffness (GPa) \n',...
               'Aggregate density = %7.4f g/cm3 \n'],drock);
fprintf(fid12,'\n\n  Self Consistent elastic stiffness (GPa) \n');

for i = 1:6
    fprintf(fid14,format320,C(i,1:6));   
    fprintf(fid12,format320,C(i,1:6));   
end

% CLOSE OUTPUT FILEs
fclose(fid14);
fclose(fid12);

% clf;

%**************************************************************************
%% Demo of Multi_seismic_plot_MTEX53_EK.m
%**************************************************************************
% Seismic velocities and anisotropy - single crystal plots
% Tested with MTEX 5.0.3 and MATLAB 2016a  10/05/2018 David Mainprice
% Tested with MTEX 5.2.beta3 and MATLAB R2016b  28/03/2019 Eunyoung Kim
%
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
% e.g. Multi_seismic_plot_MTEX5(C,Label_X,Label_Y,Label_Z,hemisphere,symm)
%
% *** Remember that the function "Multi_seismic_plot_MTEX53_EK2"
% should be in your current directory ***
%
% David Mainprice 10/05/2018
% Eunyoung Kim 28/03/2019

%**************************************************************************
% Define elastic constants
%**************************************************************************

% SC_Tensor of aggregate
%
% Crystal symmetry: triclinic -1
%
% Define density (g/cm3)
rho_SC = drock;
% ref1_SC: SC aggregate, density = ( ) g/cm3
% ref2_SC: Crystal symmetry is triclinic
ref1_SC = {['SC aggregate, density = ',num2str(drock),' g/cm3']};
ref2_SC = {'Crystal symmetry is triclinic'};

% Define Cartesian Tensor crystal symmetry and frame
cs_SC = crystalSymmetry('-1','mineral','SC aggregate');


% *Import 4th rank tensor as 6 by 6 matrix*
% elastic Cij stiffness tensor (GPa) as matrix M
M_SC = C;

% Define tensor object in MTEX
% Cij -> Cijkl - elastic stiffness tensor
% M as stiffness tensor C with MTEX tensor command
SC_Tensor = ...
    stiffnessTensor(M_SC,cs_SC,'density',rho_SC);
% % % % % [vp,vs1,vs2,pp,ps1,ps2] = velocity(C,xvector)


% Run tensor symmetry check
if(SC_Tensor.checkSymmetry == false)
    warning('MTEX:tensor','Tensor does not have the correct symmetry')
else
    disp('MTEX:Tensor symmetry is correct')
end
%**************************************************************************
% Plot Polycrystalline aggregate in specimen coordinates X,Y,Z
%**************************************************************************

% plot x to east
plotx2east
% plot a to east
plota2east

% You can change these labels to fit your application
% Label_X='a[100]'; Label_Y='b[010]';Label_Z='c[001]';
Label_X='X'; Label_Y='Y';Label_Z='Z';
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
font_size=38;
% marker_size (max min)
marker_size=30;
%
setMTEXpref('FontSize',font_size)
%
% Multi_seismic_plot_MTEX53.m is a modified version by David Mainprice and Ralf Hielscher
Multi_seismic_plot_MTEX53_EK(SC_Tensor,Label_X,Label_Y,Label_Z,hemisphere,symm,font_size)
%
drawNow(gcm,'figSize','large')
set(gcf,'renderer','zBuffer')

%*****************************************************************
%% Functions
%******************************************************************

function [ct,tc] = bond(a,e,iopt)

se = zeros(1,3);
ce = zeros(1,3);
tc = zeros(3,3);
ct = zeros(3,3);

rad = pi/180;

% Bond transformation matrices TC and CT
for i = 1:3
    se(i) = sin(e(i)*rad);
    ce(i) = cos(e(i)*rad);
end

u = (ce(3)-ce(1)*ce(2))/se(2);
v = sqrt(se(1)^2-u^2);
w = ce(2)/se(2);

ct(1,1) = a(1)*se(2);
ct(1,2) = a(2)*u;
ct(2,2) = a(2)*v;
ct(3,1) = a(1)*ce(2);
ct(3,2) = a(2)*ce(1);
ct(3,3) = a(3);
tc(1,1) = 1/a(1)/se(2);
tc(1,2) = -u/a(1)/v/se(2);
tc(2,2) = 1/a(2)/v;
tc(3,1) = -w/a(3);
tc(3,2) = (u*w-ce(1))/(a(3)*v);
tc(3,3) = 1/a(3);

% print if IOPT=1
format40   = ['\n\n',repmat(' ',1,8),'matrix cmt'];
format60   = ['\n\n',repmat(' ',1,8),'matrix tmc'];
format1000 = [repmat('%12.4f',1,3),' \n'];

if iopt == 1 
    fprintf(format40);
    for i = 1:3
        fprintf(format1000,ct(i,1:3));
    end
    fprintf(format60);    
    for i = 1:3
        fprintf(format1000,tc(i,1:3)); 
    end
end

% ct = ct';
% tc = tc';

end



function [isym,phi1max,phimax,phi2max] = eulsym(phi1max,phimax,phi2max,isym)

disp(' Euler space for crystal- *triclinic sample symmetry*');
disp(' 11 Proper point groups (Schoenflies - International)');
disp(' Crystal symmetry                  phi1 : phi : phi2');
disp(' 0 = Triclinic (C1 - 1)             360 : 180 : 360');
disp(' 1 = Monoclinic (C2 - 2)            360 : 180 : 180');
disp(' 2 = Orthorhombic (D2 - 222)        180 : 180 : 180');
disp(' 3 = Trigonal(D3 - 32) (2-fold//Y)  180 : 180 : 120');
disp(' 4 = Trigonal(C3 - 3)               360 : 180 : 120');
disp(' 5 = Tetragonal (D4 - 422)          180 : 180 :  90');
disp(' 6 = Tetragonal (C4 - 4)            360 : 180 :  90');
disp(' 7 = Hexagonal  (D6 - 622)          180 : 180 :  60');
disp(' 8 = Hexagonal  (C6 - 6)            360 : 180 :  60');
disp(' 9 = Cubic (O - 432)                180 : 180 :  90');
disp('10 = Cubic (T - 23)                 360 : 180 :  90');
disp(' ');
disp(' N.B. 0 = No crystal symmetry operations');
disp(' ');


switch isym
% Triclinic (C1 - 1) 360:180:360
    case 0
        phi1max = 360;
        phimax  = 180;
        phi2max = 360;
% Monoclinic (C2 - 2) 360:180:180
    case 1
        phi1max = 360;
        phimax  = 180;
        phi2max = 180;
% Orthorhombic (D2 - 222)  180:180:180
    case 2
        phi1max = 180;
        phimax  = 180;
        phi2max = 180;
% Trigonal(D3 - 32) (*2-fold//Y*) 180:180:120
    case 3
        phi1max = 180;
        phimax  = 180;
        phi2max = 120;
% Trigonal (C3 - 3) 360:180:120
    case 4
        phi1max = 360;
        phimax  = 180;
        phi2max = 120;
% Tetragonal (D4 - 422) 180:180:90
    case 5
        phi1max = 180;
        phimax  = 180;
        phi2max = 90;
% Tetragonal (C4 - 4) 360:180:90
    case 6
        phi1max = 360;
        phimax  = 180;
        phi2max = 90;
% Hexagonal (D6 - 622) 180:180:60
    case 7
        phi1max = 180;
        phimax  = 180;
        phi2max = 60;
% Hexagonal (C6 - 6) 360:180:60
    case 8
        phi1max = 360;
        phimax  = 180;
        phi2max = 60;
% Cubic (O - 432) 180:180:90
    case 9
        phi1max = 180;
        phimax  = 180;
        phi2max = 90;
% Cubic (T - 23) 360:180:90
    case 10
        phi1max = 360;
        phimax  = 180;
        phi2max = 90;
end 

end 




function [xr,ixopt,zr,izopt] = refustage(isym,xr,zr)
% Crystal Reference Directions(Xc,Yc,Zc) for Euler angle determination
% *** Universal stage data ***
% Cubic/Othorhombic/Tetragonal Xc=100 Yc=010 Zc=001
% Hexagonal/Trigonal Xc=m(0-1.0) Yc=a Zc=c(00.1) (*euler ODF*)
% Monoclinic Xc=a* Yc=010 Zc=001 (*clinopyroxenes*)
% Triclinic Xc=100 Yc=010 Zc=c*  (*plagioclase*)
%
% ISYM codes                       Order of point group
%  0 = Triclinic (C1 - 1)              1
%  1 = Monoclinic (C2 - 2)             2
%  2 = Orthorhombic (D2 - 222)         4
%  3 = Trigonal (D3 - 32)              6
%  4 = Trigonal (C3 - 3)               3
%  5 = Tetragonal (D4 - 422)           8
%  6 = Tetragonal (C4 - 4)             4
%  7 = Hexagonal (D6 - 622)           12
%  8 = Hexagonal (C6 - 6)              6
%  9 = Cubic (O - 432)                24
% 10 = Cubic (T - 23)                 12

switch isym
% 0 Triclinic (C1 - 1)
    case 0
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 2;
% 1 Monoclinic (C2 - 2)   (NB case of cpx X1=(100) and X3=[001])
    case 1
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 2;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 1;
% 2 Orthorhombic (D2 - 222)
    case 2
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 1;
% 3 Trigonal (D3 - 32)
    case 3
        xr(1) = 0;
        xr(2) = -1;
        xr(3) = 0;
        ixopt = 2;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 2;
% 4 Trigonal (C3 - 3)
    case 4
        xr(1) = 0;
        xr(2) = -1;
        xr(3) = 0;
        ixopt = 2;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 2;
% 5 Tetragonal (D4 - 422)
    case 5
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 1;
% 6 Tetragonal (C4 - 4)
    case 6
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 1;
% 7 Hexagonal (D6 - 622)
    case 7
        xr(1) = 00;
        xr(2) = -1;
        xr(3) = 0;
        ixopt = 2;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 2;
% 8 Hexagonal (C6 - 6)
    case 8
        xr(1) = 0;
        xr(2) = -1;
        xr(3) = 0;
        ixopt = 2;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 2;
% 9 Cubic (O - 432)
    case 9
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 00;
        zr(3) = 1;
        izopt = 1;
% 10 Cubic (T - 23)
    case 10
        xr(1) = 1;
        xr(2) = 0;
        xr(3) = 0;
        ixopt = 1;
        zr(1) = 0;
        zr(2) = 0;
        zr(3) = 1;
        izopt = 1;
end
end 



function [xr,ixopt,zr,izopt] = refebsd(xr,zr)
% Crystal Reference Directions(Xc,Yc,Zc) for EBSD
% for Euler angle in *** Channel+ ***
% Cubic/Othorhombic/Tetragonal Xc=100 Yc=010 Zc=001
% Hexagonal/Trigonal Xc=m*(10.0) Yc=a Zc=c(00.1) (*euler ODF*)
% Monoclinic Xc=a* Yc=[010] Zc=c[001]
% Triclinic Xc=a* Yc=c x a* Zc=c[001]

% Crystal Reference Directions
% for Euler angle determination in Channel+
%
% (called X2,Y2,Z2 in Channel+ manual)
% X2=a*(100) Z2=c[001] Y2= Z2 x X2

% X2=a*(100)
% IXOPT =2 pole to plane (hkl)

xr(1) = 1;
xr(2) = 0;
xr(3) = 0;
ixopt = 2;

% Z2=c[001]
% IZOPT =1 direction (uvw)
zr(1) = 0;
zr(2) = 0;
zr(3) = 1;
izopt = 1;

end 


function [vc] = ortho(cmt,tmc,iopt,v)

vc = zeros(1,3);

for i = 1:3
    vc(i) = 0;
    for j = 1:3
        if iopt == 1
            xm = cmt(i,j);  
        elseif iopt == 2
            xm = tmc(j,i); 
        end
        vc(i) = vc(i) + xm*v(j);
    end
end

end


function [c] = xprod(a,b)
% CROSS VECTOR PRODUCT C = A X B

c = zeros(1,3);
d = zeros(1,3);

d(1,1) = a(2)*b(3) - a(3)*b(2);
d(1,2) = a(3)*b(1) - a(1)*b(3);
d(1,3) = a(1)*b(2) - a(2)*b(1);

c(:,:) = d(:,:);
end


function [value] = angleAB(a,b)
% ANGLE BETWEEN DIRECTION COSINES A,B

dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
if dot > 1
    dot = 1; 
elseif dot < -1
    dot = -1; 
end

value = acos(dot)*(180/pi); % 180/pi = 57.2958

end



function [g] = gmatrix(psi1,phi,psi2)
% Bunge Euler matrix g(psi1,phi,psi2)

rad = pi/180;
g = zeros(3,3);

spsi1 = sin(rad*psi1);
sphi  = sin(rad*phi);
spsi2 = sin(rad*psi2);
cpsi1 = cos(rad*psi1);
cphi  = cos(rad*phi);
cpsi2 = cos(rad*psi2);

g(1,1) = cpsi1*cpsi2 - spsi1*spsi2*cphi;
g(1,2) = spsi1*cpsi2 + cpsi1*spsi2*cphi;
g(1,3) = spsi2*sphi;
g(2,1) = -cpsi1*spsi2 - spsi1*cpsi2*cphi;
g(2,2) = -spsi1*spsi2 + cpsi1*cpsi2*cphi;
g(2,3) = cpsi2*sphi;
g(3,1) = spsi1*sphi;
g(3,2) = -cpsi1*sphi;
g(3,3) = cphi;

end



function [C] = crot(r,ec)
% ROTATE ELASTIC CONSTANTS FROM CRYSTAL TO SPACIAL COORDINATES
% Cijkl = Rip*Rjq*Rkr*Rls*ECpqrs
% Cij   = CMik*CMjl*ECkl using method of W.L.Bond

% Preallocation
C  = zeros(6,6);
cm = zeros(6,6);

% Matrix CM
cm(1,1) = r(1,1)^2;
cm(1,2) = r(1,2)^2;
cm(1,3) = r(1,3)^2;
cm(1,4) = 2*r(1,2)*r(1,3);
cm(1,5) = 2*r(1,3)*r(1,1);
cm(1,6) = 2*r(1,1)*r(1,2);
cm(2,1) = r(2,1)^2;
cm(2,2) = r(2,2)^2;
cm(2,3) = r(2,3)^2;
cm(2,4) = 2*r(2,2)*r(2,3);
cm(2,5) = 2*r(2,3)*r(2,1);
cm(2,6) = 2*r(2,1)*r(2,2);
cm(3,1) = r(3,1)^2;
cm(3,2) = r(3,2)^2;
cm(3,3) = r(3,3)^2;
cm(3,4) = 2*r(3,2)*r(3,3);
cm(3,5) = 2*r(3,3)*r(3,1);
cm(3,6) = 2*r(3,1)*r(3,2);
cm(4,1) = r(2,1)*r(3,1);
cm(4,2) = r(2,2)*r(3,2);
cm(4,3) = r(2,3)*r(3,3);
cm(4,4) = r(2,2)*r(3,3)+r(2,3)*r(3,2);
cm(4,5) = r(2,1)*r(3,3)+r(2,3)*r(3,1);
cm(4,6) = r(2,2)*r(3,1)+r(2,1)*r(3,2);
cm(5,1) = r(3,1)*r(1,1);
cm(5,2) = r(3,2)*r(1,2);
cm(5,3) = r(3,3)*r(1,3);
cm(5,4) = r(1,2)*r(3,3)+r(1,3)*r(3,2);
cm(5,5) = r(1,1)*r(3,3)+r(1,3)*r(3,1);
cm(5,6) = r(1,1)*r(3,2)+r(1,2)*r(3,1);
cm(6,1) = r(1,1)*r(2,1);
cm(6,2) = r(1,2)*r(2,2);
cm(6,3) = r(1,3)*r(2,3);
cm(6,4) = r(2,2)*r(1,3)+r(1,2)*r(2,3);
cm(6,5) = r(1,1)*r(2,3)+r(1,3)*r(2,1);
cm(6,6) = r(2,2)*r(1,1)+r(1,2)*r(2,1);

% Cij = CMik*CMjl*ECkl

for i = 1:6
    for j = 1:6
        C(i,j) = 0;
        for k = 1:6
            for l = 1:6
                C(i,j) = C(i,j) + cm(i,k)*cm(j,l)*ec(k,l);
            end
        end 
    end 
end
end


function [S] = srot(r,es)
% ROTATE ELASTIC CONSTANTS FORM CRYSTAL TO SPACIAL COORDINATES
% Sijkl=Rip*Rjq*Rkr*Rls*ESpqrs
% Sij = SNik*SNjl*ESkl using method of W.L.Bond

% Preallocation
S  = zeros(6,6);
sn = zeros(6,6);

% Matrix SN
sn(1,1) = r(1,1)^2;
sn(1,2) = r(1,2)^2;
sn(1,3) = r(1,3)^2;
sn(1,4) = r(1,2)*r(1,3);
sn(1,5) = r(1,3)*r(1,1);
sn(1,6) = r(1,1)*r(1,2);
sn(2,1) = r(2,1)^2;
sn(2,2) = r(2,2)^2;
sn(2,3) = r(2,3)^2;
sn(2,4) = r(2,2)*r(2,3);
sn(2,5) = r(2,3)*r(2,1);
sn(2,6) = r(2,1)*r(2,2);
sn(3,1) = r(3,1)^2;
sn(3,2) = r(3,2)^2;
sn(3,3) = r(3,3)^2;
sn(3,4) = r(3,2)*r(3,3);
sn(3,5) = r(3,3)*r(3,1);
sn(3,6) = r(3,1)*r(3,2);
sn(4,1) = 2*r(2,1)*r(3,1);
sn(4,2) = 2*r(2,2)*r(3,2);
sn(4,3) = 2*r(2,3)*r(3,3);
sn(4,4) = r(2,2)*r(3,3) + r(2,3)*r(3,2);
sn(4,5) = r(2,1)*r(3,3) + r(2,3)*r(3,1);
sn(4,6) = r(2,2)*r(3,1) + r(2,1)*r(3,2);
sn(5,1) = 2*r(3,1)*r(1,1);
sn(5,2) = 2*r(3,2)*r(1,2);
sn(5,3) = 2*r(3,3)*r(1,3);
sn(5,4) = r(1,2)*r(3,3) + r(1,3)*r(3,2);
sn(5,5) = r(1,1)*r(3,3) + r(1,3)*r(3,1);
sn(5,6) = r(1,1)*r(3,2) + r(1,2)*r(3,1);
sn(6,1) = 2*r(1,1)*r(2,1);
sn(6,2) = 2*r(1,2)*r(2,2);
sn(6,3) = 2*r(1,3)*r(2,3);
sn(6,4) = r(2,2)*r(1,3) + r(1,2)*r(2,3);
sn(6,5) = r(1,1)*r(2,3) + r(1,3)*r(2,1);
sn(6,6) = r(2,2)*r(1,1) + r(1,2)*r(2,1);

% Sij = CMik*CMjl*ESkl

for i=1:6
    for j=1:6
        S(i,j) = 0;
        for k=1:6
            for l=1:6
                S(i,j)=S(i,j)+sn(i,k)*sn(j,l)*es(k,l);
            end
        end
    end
end
end 



function [tmatm] = scater(Cmat,Cinc,Gt)
% Scattering t-matrix based on method of
% J.E.Gubernatis and J.A.Krumhansl(1975) J.Appl.Phys.vol.46,pp.1875-1882.
% Cmat(6,6)  = Matrix stiffness tensor
% Cinc(6,6)  = Inclusion stiffness tensor
% Gt(6,6)    = Green's tensor
% tmatm(6,6) = t-matrix

% Preallocations
iden  = zeros(6,6); 
Am    = zeros(6,6); 
diffm = zeros(6,6); 
tmatm = zeros(6,6); 
 
% Iij - cartesian idenity matrix
for i = 1:6
    for j = 1:6
        iden(i,j) = floor(i/j)*floor(j/i);
    end
end

% transform all tensors to cartesian matrices
[cmatm] = Voigt2Kelvin(Cmat);
[cincm] = Voigt2Kelvin(Cinc);
[Gm]    = Voigt2Kelvin(Gt);

% difference matrix Ci-C*
diffm(:,:) = cincm(:,:) - cmatm(:,:);



% I+G(Ci-Cmat)  (** +ve G in our case **)
Am(:,:) = iden(:,:) + Gm(:,:)*diffm(:,:);
% Am = inv (I+G(Ci-Cmat))
% tmatm = (Ci-Cmat)*Am = (Ci-Cmat)*inv(I+G(Ci-Cmat))
tmatm(:,:) = Am(:,:)\diffm(:,:);


end 


function [Voigt] = Kelvin2Voigt(Kelvin)
% The Voigt matrix does not preserves the norm of the tensor.
% Converts Voigt to Kelvin matrix and Kelvin to Voigt matrix
% The Kelvin elastic matrix preserves the norm of the tensor
% and its eigen-values and vectors have physical meaning.
% 
% All mathematical matrix functions can be applied to Kelvin matrix 
%
% References
% W.K. Thomson (Lord Kelvin),1856 
%   Elements of a mathematical theory of elasticity, 
%   Philos. Trans. R. Soc. 146, 481-498.
% W.K. Thomson (Lord Kelvin),1856 
%   On six principal strains of an elastic solid, 
%   Phil. Trans. R. Soc. 166, 495-498.
% W.K. Thomson (Lord Kelvin),1878 
%   Mathematical theory of elasticity, 
%   Encycl. Br. 7, 819-825. 
%
% INPUT
% Kelvin matrix
% 
% OUTPUT
% Voigt matrix
%
% David Mainprice (18/07/2018)
%**************************************************************************
% Matrix conversion Voigt to Kelvin = inv(Kelvin2Voigt)
Kelvin2VoigtMatrix = [[  1      0      0     0      0       0];...
                     [   0      1      0     0      0       0];...
                     [   0      0      1     0      0       0];...
                     [   0      0      0   1/sqrt(2)  0     0];...
                     [   0      0      0     0    1/sqrt(2) 0];...
                     [   0      0      0     0      0  1/sqrt(2)]];
% Matrix Kelvin to Voigt method  
Voigt = Kelvin2VoigtMatrix*Kelvin*Kelvin2VoigtMatrix;
end


function [Kelvin] = Voigt2Kelvin(Voigt)
% The Voigt matrix does not preserves the norm of the tensor.
% Converts Voigt to Kelvin matrix and Kelvin to Voigt matrix
% The Kelvin elastic matrix preserves the norm of the tensor
% and its eigen-values and vectors have physical meaning.
% 
% All mathematical matrix functions can be applied to Kelvin matrix 
%
% References
% W.K. Thomson (Lord Kelvin),1856 
%   Elements of a mathematical theory of elasticity,
%   Philos. Trans. R. Soc. 146, 481-498.
% W.K. Thomson (Lord Kelvin),1856 
%   On six principal strains of an elastic solid, 
%   Phil. Trans. R. Soc., 166, 495-498.
% W.K. Thomson (Lord Kelvin),1878 
%   Mathematical theory of elasticity,
%   Encycl. Br. 7, 819-825. 
%
% INPUT
% Voigt matrix
% 
% OUTPUT
% Kelvin matrix
%
% David Mainprice (18/07/2018)
%**************************************************************************
% Matrix conversion Voigt to Kelvin = inv(Kelvin2Voigt)
Voigt2KelvinMatrix = [[  1      0      0     0      0       0];  ...
                     [   0      1      0     0      0       0];  ...
                     [   0      0      1     0      0       0];  ...
                     [   0      0      0   sqrt(2)  0       0];  ...
                     [   0      0      0     0    sqrt(2)   0];  ...
                     [   0      0      0     0      0  sqrt(2)]];
% Matrix Voigt to Kelvin method  
Kelvin = Voigt2KelvinMatrix*Voigt*Voigt2KelvinMatrix;
end





function [nb,x,w,ntcheb] = error4G(Cm,Errmax,eaxis)
% Error analysis of Green's Tensor for medium Cij
%
% Cm(6,6)  = Background medium
% Errmax   = Max. Error % (0.1 is a typical value) on diagonal terms
% eaxis(3) = Ellipsoid semi-axes
% nb       = number of iterations necessary to obtain
%            error < max. error
% ** revised June 1996 **

GtOld = zeros(6,6); 

% loop over possible number of integration steps
for inb = 1:100
    nb = inb*4;
    
    [ntcheb,x,w] = tcheb(nb);
    
    % Calculate Green's tensor
    [Gt] = green3(Cm,eaxis,x,w,ntcheb);
    
    diag = 0; 
    
    % upper diagonal of symmetric matrix
    for i = 1:6
        for j = i:6            
            % sum diagonal terms only to avoid rounding errors
            % on off-diagonal terms in isotropic case
            if i == j
                diag = diag + abs(Gt(i,j) - GtOld(i,j))/abs(Gt(i,j)); 
            end
            GtOld(i,j) = Gt(i,j);
        end
    end
    
    % average difference of diagonal terms as percentage
    diagav = 100*diag/6;
    
    % quit loop if convergence
    if diagav < Errmax
        [ntcheb,x,w] = tcheb(nb); 
        return;
    end
end

end


function [ntcheb,x,w] = tcheb(nb)

x  = zeros(1,nb/2);
w  = zeros(1,nb/2);

ntcheb = nb/2;
x1 = 0;
x2 = 1;

[xp,wp] = gauleg(x1,x2,nb);

for i = 1:nb/2
    x(i) = 0.5 - xp(i);
    w(i) = wp(i);
end

end



function [x,w] = gauleg(x1,x2,n)
% Given the lower and upper limits of integration X1 and X2, and given N, 
%   this routine returns 
% arrays X and W of length N, containing the abscissas and weights of the 
%   Gauss-Legendre 
% N-point quadrature formula

eps = 3*10^(-14);

x = zeros(1,n);
w = zeros(1,n);

% m = (n+1)/2;
m  = n/2;
xm = 0.5*(x1+x2);
xl = 0.5*(x2-x1);
xn = n;

for i = 1:m
    xi = i;
    z = cos(pi*(xi-0.25)/(xn+0.5));
%     fprintf('z = %g ',z);
    p1 = 1;
    p2 = 0;
    for j = 1:n
        xj = j;
        p3 = p2;
        p2 = p1;
        p1 = ((2*j-1)*z*p2-(xj-1)*p3)/xj;
    end
    pp = n*(z*p1-p2)/(z*z-1);
    z1 = z;
    z = z1-p1/pp;

    while abs(z-z1) > eps 
        p1 = 1;
        p2 = 0;

        for j = 1:n
            xj = j;
            p3 = p2;
            p2 = p1;
            p1 = ((2*j-1)*z*p2-(xj-1)*p3)/xj;
        end
        pp = n*(z*p1-p2)/(z*z-1);
        z1 = z;
        z = z1-p1/pp;  
    end      

    x(i) = xm - xl*z;
    x(n+1-i) = xm + xl*z;
    w(i) = 2*xl/((1-z*z)*pp*pp);
    w(n+1-i) = w(i);
  
end

end


  
function [G] = green3(C2,axis,x,w,ntcheb)
% Green's tensor for Eshelby inclusion problem
%   using extended analysis of Kinoshita and Mura (1971)
%   phys. stat. sol. (a) Vol.5 pp.759-768
%   for anisotropic inclusion and in anisotropic medium
%   following Appendix of Lebensohn & Tome (1993)
%   Acta metall. mater. Vol.41 pp.2611-2624.
%
% N.B. Green's Tensor for inclusion with axes parallel to
%   the reference frame of the elastic medium (X1,X2,X3).
% Based on routine by Ricardo Lebensohn & Gilles Canova
%   modified for single site and Voigt tensors
%
% C2(6,6) = Matrix (background) stiffness Voigt tensor
% G(6,6)  = Green's Voigt tensor
% axes(3) = Ellipsoidal inclusion semi-axes

ijkl = reshape([1,6,5,6,2,4,5,4,3],[3,3]);
ijv  = reshape([1,2,3,2,3,1,1,2,3,3,1,2],[6,2]);

G  = zeros(6,6);
g2 = zeros(6,6);
ck = zeros(3,3);
cp = zeros(3,3);
xk = zeros(1,3);

% POINTS & WEIGHTS FOR GAUSS INTEGRATION
nint = 2*ntcheb;

xp = zeros(1,nint);
wp = zeros(1,nint);

for i = 1:ntcheb
    xp(i)        = 0.5 + x(i);
    wp(i)        = w(i);
    xp(i+ntcheb) = 0.5 - x(i);
    wp(i+ntcheb) = w(i);
end

% INTEGRATION INTERVALS: [0,PI/2][PI/2,PI] BOTH FOR THETA AND PHI
pas = pi/2;

stheta = zeros(nint,2);
ctheta = zeros(nint,2);
sphi   = zeros(nint,2);
cphi   = zeros(nint,2);

for in = 1:2
    xx1 = pas*(in-1);
    for ig = 1:nint
        xx = xx1 + pas*xp(ig);        
        
        stheta(ig,in) = sin(xx);
        ctheta(ig,in) = cos(xx);
        sphi(ig,in)   = sin(xx);
        cphi(ig,in)   = cos(xx);
        
    end
end

% INITIALIZE SYMMETRIC Greens Tensor
G(:,:) = 0;

% BIG LOOP (OVER INTERVALS AND INTEGRATION POINTS)
for it = 1:2
    for ip = 1:2
        for jt = 1:nint
            for jp = 1:nint
                
                sth = stheta(jt,it);
                cth = ctheta(jt,it);
                sph = sphi(jp,ip);
                cph = cphi(jp,ip);
                
                % ALFA: VERSOR in K DIRECTION
                xk(1) = sth*cph;
                xk(2) = sth*sph;
                xk(3) = cth;

                % (k**2 tch_1 Tfourier(G))=f(alfa)
                % Christoffel matrix
                for i = 1:3
                    for k = 1:i
                        ckik = 0;
                        for j = 1:3
                            for l = 1:3                                
                                mm = ijkl(j,i);
                                nn = ijkl(l,k);
                                ckik = ckik+C2(nn,mm)*xk(j)*xk(l);                                
                            end
                        end
                        ck(k,i) = ckik;
                        ck(i,k) = ckik;
                    end
                end
                
                cp(:,:) = ck(:,:);

                [ck] = inv(cp);
                
                % Symmetric alfa*alfa*ck
                % /2 instead of /4, because of integration
                % in phi=[0,pi] instead of [0,2pi]
                for i = 1:6
                    i1 = ijv(i,1);
                    i2 = ijv(i,2);
                    for j = 1:6
                        j1 = ijv(j,1);
                        j2 = ijv(j,2);
                        g2(j,i) = ck(i1,j1)*xk(i2)*xk(j2) + ...
                                  ck(i2,j1)*xk(i1)*xk(j2) + ...
                                  ck(i1,j2)*xk(i2)*xk(j1) + ...
                                  ck(i2,j2)*xk(i1)*xk(j1);
                        g2(j,i) = 0.5*g2(j,i);
                    end
                end

                % F(th,ph)/ro**6 CALCULATION
                ro = sqrt((axis(1)*xk(1))^2 + ...
                          (axis(2)*xk(2))^2 + ...
                          (axis(3)*xk(3))^2);
                alfa = 2*ro;
                sgalfa = abs(alfa)/alfa;
                t = (-alfa^3 + 6*ro*alfa^2 - 6*alfa*ro^2)*sgalfa;
                f1 = 2*t/ro^6;

                % sum
                for i = 1:6
                    for j = 1:6
                        G(j,i) = G(j,i) + g2(j,i)*sth*f1*wp(jt)*wp(jp);
                    end
                end
                
                % END OF BIG LOOP
            end
        end
    end
end

fact = axis(1)*axis(2)*axis(3)/32/pi;

G(:,:) = fact*pas*pas*G(:,:);

end



function [Am,CiAm] = inter(Cmat,Cinc,Gt)
% Inclusion interaction problem based on method of
% J.R.Willis (1977) J.Mech.Phys.Solids vol25,pp.185-202.
%
% Cmat(6,6) = Matrix stiffness tensor
% Cinc(6,6) = Inclusion stiffness tensor
% Gt(6,6)   = Green's tensor
% Am(6,6)   = Interaction matrix
% CiAm(6,6) = product Ci * Am

Iden = zeros(6,6);
Am   = zeros(6,6);
CiAm = zeros(6,6);

% Iij - cartesian idenity matrix
for i = 1:6
    for j = 1:6
        Iden(i,j) = floor(i/j)*floor(j/i);
    end
end

% transform all (Voigt) tensors to cartesian (Kelvin) matrices
[Cmatm] = Voigt2Kelvin(Cmat);
[Cincm] = Voigt2Kelvin(Cinc);
[Gm]    = Voigt2Kelvin(Gt);

% I+G(Ci-Cmat)
for i = 1:6
    for j = 1:6
        Am(i,j) = Iden(i,j);
        for k = 1:6
            Am(i,j) = Am(i,j) + Gm(i,k)*(Cincm(k,j) - Cmatm(k,j));
        end
    end
end
% % % Am = iden + Gm*(Cincm - Cmatm);

% Am = inv (I+G(Ci-Cmat))
% [Am] = matinv(Am,6);
Am = inv(Am);

% CiAm = Ci*Am
for i = 1:6
    for j = 1:6
        CiAm(i,j) = 0;
        for k = 1:6
            CiAm(i,j) = CiAm(i,j) + Cincm(i,k)*Am(k,j);
        end
    end
end

end





function [r] = azr(al,xl,af,xf)
% AZ/INC Lineation & Foliation -> Rij 
%
% ellipsoid semi-axis A1
%    al: Azimuth 
%    xl: Inclination 
% ellipsoid semi-axis A3
%    af: Azimuth 
%    xf: Inclination 

r   = zeros(3,3); 
deg = 180/pi;

[x1,y1,z1] = xyz(al,xl);
[x3,y3,z3] = xyz(af,xf);

[x1,y1,z1] = dircos(x1,y1,z1);
[x3,y3,z3] = dircos(x3,y3,z3);

dot = x1*x3 + y1*y3 + z1*z3;

if dot > 1
    dot = 1; 
elseif dot < -1
    dot = -1; 
end

value = acos(dot)*deg;

fprintf(' Check 1 :- Angle A1 to A3 = %0.15g \n',value);

% repeat to round-off errors
[x2,y2,z2] = xcross(x3,y3,z3,x1,y1,z1);
[x2,y2,z2] = dircos(x2,y2,z2);
[x3,y3,z3] = xcross(x1,y1,z1,x2,y2,z2);
[x1,y1,z1] = xcross(x2,y2,z2,x3,y3,z3);
[x2,y2,z2] = xcross(x3,y3,z3,x1,y1,z1);

dot = x1*x3 + y1*y3 + z1*z3;

if dot > 1
    dot = 1; 
elseif dot < -1
    dot = -1; 
end

value = acos(dot)*deg;

fprintf(' Check 2 :- Angle A1 to A3 = %0.15g \n',value);

% stored by rows
r(1,1) = x1;
r(1,2) = y1;
r(1,3) = z1;
r(2,1) = x2;
r(2,2) = y2;
r(2,3) = z2;
r(3,1) = x3;
r(3,2) = y3;
r(3,3) = z3;

rdet = det(r);

fprintf([' Check 3 :- Det|R| should be +1 \n',...
         ' Det|R| = %0.15g \n'],rdet);

end

function [x,y,z] = xyz(az,xinc)
% az,xinc --> x,y,z (direction cosines)
%
% Right handed system 
%       x = North y = West z = Vertical

rad  = pi/180;
raz  = az*rad;
rinc = xinc*rad;

x = cos(raz)*cos(rinc);
y = -sin(raz)*cos(rinc);
z = sin(rinc);

end


function [x,y,z] = dircos(x,y,z)
% Convert to direction cosines

small = 0.0001;
xm = sqrt(x^2 + y^2 + z^2);

if abs(x) > small
    x = x/xm; 
elseif abs(x) <= small
    x = 0; 
end

if abs(y) > small
    y = y/xm; 
elseif abs(y) <= small
    y = 0; 
end

if abs(z) > small
    z = z/xm; 
elseif abs(z) <= small
    z = 0; 
end

end



function [x3,y3,z3] = xcross(x1,y1,z1,x2,y2,z2)
% Vector cross-product 1 x 2 = 3 right-handed

x3 = (y1*z2) - (z1*y2);
y3 = (z1*x2) - (x1*z2);
z3 = (x1*y2) - (y1*x2);

end


function Check_Voigt_Matrix_Symmetry_EK(M)
% Detects if M Voigt matrix is not symmetric
% It identifies which Cij are not equal Cji and prints to screen
%
% INPUT
% M = Voigt matrix (6x6)
%
% David Mainprice 1/08/2019
%
%% Detected Upper and Lower Cij not equal to Cji
% set error flag to zero
error_flag=0;
% scan Cij of upper trangle
for i=1:6
    for j=i:6
        if(i~=j)
            % check M(ij) not equal to M(ji)
            if(M(i,j)~=M(j,i))
            % error detected 
            error_flag=error_flag+1;
            fprintf('Detected Upper and Lower triangle Cij not equal to Cji \n');
            fprintf('Upper Triangle: C(%i%i) %6.2f \n',i,j,M(i,j));
            fprintf('Lower Triangle: C(%i%i) %6.2f \n',j,i,M(j,i));
            end
        end
    end
end
if error_flag >0
    disp('Matrix is not symmetric')
end 
if error_flag ==0
    disp('Matrix is symmetric')
end
end


function [  ] = Multi_seismic_plot_MTEX53_EK(C,Label_X,Label_Y,Label_Z,hemisphere,symm,font_size)
%% Modified by Ralf 13/03/2019
%
% Tested with MTEX 5.1.1 and MATLAB 2016a
%
% David Mainprice 13/03/2019
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
xlabel(['Vp Anisotropy = ',num2str(AVp,'%6.1f')],titleOpt{:})
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

xlabel(['Max Vs Anisotropy = ',num2str(maxAVs,'%6.1f')],titleOpt{:})

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
text([xvector,yvector,zvector],{Label_X,Label_Y,Label_Z},...
  'backgroundcolor','w','doNotDraw');
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

xlabel(['Vs1 Anisotropy = ',num2str(AVs1,'%6.1f')],titleOpt{:}) 

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
xlabel(['Vs2 Anisotropy = ',num2str(AVs2,'%6.1f')],titleOpt{:})

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
AVpVs1=200*(maxVpVs1-minVpVs1)./(maxVpVs1+minVpVs1);

xlabel(['Vp/Vs1 Anisotropy = ',num2str(AVpVs1,'%6.1f')],titleOpt{:})

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
AVpVs2=200*(maxVpVs2-minVpVs2)./(maxVpVs2+minVpVs2);

xlabel(['Vp/Vs2 Anisotropy = ',num2str(AVpVs2,'%6.1f')],titleOpt{:})

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
drawNow(gcm,'figSize','large')
set(gcf,'renderer','zBuffer')
saveFigure('Multi_seismic_plot_MTEX53_EK.eps')
%**************************************************************************
% reset to default MTEX colormap
%**************************************************************************
setMTEXpref('defaultColorMap',WhiteJetColorMap);
end



