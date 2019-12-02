clear all;
close all;
clc;
warning off

%% =========================
%  Directory options
%  =========================

% imageFolder='C:\Users\s115524\Dropbox\PhD Program\test';
% formulaFolder='C:\Users\s115524\Dropbox\PhD Program\Nick Verschuur Master Thesis/globalDIC - v1.0.rc - Release Candidate';
% saveoutcomeFolder='C:\Users\s115524\Drop.box\PhD Program\test';

imageFolder='/Users/nielsvonk/Dropbox/PhD Program/test';
formulaFolder='/Users/nielsvonk/Dropbox/PhD Program/Nick Verschuur Master Thesis/globalDIC - v1.0.rc - Release Candidate';
saveoutcomeFolder='/Users/nielsvonk/Dropbox/PhD Program/test';

% imageFolder='/Users/Nick/Documents/Tue/Afstuderen/Phase II/test/Climate box fiber';
% formulaFolder='/Users/Nick/Documents/Tue/Afstuderen/Phase II/test/globalDIC - v1.0.rc - Release Candidate';
% saveoutcomeFolder='/Users/Nick/Documents/Tue/Afstuderen/Phase II/test/Climate box fiber';
cd(formulaFolder);

%% =========================
%  m-file options
%  =========================

printpdf =false;                        %Save figures as PDF
printplot = false;                       %Show figures
savegdic = false;                       %Save GDIC output
RegionDet = false;                             %ROI determination
expands = 2;                              %Expandings of the mask, x rows of pixels
corr = 2;                                %microns correction for polynomial fit
% width_image = 175.7;                    %Width image Bruker 50x [um]
% height_image =132.8;                    %Height image Bruker 50x [um]
width_image = 80;                    %Width image Bruker 100x [um]
height_image =60.4;                    %Height image Bruker 100x [um]
% width_image = 39.3;                    %Width image Bruker 100x 2x[um]
% height_image =29.7;                    %Height image Bruker 100x 2x[um]
% width_image = 254.64;                    %Width image Sensofar 50x [um]
% height_image =190.90;                    %Height image Sensofar 50x [um]
% width_image = 84.83;                    %Width image Sensofar 150x [um]
% height_image =63.60;                    %Height image Sensofar 150x [um]
pixels_x = 1376;
pixels_y = 1040;
corr_x = round((corr/width_image)*pixels_x);
corr_y = round((corr/height_image)*pixels_y);

%% =========================
%  GDIC Options
%  See [help globalDIC2D] for more info about the options.
%  help globalDIC2
%  =========================

options.list        = true;                 
options.coarsesteps = [1] ; % 0,N,-N
options.convcrit    = 1e-4 ;
options.maxiter     = 500;
options.interpmeth  = 'cubic';
options.normalized  = false;
options.basis       = 'polynomial';
options.order       = 2;
options.gradg       = false;
options.logfile     = 'globalDIC2D.log' ;
options.comment     = 'Testing version 1.0.0 RC';
options.tempfile    = 'globalDIC2Dtmp' ;
options.verbose     = 2; 
options.fulloutput  = true;
options.debug       = true ;
options.plot        = true ;
options.logfile     = 'globalDIC2D.log' ;
options.datatype    = 'double';
pol                 = num2str(options.order);
%% =========================
%  Loading the images
%  =========================

% Change directory
myFolder=imageFolder;
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get information of images in the folder
filePattern = fullfile(myFolder, '*.OPD');
imagefiles = dir(filePattern);

for i=1:length(imagefiles);
ix=strfind(imagefiles(i).name,'d');
ixx=strfind(imagefiles(i).name,'_');
aa=imagefiles(i).name(ix+1:ixx(2)-1);
imagefiles(i).number=(aa);
end

imagefiless=imagefiles;
xlsfiles={imagefiless.number};
xlsfiles=xlsfiles';
sortfiles=natsortrows(xlsfiles,1);
for i=1:length(sortfiles);
    for n=1:length(imagefiless);
        if find(strcmp(imagefiless(n).number, sortfiles(i)))>=1;
        imagefiles(i)=imagefiless(n);
        else
        end
    end
end

ii=0;
q=0;
p=0;

%% =========================
%  Calculation
%  =========================

% Load the reference image(s) & deformed image(s)
i=0;
for ii=1:length(imagefiles)
    i=1+i;
    currentfilename = imagefiles(ii).name;
    ix=strfind(currentfilename,'.'); 
    storagefilename=[currentfilename(1:ix(1)-1)];
    imagenames{i} = storagefilename;
    fullFileName = fullfile(myFolder, currentfilename);
    currentimage = opdread(fullFileName);
%         currentimage = pluread(fullFileName);
    OPD.(imagenames{i}) = currentimage;                            
    A= OPD.(imagenames{i}).Z;
    images.(imagenames{i})=A;
end

n=0;
t=0;
amountcalculations=length(imagefiles)-1;

ROI_Marges(1)= 3; %left
ROI_Marges(2)= 3; %right
ROI_Marges(3)= 12; %top      
ROI_Marges(4)= 9; %bottomm

%% test
        
% Size of f (and g)
[nn mm] = size(images.(imagenames{1}));
% Pixel dimensions (in um)
pixelsize = [width_image/mm height_image/nn];

% Create pixel position columns
x = linspace(1,mm*pixelsize(1),mm);
y = linspace(1,nn*pixelsize(2),nn);

% plot ROI

if RegionDet == 1
    figure 
    images_cell = struct2cell(images);
    reference = images_cell{1};
    imagesc(x,y,reference)
    hold on
        plot(ROI_Marges(1),ROI_Marges(3),'or')
        plot(ROI_Marges(1),height_image-ROI_Marges(4),'or')
        plot(width_image-ROI_Marges(2),ROI_Marges(3),'or')
        plot(width_image-ROI_Marges(2),height_image-ROI_Marges(4),'or')
        plot([ROI_Marges(1); width_image-ROI_Marges(2)],[ROI_Marges(3); ROI_Marges(3)],'--r')
        plot([ROI_Marges(1); width_image-ROI_Marges(2)],[height_image-ROI_Marges(4); height_image-ROI_Marges(4)],'--r')
        plot([ROI_Marges(1); ROI_Marges(1)],[ROI_Marges(3); height_image-ROI_Marges(4)],'--r')
        plot([width_image-ROI_Marges(2); width_image-ROI_Marges(2)],[ROI_Marges(3); height_image-ROI_Marges(4)],'--r')
        
        plot(ROI_Marges(1)+corr,ROI_Marges(3)+corr,'ob')
        plot(ROI_Marges(1)+corr,height_image-ROI_Marges(4)-corr,'ob')
        plot(width_image-ROI_Marges(2)-corr,ROI_Marges(3)+corr,'ob')
        plot(width_image-ROI_Marges(2)-corr,height_image-ROI_Marges(4)-corr,'ob')
        plot([ROI_Marges(1)+corr; width_image-ROI_Marges(2)-corr],[ROI_Marges(3)+corr; ROI_Marges(3)+corr],'--b')
        plot([ROI_Marges(1)+corr; width_image-ROI_Marges(2)-corr],[height_image-ROI_Marges(4)-corr; height_image-ROI_Marges(4)-corr],'--b')
        plot([ROI_Marges(1)+corr; ROI_Marges(1)+corr],[ROI_Marges(3)+corr; height_image-ROI_Marges(4)-corr],'--b')
        plot([width_image-ROI_Marges(2)-corr; width_image-ROI_Marges(2)-corr],[ROI_Marges(3)+corr; height_image-ROI_Marges(4)-corr],'--b')
    set(gca, 'ydir', 'reverse');
    colorbar
    xlabel 'x-coordinate [\mum]'
    ylabel 'y-coordinate [\mum]'
    title 'ROI Determination'
    pause(30)
end

% Center the axes
x = x - mean(x);
y = y - mean(y);
[X Y] = meshgrid(x,y);
        
f=images.(imagenames{1});
valid=~isnan(f);
f = griddata(X(valid),Y(valid),f(valid),X,Y,'linear');

tilt_x=0.4;
[h w]=size(images.(imagenames{1}));
mid_x=round(w/2);
tilt_p_cell_x=tilt_x/(mid_x-1);
til_x=zeros(h,w);

tilt_y=((tilt_x/39.3)*29.7);
mid_y=round(h/2);
tilt_p_cell_y=tilt_y/(mid_y-1);
til_y=zeros(h,w);

ii=0;
for i=1:w;
    til_x(1:h,i)=tilt_x-ii*tilt_p_cell_x;
    ii=ii+1;
end

ii=0;
for i=1:h;
    til_y(i,1:w)=tilt_y-ii*tilt_p_cell_y;
    ii=ii+1;
end


%%
s=0;
ss=0;
count=0;
z = 1;
while amountcalculations ~= n
    
    n=1+n
    ss=4+ss;
    s=s+1;
    close all;
        
    % Size of f (and g)
    [nn mm] = size(images.(imagenames{1}));

    % Pixel dimensions (in um)
    pixelsize = [width_image/mm height_image/nn];
    
    x1 = round(ROI_Marges(1)/pixelsize(1));
    x2 = round(ROI_Marges(2)/pixelsize(1));
    y1 = round(ROI_Marges(3)/pixelsize(2));
    y2 = round(ROI_Marges(4)/pixelsize(2));
    
    % Create pixel position columns
    x = linspace(1,mm*pixelsize(1),mm);
    y = linspace(1,nn*pixelsize(2),nn);

    % Center the axes
    x = x - mean(x);
    y = y - mean(y);
        
    f=images.(imagenames{1});
    g=images.(imagenames{n+1});


    dis_x_max=0.5817;
    dis_x_min=0.312;
    pix_dis_max=round(dis_x_max/pixelsize(1));
    pix_dis_min=round(dis_x_min/pixelsize(1));
    [h w]=size(images.(imagenames{1}));
    pix_dif=(pix_dis_max-pix_dis_min)/h;

    ii=0;
    for i=1:h;
        X_disp(i,1:w)=pix_dis_min+ii*pix_dif;
        ii=ii+1;
    end
    
    X_disp=round(X_disp);
    g=images.(imagenames{1+1});
    for i=1:h
        g(i,:)=circshift(g(i,:),X_disp(i,1));
    end
        
    H = fspecial('gaussian', [3 3], 0.8) ;
    f = imfilter(f,H,'symmetric','conv');
    g = imfilter(g,H,'symmetric','conv');
    cd(formulaFolder);
        
    % Plot the position field figures
    fname=char(imagenames(1));
    gname=char(imagenames(n+1));
    currentfilename=[fname '  -  ' gname '  -  ' pol ' order poly'];
    logfilename=[ currentfilename '.txt'];
    cd(imageFolder)
        
    if (exist(logfilename, 'file') == 2);
        delete(logfilename);
    end
        diary (logfilename)
        cd(formulaFolder);
    if printplot == 1;
        plotop.name     = [ 'Position field: ' currentfilename] ;
        plotop.titles   = {'f','g'};
        plotop.colorlabel = {}; 
        [hh B] = globalDICplotposition(x,y,f,g,plotop); 
        hold on
    end

    % Save the position field figures
    if printpdf == 1;
        fname=char(imagenames(1));
        gname=char(imagenames(n+1));
        currentfilename=['Position field  ' fname '  -  ' gname '  -  ' pol ' order poly'];
        savepdf(currentfilename,myFolder);
        savefig(currentfilename)
    end
            
    cd(formulaFolder);
    
    % =========================
    % Start the digital image correlation
    % =========================

    % Clear the old logfile
    if isfield(options,'logfile')
        delete(options.logfile);
    end
    
    % Generate a full basis for polynomial or chebyshev functions
    % =================
    phi = dofbuild_poly(options.order);

    % Initial Guess
    % =========================
    Ndof = size(phi,1);
%   u    = zeros(Ndof,1)
    if n==1;
        u = zeros(Ndof,1)
%         cd(imageFolder)
%         uu = load('init.mat');
%         uuu = struct2cell(uu);
%         u = uuu{1}
%         cd(formulaFolder)
    else
        if convergence ==1;
            u=gdic.u
        else
            u
        end
    end
         
    options.ROI(1) = x(1)   + ROI_Marges(1);
    options.ROI(2) = x(end) - ROI_Marges(2);
    options.ROI(3) = y(1)   + ROI_Marges(3);
    options.ROI(4) = y(end) - ROI_Marges(4);

    if printplot == 1;
        subplot(1,2,1);     
        hold on
        imagesc(x,y,A);
        set(gca, 'ydir', 'reverse');
        hold on
        plot(options.ROI(1),options.ROI(3),'or')
        plot(options.ROI(1),options.ROI(4),'or')
        plot(options.ROI(2),options.ROI(3),'or')
        plot(options.ROI(2),options.ROI(4),'or')
        plot([options.ROI(1); options.ROI(2)],[options.ROI(3); options.ROI(3)],'--r')
        plot([options.ROI(1); options.ROI(2)],[options.ROI(4); options.ROI(4)],'--r')
        plot([options.ROI(1); options.ROI(1)],[options.ROI(3); options.ROI(4)],'--r')
        plot([options.ROI(2); options.ROI(2)],[options.ROI(3); options.ROI(4)],'--r')
    end

    % Save the position field figures with ROI
    if printpdf == 1;
        fname=char(imagenames(1));
        gname=char(imagenames(n+1));
        currentfilename=['Position field with ROI  ' fname '  -  ' gname '  -  ' pol ' order poly'];
        savepdf(currentfilename,myFolder);
        savefig(currentfilename)
    end                    
            
    cd(formulaFolder);

            
    % The working part
    f=images.(imagenames{1});
            
    figure(10);
    imagesc(x,y,f);
              
    g=images.(imagenames{s+1});
    g(isnan(g))=nanmean(g(:));

    [exp_mask NaN_mask]=Maskexpand(isnan(f),expands);
            
    mask_in=find(exp_mask);
    options.mask=mask_in;
    
            
    namelogfile     = [fname '  -  ' gname '  -  ' pol ' order poly' '.txt'] ;
    if n == 1
        [gdic crs convergence mask M m Phi b]  = globalDIC2D(f,g,x,y,u,phi,options);
    else 
        
        S = size(gdic.h);
        Sy = 1040-S(1)-y1-y2;
        addit_y = Sy/2;
        y11 = y1+addit_y;
        y22 = y2+addit_y;
        Sx = 1376-S(2)-x1-x2;
        addit_x = Sx/2;
        x11 = x1+addit_x;
        x22 = x2+addit_x; 
        M1 = zeros(1040,x11);
        M2 = zeros(1040,x22);
        M3 = zeros(y11,S(2));
        M4 = zeros(y22,S(2));
        H = [M3; gdic.h; M4];
        HH = [M1 H M2];
        [gdic crs convergence mask M m Phi b]  = globalDIC2D(HH,g,x,y,u,phi,options);
    end
    
    NaN_mask = NaN_mask(crs.indy,crs.indx);
    % =========================
    % Plot the results
    % =========================

    % Extract some data from the gdic structure
    Ux = gdic.Ux;
    Uy = gdic.Uy;
    Uz = gdic.Uz;
    x  = gdic.x;
    y  = gdic.y;


    % Residual
    r = gdic.r;
    r_correct=r;
    r_correct(find(NaN_mask==1))=NaN;
    gdic.r_correct=r_correct;

    % Plot residual field
    if printplot == 1;
        plotop.name  = 'Residual' ;
        plotop.titles   = {'r'};
        plotop.colorlabel = {};
        hh = globalDICplotresidual(x,y,r_correct,plotop);
    end
            
    % Save the residual field figure(s)
    if printpdf == 1;
        fname=char(imagenames(1));
        gname=char(imagenames(s+1));
        currentfilename=['Residual field  ' fname '  -  ' gname '  -  ' pol ' order poly'];
        savepdf(currentfilename,myFolder);
        savefig(currentfilename)
    end
            
    % Plot displacement fields
    if printplot == 1;
        plotop.name       = 'Displacement Fields' ;
        plotop.titles     = {'Ux','Uy','Uz'};
        plotop.colorlabel = {'U_i [\mum]'};
        cd(formulaFolder);
        hh = globalDICplotdisplacement(r_correct,x,y,Ux,Uy,Uz,plotop);
    end
            
    % Save the displacement fields figure(s)
    if printpdf == 1;
        fname=char(imagenames(1));
        gname=char(imagenames(s+1));
        currentfilename=['Displacement fields  ' fname '  -  ' gname '  -  ' pol ' order poly'];
        cd(formulaFolder);
        savepdf(currentfilename,myFolder);
        savefig(currentfilename)        
    end
                
                
                
    % Calculate strain
    dx = mean(diff(x));
    dy = mean(diff(y));

    [X Y] = meshgrid(x,y);
    [nnn mmm] = size(X);

    % Assume innitially flat surface
    Z = zeros(nnn,mmm);
    
    [dXx dXy] = gradient(X+Ux,dx,dy);
    [dYx dYy] = gradient(Y+Uy,dx,dy);
    [dZx dZy] = gradient(Z+Uz,dx,dy);

    Lx = hypot(dXx,dZx);
    Ly = hypot(dYy,dZy);

    % This strain definition works for stretched membranes, but is far from
    % universal, test, check, and verify, before using.
    Exx = Lx - 1;
    Eyy = Ly - 1;
    EX_proj(z) = mean(mean(Exx));
    EY_proj(z) = mean(mean(Eyy));
    cd(formulaFolder);
        
    if n==1;
        x_short=x;
        y_short=y;
        Z_short=gdic.f;
        [fitresult, gof, Z_surf] = create_surf_Fit(x_short, y_short, Z_short);
    end

    
    Uxx=Ux;
    Uyy=Uy;
    Uzz=Uz;
    [Exx_jan_true, Eyy_jan_true, Exy_jan] = jan_strain_initially_nonflat_true_curvature(x,y,Uxx,Uyy,Uzz,Z_surf);
    
    % Determine mean strain (Niels)
    EXX(z) = mean(mean(Exx_jan_true(corr_x:end-corr_x,corr_y:end-corr_y)));
    EYY(z) = mean(mean(Eyy_jan_true(corr_x:end-corr_x,corr_y:end-corr_y)));
    
    EXXX(z) = mean(mean(Exx_jan_true(5:end-5,5:end-5)));
    EYYY(z) = mean(mean(Eyy_jan_true(5:end-5,5:end-5)));
%     EXXX{z} = Exx_jan_true;
%     EYYY{z} = Eyy_jan_true;
    res(z) = nanmean(nanmean(gdic.r));
        
    % Plot strain fields image(s)
    if printplot == 1;
        plotop.name       = 'Strain Fields' ;
        plotop.titles     = {'\epsilon_x','\epsilon_y'};
        plotop.colorlabel = {'\epsilon [%]'};
        cd(formulaFolder);     
        hh = globalDICplotstrains(r_correct,x,y,1e2*Exx_jan_true,1e2*Eyy_jan_true,plotop);
    end  
    
    % Save strain fields image(s)
    if printpdf == 1;
        fname=char(imagenames(1));
        gname=char(imagenames(s+1));
        currentfilename=['Strain fields  ' fname '  -  ' gname '  -  ' pol ' order poly'];
        savepdf(currentfilename,myFolder);
        savefig(currentfilename)
    end

    fname=char(imagenames(1));
    gname=char(imagenames(s+1));
    currentfilename=['GDIC output  ' fname '  -  ' gname];
    GDIC_output(s).name=currentfilename;
    GDIC_output(s).gdic=gdic;
    
    if amountcalculations == n
        if savegdic == 1;
            cd(imageFolder)
            save ('GDIC_output_Pol2_p1', 'GDIC_output','-v7.3')
    %       save ('GDIC_output_Bruker_Best_VSI_50x_Y_POL_0', 'GDIC_output','-v7.3')
        end
    else
    end
    
    cd(formulaFolder);
     

    if convergence==1;
        Convergence='Yes'
    else
        Convergence='No'
    end
    diary off
    z =z+1;
end 

INIT = gdic.u;
cd(imageFolder);
save('init.mat','INIT','-v7.3')
save('strains.mat','EXX','EYY','-v7.3')
save('strains2.mat','EXXX','EYYY','-v7.3')
save('residual.mat','res','-v7.3')
save('strains_proj.mat','EX_proj','EY_proj','-v7.3')

% f=GDIC_output.gdic.f;
% f(find(NaN_mask==1))=NaN;
% g=f;
%  cd(formulaFolder);
% 
%         plotop.name     = [ 'Position field: ' currentfilename] ;
%         plotop.titles   = {'f','g'};
%         plotop.colorlabel = {}; 
%         [hh B] = globalDICplotposition(x,y,f,g,plotop); 
% % close all;    
% % cd('/Users/Nick/Documents/Tue/Afstuderen/Phase II')
