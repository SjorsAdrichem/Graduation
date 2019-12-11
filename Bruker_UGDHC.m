close all; clearvars; clc
warning off

%% =========================
%  Directory options
%  =========================
imageFolder =       'C:\Users\s130023\Desktop\Bruker data\Test UPy-PEG 10k 2\Images_processing';
formulaFolder =     'C:\Users\s130023\Documents\MATLAB';
saveoutcomeFolder = 'C:\Users\s130023\Desktop\Bruker data\Test UPy-PEG 10k 2\Outcome';

cd(formulaFolder);

%% ========================= 
%  m-file options
%  =========================

printplot = false;              %Show figures
printpdf = false;               %Save figures as (.pdf or) .fig and .png
savegdic = false;               %Save GDIC output
RegionDet = true;               %ROI determination
expands = 2;                    %Expandings of the mask, x rows of pixels
corr = 4;                       %microns correction for polynomial fit; 2 for 200x, 4 for 100x
% width_image = 175.7;            %Width image Bruker 50x [um]
% height_image = 132.8;           %Height image Bruker 50x [um]
% width_image = 147.8;            %Width image Bruker 100x 0.55x[um]
% height_image = 111.7;           %Height image Bruker 100x 0.55x[um]
width_image = 80;               %Width image Bruker 100x [um]
height_image = 60.4;            %Height image Bruker 100x [um]
% width_image = 39.3;             %Width image Bruker 100x 2x[um]
% height_image = 29.7;            %Height image Bruker 100x 2x[um]
pixels_x = 1376;
pixels_y = 1040;
corr_x = round((corr/width_image)*pixels_x);    %round number of pixels in x after correction 
corr_y = round((corr/height_image)*pixels_y);   %round number of pixels in y after correction

%% =========================
%  GDIC Options
%  =========================

options.list        = true;                         %true if coarsesteps are specified explicitly
options.coarsesteps = 1;                            %no coarsening == 1
options.convcrit    = 1e-4;                         %convergence criteria
options.maxiter     = 500;                          %max. amount of interations per step
options.interpmeth  = 'cubic';                      %interpolation method
options.normalized  = false;                        %normalized domain of polynomial base functions
options.basis       = 'polynomial';                 %phi(k) = [a b d] = x^a * y^b * e_d
options.order       = 2;                            %order of base function
options.gradg       = false;                        %use of the gradient of g
options.logfile     = 'globalDIC2D.log';            %log the correlation progress
options.comment     = 'Testing version 1.0.0 RC';   %a comment for in the logfile header
options.tempfile    = 'globalDIC2Dtmp';             %reduce the memory footprint
options.verbose     = 2;                            %sets the level of output
options.fulloutput  = true;                         %a more full output (gdic) structure
options.debug       = false;                        %creates DbugXX structures for coarse graining step containing all usefull variables
options.plot        = true;                         %generate a lot of figures for each coarse grain step
options.datatype    = 'double';                     %set floating point precision

pol                 = num2str(options.order);

%% =========================
%  Loading the images
%  =========================
if ~isfolder(imageFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', imageFolder);
    uiwait(warndlg(errorMessage));
    return
end

% Get information of images in the folder
filePattern = fullfile(imageFolder, '*.OPD');  %build full file name from parts
imagefiles = dir(filePattern);              %create struct with image information

for i = 1:length(imagefiles)                %break up name of image for information
    ix = strfind(imagefiles(i).name,'d');   
    ixx = strfind(imagefiles(i).name,'_');
    aa = imagefiles(i).name(ix+1:ixx(2)-1); %get consecutive number of image
    imagefiles(i).number = (aa);
end

% Put the images in the right order
imagefiless = imagefiles;
sortfiles = natsortrows({imagefiless.number}',1); %alphanumeric sort of the rows
for i = 1:length(sortfiles)
    for n = 1:length(imagefiless)
        if find(strcmp(imagefiless(n).number, sortfiles(i))) >= 1
            imagefiles(i) = imagefiless(n);
        else
        end
    end
end

% Load the reference image & deformed image(s)
imagenames = cell(1,length(imagefiles)); %preallocating
for i = 1:length(imagefiles)
    currentfilename = imagefiles(i).name;
    ix = strfind(currentfilename,'.');
    imagenames{i} = (currentfilename(1:ix(1)-1));   
end

OPD.(imagenames{1}) = opdread(fullfile(imageFolder, [char(imagenames(1)) '.OPD']));
A = OPD.(imagenames{1}).Z;
image_f = A;

clear ix ixx aa i n sortfiles imagefiless A filePattern currentfilename OPD

amountcalculations = length(imagefiles)-1;

%% =========================
%  ROI
%  =========================

% Determine pixelsize and plot ROI       
ROI_Marges(1)= 11; %left
ROI_Marges(2)= 12; %right
ROI_Marges(3)= 18; %top      
ROI_Marges(4)= 18; %bottom

% Size of f (and g)
[nn, mm] = size(image_f);
% Pixel dimensions (in um)
pixelsize = [width_image/mm height_image/nn];
% Create pixel position columns
x = linspace(1,mm*pixelsize(1),mm);
y = linspace(1,nn*pixelsize(2),nn);

if RegionDet == 1 %plot ROI of reference image
    figure; hold on
    imagesc(x,y,image_f)
    
    % plot ROI without corrections
    plot([ROI_Marges(1) ROI_Marges(1) width_image-ROI_Marges(2) width_image-ROI_Marges(2)], ...
        [ROI_Marges(3) height_image-ROI_Marges(4) ROI_Marges(3) height_image-ROI_Marges(4)],'or')
    plot([ROI_Marges(1) width_image-ROI_Marges(2) width_image-ROI_Marges(2) ROI_Marges(1) ROI_Marges(1)], ...
        [ROI_Marges(3) ROI_Marges(3) height_image-ROI_Marges(4) height_image-ROI_Marges(4) ROI_Marges(3)],'--r') 
    % plot ROI with corrections
    plot([ROI_Marges(1)+corr ROI_Marges(1)+corr width_image-ROI_Marges(2)-corr width_image-ROI_Marges(2)-corr], ...
        [ROI_Marges(3)+corr height_image-ROI_Marges(4)-corr ROI_Marges(3)+corr height_image-ROI_Marges(4)-corr],'ob')
    plot([ROI_Marges(1)+corr width_image-ROI_Marges(2)-corr width_image-ROI_Marges(2)-corr ROI_Marges(1)+corr ROI_Marges(1)+corr], ...
        [ROI_Marges(3)+corr ROI_Marges(3)+corr height_image-ROI_Marges(4)-corr height_image-ROI_Marges(4)-corr ROI_Marges(3)+corr],'--b')

    set(gca,'YDir','Reverse');
    xlabel('x-coordinate [\mum]'); ylabel('y-coordinate [\mum]'); title('ROI determination')
    colorbar

    strainPoint = [40, 30]; %point [in micron] where strain is determined in ROI 
    plot(strainPoint(1),strainPoint(2),'ko','MarkerSize',10);
    strainPoint(1) = round((strainPoint(1)-ROI_Marges(1))/pixelsize(1));
    strainPoint(2) = round((strainPoint(2)-ROI_Marges(3))/pixelsize(2));
    
    pause(3)
end
clear x y

% return 

%% Calculation loop
[EXX, EYY, EXY, EXX_proj, EYY_proj, EXY_proj, EXXp, EYYp, EXYp, res] = deal(zeros(1,length(imagefiles)-1)); %preallocating

fname = char(imagenames(1));

x1 = ROI_Marges(1)/pixelsize(1);
x2 = ROI_Marges(2)/pixelsize(1);
y1 = ROI_Marges(3)/pixelsize(2);
y2 = ROI_Marges(4)/pixelsize(2);

n = 1; %image number, and index for saving
while amountcalculations ~= n-1
    fprintf('n = %d \n',n);
    close all;

    OPD.(imagenames{n+1}) = opdread(fullfile(imageFolder, [char(imagenames(n+1)) '.OPD']));
    A = OPD.(imagenames{n+1}).Z;
    image_g = A;
    clear OPD A
    
    % Create pixel position columns
    x = linspace(1,mm*pixelsize(1),mm);
    y = linspace(1,nn*pixelsize(2),nn);
    % Center the axes
    x = x - mean(x);
    y = y - mean(y);
        
%     f = image_f;
%     g = image_g;
%         
%     H = fspecial('gaussian', [3 3], 0.8) ;
%     f = imfilter(f,H,'symmetric','conv');
%     g = imfilter(g,H,'symmetric','conv');
%     cd(formulaFolder);
        
    % Plot the position field figures
    gname = char(imagenames(n+1));
    currentfilename = [fname '  -  ' gname '  -  ' pol ' order poly'];
    logfilename = [ currentfilename '.txt'];
    cd(imageFolder)
        
    if (exist(logfilename, 'file') == 2)
        delete(logfilename);
    end
        diary (logfilename)
        cd(formulaFolder);
        
%     if printplot == 1
%         plotop.name     = [ 'Position field: ' currentfilename] ;
%         plotop.titles   = {'f','g'};
%         plotop.colorlabel = {}; 
%         [HH, B] = globalDICplotposition(x,y,f,g,plotop); 
%         hold on
%     end

%     % Save the position field figures
%     if printpdf == 1
%         currentfilename = ['Position field  ' fname '  -  ' gname '  -  ' pol ' order poly'];
%         cd(saveoutcomeFolder);
% %         savepdf(currentfilename,imageFolder);
%         savefig(currentfilename)
%         saveas(gcf,[currentfilename '.png'])
%         cd(formulaFolder);
%     end
            
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

    if n == 1
        u = zeros(Ndof,1)
%         cd(imageFolder)
%         uu = load('init.mat');
%         uuu = struct2cell(uu);
%         u = uuu{1}
%         load('h_UGDHC.mat');
%         cd(formulaFolder)
    else
        if convergence == 1
            u = gdic.u
        else
            u
        end
    end
         
    options.ROI(1) = x(1)   + ROI_Marges(1);
    options.ROI(2) = x(end) - ROI_Marges(2);
    options.ROI(3) = y(1)   + ROI_Marges(3);
    options.ROI(4) = y(end) - ROI_Marges(4);

%     if printplot == 1
%         figure; hold on
%         imagesc(x,y,f);        
%         plot([options.ROI(1) options.ROI(1) options.ROI(2) options.ROI(2)], ...
%             [options.ROI(3) options.ROI(4) options.ROI(3) options.ROI(4)],'or')
%         plot([options.ROI(1) options.ROI(2) options.ROI(2) options.ROI(1) options.ROI(1)], ...
%             [options.ROI(3) options.ROI(3) options.ROI(4) options.ROI(4) options.ROI(3)],'--r')
%         set(gca, 'ydir', 'reverse');
%    end

%     % Save the position field figures with ROI
%     if printpdf == 1
%         currentfilename=['Position field with ROI  ' fname '  -  ' gname '  -  ' pol ' order poly'];
%         cd(saveoutcomeFolder)    
% %         savepdf(currentfilename,imageFolder);
%         savefig(currentfilename)
%         saveas(gcf,[currentfilename '.png'])
%         cd(formulaFolder);
%     end                                
            
    cd(formulaFolder);
            
    % The working part
    f = image_f;            
    figure(10);
    imagesc(x,y,f);              
    g = image_g;
    g(isnan(g)) = nanmean(g(:));

    [exp_mask, NaN_mask] = Maskexpand(isnan(f),expands);           
    mask_in = find(exp_mask);
    
    % comment for mask 0 expands Niels ----
    options.mask = mask_in;
    %--------------------------------------
       
    namelogfile     = [fname '  -  ' gname '  -  ' pol ' order poly' '.txt'];
    
    if n == 1 && str2double(imagefiles(2).number) == 1 + str2double(imagefiles(1).number)
        [gdic, crs, convergence, mask, M, m, Phi, b] = globalDIC2D(f,g,x,y,u,phi,options);
    else
        add_left = crs.indx(1)-1;
        add_right = pixels_x-crs.indx(end);
        add_top = crs.indy(1)-1;
        add_bottom = pixels_y-crs.indy(end);

        M3 = zeros(add_top,length(crs.indx));    	%zeros matrix left of ROI
        M4 = zeros(add_bottom,length(crs.indx)); 	%zeros matrix right of ROI
        M1 = zeros(pixels_y,add_left);              %zeros matrix on top of ROI
        M2 = zeros(pixels_y,add_right);             %zeros matrix on bottom of ROI
        hh = [M3; gdic.h; M4];                      %H with right # of pixels in x-dir
        hhh = [M1 hh M2];                           %H with right # of pixels in y-dir
        [gdic, crs, convergence, mask, M, m, Phi, b] = globalDIC2D(hhh,g,x,y,u,phi,options);
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
    r_correct = r;
    r_correct(find(NaN_mask == 1)) = NaN;
    gdic.r_correct = r_correct;

%     % Plot residual field
%     if printplot == 1
%         plotop.name  = 'Residual' ;
%         plotop.titles   = {'r'};
%         plotop.colorlabel = {};
%         HH = globalDICplotresidual(x,y,r_correct,plotop);
%     end
            
%     % Save the residual field figure(s)
%     if printpdf == 1
%         currentfilename = ['Residual field  ' fname '  -  ' gname '  -  ' pol ' order poly'];
%     	cd(saveoutcomeFolder)
% %         savepdf(currentfilename,imageFolder);
%         savefig(currentfilename)
%         saveas(gcf,[currentfilename '.png'])
%         cd(formulaFolder)
%     end
            
%     % Plot displacement fields
%     if printplot == 1
%         plotop.name = 'Displacement Fields';
%         plotop.titles = {'Ux','Uy','Uz'};
%         plotop.colorlabel = {'U_i [\mum]'};
%         HH = globalDICplotdisplacement(r_correct,x,y,Ux,Uy,Uz,plotop);
%     end
            
%     % Save the displacement fields figure(s)
%     if printpdf == 1
%         currentfilename = ['Displacement fields  ' fname '  -  ' gname '  -  ' pol ' order poly'];
%         cd(saveoutcomeFolder);
% %         savepdf(currentfilename,imageFolder);
%         savefig(currentfilename)        
%         saveas(gcf,[currentfilename '.png'])
%         cd(formulaFolder)
%     end                
                               
    % ========================
    % Calculate strain
    % ========================
%     dx = mean(diff(x));
%     dy = mean(diff(y));
% 
%     [X, Y] = meshgrid(x,y);
%     [nnn, mmm] = size(X);
% 
%     % Assume innitially flat surface
%     Z = zeros(nnn,mmm);
    
    if n == 1
        x_short = x;
        y_short = y;
        Z_short = gdic.f;
        [fitresult, gof, Z_surf] = create_surf_Fit(x_short, y_short, Z_short);
    end
    
    Uxx = Ux;
    Uyy = Uy;
    Uzz = Uz;
    [Exx_jan_true, Eyy_jan_true, Exy, Exx, Eyy] = jan_strain_initially_nonflat_true_curvature(x,y,Uxx,Uyy,Uzz,Z_surf);
    
    % Determine mean strain (Niels)
    EXX(n) = mean(mean(Exx_jan_true(corr_x:end-corr_x,corr_y:end-corr_y)));
    EYY(n) = mean(mean(Eyy_jan_true(corr_x:end-corr_x,corr_y:end-corr_y)));
    EXY(n) = mean(mean(Exy(corr_x:end-corr_x,corr_y:end-corr_y)));
    
    EXX_proj(n) = mean(mean(Exx(corr_x:end-corr_x,corr_y:end-corr_y)))-1;
    EYY_proj(n) = mean(mean(Eyy(corr_x:end-corr_x,corr_y:end-corr_y)))-1;
    EXY_proj(n) = mean(mean(Exy(corr_x:end-corr_x,corr_y:end-corr_y)));
        
    EXXp(n) = Exx_jan_true(strainPoint(2),strainPoint(1));
    EYYp(n) = Eyy_jan_true(strainPoint(2),strainPoint(1));
    EXYp(n) = Exy(strainPoint(2),strainPoint(1));
       
    res(n) = nanmean(nanmean(gdic.r));
        
%     % Plot strain fields image(s)
%     if printplot == 1
%         plotop.name = 'Strain Fields';
%         plotop.titles = {'\epsilon_x','\epsilon_y'};
%         plotop.colorlabel = {'\epsilon [%]'};
%         cd(formulaFolder);     
%         HH = globalDICplotstrains(r_correct,x,y,Exx_jan_true,Eyy_jan_true,plotop);
%     end  
    
%     % Save strain fields image(s)
%     if printpdf == 1
%         currentfilename = ['Strain fields  ' fname '  -  ' gname '  -  ' pol ' order poly'];
%         cd(saveoutcomeFolder)
% %         savepdf(currentfilename,imageFolder);
%         savefig(currentfilename)
%         saveas(gcf,[currentfilename '.png'])
%         cd(formulaFolder);     
%     end

%     currentfilename=['GDIC output  ' fname '  -  ' gname];
%     GDIC_output(s).name=currentfilename;
%     GDIC_output(s).gdic=gdic;
% 
%     if amountcalculations == n
%         if savegdic == 1
%             cd(saveoutcomeFolder)
%             save ('GDIC_output_Pol2_p1', 'GDIC_output','-v7.3')
% %             save ('GDIC_output_Bruker_Best_VSI_50x_Y_POL_0', 'GDIC_output','-v7.3')
%         end
%     end
    
    cd(formulaFolder);
     
    if convergence == 1
        Convergence = 'Yes';
    else
        Convergence = 'No';
    end
    fprintf('\nConvergence = %s \n \n',Convergence);
    
    diary off
    n = n+1;        
end

%% saving data in .mat and plotting strains
INIT = gdic.u;
cd(imageFolder);
save('init.mat','INIT','-v7.3')
% save('h_UGDHC.mat','gdic','crs','-v7.3')
save('strains.mat','EXX','EXY','EYY','EXXp','EXYp','EYYp','EXX_proj','EXY_proj','EYY_proj','-v7.3') 
save('residual.mat','res','-v7.3')
imagenames_set = {imagenames{1}; imagenames{2}; imagenames{end}};
save('set_data.mat','ROI_Marges','imagenames_set','-v7.3')

figure; hold on
plot(EXX,'b.')
plot(EYY,'r.')

% close all
