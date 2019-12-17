close all; clearvars; clc

% width_image = 80;       %Width image Bruker 100x [um]
% height_image = 60.4;	%Height image Bruker 100x [um]
width_image = 147.8;   	%Width image Bruker 100x 0.55x[um]
height_image = 111.7; 	%Height image Bruker 100x 0.55x[um]

%Goede samples batch1:   5 2; 6 2LMR; 7 2MR; 8 1MR; 9 2LMR; 11 2 LMR
photoFolder = 'C:\Users\s130023\Desktop\Bruker data\Nieuwe fibers\Batch1\';

cd(photoFolder)
topography = opdread(fullfile([photoFolder 'Sample8 2M.OPD']));
x = topography.x ;
y = topography.y ;
z = topography.Z ;

figure
imagesc(x,y,z)
set(gca,'YDir','Reverse')
axis([-1 (width_image+1) -1 (height_image+1)])
daspect([1 1 1])
xlabel('x-coordinate [\mum]'); ylabel('y-coordinate [\mum]')
h = colorbar;
ylabel(h,'relative height [\mum]','FontSize',12)

figure
surf(x,y,z)
set(gca,'YDir','Reverse');
daspect([1 1 1])
xlabel('x-coordinate [\mum]'); ylabel('y-coordinate [\mum]'); zlabel('relative height [\mum]')


