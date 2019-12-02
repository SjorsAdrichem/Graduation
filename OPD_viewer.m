close all; clearvars; clc

width_image = 80;     	%Width image Bruker 100x [um]
height_image = 60.4;	%Height image Bruker 100x [um]

photoFolder = 'C:\Users\s130023\Desktop\Bruker data\Test UPy-PEG 10k\Images all\';
topography = opdread(fullfile([photoFolder 'M100x_d516_T12_WL_13_32_10.OPD']));
x = topography.x ;
y = topography.y ;
z = topography.Z ;

figure
imagesc(x,y,z)
set(gca,'YDir','Reverse');
axis([-1 (width_image+1) -1 (height_image+1)])
xlabel('x-coordinate [\mum]'); ylabel('y-coordinate [\mum]')
h = colorbar;
ylabel(h,'relative height [\mum]','FontSize',11)


