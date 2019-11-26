close all; clearvars; clc

imageFolder = 'C:\Users\s130023\Desktop\Bruker data\Test UPy-PEG 10k 2\Images all';
filePattern = fullfile(imageFolder, '*.OPD');
imagefiles = dir(filePattern);

% load('movie'); %RH and strain data = same length as amount of images
load('RH_movie'); %RH

for i = 1:length(imagefiles)
    names{i} = imagefiles(i).name;
end

names_sort = natsortfiles(names);

for i = 1:length(names)
    FileName = names_sort{i};
    A = opdread(FileName);
    x = A.x;
    y = A.y;
    z = A.Z;
    
    fig1 = figure('Renderer', 'painters', 'Position', [100 100 1500 600]);
    hold on
    subplot(1,2,1)
    imagesc(x,y,z);
    xlabel('x [\mum]'); ylabel('y [\mum]')
    title(strcat(num2str(i*round(12.8)),' seconds'),'FontSize',24)
    caxis([-30 0]);
    h = colorbar;
    set(get(h,'title'),'string','relative z [\mum]','Fontsize',24);
    set(gca,'FontSize',24)
    subplot(1,2,2)
    hold on
%     plot(RH_plot(1:i),eps_long(1:i),'b.','linewidth',15)
%     plot(RH_plot(1:i),eps_trans(1:i),'r.','linewidth',15)
%     axis([28 90 -0.01 0.35]);
%     xlabel('Relative Humidity [%]'); ylabel('Strain [-]');
%     legend('\epsilon_{ll}', '\epsilon_{tt}','location','NorthWest','FontSize',24);
%     set(gca,'FontSize',24)
% %     text(35,2,strcat(num2str(i*round(12.8)),' seconds'),'FontSize',24);
    plot(RH_plot(i,2),RH_plot(i,3),'k.','linewidth',15)
    axis([0 12 25 95])
    xlabel('Time [h]'); ylabel('RH [5]')
    set(gca,'FontSize',24)
    set(gcf,'Color','w');
    M(i) = getframe(fig1);
    close all;
    i
    
end
 
V = VideoWriter('Expansion without strain','MPEG-4');
V.FrameRate = 60;
V.Quality = 100;
open(V)
writeVideo(V,M)
close(V)
