close all; clearvars; clc

imageFolder = 'C:\Users\s130023\Desktop\Bruker data\Test viscose\Images_all';
filePattern = fullfile(imageFolder, '*.OPD');
imagefiles = dir(filePattern);

load('RH_movie'); %RH and strain data = same length as amount of images
RH_plot = [0; RH_plot];
eps_long = [0, eps_long];
eps_trans = [0, eps_trans];

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
    caxis([-30 0]);
    h = colorbar;
    set(get(h,'title'),'string','relative z [\mum]','Fontsize',24);

    subplot(1,2,2)
    hold on
%     plot(eps_long_final(1:i),'b.')
%     plot(eps_trans_final(1:i),'r.')
%     xlabel('image #'); ylabel('\epsilon [-]')   
%     axis([0 1700 -.1 0.25])
%     yyaxis right
%     plot(RH_plot(1:i),'k.')
%     ylabel('RH [%]')
%     set(gca,'YColor','k','YLim',[30 95],'FontSize',24)
%     legend('\color{blue} \epsilon_{long}','\color{red} \epsilon_{trans}', ...
%         '\color{black} RH')
%     set(legend,'Location','South','Box','On','FontSize',16)
    
    plot(RH_plot(1:i),eps_long(1:i),'b.','linewidth',15)
    plot(RH_plot(1:i),eps_trans(1:i),'r.','linewidth',15)
    axis([45 95 -0.005 0.09])
    xlabel('Relative Humidity [%]'); ylabel('Strain [-]');
    legend('\epsilon_{ll}', '\epsilon_{tt}','location','NorthWest','FontSize',24);
    set(gca,'FontSize',24)

    set(gcf,'Color','w');
    title(strcat(num2str(i*round(12.8)),' seconds'),'FontSize',24)
    M(i) = getframe(fig1);
    close all;
    i
    
end
clear imagefiles A x y z eps_long_final eps_trans_final fig1 FileName filePattern h i names names_sort RH_plot

%%
V = VideoWriter('Viscose swelling','MPEG-4');
V.FrameRate = 60;
V.Quality = 100;
open(V)
writeVideo(V,M)
close(V)
