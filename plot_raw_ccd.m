clc
close all
clear


%% ############################ READ HDF DATA ########################### %%
type = 'IAW';
shotnos = [86801,92522,92524,92525,92527,92528,92530,92531,92532,92533,92534,92535,92536,92537,92538];
fname_epw = "epw_ccd_s";
fname_iaw = "iaw_ccd_s";
EPW_data = zeros(length(shotnos),1024,1024);
IAW_data = zeros(length(shotnos),1024,1024);

fp = '/Users/soham/Documents/tsadar_4-14-25/tsadar/tsadar/external/data';

for i = 1:length(shotnos)
    shotno = string(shotnos(i));

    if shotno == "86801"
        EPW = hdfinfo(fullfile(fp,"EPW_CCD-s"+shotno+".hdf"));
        IAW = hdfinfo(fullfile(fp,"IAW_CCD-s"+shotno+".hdf"));
    else
        EPW = hdfinfo(fullfile(fp,fname_epw+shotno+".hdf"));
        IAW = hdfinfo(fullfile(fp,fname_iaw+shotno+".hdf"));
    end

    for j = 1:length(EPW.SDS.Attributes)
            AttName(i) = string(EPW.SDS.Attributes(i).Name);
            AttValue{i} = EPW.SDS.Attributes(i).Value;
    end

    data = hdfread(EPW.SDS);
    data1 = rot90(double(squeeze(data(1,:,:))),-1);
    %data1 = flip(data1,2);
    EPW_data(i,:,:) = data1;
    data = hdfread(IAW.SDS);
    data1 = rot90(double(squeeze(data(1,:,:))),-1);
    data1 = flip(data1,1);
    %data1 = flip(data1,2);
    IAW_data(i,:,:) = data1;
end

%% ##################### Axes and Calibration ######################### %%
cal.magE = 5.13;
cal.magI = 2.87; %3.8; (3.8 is the modified one to match the length of both features)
cal.EPWtcc = 565;%1024-456.1;
cal.IAWtcc = 625;%519; %625; (625 is the modified one to match the length of both features)
cal.EPWDisp = 0.27093;
cal.IAWDisp = 0.0057;
cal.EPWoff = 384.80361;
%cal.IAWoff = 529.49; Joe Katz
%cal.IAWoff = 529.62;  % 92534
%cal.IAWoff = 523.6899;  % 92535 with IAWDisp = +0.0057
cal.IAWoff = 523.74;
%cal.IAWoff = 529.58;  % 92527 
%cal.IAWoff = 529.75;  % 92528 bottom offset = 523.92
cal.cwave = 526.5;
xepw = (linspace(1,1024,1024)-cal.EPWtcc)*cal.magE;
xiaw = (linspace(1,1024,1024)-cal.IAWtcc)*cal.magI;
xiaw = 4 - xiaw/1000;
xepw = 4 - xepw/1000;
yepw = (linspace(1,1024,1024)*cal.EPWDisp)+cal.EPWoff;
yiaw = (linspace(1,1024,1024)*cal.IAWDisp)+cal.IAWoff;
%% ################### Plot ##################### %%
% 
 shot = 92537;
 n = find(shotnos == shot);  % file index(shotno) to plot

% ################### Subtract Background ########################### %
% line = 950;
% bg = [mean(squeeze(EPW_data(n,:,line))),mean(squeeze(IAW_data(n,:,line)))];
% EPW_data(n,:,:) = EPW_data(n,:,:)-bg(1);
% IAW_data(n,:,:) = IAW_data(n,:,:)-bg(2);
% ################################################################### %

elim = [-800,500];
ilim = [-800,500];
%% #################### Plot EPW ccd image ########################### %%
% subplot(2,1,1)
% box on
% hold on
% imagesc(xepw,yepw,squeeze(EPW_data(n,:,:)))
% %pcolor(squeeze(EPW_data(n,:,:)))
% shading flat
% colormap turbo
% xlim([3.55,4.75])
% ylim([min(yepw),max(yepw)])
% % % plot epw boundary lines
% %plot(linspace(ilim(1),ilim(1),length(xepw)),yepw,"w:",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xepw)),yepw,"w:",LineWidth=2)
% % plot central wavelength
% plot(xepw,linspace(cal.cwave,cal.cwave,length(yepw)),"w--",LineWidth=2)
% xlabel('x(mm)')
% ylabel('\lambda(nm)')
% fontsize(16,"points")
% ax1 = gca;
% ax1.LineWidth = 1;
% % ax1.XTick = [];
% % ax1.YTick = [];
% ax1.TickDir = 'out';
% set(gca, 'XDir', 'reverse')
% %axis("tight")
% hold off
% %title("Shot no."+string(shot))
% title("EPW")
% set(gcf,'units','inches','position',[500,500,7.5,5.7])

% %% #################### Plot IAW ccd image ########################### %%
%subplot(2,1,2)
box on
hold on
%imagesc(xiaw,yiaw,(squeeze(IAW_data(n,:,:))))
imagesc((squeeze(IAW_data(n,:,:))))
%pcolor(squeeze(IAW_data(n,:,:)))
% plot central wavelength %
%plot(xiaw,linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
shading flat
colormap turbo
% % plot iaw boundary lines %
%plot(linspace(ilim(1),ilim(1),length(xiaw)),yiaw,"w",LineWidth=2)
%plot(linspace(ilim(2),ilim(2),length(xiaw)),yiaw,"w",LineWidth=2)
%plot(linspace(elim(1),elim(1),length(xiaw)),yiaw,"w:",LineWidth=2)
%plot(linspace(elim(2),elim(2),length(xiaw)),yiaw,"w:",LineWidth=2)
%xlim([3.55,4.75])
xlim([300,800])
%ylim([524,528])
xlabel('x(mm)')
ylabel('\lambda(nm)')
fontsize(16,"points")
ax2 = gca;
ax2.LineWidth = 1;
set(gca, 'XDir', 'reverse')
% ax2.XTick = [];
% ax2.YTick = [];
ax2.TickDir = 'out';
title("IAW")
%title("Gap = [600,100]")
hold off
fontname("SansSerif")
%axis("tight")
set(gcf,'units','inches','position',[500,500,7.5,5.7])

%% #################### Plot IAW lineouts ########################### %%
% subplot(2,1,2)
% box on
% hold on
% pos = [-400,-250,-100];
% for j=1:length(pos)
%     x(j) = find(abs(xiaw - pos(j))<cal.magI,1);
%     plot(squeeze(IAW_data(n,:,x(j))),'DisplayName','x = '+string(pos(j))+'\mum',LineWidth=2)
% end
% % 92528 ~ 572 pixel
% % 92534 ~ 548 pixel
% % 92527 ~ 540 pixel
% plot(linspace(548,548,200),linspace(1,200,200),"k:",'DisplayName','pixel = 548',LineWidth=2)
% 
% legend 
% ylabel('Amplitude(arb)')
% xlabel('\lambda(pixels)')
% title("Lineouts")
% ax3 = gca;
% ax3.LineWidth = 1;
% fontsize(16,"points")
% fontname("SansSerif")
% axis("tight")


%% #################### Plot IAW stitch image ########################### %%

%% 92527 - 92530 %%
% shot = [92527,92530];
% for i = 1:length(shot)
%     n(i) = find(shotnos == shot(i));  % file index(shotno) to plot
% end
% pnt = [6,5.2];
% tcc = [632,620];
% x = zeros(length(shot),1024);
% for i=1:length(shot)
%     x(i,:) = ((linspace(1,1024,1024)-tcc(i))*cal.magI)/1000;
%     x(i,:) = pnt(i) - x(i,:);
% end
% % ################### Subtract Background ########################### %
% % line = 950;
% % bg = [mean(squeeze(EPW_data(n,:,line))),mean(squeeze(IAW_data(n,:,line)))];
% % EPW_data(n,:,:) = EPW_data(n,:,:)-bg(1);
% % IAW_data(n,:,:) = IAW_data(n,:,:)-bg(2);
% % ################################################################### %
% 
% elim = [-300,750];
% ilim = [-300,750];
% tiledlayout(1,2,"TileSpacing", "tight")
% %subplot(1,2,1)
% nexttile
% box on
% hold on
% imagesc(x(2,:),yiaw,squeeze(IAW_data(n(2),:,:)))
% % plot central wavelength %
% plot(x(2,:),linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
% shading flat
% colormap turbo
% xlim([4.5,5.5])
% ylim([min(yiaw),max(yiaw)])
% % % plot epw boundary lines
% %plot(linspace(ilim(1),ilim(1),length(xiaw)),yepw,"w:",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xiaw)),yepw,"w:",LineWidth=2)
% % plot central wavelength
% xlabel('x(mm)')
% ylabel('\lambda(nm)')
% fontsize(16,"points")
% ax1 = gca;
% ax1.LineWidth = 1;
% 
% % ax1.XTick = [];
% % ax1.YTick = [];
% ax1.TickDir = 'out';
% %axis("tight")
% hold off
% %title("Shot no."+string(shot))
% title("TS image 1")
% 
% 
% %subplot(1,2,2)
% nexttile
% box on
% hold on
% imagesc(x(1,:),yiaw,(squeeze(IAW_data(n(1),:,:))))
% % plot central wavelength %
% plot(x(1,:),linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
% shading flat
% colormap turbo
% xlim([5.5,6.7])
% ylim([min(yiaw),max(yiaw)])
% % % plot iaw boundary lines %
% %plot(linspace(ilim(1),ilim(1),length(xiaw)),yiaw,"w",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xiaw)),yiaw,"w",LineWidth=2)
% %plot(linspace(elim(1),elim(1),length(xiaw)),yiaw,"w:",LineWidth=2)
% %plot(linspace(elim(2),elim(2),length(xiaw)),yiaw,"w:",LineWidth=2)
% 
% xlabel('x(mm)')
% %ylabel('\lambda(nm)')
% fontsize(16,"points")
% ax2 = gca;
% ax2.LineWidth = 1;
% % ax2.XTick = [];
% % ax2.YTick = [];
% ax2.TickDir = 'out';
% yticklabels({})
% title("TS image 2")
% %title("Gap = [600,100]")
% hold off
% fontname("SansSerif")
% %set(gcf,'position',[500,500,500,800])

%% 92537 - 92536 - 92535 %%

% shot = [92534,92533,92532];
% for i = 1:length(shot)
%     n(i) = find(shotnos == shot(i));  % file index(shotno) to plot
% end
% pnt = [3, 5, 4];
% tcc = [625, 630, 625];
% x = zeros(length(shot),1024);
% for i=1:length(shot)
%     x(i,:) = ((linspace(1,1024,1024)-tcc(i))*cal.magI)/1000;
%     x(i,:) = pnt(i) - x(i,:);
% end
% % ################### Subtract Background ########################### %
% % line = 950;
% % bg = [mean(squeeze(EPW_data(n,:,line))),mean(squeeze(IAW_data(n,:,line)))];
% % EPW_data(n,:,:) = EPW_data(n,:,:)-bg(1);
% % IAW_data(n,:,:) = IAW_data(n,:,:)-bg(2);
% % ################################################################### %
% 
% elim = [-300,750];
% ilim = [-300,750];
% tiledlayout(1,3,"TileSpacing", "tight")
% 
% nexttile
% box on
% hold on
% imagesc(x(1,:),yiaw,squeeze(IAW_data(n(1),:,:)))
% % plot central wavelength %
% plot(x(1,:),linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
% shading flat
% colormap turbo
% xlim([2.4,4.0])
% ylim([min(yiaw),max(yiaw)])
% % % plot epw boundary lines
% %plot(linspace(ilim(1),ilim(1),length(xiaw)),yepw,"w:",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xiaw)),yepw,"w:",LineWidth=2)
% % plot central wavelength
% xlabel('x(mm)')
% ylabel('\lambda(nm)')
% fontsize(16,"points")
% ax1 = gca;
% ax1.LineWidth = 1;
% 
% % ax1.XTick = [];
% % ax1.YTick = [];
% ax1.TickDir = 'out';
% %axis("tight")
% hold off
% %title("Shot no."+string(shot))
% title("92534")
% 
% 
% 
% nexttile
% box on
% hold on
% imagesc(x(3,:),yiaw,(squeeze(IAW_data(n(3),:,:))))
% % plot central wavelength %
% plot(x(3,:),linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
% shading flat
% colormap turbo
% xlim([3.4,5])
% ylim([min(yiaw),max(yiaw)])
% % % plot iaw boundary lines %
% %plot(linspace(ilim(1),ilim(1),length(xiaw)),yiaw,"w",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xiaw)),yiaw,"w",LineWidth=2)
% %plot(linspace(elim(1),elim(1),length(xiaw)),yiaw,"w:",LineWidth=2)
% %plot(linspace(elim(2),elim(2),length(xiaw)),yiaw,"w:",LineWidth=2)
% 
% xlabel('x(mm)')
% %ylabel('\lambda(nm)')
% fontsize(16,"points")
% ax2 = gca;
% ax2.LineWidth = 1;
% % ax2.XTick = [];
% % ax2.YTick = [];
% ax2.TickDir = 'out';
% yticklabels({})
% title("92532")
% %title("Gap = [600,100]")
% hold off
% fontname("SansSerif")
% %set(gcf,'position',[500,500,500,800])
% 
% 
% nexttile
% box on
% hold on
% imagesc(x(2,:),yiaw,squeeze(IAW_data(n(2),:,:)))
% % plot central wavelength %
% plot(x(2,:),linspace(cal.cwave,cal.cwave,length(yiaw)),"w--",LineWidth=2)
% shading flat
% colormap turbo
% xlim([4.4,6])
% ylim([min(yiaw),max(yiaw)])
% ax3 = gca;
% ax3.LineWidth = 1;
% ax3.TickDir = 'out';
% yticklabels({})
% % % plot epw boundary lines
% %plot(linspace(ilim(1),ilim(1),length(xiaw)),yepw,"w:",LineWidth=2)
% %plot(linspace(ilim(2),ilim(2),length(xiaw)),yepw,"w:",LineWidth=2)
% % plot central wavelength
% xlabel('x(mm)')
% fontsize(16,"points")
% title("92533")
% ax3 = gca;
% ax3.LineWidth = 1;
