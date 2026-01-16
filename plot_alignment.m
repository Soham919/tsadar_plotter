clc
close all
clear


%% ############################ READ HDF DATA ########################### %%
fname_epw = ["EPW_526p5cw_1500um.hdf", "EPWCCD_0cw__50um_slit_500um_off_tcc.hdf", "EPWCCD_0cw__50um_slit_2500um_off_tcc.hdf", "EPWCCD_0cw__50um_slit_4500um_off_tcc.hdf"];
fname_iaw = ["IAWCCD_0nmcw_1500toP4hdf.hdf","IAWCCD_532cw_500um_off_tcc.hdf", "IAWCCD_532cw_1300um_off_tcc.hdf", "IAWCCD_532cw_1500um_off_tcc.hdf", "IAWCCD_532cw_2500um_off_tcc.hdf", "IAWCCD_532cw_4500um_off_tcc.hdf"];
fname_TPDI = "TPDI_1500um.hdf";
epw_algn = zeros(length(fname_epw),2,1024,1024);
iaw_algn = zeros(length(fname_iaw),2,1024,1024);

for i = 1:length(fname_iaw)
    buf = hdfinfo(fname_iaw(i));
    iaw_algn(i,:,:,:) = hdfread(buf.SDS);
end

for i = 1:length(fname_epw)
    buf = hdfinfo(fname_epw(i));
    epw_algn(i,:,:,:) = hdfread(buf.SDS);
end

buf = hdfinfo(fname_TPDI);
tpdi = hdfread(buf.SDS);

%% ##################### Axes and Calibration ######################### %%
cal.cwave = 526.5;
x = linspace(1,1024,1024);
y = linspace(1,1024,1024);
%% ################### Plot ##################### %%
n = 5;
xlim = [500,570];
ylim = [350,450];
%% #################### Plot EPW/IAW algn image ########################### %%
subplot(2,2,1)
box on
hold on
%pcolor(squeeze(iaw_algn(n,1,:,:)))
pcolor(x(xlim(1):xlim(2)),y(ylim(1):ylim(2)),squeeze(iaw_algn(n,1,ylim(1):ylim(2),xlim(1):xlim(2))))
plot(linspace(538,538,ylim(2)-ylim(1)+1),y(ylim(1):ylim(2)),"w:",LineWidth=2)    % input best eyeball value here
plot(x(xlim(1):xlim(2)),linspace(400,400,xlim(2)-xlim(1)+1),"w:",LineWidth=2)    % input best eyeball value here
shading flat
colormap turbo
xlabel('x-px')
ylabel('y-px')
fontsize(16,"points")
ax1 = gca;
ax1.LineWidth = 1;
% ax1.XTick = [1,256,512,768,1024];
% ax1.YTick = [1,256,512,768,1024];
ax1.TickDir = 'out';
title(fname_iaw(n))
ax1.XLim = [xlim(1),xlim(2)];
ax1.YLim = [ylim(1),ylim(2)];
hold off

%% #################### Plot lam lineouts ########################### %%
subplot(2,2,2)
box on
hold on
pos = linspace(xlim(1),xlim(2),10);

for j=1:length(pos)
    lo = squeeze(iaw_algn(n,1,:,ceil(pos(j))));
    plot(lo(ylim(1):ylim(2)),y(ylim(1):ylim(2)),'DisplayName','x-px = '+string(ceil(pos(j))),LineWidth=2)
end

legend 
ylabel('y-px')
xlabel('Amplitude(arb)')
title("x-lineouts")
ax2 = gca;
ax2.LineWidth = 1;
ax2.XLim = [ylim(1),ylim(2)];
%ax2.XTick = [1,256,512,768,1024];
fontsize(16,"points")
fontname("SansSerif")
axis("tight")

%% #################### Plot X lineouts ########################### %%
subplot(2,2,3)
box on
hold on
posx = linspace(ylim(1),ylim(2),10);

for j=1:length(posx)
    lo = squeeze(iaw_algn(n,1,ceil(posx(j)),:));
    plot(x(xlim(1):xlim(2)),-lo(xlim(1):xlim(2)),'DisplayName','y-px = '+string(ceil(posx(j))),LineWidth=2)
end

legend 
ylabel('Amplitude(arb)')
xlabel('x-px')
title("y-lineouts")
ax3 = gca;
ax3.LineWidth = 1;
ax3.XLim = [xlim(1),xlim(2)];
%ax3.XTick = [1,256,512,768,1024];
fontsize(16,"points")
fontname("SansSerif")
%axis("tight")
set(gcf,'position',[500,500,1000,800])
