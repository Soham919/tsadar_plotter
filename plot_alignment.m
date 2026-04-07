clc
close all
clear


%% ############################ READ HDF DATA ########################### %%
% fname.epw = ["EPW_526p5cw_1500um.hdf", "EPWCCD_0cw__50um_slit_500um_off_tcc.hdf", "EPWCCD_0cw__50um_slit_2500um_off_tcc.hdf", "EPWCCD_0cw__50um_slit_4500um_off_tcc.hdf"];
% fname.iaw = ["IAWCCD_0nmcw_1500toP4hdf.hdf","IAWCCD_532cw_500um_off_tcc.hdf", "IAWCCD_532cw_1300um_off_tcc.hdf", "IAWCCD_532cw_1500um_off_tcc.hdf", "IAWCCD_532cw_2500um_off_tcc.hdf", "IAWCCD_532cw_4500um_off_tcc.hdf"];
% fname_TPDI = "TPDI_1500um.hdf";
data_fp = '\\profiles\Users$\sban\Documents\Kinshock\Kinshock-26A\Data\Alignment\';
cd(data_fp)
fname.epw = ["EPWROSS_2w_532cw_100ph.hdf","EPWCCD_2w_0cw_50slit.hdf", "EPWCCD_2w_0cw_100ph.hdf", "EPWCCD_2w_526p5cw_50slit.hdf", "EPWCCD_2w_532cw_50slit.hdf", "EPWCCD_2w_532cw_100ph.hdf"];
fname.iaw = ["IAWROSS_2w_532cw_100ph.hdf","IAWCCD_2w_0cw_100ph.hdf","IAWCCD_2w_0cw_100slit.hdf","IAWCCD_2w_532cw_100ph.hdf", "IAWCCD_2w_532cw_100slit.hdf"];
%fname_TPDI = "TPDI_1500um.hdf";

epw_algn = zeros(length(fname.epw),2,1024,1024);
iaw_algn = zeros(length(fname.iaw),2,1024,1024);

for i = 1:length(fname.iaw)
    buf = hdfinfo(fname.iaw(i));
    iaw_algn(i,:,:,:) = hdfread(buf.SDS);
end

for i = 1:length(fname.epw)
    buf = hdfinfo(fname.epw(i));
    epw_algn(i,:,:,:) = hdfread(buf.SDS);
end

%buf = hdfinfo(fname_TPDI);
%tpdi = hdfread(buf.SDS);

%% ##################### Axes and Calibration ######################### %%
cal.cwave = 526.5;
x = linspace(1,1024,1024);
y = linspace(1,1024,1024);
%% ################### Plot ##################### %%
while true
    type = input("epw or iaw?\n", 's');
    if type == "epw"
        img = epw_algn;
    elseif type == "iaw"
        img = iaw_algn;
    else
        fprintf("You gotta type either epw or iaw(capsensitive)")
        continue
    end
    
    fprintf('Files available: \n');
    for j = 1:length(fname.(type))
        fprintf('  [%g] %s\n', j, fname.(type)(j));
    end
    
    n = input("Which one would you like to plot: \n");
    if n <= length(fname.(type))
        break
    end
end
%xlim = [500,570];
%ylim = [350,450];

bool1 = true;
bool2 = true;

figure;
fig = gcf;
set(fig,'Position',[100,100,600,600])
if type == "epw"
    img = rot90(squeeze(img(n,1,:,:)),-1);
    imagesc(img)
else 
    img = flip((rot90(squeeze(img(n,1,:,:)),-1)),1);
    imagesc(img)
end
ax = gca;
while or(bool1,bool2)
    flag = input("Change x-limits? : ","s");
    if or(isempty(flag), ismember(flag, ["false","False","n","N","no","No"]))
        bool1 = false;
        xlim = [ceil(ax.XLim(1)),floor(ax.XLim(end))];
    else
        temp = string(flag).split([",", " ", "[", "]", "(", ")"]);
        temp = str2double(temp(temp.strlength > 0));
        temp = temp(~isnan(temp));
        xlim = temp;
        
    end

    flag = input("Change y-limits? : ","s");
    if or(isempty(flag), ismember(flag, ["false","False","n","N","no","No"]))
        bool2 = false;
        ylim = [ceil(ax.YLim(1)),floor(ax.YLim(end))];
    else
        temp = string(flag).split([",", " ", "[", "]", "(", ")"]);
        temp = str2double(temp(temp.strlength > 0));
        temp = temp(~isnan(temp));
        ylim = temp;
    end
    % delete(fig)
    % figure;
    % fig = gcf;
    % ax = gca;
    % imagesc(squeeze(img(n,1,:,:)))
    ax.XLim = xlim;
    ax.YLim = ylim;
end

hold on
bool1 = true;
bool2 = true;
xval = ceil((xlim(2)-xlim(1))/2);
yval = ceil((ylim(2)-ylim(1))/2);
plot(linspace(xval,xval,ylim(2)-ylim(1)+1),y(ylim(1):ylim(2)),"w:",LineWidth=2)    % input best eyeball value here
plot(x(xlim(1):xlim(2)),linspace(yval,yval,xlim(2)-xlim(1)+1),"w:",LineWidth=2)    % input best eyeball value here
while or(bool1,bool2)
    newc = input(['\nAdjust x, y allignment? [', xval, ', ', yval, ']: '], 's');
    if or(isempty(newc), ismember(newc, ["false", "False", "n", "N", "no", "No"]))
        bool1 = false;
        bool2 = false;
    else
        temp = string(newc).split([",", " ", "[", "]", "(", ")"]);
        temp = str2double(temp(temp.strlength > 0));

        if or(sum(~isfinite(temp)), length(temp) ~= 2)
            warning('Input must be two numbers');
        else
            xval = temp(1);
            yval = temp(2);
        end
    end
    delete(fig)
    figure;
    fig = gcf;
    set(fig,'Position',[100,100,600,600])
    ax = gca;
    hold on
    imagesc(img)
    ax.XLim = xlim;
    ax.YLim = ylim;
    plot(linspace(xval,xval,ylim(2)-ylim(1)+1),y(ylim(1):ylim(2)),"w:",LineWidth=2)    % input best eyeball value here
    plot(x(xlim(1):xlim(2)),linspace(yval,yval,xlim(2)-xlim(1)+1),"w:",LineWidth=2)    % input best eyeball value here
end

%% #################### Plot EPW/IAW algn image ########################### %%
t = tiledlayout(2,2);
ax1 = nexttile;
box on
hold on
%pcolor(squeeze(iaw_algn(n,1,:,:)))
imagesc(ax1,x(xlim(1):xlim(2)),y(ylim(1):ylim(2)),(img(ylim(1):ylim(2),xlim(1):xlim(2))))
plot(ax1,linspace(xval,xval,ylim(2)-ylim(1)+1),y(ylim(1):ylim(2)),"w:",LineWidth=2)    % input best eyeball value here
plot(ax1,x(xlim(1):xlim(2)),linspace(yval,yval,xlim(2)-xlim(1)+1),"w:",LineWidth=2)    % input best eyeball value here
shading flat
colormap turbo
xlabel('x-px')
ylabel('y-px')
fontsize(16,"points")
%ax1 = gca;
ax1.LineWidth = 1;
% ax1.XTick = [1,256,512,768,1024];
% ax1.YTick = [1,256,512,768,1024];
ax1.TickDir = 'out';
title(fname.(type)(n))
ax1.XLim = [xlim(1),xlim(2)];
ax1.YLim = [ylim(1),ylim(2)];
hold off

%% #################### Plot lam lineouts ########################### %%
ax2 = nexttile;
box on
hold on
pos = linspace(xlim(1),xlim(2),10);

for j=1:length(pos)
    lo = img(:,ceil(pos(j)));
    plot(ax2,lo(ylim(1):ylim(2)),y(ylim(1):ylim(2)),'DisplayName','x-px = '+string(ceil(pos(j))),LineWidth=2)
end
plot(ax2,ax2.XLim(1):ax2.XLim(2),linspace(yval,yval,ax2.XLim(2)-ax2.XLim(1)+1),"b--",'DisplayName','y-px = '+string(yval),LineWidth=2)    % input best eyeball value here
legend 
ylabel('y-px')
xlabel('Amplitude(arb)')
title("x-lineouts")
%ax2 = gca;
ax2.LineWidth = 1;
ax2.XLim = [ylim(1),ylim(2)];
%ax2.XTick = [1,256,512,768,1024];
fontsize(16,"points")
fontname("SansSerif")
axis("tight")

%% #################### Plot X lineouts ########################### %%
ax3 = nexttile;
box on
hold on
posx = linspace(ylim(1),ylim(2),10);

for j=1:length(posx)
    lo = img(ceil(posx(j)),:);
    plot(ax3,x(xlim(1):xlim(2)),-lo(xlim(1):xlim(2)),'DisplayName','y-px = '+string(ceil(posx(j))),LineWidth=2)
end
plot(ax3,linspace(xval,xval,ax3.YLim(2)-ax3.YLim(1)+1),ax3.YLim(1):ax3.YLim(2),"b--",'DisplayName','x-px = '+string(xval),LineWidth=2)    % input best eyeball value here
legend 
ylabel('Amplitude(arb)')
xlabel('x-px')
title("y-lineouts")
%ax3 = gca;
ax3.LineWidth = 1;
ax3.XLim = [xlim(1),xlim(2)];
%ax3.XTick = [1,256,512,768,1024];
fontsize(16,"points")
fontname("SansSerif")
%axis("tight")
%figure('Position', [100, 100, 800, 600]);
