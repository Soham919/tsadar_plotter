clc
close all
clear

fp = "\\profiles\Users$\sban\Documents\Kinshock\Kinshock-26A\Gasjet\";
file = "D-GJ-C-232_H2He_1500psi.h5";
fp = fullfile(fp,file);
info = h5info(fp);

dat = h5read(fp,'/pressure');

imagesc(dat)
colormap('jet')
colorbar