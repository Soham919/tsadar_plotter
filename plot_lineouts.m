clc 
clear
close all

fp = "/Users/soham/Documents/Quick Access/TSfit_Summer_2024/csv";
files = dir(fp);
files = files([files.isdir]);
disp("Shots present in folder:\n")
for k = 1:length(files)
    disp(files(k).name)
end

shot = input("Load lineouts for shots :\n",'s');
shots = strsplit(shot, {',',' '});

for i = 1:length(shots)
    fprintf("files for shot %s :",string(shots(i)))
    fp2 = fullfile(fp,shots(i));
    files = dir(fp2);
    for k = 1:length(files)
        fprintf("%s, ",files(k).name)
    end

    files2 = input("load ")
    
end
