close all
clc
clear
addpath '\\profiles\users$\sban\Documents\MATLAB\Lilac interpreters'
%% Import Lilac.hdf5 files

baseDir = ".";
choosing = false;
baseDir = fullfile(baseDir,"..","Kinshock","Kinshock-26A","Data");
while ~choosing
    fl = lss(baseDir);
    fl = fl(3:end);
    fprintf("Folders present:\n")
    for i = 1:length(fl)
        fprintf("\n %d. %s", i, fl(i))
    end

    which = input("\n Which one do you want to see?:\n",'s');
    which = round(str2num(which));
    if or(which < 0, which > length(fl))
        fprintf('***Error: %g is not an option.***\n', which);
    else
        baseDir = fullfile(baseDir,fl(which));
        choosing = ~isfile(baseDir);
    end
end

fl = lss(baseDir);
data = hdf


