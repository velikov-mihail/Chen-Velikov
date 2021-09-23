
clear
clc

restoredefaultpath; % Start with the default path
matlabPackagePath='D:/MATLAB AP Package for Chen and Velikov (JFQA)/'; % Path to the MATLAB asset pricing package
paperCodePath='D:/MATLAB AP Package for Chen and Velikov (JFQA)/Chen and Velikov'; % Path to the code for Chen and Velikov
inputsPath='D:/MATLAB AP Package for Chen and Velikov (JFQA)/Chen and Velikov/Inputs/'; % Path to the folder with inputs that should contain 

% Add the relevant folders (with subfolders) to the path
addpath(genpath([matlabPackagePath,'Data']))
addpath(genpath([matlabPackagePath,'Functions']))
addpath(genpath([matlabPackagePath,'Library Update']))
addpath(genpath([paperCodePath]))
addpath(genpath([inputsPath]))
cd(paperCodePath)

%% Add the directories

if ~exist([pwd,'Data'], 'dir')
    mkdir(['Data'])
end
if ~exist([pwd,'Results'], 'dir')
    mkdir(['Results'])
end

addpath(genpath(pwd));

%% Make trading costs measures

run('make_tcost_measure.m');

%% Run unmitigated strategies

run('run_unmitigated_strategies.m');

%% Run mitigated strategies

run('run_mitigated_strategies.m'); 

%% Run Novy-Marx and Velikov (RFS, 2016) reconciliation results

run('run_nmv_reconciliation_results.m'); 

%% Run combination strategy results

run('run_combination_strategies.m');

%% Organize results

run('organize_results.m');

%% Print tables

run('make_tables.m');

%% Make and store the figures

run('make_figures.m');

%% Additional code: replicate tcost figures from original papers

open replicate_tcost_figures

