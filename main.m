clear;

if ~exist(fullfile(pwd, 'function'), 'dir')
    mkdir('function');
end

if ~exist(fullfile(pwd, 'data'), 'dir')
    mkdir('data');
end

addpath("function;data;leg_calculator;model;simulator");

if ~exist(fullfile(pwd, 'data', 'model_sjtu.mat'), 'file')
    disp('run model_sjtu.m');
    run("model/model_sjtu.m");
else
    load('model_sjtu.mat');
end

% if ~exist(fullfile(pwd, 'data', 'leg_calculator.mat'), 'file')
%     disp('run leg_calculator.m');
%     run("leg_calculator/leg_calculator.m");
% else
%     load('leg_calculator.mat');
% end

if ~exist(fullfile(pwd, 'data', 'five_link.mat'), 'file')
    disp('run five_link_calculator.m');
    run("leg_calculator/five_link_calculator.m");
else
    load('five_link.mat');
end
