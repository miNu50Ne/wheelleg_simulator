clear;

if ~exist(fullfile(pwd, 'function'), 'dir')
    mkdir('function');
end

addpath("function");

run("model_sjtu.m");
run("leg_calculator.m");
run("five_link_calculator.m");

