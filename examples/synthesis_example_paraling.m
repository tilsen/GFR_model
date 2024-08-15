
clear; close all
[scriptdir] = fileparts(mfilename('fullpath'));
cd(scriptdir);

addpath(genpath(['..' filesep])) % run from directory that contains this script

initial_values = [105, 135, 150, 157.5, 161.25];
span_values = [20, 60, 80, 90, 95];

for i = 1:length(initial_values) % make 5 steps of paralinguistic variations

    gests = {{'H' 'L'}}; 
    data.t = {0:0.001:0.5};
    data.y = {[initial_values(i), initial_values(i)-0.01]}; % initial values change across conditions
    M = gen_structure(gests,data);
    PAR = gen_params(M,data);

    % specify parameter values (same for all conditions)
    PAR = set_parameter(PAR,'targ_01',1); PAR = set_parameter(PAR,'targ_02',0); 

    PAR = set_parameter(PAR,'ons_01',0); PAR = set_parameter(PAR,'ons_02',0.2); 
    PAR = set_parameter(PAR,'dur_01',0.4); PAR = set_parameter(PAR,'dur_02',0.3); 
    PAR = set_parameter(PAR,'ramp',0.1);
    
    % specify register variables (span differs across conditions)
    PAR = set_parameter(PAR,'floor_01',95); PAR = set_parameter(PAR,'decl',-10);
    PAR = set_parameter(PAR,'span_01',span_values(i)); 


    M = assign_params(M,PAR);

    [y,M] = f0mod(M.BETA.value,M);

    % plot
    h = plot_model(M,y(:,1));

end