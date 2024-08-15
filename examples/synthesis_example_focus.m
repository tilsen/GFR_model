
close all
[scriptdir] = fileparts(mfilename('fullpath'));
cd(scriptdir);

addpath(genpath(['..' filesep])) % run from directory that contains this script

title_string = {'neutral','initial focus','final focus'};

for i = 1:length(title_string)

    gests = {{'H' 'L'}, {'H' 'L'}, {'H' 'L'}};  % assume three prosodic units, each with H and L gestures
    data.t = {0:0.001:0.5};
    data.pw_t0 = {[0 0.2 0.3]}; data.pw_t1 = {[0.2 0.3 0.5]};
    data.y = {[140, 139.99]};
    M = gen_structure(gests,data);
    PAR = gen_params(M,data);

    % specify parameter values (same for all focus conditions)
    PAR = set_parameter(PAR,'targ_01',1); PAR = set_parameter(PAR,'targ_02',0); 

    PAR = set_parameter(PAR,'ons_01',0); PAR = set_parameter(PAR,'ons_02',0.075); 
    PAR = set_parameter(PAR,'ons_03',0); PAR = set_parameter(PAR,'ons_04',0.025);
    PAR = set_parameter(PAR,'ons_05',0); PAR = set_parameter(PAR,'ons_06',0.075);

    PAR = set_parameter(PAR,'dur_01',0.125); PAR = set_parameter(PAR,'dur_02',0.125); 
    PAR = set_parameter(PAR,'dur_03',0.075); PAR = set_parameter(PAR,'dur_04',0.075);
    PAR = set_parameter(PAR,'dur_05',0.125); PAR = set_parameter(PAR,'dur_06',0.125);

    % specify register variables (span differs across conditions)
    PAR = set_parameter(PAR,'decl',-5);
    PAR = set_parameter(PAR,'floor_01',130); PAR = set_parameter(PAR,'floor_02',120); PAR = set_parameter(PAR,'floor_03',110); 
    if i == 1 
        PAR = set_parameter(PAR,'span_01',30); PAR = set_parameter(PAR,'span_02',30); PAR = set_parameter(PAR,'span_03',30); 
    elseif i == 2
        PAR = set_parameter(PAR,'span_01',40); PAR = set_parameter(PAR,'span_02',15); PAR = set_parameter(PAR,'span_03',15); 
    else
        PAR = set_parameter(PAR,'span_01',30); PAR = set_parameter(PAR,'span_02',30); PAR = set_parameter(PAR,'span_03',40); 
    end

    M.register = 'byunit';
    M = assign_params(M,PAR);

    [y,M] = f0mod(M.BETA.value,M);

    % plot
    h = plot_model(M,y(:,1));
    axes(h.ax_f0); title(['(' num2str(i) ') ' title_string{i}],'Fontsize',30);

end