
close all
[scriptdir] = fileparts(mfilename('fullpath'));
cd(scriptdir);

addpath(genpath(['..' filesep])) % run from directory that contains this script

tunes = {{'H' 'H' 'H'}, {'H' 'H' 'L'}, {'H' 'L' 'H'}, {'H' 'L' 'L'}, ...
        {'L' 'H' 'H'}, {'L' 'H' 'L'}, {'L' 'L' 'H'}, {'L' 'L' 'L'}};

for i = 1:8 % for eight tunes
    
    gests = {tunes{i}};
    data.t = {0:0.001:0.45};
    if (i == 2) | (i == 6)
        data.t = {0:0.001:0.37}; % truncate final gesture for HHL and LHL to model upstep
    end
    data.y = {[105, 104.99]};
    M = gen_structure(gests,data);
    PAR = gen_params(M,data);

    % specify parameter values (same for all tunes)
    PAR = set_parameter(PAR,'ons_01',0); PAR = set_parameter(PAR,'ons_02',0.2); PAR = set_parameter(PAR,'ons_03',0.3);
    PAR = set_parameter(PAR,'dur_01',0.25); PAR = set_parameter(PAR,'dur_02',0.15); PAR = set_parameter(PAR,'dur_03',0.15);
    PAR = set_parameter(PAR,'gain',30); 
    PAR = set_parameter(PAR,'floor_01',95); PAR = set_parameter(PAR,'span_01',35); PAR = set_parameter(PAR,'decl',-10);

    % specify gestural targets (H: 1, L: 0)
    if i <= 4
        PAR = set_parameter(PAR,'targ_01',1); PAR = set_parameter(PAR,'targ_02',0); 
    else
        PAR = set_parameter(PAR,'targ_01',0); PAR = set_parameter(PAR,'targ_02',1); 
    end

    M = assign_params(M,PAR);

    [y,M] = f0mod(M.BETA.value,M);

    
    % plot
    h = plot_model(M,y(:,1));
    axes(h.ax_f0); title(['(' num2str(i) ') ' [tunes{i}{1} tunes{i}{2} tunes{i}{3}]],'Fontsize',30);
   
end

