
[scriptdir] = fileparts(mfilename('fullpath'));
cd(scriptdir);

addpath(genpath(['..' filesep])) %run from directory that contains this script

load(['..' filesep 'data' filesep 'data_optimization_example.mat'],'data');

gests = {{'H' 'L'},{'H' 'L'}};
M = gen_structure(gests,data); % generates gestures and gestural activation events tables
PAR = gen_params(M,data);

%impose different target bounds for H vs. L gestures:
hids = arrayfun(@(c){sprintf('%02d',c)},M.G.id(contains(M2.G.name,'H')));
lids = arrayfun(@(c){sprintf('%02d',c)},M.G.id(contains(M2.G.name,'L')));
hix = ismember(PAR.name,append('targ_',hids));
lix = contains(PAR.name,append('targ_',lids));
PAR.lb(hix) = 0.5;  %the lower bound of the H gesture
PAR.ub(lix) = 0.5;  %the upper bound of the L gesture

M = assign_params(M,PAR);
disp(M)

%gestures in the same pwrd cannot overlap more than this value:
M = ineq_constraints(M,'gest_overlap',0.100); 

%gestures cannot extend past pwrd by more than this value
M = ineq_constraints(M,'gest_endpw',0);

options = optimoptions('patternsearch', ...
    'Display','none', ...
    'InitialMeshSize',0.1, ...
    'UseParallel',true, ...
    'MaxIterations',1000);

cost_fcn = @(b)M.cost_fcn(b,M);

%run the optimization (will take a while):
[betahat,cost,exitflag,output] = patternsearch(cost_fcn,M.BETA.b0,M.CON.A,M.CON.b,[],[],M.BETA.lb,M.BETA.ub,options);

if exitflag<0
    disp(output.message);
    return;
end

M.BETA.value = betahat;
[y,M] = f0mod(betahat,M);

figure; plot_model(M,data.y{1},gcf);
figure; plot_params(M,gcf);

%save("../data/params_optimization_example.mat",'betahat','cost');