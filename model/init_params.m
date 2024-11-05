function [PAR] = init_params(PAR,D,A,G,P,varargin)

dbstop if error;

p = inputParser;
p.KeepUnmatched = true;

def_register = 'shared';

def_targ_lb_method = 0;
def_targ_ub_method = 1;

def_ons_lb_method = 0;
def_ons_ub_method = 'maxdur';

def_dur_lb_method = 0.200;
def_dur_ub_method = 'maxdur';

def_ramp_lb_method = 0.02; % 0
def_ramp_ub_method = 0.1;  % 0.05;

def_gain_lb_method = 3; %0;
def_gain_ub_method = 200; % inf;

def_floor_lb_method = 'min';
def_floor_ub_method = 'initial'; 
def_span_lb_method = 'range'; 
def_span_ub_method = {'range-n',2}; % range*2;

def_decl_lb_method = 'range/time'; %-inf
def_decl_ub_method = 0;

def_targ_b0_method = 'H/L_max/min';
def_ons_b0_method = 'sequential';
def_dur_b0_method = 0.200;

def_ramp_b0_method = 'bounds_mid'; % 'ub'
def_gain_b0_method = 'bounds_mid'; % 100

def_floor_b0_method = 'lb';
def_span_b0_method = 'range'; 
def_decl_b0_method = 'ub';

addRequired(p,'PAR');
addRequired(p,'D');
addRequired(p,'A');
addRequired(p,'G');
addRequired(p,'P');
addParameter(p,'register',def_register);
addParameter(p,'targ_lb_method',def_targ_lb_method);
addParameter(p,'targ_ub_method',def_targ_ub_method);
addParameter(p,'ramp_lb_method',def_ramp_lb_method);
addParameter(p,'ramp_ub_method',def_ramp_ub_method);
addParameter(p,'gain_lb_method',def_gain_lb_method);
addParameter(p,'gain_ub_method',def_gain_ub_method);
addParameter(p,'ons_lb_method',def_ons_lb_method);
addParameter(p,'ons_ub_method',def_ons_ub_method);
addParameter(p,'dur_lb_method',def_dur_lb_method);
addParameter(p,'dur_ub_method',def_dur_ub_method);
addParameter(p,'floor_lb_method',def_floor_lb_method);
addParameter(p,'floor_ub_method',def_floor_ub_method);
addParameter(p,'span_lb_method',def_span_lb_method);
addParameter(p,'span_ub_method',def_span_ub_method);
addParameter(p,'decl_lb_method',def_decl_lb_method);
addParameter(p,'decl_ub_method',def_decl_ub_method);
addParameter(p,'targ_b0_method',def_targ_b0_method);
addParameter(p,'ramp_b0_method',def_ramp_b0_method);
addParameter(p,'gain_b0_method',def_gain_b0_method);
addParameter(p,'ons_b0_method',def_ons_b0_method);
addParameter(p,'dur_b0_method',def_dur_b0_method);
addParameter(p,'floor_b0_method',def_floor_b0_method);
addParameter(p,'span_b0_method',def_span_b0_method);
addParameter(p,'decl_b0_method',def_decl_b0_method);

parse(p,PAR,D,A,G,P,varargin{:});

r = p.Results;

if isstruct(D), D = struct2table(D); end

if ~ismember('pw_t0',D.Properties.VariableNames)
    D.pw_t0 = {D.t{1}(1)};
    D.pw_t1 = {D.t{1}(end)};
end

if ~iscell(D.t)
    D.t = {D.t};
end
if ~iscell(D.y)
    D.y = {D.y};
end

%add f0 trajectories by prosodic unit
for i=1:length(D.pw_t0{1})
    if i<length(D.pw_t0{1})
        D.y_pw{1,i} = D.y{1}(D.t{1}>=D.pw_t0{1}(i) & D.t{1}<D.pw_t1{1}(i));
    else
        D.y_pw{1,i} = D.y{1}(D.t{1}>=D.pw_t0{1}(i) & D.t{1}<=D.pw_t1{1}(i));
    end
end

inputs = {D,A,G,P};

%
pars = {'targ','ramp','gain','ons','dur',...
    'floor','span','decl'};

for i=1:length(pars)
    PAR = apply_bound(PAR,inputs,'lb',pars{i},r.([pars{i} '_lb_method']),r.register);
    PAR = apply_bound(PAR,inputs,'ub',pars{i},r.([pars{i} '_ub_method']),r.register);
    PAR = apply_init(PAR,inputs,pars{i},r.([pars{i} '_b0_method']),r.register);
end

end

%% initial values
function [PAR] = apply_init(PAR,inputs,name,method,register)

[get_values_inputs,ixs] = collect_inputs(method,name,inputs,PAR);

if isnumeric(method)
    PAR.b0(ixs) = method;
else    
    x = get_values(get_values_inputs,ixs);

    %handle case where multiple units are specified but register is shared
    if (strcmp(register,'shared') && ismember(name,{'floor','span'}))
        if iscell(method)
            meth = method{1};
            rangen = method{2};
        else 
            meth = method;
        end
        switch(meth)
            case 'range'
                x = max(inputs{1}.y{1})-min(inputs{1}.y{1});
            case 'range-n'
                x = rangen*max(inputs{1}.y{1})-min(inputs{1}.y{1});                
            case 'initial'
                x = x(1);
            case 'min'
                x = min(x);
            case 'max'
                x = max(x);
        end
    end
    
    PAR.b0(ixs) = x;    
end

end

%% bounds
function [PAR] = apply_bound(PAR,inputs,bound,name,method,register)

[get_values_inputs,ixs] = collect_inputs(method,name,inputs,PAR);

if isnumeric(method)
    PAR.(bound)(ixs) = method;
else    
    x = get_values(get_values_inputs,ixs);

    %handle case where multiple units are specified but register is shared
    if (strcmp(register,'shared') && ismember(name,{'floor','span'}))
        if iscell(method)
            meth = method{1};
            rangen = method{2};
        else 
            meth = method;
        end
        switch(meth)
            case 'range'
                x = max(inputs{1}.y{1})-min(inputs{1}.y{1});
            case 'range-n'
                x = rangen*max(inputs{1}.y{1})-min(inputs{1}.y{1});                
            case 'initial'
                x = x(1);
            case 'min'
                x = min(x);
            case 'max'
                x = max(x);
        end
    end

    PAR.(bound)(ixs) = x;    
end

end

%%
function [outputs,ixs] = collect_inputs(method,name,inputs,PAR)

method_pars = {};
if iscell(method)
    method_pars = method(2:end);
    method = method{1};
end

D = inputs{1};
A = inputs{2};
G = inputs{3};
P = inputs{4};
P.dur = P.t1-P.t0;

ixs = find(contains(PAR.name,string(name) + '_' + digitsPattern(1,2)));

if isempty(ixs)
    ixs = find(strcmp(PAR.name,string(name)));
end

%prosodic unit ids:
uIds = [];
gnames = {};

switch(PAR.replicates{ixs(1)})
    case 'gestures'
        gnames = G.name(ixs);

    case 'activationEvents'
        uIds = A.pw;

    case 'units'
        uIds = P.id;

    otherwise %utterance
        uIds = 1;
end

outputs = {method,D,PAR,A,G,P,uIds,gnames,method_pars};

end

%% parameter value initialization methods
function [x] = get_values(inputs,ixs)
method = inputs{1};
D = inputs{2};
PAR = inputs{3};
A = inputs{4};
G = inputs{5};
P = inputs{6};
uIds = inputs{7};
gnames = inputs{8};
method_pars = inputs{9};
switch(method)
    case 'min'
        x = cellfun(@(c)min(c),D.y_pw(1,uIds));
    case 'range'
        x = cellfun(@(c)range(c),D.y_pw(1,uIds));
    case 'range-n'
        x = method_pars{1}*cellfun(@(c)range(c),D.y_pw(1,uIds));
    case {'init','initial'}
        x = cellfun(@(c)c(1),D.y_pw(1,uIds));
    case 'end'
        x = cellfun(@(c)c(end),D.y_pw(1,uIds));
    case 'endt'
        x = P.t1(uIds);
    case 'startt'
        x = P.t0(uIds);
    case 'maxdur'
        x = P.dur(uIds);
    case 'range/time'
        x = cellfun(@(c)-range(c),D.y_pw(1,uIds))./P.dur(uIds);
    case {'H/L_max/min'}
        x = nan(1,numel(gnames));
        x(contains(gnames,'H')) = 1;
        x(contains(gnames,'L')) = 0;
    case 'bounds_mid'
        x = (PAR.ub(ixs)-PAR.lb(ixs))/2;
    case 'lb'
        x = PAR.lb(ixs);
    case 'ub'
        x = PAR.ub(ixs);
    case 'sequential'            
        B0 = [];
        for i=1:height(P)
            ixs_A = find(A.pw==i);            
            b0 = linspace(0,P.dur(i),numel(ixs_A)+1);
            B0 = [B0 b0(1:end-1)];
        end
        x = B0;
end
end
