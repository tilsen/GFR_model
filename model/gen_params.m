function [PAR] = gen_params(M,D,varargin)

%{
    M: model structure
    D: data structure    
%}

dbstop if error;
p = inputParser;
p.KeepUnmatched = true;

def_ramping = 'shared';             % share ramping param across gests
def_register = 'byunit';            % separate register params by prosodic unit
def_declination = 'shared';         % value or shared or by prosodic unit

valid_ramping = {'none' 'shared' 'bygesture'};
valid_register = {'none' 'shared' 'byunit'};
valid_declination = {'shared' 'byunit'};

addRequired(p,'M',@(x)isstruct(x) & isfield(x,'A') & isfield(x,'G') & isfield(x,'P'));
addRequired(p,'D',@(x)isstruct(x));

addParameter(p,'ramping',def_ramping, @(x)ismember(x,valid_ramping));
addParameter(p,'register',def_register, @(x)ismember(x,valid_register));
addParameter(p,'declination',def_declination, @(x)isnumeric(x) | ismember(x,valid_declination));

parse(p,M,D,varargin{:});

A = M.A;
G = M.G;
P = M.P;

if ~iscell(D.y)
    D.y = {D.y};
end

empirical_trajectory = length(D.y{1})==length(M.t);

r = p.Results;

par = readtable('parameters.csv');

PAR=[];
for i=1:height(par)
    switch(par.replicates{i})
        case 'gestures'
            X = G;
            parx = add_pars([],...
                cellfun(@(c){sprintf('%s_%02d',par.parameterName{i},c)},num2cell(X.id)),...
                cellfun(@(c){sprintf('%s: %02d',par.description{i},c)},num2cell(X.id)));

        case 'activationEvents'
            X = A;
            parx = add_pars([],...
                cellfun(@(c){sprintf('%s_%02d',par.parameterName{i},c)},num2cell(X.id)),...
                cellfun(@(c){sprintf('%s: %02d',par.description{i},c)},num2cell(X.id)));

        case 'units'
            X = P;
            
            switch(par.parameterName{i})
                case {'floor','span','decl'}
                    switch(r.register)
                        case 'shared'
                            switch(par.parameterName{i})
                                case {'floor','span'}
                                    parx = add_pars([],sprintf('%s_%02d',par.parameterName{i},1),par.description{i});                            
                                otherwise
                                    parx = add_pars([],par.parameterName{i},par.description{i});                                                           
                            end

                        case 'byunit'
                            parx = add_pars([],...
                                cellfun(@(c){sprintf('%s_%02d',par.parameterName{i},c)},num2cell(X.id)),...
                                cellfun(@(c){sprintf('%s: %02d (Hz)',par.description{i},c)},num2cell(X.id)));
                    end

                otherwise
                    parx = add_pars([],...
                        cellfun(@(c){sprintf('%s_%02d',par.parameterName{i},c)},num2cell(X.id)),...
                        cellfun(@(c){sprintf('%s: %02d',par.description{i},c)},num2cell(X.id)));
            end

        case 'utterance'
            parx = add_pars([],par.parameterName{i},par.description{i});

        otherwise
            switch(par.parameterName{i})
                case 'ramp'
                    switch(r.ramping)
                        case 'none'

                        case 'shared'
                            parx = add_pars([],par.parameterName{i},par.description{i});

                        case 'bygesture'
                            parx = add_pars([],...
                                arrayfun(@(c){sprintf('ramp_%02i',c)},G.id),...
                                arrayfun(@(c){sprintf('ramping for gesture: %s',c)},G.name));
                    end  
                otherwise
                    parx = add_pars([],par.parameterName{i},par.description{i});
            end
    end

    parx = set_par_info_from_defaults(parx,par(i,:));
    PAR = [PAR; parx];
end

%determine parameter bounds/initial guesses if an empirical trajectory was provided
if empirical_trajectory
    PAR = init_params(PAR,D,A,G,P,varargin{:}); 
end

end

%%
function [P] = add_pars(PAR,names,descr,fixed,values)

if ~exist('values','var'), values = nan; end   %default is no value
if ~exist('fixed','var'), fixed = false; end  %default is varied parameter
if ~iscell(names), names = {names}; end
if ~iscell(descr), descr = {descr}; end

if numel(fixed)~=numel(names)
    fixed = repmat(fixed,numel(names),1);
end
if numel(values)~=numel(names)
    values = repmat(values,numel(names),1);
end

for i=1:length(names)
    P(i).name = names{i};
    P(i).fixed = logical(fixed(i));
    P(i).value = values(i);
    P(i).lb = nan;
    P(i).b0 = nan;
    P(i).ub = nan;
    P(i).descr = descr{i};      
end
P = struct2table(P);

if ~isempty(PAR)
    P = [PAR; P];
end

end

%%
function [parx] = set_par_info_from_defaults(parx,par)
parx.fixed = logical(par.deffixed(1)*ones(height(parx),1));
parx.lb = par.deflb(1)*ones(height(parx),1);
parx.b0 = par.defb0(1)*ones(height(parx),1);
parx.ub = par.defub(1)*ones(height(parx),1);
parx.replicates = repmat(par.replicates(1),height(parx),1);

% need to convert char to cell when there is single gesture
if size(parx,1)==1 & strcmp(parx.replicates,'gestures')
    parx.name = {parx.name};
    parx.descr = {parx.descr};
end
end


