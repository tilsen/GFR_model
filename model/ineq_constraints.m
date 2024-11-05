function [M] = ineq_constraints(M,constraint,value,varargin)

p = inputParser();
p.KeepUnmatched = true;

def_value = 0;

addRequired(p,'M');
addRequired(p,'constraint');
addOptional(p,'value',def_value);

parse(p,M,constraint,value);
r = p.Results;

if ~isfield(M,'CON')
    M.CON.A = zeros(0,height(M.BETA));
    M.CON.b = zeros(0,1);
    M.CON.constraints = {};
end

%%
switch(constraint)
    case {'gest_overlap','gestural_overlap'}
        M.A = sortrows(M.A,{'pw','order'}); 
        pw = unique(M.A.pw);
        for i=1:length(pw)
            gixs = find(M.A.pw==pw(i));
            if numel(gixs)>1
                for j=1:length(gixs)-1
                    ix0_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j)))));
                    ix0_dur = find(ismember(M.BETA.name,sprintf('dur_%02i',M.A.id(gixs(j)))));
                    ix1_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j+1)))));
                    M.CON.A(end+1,[ix0_ons ix0_dur ix1_ons]) = [1 1 -1];
                    M.CON.b(end+1,1) = r.value;
                    M.CON.constraints{end+1} = ...
                        {constraint,M.BETA.name{ix0_ons},'+',M.BETA.name{ix0_dur},'-',M.BETA.name{ix1_ons},'<',M.CON.b(end)};
                    test_constraint(M.CON,M.BETA.b0);
                end
            end
        end
    case {'gest_endoverlap','gestural_endoverlap'}
        M.A = sortrows(M.A,{'pw','order'}); 
        pw = unique(M.A.pw);
        for i=1:length(pw)
            gixs = find(M.A.pw==pw(i));
            if numel(gixs)>1
                for j=1:length(gixs)-1
                    ix0_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j)))));
                    ix0_dur = find(ismember(M.BETA.name,sprintf('dur_%02i',M.A.id(gixs(j)))));
                    ix1_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j+1)))));
                    ix1_dur = find(ismember(M.BETA.name,sprintf('dur_%02i',M.A.id(gixs(j+1)))));
                    M.CON.A(end+1,[ix0_ons ix0_dur ix1_ons ix1_dur]) = [1 1 -1 -1];
                    M.CON.b(end+1,1) = r.value;
                    M.CON.constraints{end+1} = ...
                        {constraint,M.BETA.name{ix0_ons},'+',...
                        M.BETA.name{ix0_dur},'-',...
                        M.BETA.name{ix1_ons},'-',...
                        M.BETA.name{ix1_dur},'<',M.CON.b(end)};
                    test_constraint(M.CON,M.BETA.b0);
                end
            end
        end        
    case {'gest_underlap','gestural_underlap'}
        M.A = sortrows(M.A,{'pw','order'}); 
        pw = unique(M.A.pw);
        for i=1:length(pw)
            gixs = find(M.A.pw==pw(i));
            if numel(gixs)>1
                for j=1:length(gixs)-1
                    ix0_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j)))));
                    ix0_dur = find(ismember(M.BETA.name,sprintf('dur_%02i',M.A.id(gixs(j)))));
                    ix1_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j+1)))));
                    M.CON.A(end+1,[ix0_ons ix0_dur ix1_ons]) = [-1 -1 1];
                    M.CON.b(end+1,1) = r.value;
                    M.CON.constraints{end+1} = ...
                        {constraint,M.BETA.name{ix1_ons},'-',M.BETA.name{ix0_dur},'-',M.BETA.name{ix0_ons},'<',M.CON.b(end)};
                    test_constraint(M.CON,M.BETA.b0);
                end
            end
        end        
    case {'gest_endpw'}
        M.A = sortrows(M.A,{'pw','order'}); 
        pw = unique(M.A.pw);
        for i=1:length(pw)
            gixs = find(M.A.pw==pw(i));
            for j=1:length(gixs)
                ix_ons = find(ismember(M.BETA.name,sprintf('ons_%02i',M.A.id(gixs(j)))));
                ix_dur = find(ismember(M.BETA.name,sprintf('dur_%02i',M.A.id(gixs(j)))));
                M.CON.A(end+1,[ix_ons ix_dur]) = [1 1];
                M.CON.b(end+1,1) = r.value+(M.pw_t1(i)-M.pw_t0(i));
                M.CON.constraints{end+1} = ...
                    {constraint,M.BETA.name{ix_ons},'+',M.BETA.name{ix_dur},'<',M.CON.b(end)};
                test_constraint(M.CON,M.BETA.b0);
            end
        end       
    case {'floor_end'}
        pw = varargin{1};
        for i=1:length(pw)

            floorpar = sprintf('floor_%02i',pw(i));
            declpar = 'decl'; %ToDo: handle pwrd-specific declination
            ix_floor = find(ismember(M.BETA.name,floorpar));
            ix_decl = find(ismember(M.BETA.name,declpar));
            pwrd_dur = M.P.t1(pw(i))-M.P.t0(pw(i));

            M.CON.A(end+1,[ix_floor ix_decl]) = [1 pwrd_dur];
            M.CON.b(end+1,1) = r.value;
            M.CON.constraints{end+1} = ...
                {constraint,M.BETA.name{ix_floor},'+',[num2str(pwrd_dur) '*' M.BETA.name{ix_decl}],'<',M.CON.b(end)};
            test_constraint(M.CON,M.BETA.b0);

        end

end

end

%%
function [] = test_constraint(CON,b0)

A = CON.A(end,:);
b = CON.b(end,:);
constr = CON.constraints{end};

if ~(A*b0 < b)
    fprintf('init param violates constraint: ');
    for j=1:length(constr)
        if isnumeric(constr{j})
            fprintf('%1.3f\t',constr{j});
        else
            fprintf('%s\t',constr{j});
        end
        
    end
    fprintf('\n');
end

end