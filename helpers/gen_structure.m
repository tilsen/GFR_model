function [M] = gen_structure(gests,data)

%convert each set of gests to cell
for i=1:length(gests)
    if ~iscell(gests{i}), gests{i} = gests(i); end
end

%parse gests input
A = [];
for i=1:length(gests)
    for j=1:length(gests{i})
        A(end+1).name = gests{i}{j};
        A(end).pw = i;
    end
end
A = struct2table(A);
A.order = (1:height(A))';
A.id = (1:height(A))';

%collect non-unique:
pws = cellfun(@(c){A.pw(ismember(A.name,c))'},unique(A.name,'stable'));
orders = cellfun(@(c){A.order(ismember(A.name,c))'},unique(A.name,'stable'));
G = table(unique(A.name,'stable'),pws,orders,'VariableNames',{'name' 'pw' 'order'});
G.id = (1:height(G))';


%gesture ids
for i=1:height(A)
    A.gestid(i) = find(cellfun(@(c)ismember(A.id(i),c),G.order));
end

%default color
G.color = repmat([.5 .5 .5],height(G),1);

%H: colors
G.color(contains(G.name,'H'),:) = repmat([.8 .8 .8],sum(contains(G.name,'H')),1);

A.color = G.color(A.gestid,:);

P.id = unique(A.pw);
P = struct2table(P);

cols = fieldnames(data);

if ~iscell(data.t)
    data.t = {data.t};
end

if ~iscell(data.y)
    data.y = {data.y};
end

%enforce first time=0
data.t{1} = data.t{1}-data.t{1}(1);

%prosodic units not specified
if any(~ismember('pw_t0',cols) | ~ismember('pw_t1',cols))
    fprintf('No prosodic units specified, assuming one unit\n.');
    if height(P)>1, fprintf('ERROR: mismatch in prosodic units between gestures and data\n'); return; end
    
    P.t0 = data.t{1}(1);
    P.t1 = data.t{1}(end);
else
    if height(P)~=length(data.pw_t0{1})
        fprintf('ERROR: mismatch in prosodic units between gestures and data\n'); return; end

    P.t0 = data.pw_t0{1}(:);
    P.t1 = data.pw_t1{1}(:);
end

M = init_model(data,P,G,A);

M.P = P;
M.G = G;
M.A = A;

end