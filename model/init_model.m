function [M] = init_model(D,P,G,A)

%ToDo: add input parser allowing for control of the fixed parameters below

%to avoid numeric precision issues in indexing time frames
floorround = @(x)floor(1000*x)/1000;

% add time vector to model;
M.t = floorround(D.t{1}(:));

% add prosodic word timepoint indices

for i=1:height(P)
    M.pw_ixs(M.t>=floorround(P.t0(i)) & M.t<floorround(P.t1(i)),i) = true;    
end
M.pw_ixs(find(M.pw_ixs(:,end),1,'last')+1,end) = true;

% simulation parameters
M.dt = M.t(2)-M.t(1);
M.Nf = 100; 
M.f0 = zeros(1,M.Nf);           % number of field coordinates   
M.f = linspace(0,1,M.Nf);       % normalized F0 field coordinate values
M.Nt = length(M.t);             % number of time points

for i=1:size(M.pw_ixs,2)
    M.pw_t0(i) = M.t(find(M.pw_ixs(:,i),1,'first'));
    M.pw_t1(i) = M.t(find(M.pw_ixs(:,i),1,'last'));
end

M.A = A;
M.G = G;

if height(D)>1
    D = D(1,:);
end

if ~iscell(D.y)
    D.y = {D.y};
end

%initial values
M.x0 = D.y{1}(1);

if length(D.y{1})>1
    M.v0 = diff(D.y{1}(1:2))/M.dt;
else
    M.v0 = 0;
end

M.model_fcn = @(b,m)f0mod(b,m);  %need to add model function handle
M.cost_fcn = @(b,m)f0mod_cost(b,m,D.y{1});

end

