function [varargout] = f0mod(beta,M)

% setup
Nf = M.Nf;         % number of field coordinates
f = M.f;           % normalized F0 field coordinates
dt = M.dt;         % timestep (=1/X.TV_Fs)
t = M.t;           % time points
Nt = M.Nt;         % number of time points
max_act = 1;
Nx = 3;            % number of vocal tract system variables = 3 (f0, flr, span)
pre_pw = -0.05;    % implement register change before unit onset (ToDo: parameterize)            

% combine free and fixed parameters:
B = M.BETA;
if size(beta,1) < size(beta,2) % if beta is a row vector
    B.value = beta';
else
    B.value = beta;
end

B = [B(:,{'name' 'value'}); M.FIXED(:,{'name' 'value'})];

% parameters
beta        = B.value'; %transpose so params are columns
Tg          = beta(contains(B.name,'targ_' + digitsPattern(1,2)));        % gestural targets
Tneut       = beta(contains(B.name,'targ_neut'));   % neutral attractor target
Fgain       = beta(strcmp(B.name,'gain'));          % field gain
Fgain_neut  = beta(strcmp(B.name,'gain_neut'));          % field gain
K_floor     = beta(strcmp(B.name,'stiff_floor'));       % floor stiffness
K_span      = beta(strcmp(B.name,'stiff_span'));       % span stiffness
ramp        = beta(contains(B.name,'ramp'));        % activation ramping
decl        = beta(contains(B.name,'decl'));        % declination rates
flr         = beta(contains(B.name,'floor_' + digitsPattern(1,2)));       % pwrd-initial floors
span        = beta(contains(B.name,'span_' + digitsPattern(1,2)));        % span
ons         = beta(contains(B.name,'ons_' + digitsPattern(1,2)));         % onset times (rel. to pword)
dur         = beta(contains(B.name,'dur_' + digitsPattern(1,2)));         % durations
sigma       = beta(contains(B.name,'sigma_' + digitsPattern(1,2)));
sigma_neut  = beta(contains(B.name,'sigma_neut'));

gaussfcn = @(mu)exp(-((mu(:)-f)./sigma(:)).^2); %for gestural systems
gaussfcn_sd = @(mu,sd)exp(-((mu(:)-f)/sd).^2);  %for neutral system

% convert to absolute onset times
if isscalar(M.pw_t0)
    ons = ons+M.pw_t0(M.A.pw)';
else
    ons = ons+M.pw_t0(M.A.pw);
end

% prevents error when ramp is 0
ramp = max(ramp,1e-5); 

% expand params
if isscalar(ramp), ramp = repmat(ramp,1,length(ons)); end

% gestural activation function:
actfcn = @(ons,dur,ramp_per)(t>=ons & t<=(ons+dur)).*...
    min((t-ons)      /ramp_per,max_act).*...
    min(((ons+dur)-t)/ramp_per,max_act); 

% gestural activation timeseries (row for each activation event):
A = cell2mat(arrayfun(@(ons,dur,ramp){actfcn(ons,dur,ramp)},ons,dur,ramp));

% gestural activation timeseries (row for each gesture):
G = cell2mat(cellfun(@(c){sum(A(:,c),2)},M.G.order'));

% initialize variables
x = zeros(Nt,Nx);       %positions
v = zeros(Nt,Nx);       %velocities
rflr = zeros(Nt,1);     %dynamic register floor
rspan = zeros(Nt,1);    %dynamic register span

% dynamic register changes
% 1) times of register changes relative to acoustic pwrd onsets
d_pw_ons = pre_pw*ones(1,length(M.pw_t0));
d_pw_ons(1) = 0;    

% 2) indices of register changes
if ~isfield(M,'register') || strcmp(M.register,'shared')
    reg_ons_ixs = 1;
else
    reg_ons_ixs = 1 + floor(M.pw_t0/dt) + floor(d_pw_ons/dt);
end
reg_ons = false(Nt,size(M.pw_ixs,2));
for i=1:length(reg_ons_ixs)
    reg_ons(reg_ons_ixs(i),i) = true;
end

% fixed neutral attractor amplitude, target, and sigma
Fneut = Fgain_neut*gaussfcn_sd(Tneut,sigma_neut);

% field activation function
Fact_fcn = @(gact)sum([gact(:).*gaussfcn(Tg); Fneut],1);

% field centroid function
Fcent_fcn = @(F)sum(F.*f,2)./sum(F,2);

% registerize (convert normalized F0 value to actual value)
reg_fcn = @(nT,fl,sp)fl + nT.*sp;

% de-registerize (convert actual F0 value to normalized tract variable value)
dereg_fcn = @(T,fl,sp)(T-fl)./sp;

% time series of normalized targets and stiffnesses
% (can be calculated outside of loop since field dynamics d/n depend on previous field state)
F = cell2mat(arrayfun(@(c){Fact_fcn(G(c,:))},(1:Nt)'));
nT = Fcent_fcn(F);

K = [Fgain*sum(F,2) K_floor*ones(Nt,1) K_span*ones(Nt,1)];
d = 2*sqrt(K);      %critical damping

% initial values
rflr(1) = flr(1);                                                   
rspan(1) = span(1);
%[target floor span]
x(1,:) = [dereg_fcn(M.x0,rflr(1),rspan(1))     rflr(1) rspan(1)];   % initial tract variable / register positions
v(1,:) = [M.v0/rspan(1)                          0      0];         % initial tract variable / register velocities

%%
for i=2:Nt
    
    %register floor/span targets:
    if any(reg_ons(i,:))
        rflr(i) = flr(reg_ons(i,:));
        rspan(i) = span(reg_ons(i,:));
    else
        dflr = decl;
        rflr(i) = rflr(i-1) + dflr*dt;
        rspan(i) = rspan(i-1);
    end
    
    % update f0 and register variable positions and velocities
    dx = v(i-1,:);
    dv = -d(i-1,:).*v(i-1,:) - K(i-1,:).*(x(i-1,:)-[nT(i-1,:) rflr(i-1) rspan(i-1)]);
    
    x(i,:) = x(i-1,:) + dx*dt;
    v(i,:) = v(i-1,:) + dv*dt;
 
end

%apply register
T = reg_fcn(nT,x(:,2),x(:,3));

%%

%scale v to Hz:
v(:,1) = v(:,1).*x(:,3);

%output position and velocity of f0
y = [reg_fcn(x(:,1),x(:,2),x(:,3)) v(:,1)];

varargout{1} = y;

if nargout==2

    M.x = x;
    M.v = v;
    M.rflr = rflr;
    M.rspan = rspan;
    M.F = F;
    M.nT = nT;
    M.T = T;
    M.a = A;
    M.g = G;
    M.reg_fcn = reg_fcn;
    M.dereg_fcn = dereg_fcn;
    M.reg_ons = reg_ons;
    M.Y = y;

    varargout{2}= M;

end
%confirm that velocity is change of position:
%vel_correct = all(abs(v(1:end-1,1)-(diff(x(:,1)).*x(1:end-1,3)/M.dt))<1e-3);

end


