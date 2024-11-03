function [h] = plot_model(M,varargin)

p = inputParser();

def_empirical = true;
def_initial = false; %plot initial f0 fit
def_optimized = true; %plot optimized fit
def_pword = true; 
def_score_opt = true; 
def_score_init = false;
def_register_opt = true;
def_register_init = false;
def_field_opt = true;
def_field_init = false;
def_score_panel_prop = 0.25; %proportion of height for gestural score
def_field_panel_prop = 0.25; %proportion of height for dynamic field
def_panel_order = {'f0','score_init','score_opt','field_init','field_opt'};
def_linewidth = 2;
def_linecolors = lines(3); %empirical, initial, final
def_fontsize = 16;

addRequired(p,'M',@(x)isstruct(x) & isfield(x,'BETA') & isfield(x,'model_fcn'));
addOptional(p,'y_emp',[],@(x)isempty(x) || (isvector(x) & isnumeric(x)));
addOptional(p,'parent',[],@(x)ishandle(x));

%data/model fits
addParameter(p,'empirical',def_empirical);
addParameter(p,'initial',def_initial);
addParameter(p,'optimized',def_optimized);
addParameter(p,'pword',def_pword);
addParameter(p,'register_opt',def_register_opt);
addParameter(p,'register_init',def_register_init);

%additional panels
addParameter(p,'score_opt',def_score_opt);
addParameter(p,'score_init',def_score_init);
addParameter(p,'field_opt',def_field_opt);
addParameter(p,'field_init',def_field_init);
addParameter(p,'score_panel_prop',def_score_panel_prop);
addParameter(p,'field_panel_prop',def_field_panel_prop);
addParameter(p,'panel_order',def_panel_order);

%graphics
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'linecolors',def_linecolors);
addParameter(p,'fontsize',def_fontsize);

parse(p,M,varargin{:});
r = p.Results;

if isempty(r.y_emp)
    r.empirical = false;
end

BETA = M.BETA;
modelfcn = M.model_fcn;

axmarg = [0.10 0.10 0.01 0.05];
if isempty(r.parent)
    h.init_ax = stf([1 1],axmarg);
    h.figh = h.init_ax.Parent;
else
    switch(r.parent.Type)
        case 'figure'
            h.figh = r.parent;
            h.init_ax = stf([1 1],axmarg,'parent',r.parent);
        case 'axes'
            h.figh = r.parent.Parent;
            h.init_ax = r.parent;
        otherwise
            fprintf('ERROR: parent must be a figure or axes handle\n'); return;
    end
end

%check for optimized parameters
if any(isnan(M.BETA.value))
    r.optimized = false;
    r.score_opt = false;
    r.field_opt = false;
    r.initial = true;
end

%get initial fit
if any([r.initial r.score_init r.field_init r.register_init])
    [y_init,M0] = modelfcn(BETA.b0,M);
end

if r.optimized && ~any(isnan(BETA.value))
    [y_opt,M] = modelfcn(BETA.value,M);
end

%% axes setup
%determine which panels to plot
panels = r.panel_order;

r.f0 = true; %always plot f0

keep_panels = false(size(panels));
for i=1:length(panels)
    if (~isfield(r,panels{i}))
        fprintf('Error: unknown panel type: %s\n',panels{i});
        return;
    end
    keep_panels(i) = r.(panels{i});
end
panels = panels(keep_panels);
npanels = numel(panels);

if npanels>1

    %calculate panel dimensions
    yp = [];
    r.f0 = true; %always plot f0
    h.ix_f0 = 1;
    ix = 1;
    for i=1:length(panels)
        if (r.(panels{i}))
            switch(panels{i})
                case 'score_init'
                    yp(i) = r.score_panel_prop;                    
                case 'score_opt'
                    yp(i) = r.score_panel_prop;                    
                case 'field_init'
                    yp(i) = r.field_panel_prop;                    
                case 'field_opt'
                    yp(i) = r.field_panel_prop;  
                case 'f0'
                    yp(i) = nan;
            end            
            h.(['ix_' panels{i}]) = ix; ix=ix+1;
        end
    end
    
    %give f0 panel the remainder (and a mininum of 25%):
    yp(1) = max(0.25*nansum(yp),1-nansum(yp)); %#ok<*NANSUM>

    %renormalize in case proportions were over-allocated
    yp = yp/sum(yp);

    panix = repelem((1:npanels)',ceil(yp*100),1);

    ax = stfig_subaxpos(h.init_ax,panix,[0 0 0 0 0 0.01]);
    h.ax = ax;
    
    for i=1:length(panels)
        h.(['ax_' panels{i}]) = ax(i);
    end    
    delete(h.init_ax);
else
    h.ax = h.init_ax;
    h.ax_f0 = h.init_ax;
    
end

%% text and graphics scaling

switch(get(h.figh,'Units'))
    case 'normalized'
        def_fig_dims = [1 1];
    case 'pixels'
        def_fig_dims = [1920 1080];
end

fig_dims = get(h.figh,"Position");
rescalexy = fig_dims(3:4)./def_fig_dims;
rescale = max(rescalexy);

r.fontsize = ceil(r.fontsize*rescale);

%%
% empirical data
if r.empirical && ~isempty(r.y_emp)  
    h.y_emp = plot(M.t,r.y_emp,'parent',h.ax_f0,...
        'linewidth',r.linewidth,'color',r.linecolors(1,:)); 
    hold(h.ax_f0,'on');
    axis(h.ax_f0,'tight');
end

%initial fit
if r.initial
    h.y_init = plot(M.t,y_init(:,1),'parent',h.ax_f0,...
        'linewidth',r.linewidth,'color',r.linecolors(2,:)); 
    hold(h.ax_f0,'on');
    axis(h.ax_f0,'tight');
end

%optimized fit
if r.optimized
    h.y_opt = plot(M.t,y_opt(:,1),'parent',h.ax_f0,...
        'linewidth',r.linewidth,'color',r.linecolors(3,:)); 
    hold(h.ax_f0,'on');
    axis(h.ax_f0,'tight');
end

%initial register
if r.register_init && exist('M0','var')
    h.reg_init(1) = plot(M0.t,M0.x(:,2),'parent',h.ax_f0,...
        'linestyle','--','linewidth',r.linewidth,'color',r.linecolors(2,:)); 
    h.reg_init(2) = plot(M0.t,sum(M0.x(:,[2 3]),2),'parent',h.ax_f0,...
        'linestyle','--','linewidth',r.linewidth,'color',r.linecolors(2,:));       
end

%optimized register
if r.optimized && r.register_opt
    h.reg_opt(1) = plot(M.t,M.x(:,2),'parent',h.ax_f0,...
        'linestyle','--','linewidth',r.linewidth,'color',r.linecolors(3,:)); 
    h.reg_opt(2) = plot(M.t,sum(M.x(:,[2 3]),2),'parent',h.ax_f0,...
        'linestyle','--','linewidth',r.linewidth,'color',r.linecolors(3,:));      
end

%pword
if r.pword && isfield(M,'pw_t0')
    for i=1:length(M.pw_t0)
        h.pwrd(i,1) = plot(M.pw_t0(i)*[1 1],h.ax_f0.YLim,'k:','parent',h.ax_f0);
        %h.pwrd(i,2) = plot(M.pw_t1(i)*[1 1],h.ax_f0.YLim,'k:','parent',h.ax_f0);
    end
end

%gestural scores
if r.score_init 
    score_init = draw_score(M0,h.ax_score_init,r,'init');
    h.score_init = score_init;
end

if r.optimized && r.score_opt 
    score_opt = draw_score(M,h.ax_score_opt,r,'opt');
    h.score_opt = score_opt;
end

%target planning fields
if r.field_init 
    h.field_init = imagesc(M0.t,M0.f,M0.F','parent',h.ax_field_init);
    hold(h.ax_field_init,'on');
    h.field_centroid_init = plot(M0.t,M0.nT,'Parent',h.ax_field_init,'color','g','LineWidth',2);
    
    colormap(h.ax_field_init,hot(256));
    set(h.ax_field_init,'YDir','normal');
end

if r.field_opt 
    h.field_opt = imagesc(M.t,M.f,M.F','parent',h.ax_field_opt);
    hold(h.ax_field_opt,'on');
    h.field_centroid_opt = plot(M.t,M.nT,'Parent',h.ax_field_opt,'color','g','LineWidth',2);        
    colormap(h.ax_field_opt,hot(256));
    set(h.ax_field_opt,'YDir','normal');
end

%----------------
set(h.ax,'TickDir','out','Ticklen',0.003*[1 1],'Box','on','fontsize',14,'XLim',minmax(M.t(:)'));
set(h.ax_f0,'Box','off','YGrid','on');

set(h.ax(1:end-1),'XTickLabel',[]);
set(h.ax(2:end-1),'XTick',[]);
set(h.ax(2:end),'YTick',[]);

h.ax_f0.YLim = h.ax_f0.YLim + 0.01*diff(h.ax_f0.YLim)*[-1 1];
if r.score_init
    h.ax_score_init.YLim = h.ax_score_init.YLim + 0.05*diff(h.ax_score_init.YLim)*[-1 1];
end
if r.score_opt
    h.ax_score_opt.YLim = h.ax_score_opt.YLim + 0.05*diff(h.ax_score_opt.YLim)*[-1 1];
end

score_opt_label = 'optim.';
field_opt_label = 'optim.';
if (~r.score_init)
    score_opt_label = 'gestures';
end
if (~r.field_init)
    field_opt_label = 'field';
end

h.ylab(1) = ylabel(h.ax_f0,'F0 (Hz)','fontsize',r.fontsize);
if r.score_init, h.ylab(end+1) = ylabel(h.ax_score_init,'init.','fontsize',r.fontsize); end
if r.score_opt,  h.ylab(end+1) = ylabel(h.ax_score_opt,score_opt_label,'fontsize',r.fontsize); end
if r.field_init, h.ylab(end+1) = ylabel(h.ax_field_init,'init.','fontsize',r.fontsize); end
if r.field_opt,  h.ylab(end+1) = ylabel(h.ax_field_opt,field_opt_label,'fontsize',r.fontsize); end

h.xlab(1) = xlabel(h.ax(end),'time (s)','fontsize',def_fontsize); 

drawnow;

end
