function [] = plot_params(M,figh)

%{
ToDo: 
plot contraints?
group parameters by common bounds?
implement input parser and parameterize inputs
add gesture names when relevant?
%}

colors = [0 0 0; lines(2)];

color_bound = colors(1,:);
color_init = colors(2,:);
color_opt = colors(3,:);
color_pinned = [1 0 0];

ms_init = 16;
ms_opt = 12;
marker_init = 'o';
marker_opt = 'o';

fontsize_names = 16;
fontsize_axes = 14;

B = [M.BETA; M.FIXED];

B.lb(isnan(B.lb)) = -inf;
B.ub(isnan(B.ub)) = inf;

%xlims
B.xlim = [B.lb B.ub];

lbinf = isinf(B.lb);
B.xlim(lbinf,1) = B.value(lbinf)-max(1,2*abs(B.value(lbinf) - B.ub(lbinf)));

ubinf = isinf(B.ub);
B.xlim(ubinf,2) = B.value(ubinf)+max(1,2*abs(B.value(ubinf) - B.lb(ubinf)));

B.xlim(B.fixed,:) = B.value(B.fixed) + [-1 1];

B.xlim(isnan(B.xlim(:,1)),1) = B.xlim(isnan(B.xlim(:,1)),2)-1;
B.xlim(isnan(B.xlim(:,2)),2) = B.xlim(isnan(B.xlim(:,2)),1)+1;

B.xlim = B.xlim + 0.05*[-1 1].*diff(B.xlim,[],2);

B.xlim(diff(B.xlim,[],2)==0,:) = B.xlim(diff(B.xlim,[],2)==0,:) + [0 1e12];

tol = 1e-3;

B.pinned = abs(B.value-[B.lb B.ub])<tol;

%% axes setup
ncol = 3;
if height(B)<20
    ncol = 2;
end
nrow = ceil(height(B)/ncol);
axpan = reshape(1:(nrow*ncol),[],ncol);

if nargin==2
    ax = stf(axpan,[0.025 0.05 0.025 0.01],[0.05 0.05],'parent',figh);
else
    ax = stf(axpan,[0.025 0.05 0.025 0.01],[0.05 0.05]);
    figh = ax(1).Parent;
end

%% text and graphics scaling

def_fig_dims = [1920 1080];
fig_dims = get(figh,"Position");

if ~strcmp(figh.Units,'Pixels') && max(fig_dims)<=1
    orig_units = figh.Units;
    figh.Units = 'pixels';
    fig_dims = figh.Position;
    figh.Units = orig_units;
end

rescalexy = fig_dims(3:4)./def_fig_dims;
rescale = max(rescalexy);

ms_opt = round(ms_opt*rescale);
ms_init = round(ms_init*rescale);
fontsize_names = round(fontsize_names*rescale);
fontsize_axes = round(fontsize_axes*rescale);

%%
yo = 0.5;
yy = [0 1];
for i=1:height(B)
    axes(ax(i));
    
    hold(ax(i),'on');
    set(ax(i),'Ylim',yy,'XLim',B.xlim(i,:));

    col = color_bound; lw = 2;
    if B.pinned(i,1), col = color_pinned; lw = 4; end
    switch(isfinite(B.lb(i)))
        case true
            h.lb(i) = plot(B.lb(i)*[1 1],yy,'color',col,'linew',lw);
        case false

    end

    col = color_bound; lw = 2;
    if B.pinned(i,2), col = color_pinned; lw = 4; end    
    switch(isfinite(B.ub(i)))
        case true
            h.ub(i) = plot(B.ub(i)*[1 1],yy,'color',col,'linew',lw);
        case false

    end    
    
    h.init(i) = plot(B.b0(i),yo,'k',...
        'marker',marker_init,'markerfacecolor',color_init,'markersize',ms_init);
    h.opt(i) = plot(B.value(i),yo,'k',...
        'marker',marker_opt,'markerfacecolor',color_opt,'markersize',ms_opt);

    h.th(i) = text(mean(xlim),1,B.name{i},...
        'verti','top','hori','center','fontsize',fontsize_names,'interpreter','none');

    switch(B.fixed(i))
        case true
            set(gca,'Box','on','XTick',B.value(i));           

        case false
            set(gca,'Box','off','Ycolor','none');
    end

end

%-------
set(ax,'fontsize',fontsize_axes,'Ytick',[],'XGrid','on');

delete(ax(height(B)+1:end));
drawnow;

end
