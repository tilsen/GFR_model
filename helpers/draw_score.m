function [h_sc] = draw_score(M,ax,r,source,varargin)

p = inputParser;

addRequired(p,'M');
addRequired(p,'ax');
addRequired(p,'r');
addRequired(p,'source');
addParameter(p,'colors',[]);
addParameter(p,'tierMethod','auto',@(x)ismember(x,{'auto' 'gestureid'}));

parse(p,M,ax,r,source,varargin{:});

res = p.Results;

A = M.A;
G = M.G;
a = M.a;
t = M.t;

% default colors
if (isempty(res.colors) && ~ismember('color',A.Properties.VariableNames))
    gestures = unique(A.name);
    colors = lines(numel(gestures));
    for i=1:length(gestures)
        ixs = ismember(A.name,gestures{i});
        A.color(ixs,:) = repmat(colors(i,:),sum(ixs),1);
    end
elseif ~ismember('color',A.Properties.VariableNames)
    for i=1:height(A)
        A.color(i,:) = res.colors(mod(i-1,size(res.colors,1))+1,:);
    end
end

switch(source)
    case 'init'
        beta = M.BETA.b0;
    otherwise
        beta = M.BETA.value;
end
M.pw_t0 = M.pw_t0(:);
A.ons = beta(contains(M.BETA.name,'ons')) + M.pw_t0(A.pw);
A.dur = beta(contains(M.BETA.name,'dur'));
A.off = A.ons+A.dur;

da = @(c)abs(t-c)==min(abs(t-c));
A.onsix = arrayfun(@(c)find(da(c),1,'first'),A.ons);
A.offix = arrayfun(@(c)find(da(c),1,'first'),A.off);

%calculate y-positions
A = sortrows(A,'ons');
A.overlaps = A.offix > circshift(A.onsix,1);

switch(res.tierMethod)
    case 'auto'
        ypos = 0;
        for i=1:height(A)
            if A.overlaps(i) || ~ismember(A.name{i},A.name(1:i-1))
                ypos = ypos+1;
                A.ypos(i) = ypos;
            else
                A.ypos(i) = ypos;
            end
        end               
        A.ypos = max(A.ypos)-A.ypos;
    case 'gestureid'
        A.ypos = A.gestid-min(A.gestid);        
end

dy = 1;
dyy = dy-0.05;

for i=1:height(A)
    tt = t([A.onsix(i)*[1 1] A.offix(i)*[1 1]]);

    yo = A.ypos(i);
    yy = yo+[0 dyy dyy 0];
    h_sc.fh(i) = fill(tt,yy,A.color(i,:),...
        'FaceAlpha',0.5,'edgecolor','none','parent',ax); hold(ax,'on');

    tixs = A.onsix(i):A.offix(i);
    aa = yo + a(tixs,i)*dyy;
    h_sc.lh(i) = plot(t(tixs),aa,'-','color',A.color(i,:),'linew',2,'parent',ax);

    h_sc.th(i) = text(mean(t(tixs)),yo,G.name(A.gestid(i)),...
        'hori','center','verti','bot','parent',ax,'fontsize',r.fontsize,'fontname','calibri');
end

end
