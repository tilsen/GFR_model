function [cost] = f0mod_cost(b,M,y)
    yf = f0mod(b,M);
    cost = nanmean((yf(:,1)'-y).^2); %#ok<NANMEAN> 
end