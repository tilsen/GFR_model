function [cost] = f0mod_cost(b,M,y)
    yf = f0mod(b,M);
    if isfield(M,'weights')
        cost = nanmean((M.weights.*(yf(:,1)'-y).^2)); %#ok<NANMEAN> 
    else
        cost = nanmean((yf(:,1)'-y).^2); %#ok<NANMEAN> 
    end
end