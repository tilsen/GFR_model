function [M] = add_cost_weights(M,weights)

if length(weights)~=length(M.f0)
    fprintf('ERROR: weights must have same length as data')
    return
end

if ~all(size(weights)==size(M.f0))
    weights = weights';
end

M.weights = weights;

end