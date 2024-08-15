function [b] = random_perturb_params(B,b0,perturb_mag)

if nargin==2
    perturb_mag = 1;
end

b = b0 + perturb_mag*(rand(height(B),1)-0.5).*(B.ub-B.lb);
b = max(b,B.lb);
b = min(b,B.ub);

end