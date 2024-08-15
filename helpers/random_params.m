function [b] = random_params(B)

b = B.lb + rand(height(B),1).*(B.ub-B.lb);

end