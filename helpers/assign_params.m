function [M] = assign_params(M,PAR)

M.FIXED = PAR(PAR.fixed,:);
M.BETA = PAR(~PAR.fixed,:);

M.FIXED.value = M.FIXED.b0;
M.BETA.value(isnan(M.BETA.value)) = M.BETA.b0(isnan(M.BETA.value)); %uses initial values wherever value is not specified

bad_BETA = M.BETA(M.BETA.ub<M.BETA.lb,:);
for i=1:height(bad_BETA)
    fprintf('warning: non-sane parameter bounds: %s [%1.3f, %1.3f]\n',bad_BETA.name{i},bad_BETA.lb(i),bad_BETA.ub(i));
end

bad_BETA = M.BETA((M.BETA.b0>M.BETA.ub) | (M.BETA.b0<M.BETA.lb),:);
for i=1:height(bad_BETA)
    fprintf('warning: non-sane initial parameter value: %s [%1.3f] %1.3f [%1.3f]\n', ...
        bad_BETA.name{i},bad_BETA.lb(i),bad_BETA.b0(i),bad_BETA.ub(i));
end

end