function [BETA,FIXED] = f0mod_separate_params(PAR)

FIXED = PAR(PAR.fixed,:);
BETA = PAR(~PAR.fixed,:);

end