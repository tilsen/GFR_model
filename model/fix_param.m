function [PAR] = fix_param(PAR,name,value)

PAR.value(contains(PAR.name,name)) = value;
PAR.fixed(contains(PAR.name,name)) = true;

end

