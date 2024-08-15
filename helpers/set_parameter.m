function [PAR] = set_parameter(PAR,name,value)

PAR.value(strcmp(PAR.name,name)) = value;

end

