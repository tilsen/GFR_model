function [value] = get_parameter(PAR,name)

value = PAR.value(strcmp(PAR.name,name));

end

