function [out] = get_parameter(PAR,name,field)

if nargin==2
    field = 'value';
end

out = PAR.(field)(strcmp(PAR.name,name));

end

