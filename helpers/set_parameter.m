function [PAR] = set_parameter(PAR,name,value,field)

if nargin==3
    field = 'value';
end

PAR.(field)(strcmp(PAR.name,name)) = value;

end

