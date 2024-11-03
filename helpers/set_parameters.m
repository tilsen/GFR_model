function [PAR] = set_parameters(PAR,name,varargin)

for i=1:2:length(varargin)
    PAR.(varargin{i})(strcmp(PAR.name,name)) = varargin{i+1};
end

end

