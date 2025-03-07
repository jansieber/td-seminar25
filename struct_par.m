function ip=struct_par(varargin)
parnames=varargin;
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
end
