function obj = subsasgn(obj,s,varargin)
obj.data = builtin('subsasgn',obj.data,s,varargin{:});
end
