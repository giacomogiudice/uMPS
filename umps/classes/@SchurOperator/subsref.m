function varargout = subsref(obj,s)
[varargout{1:nargout}] = builtin('subsref',obj.data,s);
end
