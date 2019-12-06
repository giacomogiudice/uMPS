classdef SchurOperator
    properties (Access = protected)
        data = {};
    end
    methods
        function obj = SchurOperator(varargin)
        	if ~nargin
        		obj.data = {};
        	else
	            obj.data = cell(varargin{:});
	        end
        end
    end
end

