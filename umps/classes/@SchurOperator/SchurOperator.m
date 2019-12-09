classdef SchurOperator
    properties (Access = protected)
        data = {};
    end
    methods
        function obj = SchurOperator(varargin)
        	if (nargin == 1 && diff(size(varargin{1}))) || (nargin > 1 && all(cellfun(@iscell,varargin)))
                obj.data = buildFromTerms(varargin{:});
            elseif nargin
	            obj.data = cell(varargin{:});
	        end
        end
    end
end

function C = buildFromTerms(varargin)
% Preprocessing: extract scalars
coeff = ones(1,nargin);
for j = 1:nargin
    t = varargin{j};
    if isscalar(t{1})
        % The first term is a coefficient
        coeff(j) = t{1};
        varargin{j}(1) = [];
    end
end

% Get size of resulting Schur operator
chi = sum(cellfun(@(t) length(t) - 1,varargin)) + 2;
C = cell(chi,chi);
C{1,1} = 1;
C{chi,chi} = 1;
C{chi,1} = 0;

% Now populate C with the corresponding terms
ind = 1;
for j = 1:nargin
    t = varargin{j};
    if length(t) == 1
        % Single-body term
        O = t{1};
        assert(ismatrix(O) && ~diff(size(O)),'Expected square matrix inside term cell.');
        C{chi,1} = C{chi,1} + coeff(j)*O;
    else
        for k = 1:length(t)
            O = t{k};
            assert(ismatrix(O) && ~diff(size(O)),'Expected square matrix inside term cell.');
            if k == 1
                C{ind+1,1} = coeff(j)*O;
            elseif k == length(t)
                C{chi,ind+k-1} = O;
            else
                C{ind+k,ind+k-1} = O;
            end
        end
        ind = ind + (k-1);
    end
end
if C{chi,chi} == 0
    C{chi,chi} = [];
end
end

