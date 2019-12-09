classdef TwoSiteOperator < double
    methods
        function obj = TwoSiteOperator(varargin)
            if nargin && all(cellfun(@iscell,varargin))
                varargin = buildFromTerms(varargin{:});
            end
            obj@double(varargin);
        end
    end
end

function g = buildFromTerms(varargin)
coeff = 1;
g = 0;
for j = 1:nargin
    term = varargin{j};
    assert(iscell(term),'Expected cell containing single term.');
    if isscalar(term{1})
        % The first term is a coefficient
        coeff = term{1};
        term(1) = [];
    else
        coeff = 1;
    end
    if length(term) == 1
        % Single-body term
        O = term{1};
        assert(ismatrix(O) && ~diff(size(O)),'Expected square matrix inside term cell.');
        I = eye(size(O,1));
        g = g + coeff/2*ncon({O,I},{[-1,-3],[-2,-4]});
        g = g + coeff/2*ncon({I,O},{[-1,-3],[-2,-4]});
    elseif length(term) == 2;
        % Two-body term
        O1 = term{1};
        O2 = term{2};
        assert(ismatrix(O1) && ~diff(size(O1)),'Expected square matrix inside term cell.');
        assert(ismatrix(O2) && ~diff(size(O2)),'Expected square matrix inside term cell.');
        g = g + coeff*ncon({O1,O2},{[-1,-3],[-2,-4]});
    end
end
end
