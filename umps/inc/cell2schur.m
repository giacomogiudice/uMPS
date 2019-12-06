function W = combineSchur(M)
% COMBINESCHUR Combines a cell array of Schur operators to a single operator.
%
%  W = combineSchur(M)
% The (N x 1) cell array M of Schur operators (cell array of (d x d)
% matrices or scalars) is combined into a single operator in Schur form.
% Each element of the resulting Schur operator is a (Nx1) cell array,
% corresponding to the operator on each site.

W = M{1};
for n = 2:length(M)
    W = combine_twoelements(W,M{n});
end
end

function Wout = combine_twoelements(W1,W2)
chi = size(W1,1);
assert(isequal(size(W1),size(W2),[chi,chi]),'Size mismatch between input operators.');
Wout = SchurOperator(chi,chi);
for a = 1:chi
    for c = 1:a
        for b = 1:c
            if ~isempty(W1{a,c}) && ~isempty(W2{c,b})
                O = {};
                O1 = W1{a,c};
                O2 = W2{c,b};
                if iscell(O1)
                    O = [O1,repelem({O2},size(O1,1),1)];
                elseif iscell(O2)
                    O = [repelem({O1},size(O2,1),1),O2];
                else
                    O = {O1,O2};
                end

                Wout{a,b} = [Wout{a,b};O];
            end
        end
    end
end

end
