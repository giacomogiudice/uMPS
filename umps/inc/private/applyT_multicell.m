function M = applyT_multicell(M,A_1,O,A_2,direction)
N = length(A_1);
if ~iscell(O)
    O = repmat({O},1,N);
end
if direction == 'l'
    for n = 1:N
        M = applyT(M,A_1{n},O{n},A_2{n},'l');
    end
elseif direction == 'r'
    for n = N:(-1):1
        M = applyT(M,A_1{n},O{n},A_2{n},'r');
    end
else
    error(['Unrecognized direction ' direction '.']);
end
end
