tensorfun(f,g,s)
% TENSORFUN returns handle to f(g((t)) when f accepts vectors only
if prod(size(s)) == 1
    s = [s,1];
end
