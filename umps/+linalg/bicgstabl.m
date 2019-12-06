function varargout = bicgstabl(A,B,options)
% BICGSTABL Wrapper for bicstab supporting tensors.
%
% x = bicgstabl(A,b)
% Solves the linear system A*x = b using the default settings.
%
% [x,...] = bicgstabl(A,b,options)
% Solves the linear system A*x = b using the settings defined in the
% structure options. The options are specified in the VUMPS settings, and
% the outputs are those of the defined solver.
%
% [X,...] = bicgstabl(Afun,B,...)
% Solves the linear system A*x = b when Afun is a linear operator, and B
% is a tensor. The returned object X has the same shape as B.
%
% See also: BICGSTABL, ALGORITHMS/VUMPS.

if isnumeric(A)
    Afun = A;
    b = B;
else
    s = size(B);
    n = prod(s);
    Afun = @(x) reshape(A(reshape(x,s)),[n,1]);
    b = reshape(B,[n,1]);
end

if exist('options','var') && ~isempty(options) && ~isempty(fieldnames(options))
    if isfield(options,'v0') && ~isempty(options.v0)
        options.v0 = reshape(options.v0,[n,1]);
    end
    [varargout{1:nargout}] = bicgstabl(Afun,b,options.tol,options.maxit,[],[],options.v0);
else
    [varargout{1:nargout}] = bicgstabl(Afun,b);
end
varargout{1} = reshape(varargout{1},s);
end
