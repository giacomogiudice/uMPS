function varargout = gmres_(Afun,b,options)
% GMRES_ Wrapper for gmres.
%
% x = gmres_(A,b)
% Solves the linear system A*x = b using the default settings.
%
% [x,...] = gmres_(A,b,options)
% Solves the linear system A*x = b using the settings defined in the
% structure options. The options are specified in the VUMPS settings, and
% the outputs are those of the defined solver.
%
% See also: GMRES, VUMPS, VUMPS_SETTINGS.

if exist('options','var') && ~isempty(options) && ~isempty(fieldnames(options))
    [varargout{1:nargout}] = gmres(Afun,b,numel(b),options.tol,options.maxit,[],[],options.v0);
else
    [varargout{1:nargout}] = gmres(Afun,b,numel(b));
end
end
