function varargout = bicgstabl_(Afun,b,options)
% BICGSTAB_ Wrapper for bicstabl.
%
% x = bicgstabl_(A,b)
% Solves the linear system A*x = b using the default settings.
%
% [x,...] = bicgstabl_(A,b,options)
% Solves the linear system A*x = b using the settings defined in the
% structure options. The options are specified in the VUMPS settings, and
% the outputs are those of the defined solver.
%
% See also: BICGSTABL, VUMPS, VUMPS_SETTINGS.

if exist('options','var') && ~isempty(options) && ~isempty(fieldnames(options))
    [varargout{1:nargout}] = bicgstabl(Afun,b,options.tol,options.maxit,[],[],options.v0);
else
    [varargout{1:nargout}] = bicgstabl(Afun,b);
end
end
