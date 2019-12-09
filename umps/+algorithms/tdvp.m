function varargout = tdvp(H,A_left,A_right,C,A,settings)
% TDVP Perform time-evolution in the MPS manifold.

% Parse settings
if ~exist('settings','var') || isempty(settings)
	settings = tdvp_settings();
else
	settings = tdvp_settings(settings);
end

% Perform size checks on input tensors
assert(isequal(size(A_left),size(A_right),size(A)),'Size mismatch between input tensors.');
assert(isequal([size(A_left,1),size(A_left,2)],size(C)),'Size mismatch between canonical forms and central tensor.');
[D,~,d] = size(A);

% Perform checks on input operator
switch class(H)
    case 'double'
        chi = size(H,1);
        assert(isequal(size(H),[chi,chi,d,d]),'Input operator should be a rank-4 tensor.');
    case 'SchurOperator'
        chi = size(H,1);
        assert(isequal(size(H),[chi,chi]),'Size mismatch in input operator.');
    case 'TwoSiteOperator'
        assert(isequal(size(H),[d,d,d,d]),'Input operator should be a rank-4 tensor.');
    case 'cell'
        assert(size(H,1) == 1,'Input operator must be an (1 x n) cell array.');
        % Todo implement multicell mode
        error('Multicell mode not yet supported.')
        return
    otherwise
        error('Unrecognized class of input operator.');
end

% Initialize stats log
stats = struct;
savestats = false;
customeval = false;
if nargout >= 6
	savestats = true;
	stats.err = zeros(1,settings.niter);
	stats.energy = zeros(1,settings.niter);
	stats.energydiff = zeros(1,settings.niter);
	if isfield(settings,'customfun') && isa(settings.customfun,'function_handle')
		customeval = true;
		stats.customvalue = cell(1,settings.niter);
	end
end

% Error in gauge transformation relationships
err = error_gauge(A_left,A_right,C,A);
% Generate the environment blocks
[B_left,B_right,energy_prev] = update_environments(A_left,A_right,C,H,[],[],settings);
% Main time-stepping loop
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\t    Bond Dim\tLap Time [s]\n')
	fprintf('   0\t%12g\n',energy_prev);
end
output.flag = 0;
for iter = 1:settings.niter
	tic
	if ~all([isreal(B_left),isreal(B_right)])
		settings.isreal = false;
	end
	% Perform time-integration step for A and C
	[A,C] = tdvp_step(settings.timestep,A_left,A_right,C,[],A,H,B_left,B_right,settings);
	% Update the canonical forms
	[A_left,A_right] = update_canonical(A,C);
	% Update the environment blocks
	[B_left,B_right,energy] = update_environments(A_left,A_right,C,H,B_left,B_right,settings);
	laptime = toc;
	% Get error
	err = error_gauge(A_left,A_right,C,A);
	energydiff = energy - energy_prev;
	% Print results of iteration
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12g%12d%12.1f\n',iter,energy,energydiff,err,D,laptime);
	end
	if savestats
		stats.err(iter) = err;
		stats.energy(iter) = energy;
		stats.energydiff(iter) = energydiff;
		if customeval
			try
				stats.customvalue{iter} = settings.customfun(A_left,A_right,C,A,H,B_left,B_right,settings);
			catch me
				warning('Failed to evaluate custom function.');
				disp(getReport(me,'extended','hyperlinks','on'));
				stats.customfun{iter} = [];
			end
		end
	end
	energy_prev = energy;
end
% Set output information
output.iter = iter;
output.energy = energy;
if nargout >= 5
	blocks.left =  B_left;
	blocks.right = B_right;
end
% Choose gauge in which C is diagonal
[U,S,V] = svd(C,'econ');
A_left = ncon({U',A_left,U},{[-1,1],[1,2,-3],[2,-2]});
A_right = ncon({V',A_right,V},{[-1,1],[1,2,-3],[2,-2]});
A = ncon({U',A,V},{[-1,1],[1,2,-3],[2,-2]});
C = S/norm(S,'fro');
% Place results in varargout
varargout = {A_left,A_right,C,A,output};
if nargout >= 5
    [B_left,B_right] = update_environments(A_left,A_right,C,H,B_left,B_right,settings);
    blocks.left =  B_left;
    blocks.right = B_right;
    varargout{6} = blocks;
end
if savestats
    varargout{7} = stats;
end
end
