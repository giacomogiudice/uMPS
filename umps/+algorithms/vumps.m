function varargout = vumps(H,D_list,d,settings)
% VUMPS Optimize an MPS to find the  extremal eigenvector of an operator.
%
% The VUMPS algorithm attempts to find the smallest (or largest)
% eigenvector of an operator as a uniform MPS. This is achieved by
% simultaneously optimizing for the central tensor and the gauging matrix,
% and subsequently finding the canonical forms that best match these
% optimizations. These procedure is iterated until convergence. For a
% detailed description of the algorithm, see the original paper on the <a
% href="matlab:web('https://arxiv.org/abs/1701.07035')">arXiv</a>.
%
% [A_left,A_right,C,A] = vumps(H,D_list,d)
% Given an operator H acting on a space of physical dimetion d, a list of
% bond dimensions D_list, find the MPS approximation of its leading
% eigenvector as an MPS with canonical tensor A_left and A_right, with a
% gauging matrix C and central tensor A.
% [A_left,A_right,C,A] = vumps(H,D_list,d,tol)
% Provide a halting condition based on the gauge error. This is equivalent
% to choosing the settings.tol field.
%
% [A_left,A_right,C,A] = vumps(H,D_list,d,settings)
% Pass specific options through a settings field. A detailed list of the
% possible options is provided below.
%
% [A_left,A_right,C,A,output] = vumps(H,D_list,d,...)
% The output structure contains information on the run, and energy and
% variance of the output state.
%
% [A_left,A_right,C,A,output,blocks] = vumps(H,D_list,d,...)
% The blocks structure contains the environment tensors of the final state.
%
% [A_left,A_right,C,A,output,blocks,stats] = vumps(H,D_list,d,...)
% Additionally returns observables
%
% INDEXING CONVENTION
% The indexing convention for the MPS tensors is (left,right,top), while
% MPO tensors are ordered as (left,right,top,bottom). In the specific case
% of two-site operators, the index ordering is
% (top_left,top_right,bottom_left,bottom_right).
%
% SUPPORTED OPERATOR TYPES
% The algorithm supports different classes of input operators, specifically
%   double          - a rank-4 tensor corresponding to the MPO with the
%                   indexing convention defined above.
%   SchurOperator   - A wrapper for a lower-triangular square cell array
%                   of (d x d) matrices. cell array dimension is to be
%                   understood as the virtual MPO dimension, while d is to
%                   be understood as the physical dimension. The single
%                   matrices can be replaced with scalars if they act
%                   proportional to the identity. Remember to cast the
%                   cell array to a 'SchurOperator' otherwise it will be
%                   interpreted as 'multicell' mode.
%   TwoSiteOperator - a two-body operator. The Hamiltonian is then
%                   understood to be the sum of these gates on each link.
%   cell            - a cell array of operators, either double or
%                   SchurOperators. This is used to represent a non-
%                   trivial unit cell. See below for specific information.
%
% MULTICELL MODE
% For non-trivial unit cells, one can provide a cell array of operators,
% either double or SchurOperators. Each cell acts on one site of the unit
% cell. The behavior is roughly the same, except that the optimized
% tensors A_left, A_right, C  and A will be corresponding cell arrays.
% The cost of the algorithm is slightly higher than N times VUMPS, since
% one must solve for an additional gauging matrix. Environments are
% recomputed at each time for fastest convergence.
%
% INPUT
% H         - Input operator which can be of the types defined above.
%   double      - a rank-4 MPO tensor with the convention defined above
%   SchurOperator   - a lower-triangular cell array of (d x d) matrices. The
%               cell array dimension is to be understood as the virtual MPO
%               dimension, while d is to be understood as the physical
%               dimension. The single matrices can be replaced with scalars
%               if they act proportional to the identity.
%   TwoSiteOperator     - an operator acting on two sites.
%   cell   - a cell array of operators, either in Schur form or
%               generic. See vumps_multicell for specific options about
%               this mode.
% D_list    - Vector of bond dimensions to use for the optimization. By
%           default, the bond dimension is increased in logarithmically
%           spaces intervals between unity and the desired tolerance.
% d         - Physical dimension of the desired MPS.
% settings  - (optional) settings structure for VUMPS (see
%           vumps_settings for detailed options).
% OUTPUT
% A_left    - Left-canonical tensor of the optimized MPS.
% A_right   - Right-canonical tensor of the optimized MPS.
% C         - Gauging matrix of the optimized MPS.
% A         - Center-site tensor of the optimized MPS/
% output    - a structure containing information on the run. Its fields are
%   flag            - termination flag, which indicates tolerance met (0),
%                   stagnation in energy (1), or reached maximum number of
%                   iterations (2).
%   iter            - number of iterations.
%   err             - final tolerance at termination.
%   energy          - corresponding eigenvalue per site of the MPS.
%   energyvariance  - approximate energy variance of the optimized MPS.
% blocks    - a structure containing the fields left and right,
%           corresponding to the (infinite) environment tensors.
% stats     - a structure with detainled information on the run. These are
%   err         - error at each iteration.
%   energy      - energy at each iteration
%   energydiff  - difference of energy between current and previous
%               iteration
%   bond        - bond dimension at each iteration.
%
% EXAMPLE
% Solve the Ising model at transverse field h = 1/2
% X = [0,1;1,0]; Z = [1,0;0,-1];
% H = SchurOperator({1,[],[];-Z,[],[];X/2,Z,1});
% D = 16;    % The bond dimension
% d = 2 ;    % The physical dimension
% [A_left,A_right,C,A,output,~,stats] = vumps(H,D,d);
%
% SUPPORTED OPTIONS
% The settings object has a tree-like structure. The following fields are
% valid:
% tol       - halting condition on the gauge error (default: 1e-12).
% maxit     - maximum number of iterations (default: 100).
% isreal    - tries to solve only for real tensors. This works if the
%           Hamiltonian is real and the state is invariant under
%           reflection.
% verbose   - if true, prints results of each iteration to stdout.
% eigsolver - structure that defines the iterative eigensolver to use and
%           its internal settings. These are:
%   handle  - function handle to the eigensolver (default: @linalg.eigs).
%   mode    - Specifies which eigenvalue to find (default: 'sr').
%   options - Internal options to provide to the eigensolver. These are
%   (more can be provided depending on the specific solver):
%       issym           - flag for symmetric problems (default: true).
%       isreal          - flag for real problems (default: false).
%       tol             - halting tolerance for the eigensolver
%                       (default: eps).
%       dynamictol      - adjust the tolerance dynamically to the global
%                       error, according to the formula
%                       tol = min(maxtol,max(dynamicfactor*err,mintol))
%                       (default: true).
%       dynamicfactor   - factor by which the tolerance is adjusted
%                       dynamically (see dynamictol) .
%       mintol          - lower bound for the tolerance (see dynamictol)
%                       (default: eps).
%       maxtol          - upper bound for the tolerance (see dynamictol)
%                       (default: 1e-3).
% linsolver
%   handle  - function handle to the linear solver. Because MATLAB does not
%           have a consistent way of passing options, the following
%           wrappers are provided in linalg: bicgstab, @bicgstabl, @gmres
%           (default: @linalg.bicgstab_).
%   options - Internal options to provide to the linear equation solver.
%   These are (more can be provided depending on the specific solver):
%       maxit           - maximum number of iterations of the solver
%                       (default: 100).
%       tol             - halting tolerance for the linear solver
%                       (default: eps).
%       dynamictol      - adjust the tolerance dynamically to the global
%                       error, according to the formula
%                       tol = min(maxtol,max(dynamicfactor*err,mintol))
%                       (default: true).
%       dynamicfactor   - factor by which the tolerance is adjusted
%                       dynamically (see dynamictol) .
%       mintol          - lower bound for the tolerance (see dynamictol)
%                       (default: eps).
%       maxtol          - upper bound for the tolerance (see dynamictol)
%                       (default: 1e-3).
% initial   - structure to provide as initial guess. For a valid guess, all
% the following fields muste ge provided
%   A_left  - MPS tensor(s) in left-canonical form
%   A_right - MPS tensor(s) in right-canonical form
%   C       - central tensor(s) such that A_left*C = C*A_right
% advice    - structure used internally to pass guesses to the iterative
%           solvers
%
% SETTINGS EXAMPLE
%   settings.tol = 1e-8;
%   settings.maxit = 15;
%   settings = vumps_settings(settings);
%
% See also: error_gauge, error_variance.


% Parse settings
if ~exist('settings','var') || isempty(settings)
    settings = vumps_settings();
elseif isnumeric(settings) && isscalar(settings)
    tol = settings;
    settings = vumps_settings();
    settings.tol = tol;
else
    settings = vumps_settings(settings);
end
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
        [varargout{1:nargout}] = vumps_multicell(H,D_list,d,settings);
        return
    otherwise
        error('Unrecognized class of input operator.');
end

% Do checks on D_list
assert(all(diff(D_list) > 0),'D_list must be a vector of positive integers in ascending order.');
D = D_list(1);
bond_ind = 1;
growtol = logspace(0,log10(settings.tol),length(D_list)+1);
growtol(1) = [];
growtol(end) = 0;

if all(isfield(settings.initial,{'A_left','A_right','C'}))
    % The initial conditions are provided
    A_left = settings.initial.A_left;
    A_right = settings.initial.A_right;
    C = settings.initial.C;
    assert(isequal(size(A_left),size(A_right)),'Size mismatch between left and right canonical forms.');
    assert(isequal([size(A_left,1),size(A_left,2)],size(C)),'Size mismatch between canonical forms and central tensor.');
    A = ncon({A_left,C},{[-1,1,-3],[1,-2]});
    D = size(C,1);
    assert(D <= D_list(1),'Bond dimension of initial MPS must be smaller or equal to first element of D_list.');
    if settings.isreal && ~all([isreal(A_left),isreal(A_right),isreal(C)])
        settings.isreal = false;
        settings.eigsolver.options.isreal = false;
    end
else
    % Nothing is provided, generate at random
    A = randn([D,D,d]);
    if ~settings.isreal
        A = A + 1i*randn([D,D,d]);
    end
    C = diag(rand(D,1));
    [A_left,A_right] = update_canonical(A,C);
end
% Initialize stats log
stats = struct;
savestats = false;
if nargout >= 6
    savestats = true;
    stats.err = zeros(1,settings.maxit);
    stats.energy = zeros(1,settings.maxit);
    stats.energydiff = zeros(1,settings.maxit);
    stats.bond = zeros(1,settings.maxit);
end
% Error in gauge transformation relationships
err = error_gauge(A_left,A_right,C,A);
% Update tolerances
settings = update_tol(settings,err);
% Generate the environment blocks
[B_left,B_right,energy_prev] = update_environments(A_left,A_right,C,H,[],[],settings);
% Main VUMPS loop
output.flag = 2;
if settings.verbose
    fprintf('Iter\t      Energy\t Energy Diff\t Gauge Error\t    Bond Dim\tLap Time [s]\n')
    fprintf('   0\t%12g\n',energy_prev);
end
for iter = 1:settings.maxit
    tic
    % Solve effective problems for A and C
    if ~all([isreal(B_left),isreal(B_right)])
        settings.isreal = false;
        settings.eigsolver.options.isreal = false;
    end
    [A,C] = solve_local(A_left,A_right,C,[],A,H,B_left,B_right,settings);
    % Update the canonical forms
    [A_left,A_right] = update_canonical(A,C);
    % Update the environment blocks
    [B_left,B_right,energy] = update_environments(A_left,A_right,C,H,B_left,B_right,settings);
    % Update tolerances
    settings = update_tol(settings,err);
    laptime = toc;
    % Get error
    err = error_gauge(A_left,A_right,C,A);
    energydiff = energy - energy_prev;
    % Print results of interation
    if settings.verbose
        fprintf('%4d\t%12g\t%12g\t%12g%12d%12.1f\n',iter,energy,energydiff,err,D_list(bond_ind),laptime);
    end
    if savestats
        stats.err(iter) = err;
        stats.energy(iter) = energy;
        stats.energydiff(iter) = energydiff;
        stats.bond(iter) = D;
    end
    if bond_ind == length(D_list)
        if err < settings.tol
            % Stopping condition on the gauge error
            output.flag = 0;
            break
        elseif abs(energy_prev - energy) < eps
            % Stopping condition on stagnation
            output.flag = 1;
            break
        end
    elseif err < growtol(bond_ind)
        % Increase bond dimension
        bond_ind = bond_ind + 1;
        D = D_list(bond_ind);
        [A_left,A_right,C,A,B_left,B_right] = increasebond(D,A_left,A_right,C,H,B_left,B_right);
        [B_left,B_right] = update_environments(A_left,A_right,C,H,B_left,B_right,settings);
    end
    energy_prev = energy;
end
% Set output information
output.iter = iter;
output.err = err;
output.energy = energy;
output.energyvariance = error_variance(A_left,A_right,C,H,B_left,B_right);
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
    stats.err = stats.err(1:iter);
    stats.energy = stats.energy(1:iter);
    stats.energydiff = stats.energydiff(1:iter);
    stats.bond = stats.bond(1:iter);
    varargout{7} = stats;
end
end
