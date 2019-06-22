function [B,E] = fixedblock(H,A,direction,settings)
if strcmp(settings.mode,'generic')
	[B,E] = fixedblock_generic(H,A,direction,settings);
elseif strcmp(settings.mode,'schur')
	[B,E] = fixedblock_schur(H,A,direction,settings);
elseif strcmp(settings.mode,'multicell')
	[B,E] = fixedblock_multicell(H,A,direction,settings);
elseif strcmp(settings.mode,'twosite')
	[B,E] = fixedblock_twosite(H,A,direction,settings);
else
	error(['Unrecognized mode ' settings.mode '.']);
end
end

function [B,E] = fixedblock_generic(H,A,direction,settings)
if exist('settings','var') && isfield(settings,'advice') && isfield(settings.advice,'B')
	settings.eigsolver.options.v0 = reshape(settings.advice.B,[],1);
end
[B,E] = fixedpoint(A,H,A,direction,settings);
end

function [B,E] = fixedblock_schur(H,A,direction,settings)
chi = size(H,1);
[D,~,d] = size(A);

% Compute dominant eigenvectors of A
[rho_left,rho_right] = dominant_eigenvector_canonical(A,direction,settings);

% Start from corner and solver for each B_a
if direction == 'l'
	corner = chi;
	ind = (chi-1):-1:1;
	rho = rho_right;
else
	corner = 1;
	ind = 2:1:chi;
	rho = rho_left;
end
assert(H{corner,corner} == 1,'Element H{%d,%d} must be 1.',corner,corner);
B(:,corner) = reshape(eye(D,D),[D*D,1]);
for a = ind
	% Compute new Y_a
	Y = zeros([D^2,1]);
	if direction == 'l'
		for b = chi:-1:(a+1)
			if ~isempty(H{b,a})
				Y = Y + applyTv(B(:,b),A,H{b,a},A,'l');
			end
		end
	else
		for b = 1:(a-1)
			if ~isempty(H{a,b})
				Y = Y + applyTv(B(:,b),A,H{a,b},A,'r');
			end
		end
	end
	% Compute new B_a
	if isempty(H{a,a})
		B(:,a) = Y;
	else
		assert(isscalar(H{a,a}),'Element H{%d,%d} is not a scalar.',a,a);
		% Load initial guess
		if isstruct(settings) & isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		if H{a,a} == 1
			v_left = rho_left;
			v_right = rho_right;
		elseif abs(H{a,a}) < 1
			v_left = 0;
			v_right = 0;
		else
			error('Diagonal element in H has absolute value larger than 1.');
		end
		Tfun = @(x) applyTv(x,A,H{a,a},A,direction);
		B(:,a) = geometric_series(Y,Tfun,v_left,v_right,direction,settings);
	end
end
B = reshape(B,[D,D,chi]); 
E = real(Y.'*rho);
end

function [B,E] = fixedblock_twosite(H,A,direction,settings)
[D,~,d] = size(A);

% Compute dominant eigenvectors of A
[rho_left,rho_right] = dominant_eigenvector_canonical(A,direction,settings);

% Compute relevant boundary
block = ncon({A,A},{[-1,1,-2],[1,-4,-3]});
if direction == 'l'
	h = ncon({conj(block),H,block},{[5,1,2,-1],[1,2,3,4],[5,3,4,-2]});
	rho = rho_right;
elseif direction == 'r'
	h = ncon({conj(block),H,block},{[-1,1,2,5],[1,2,3,4],[-2,3,4,5]});
	rho = rho_left;
else
	error(['Unrecognized direction ' direction '.']);
end

% Solve linear system
if isstruct(settings) & isfield(settings.advice,'B')
	settings.linsolver.options.v0 = reshape(settings.advice.B,[D*D,1]);
end
h = reshape((h+h')/2,[D*D,1]);
Tfun = @(x) applyTv(x,A,1,A,direction);
B = geometric_series(h,Tfun,rho_left,rho_right,direction,settings);
B = reshape(B,[D,D]);
E = real(h.'*rho);
end

function [B,E] = fixedblock_multicell(H,A,direction,settings)
if ~iscell(H{1})
	[B,E] = fixedblock_multicell_generic(H,A,direction,settings);
	return
end
N = length(H);
assert(isequal(size(H),size(A)),'Size mismatch in input cell arrays.');
chi = size(H{1},1);
[D,~,d] = size(A{1});
id = reshape(eye([D,D]),[D*D,1]);
B = zeros([D*D,chi]);
Hall = combineMPO(H);

[rho_left,rho_right] = dominant_eigenvector_canonical(A,direction,settings);

% Start from corner and solver for each B_a
if direction == 'l'
	corner = chi;
	ind = (chi-1):-1:1;
	rho = rho_right;
else
	corner = 1;
	ind = 2:1:chi;
	rho = rho_left;
end
assert(all(cellfun(@(h) h == 1,Hall{corner,corner})),'All {%d,%d} H must be 1.',corner,corner);
B(:,corner) = reshape(eye(D,D),[D*D,1]);
for a = ind
	% Compute new Y_a
	Y = zeros([D^2,1]);
	if direction == 'l'
		for b = chi:-1:a
			if ~isempty(Hall{b,a})
				for s = 1:size(Hall{b,a},1)
					Y = Y + applyTv(B(:,b),A,Hall{b,a}(s,:),A,'l');
				end
			end
		end
	else
		for b = 1:(a-1)
			if ~isempty(Hall{a,b})
				for s = 1:size(Hall{a,b},1)
					Y = Y + applyTv(B(:,b),A,Hall{a,b}(s,:),A,'r');
				end
			end
		end
	end
	% Compute new B_a
	Hdiag = Hall{a,a}; 
	if isempty(Hdiag)
		B(:,a) = Y;
	else
		assert(all(cellfun(@(h) isscalar(h),Hdiag)),'One or more diagonal elements in H are not scalar.');
		% Load initial guess
		if isstruct(settings) & isfield(settings.advice,'B')
			settings.linsolver.options.v0 = reshape(settings.advice.B(:,:,a),[D*D,1]);
		end
		if all(cellfun(@(h) h == 1,Hdiag))
			v_left = rho_left;
			v_right = rho_right;
		elseif abs(H{a,a}) < 1
			v_left = 0;
			v_right = 0;
		else
			error('Diagonal element in H has absolute value larger than 1.');
		end
		Tfun = @(x) prod(cell2mat(Hdiag))*applyTv(x,A,1,A,direction);
		B(:,a) = geometric_series(Y,Tfun,v_left,v_right,direction,settings);
	end
end
B = reshape(B,[D,D,chi]);
E = real(Y.'*rho)/N;
end


function [B,E] = fixedblock_multicell_generic(H,A,direction,settings)
N = length(H);
chi = size(H{1},1);
[D,~,d] = size(A{1});
if exist('settings','var') && isfield(settings,'advice') && isfield(settings.advice,'B')
	settings.eigsolver.options.v0 = reshape(settings.advice.B,[],1);
end
[B,E] = fixedpoint(A,H,A,direction,settings);
E = E/N;
end

function [v_left,v_right] = dominant_eigenvector_canonical(A,direction,settings)
if iscell(A)
	[D,~,d] = size(A{1});
else
	[D,~,d] = size(A);
end
id = reshape(eye(D,D),[D*D,1]);
if direction == 'l'
	% Left block
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		settings.eigsolver.options.v0 = reshape((C*C').',[D*D,1]);
	end
	rho = fixedpoint(A,1,A,'r');
	rho = reshape(rho/trace(rho),[D*D,1]);
	v_left = id;
	v_right = rho;
elseif direction == 'r'
	% Right block
	if isstruct(settings) & isfield(settings.advice,'C')
		C = settings.advice.C;
		settings.eigsolver.options.v0 = reshape(C'*C,[D*D,1]);
	end
	rho = fixedpoint(A,1,A,'l');
	rho = reshape(rho/trace(rho),[D*D,1]);
	v_left = rho;
	v_right = id;
else
	error(['Unrecognized direction ' direction '.']);
end
end


function v = applyTvNtimes(v,A1,H,A2,direction)
N = length(A1);
if direction == 'l'
	ind = 1:N;
elseif direction == 'r'
	ind = N:(-1):1;
end
for n = ind
	v = applyTv(v,A1{n},H{n},A2{n},direction);
end
end