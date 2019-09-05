% UMPS
%
% Files
%   applyH2s         - Apply a 2-site MPS tensor to an effective Hamiltonian.
%   applyHA          - Apply an MPS tensor to an effective Hamiltonian.
%   applyHC          - Apply a central tensor to an effective Hamiltonian.
%   applyT           - Computes the application of a tensor to the transfer matrix.
%   applyTv          - Wrapper for applyT, but accepts a vector and returns a vector.
%   bicgstab_        - Wrapper for bicstab.
%   bicgstabl_       - BICGSTAB_ Wrapper for bicstabl.
%   canonical        - Generate the canonical forms for an MPS tensor.
%   cell2tensor      - Convert a Schur operator into an MPO tensor.
%   combineMPO       - Combines a cell array of Schur operators to a single operator.
%   error_gauge      - Return error in the gauge between the MPS tensors.
%   error_variance   - Approximates the error in the variance.
%   fixedblock       - Compute environment tensor for an MPS tensor.
%   fixedpoint       - Computes the leading eigenvector of a transfer matrix.
%   geometric_series - Compute the geometric sum of a matrix applied to a vector.
%   gmres_           - Wrapper for gmres.
%   increasebond     - Optimally increase the bond dimension of the MPS tensors.
%   nullspace        - Compute the nullspace of an MPS tensor.
%   reduceMPO        - Reduces the dimension of an MPO tensor by removing linear dependencies.
%   update_canonical - Update the canonical forms of an MPS tensor.
%   update_tol       - Update tolerance for the internal solvers given the error.
%   vumps            - Find extremal eigenvalues of an MPO or a Hamiltonian.
%   vumps_multicell  - Version of VUMPS for non-trivial unit cells.
%   vumps_settings   - Construct settings structure for vumps problem.
