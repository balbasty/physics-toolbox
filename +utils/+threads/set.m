function n = set(n)
% Set the number of OMP threads.
%
% FORMAT n = utils.threads.set(n)
%
% If n == -1, use the machine's max number of processors (omp_num_procs())
% Returns the actual number of threads used.
if nargin > 0
    setenv('SPM_NUM_THREADS',sprintf('%d',n));
    try, spm_diffeo; end
end
n = sscanf(getenv('SPM_NUM_THREADS'),'%d');
