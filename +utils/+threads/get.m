function n = get
% Get the current number of OMP threads.
%
% FORMAT n = utils.threads.get()
n = sscanf(getenv('SPM_NUM_THREADS'),'%d');
