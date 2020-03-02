function mkpath(path)
% Create a path (= recursively create all parent directories)
%
% FORMAT utils.mkpath(path)
path = strsplit(path, filesep);
if isempty(path{1}), path{1} = filesep; end
% TODO: ^ not sure about this, especially on windows...
for i=1:numel(path)
    folder = fullfile(path{1:i});
    if ~exist(folder, 'dir')
        st = mkdir(folder);
        if ~st
            error('Cannot create output folder %s', folder);
        end
    end
end