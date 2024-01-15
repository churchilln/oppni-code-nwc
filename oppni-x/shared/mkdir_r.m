function mkdir_r(pathstr)
% 
% mkdir_r(pathstr) --> populates full directory structure as specified
% 

if(~iscell(pathstr))
    pathstr = {pathstr};
end

for(n=1:numel(pathstr))

    ind = strfind(pathstr{n},'/');
    for i = 1:length(ind)
        current_path = pathstr{n}(1:ind(i)-1);
        if ~isempty(current_path)
            if ~exist(current_path,'dir')
                mkdir(current_path);
            end
        end
    end
    if ~exist(pathstr{n},'dir')
        mkdir(pathstr{n});
    end
end