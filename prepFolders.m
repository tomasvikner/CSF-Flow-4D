
% dcp = dir(base_dir);
% dcp = dcp(~ismember({dcp.name}, {'.','..'}));
% 
% outdir = 'TEMP';
% if ~exist(outdir, 'dir')
%     mkdir(outdir)
% end

% for i = 1:numel(dcp)
%     dn = dcp(i).name;
%     outfolder = fullfile(outdir, dn);
%     if ~exist(outfolder, 'dir')
%         mkdir(outfolder)
%     end
% end