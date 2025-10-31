vold = '/Volumes/radiology';
groups = fullfile(vold, 'Groups');
group = fullfile(groups, 'CVMRIGroup');
users = fullfile(group, 'Users');
user = fullfile(users, 'txv016');
wrap2 = fullfile(user, 'WRAP2');
niid = fullfile(wrap2, 'niis', 'niis');
base_dir = [];
base_dir = niid;
dcp = dir(base_dir);
dcp = dcp(~ismember({dcp.name}, {'.','..'}));

outdir = '/Users/TXV016/SaveToRad/WRAP2PROC25/';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
for i = 1:numel(dcp)
    dn = dcp(i).name;
    outfolder = fullfile(outdir, dn);
    if ~exist(outfolder, 'dir')
        mkdir(outfolder)
    end
end