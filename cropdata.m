function data = cropdata(data) 
zrange = 11:70;
xrange = 71:442;
yrange = 111:482;
fns = fields(data);
for j = 1:numel(fns)
    fn = fns{j};
    ddims = ndims(data.(fn));
    if ddims == 4
        data.(fn) = data.(fn)(xrange, yrange, zrange, :);
    elseif ddims == 3
        data.(fn) = data.(fn)(xrange, yrange, zrange);
    end
end
end