function s = formatRegShift(v)
    parts = cell(1, 3);
    for k = 1:3
        vk = round(v(k));
        if vk > 0
            parts{k} = sprintf('+%d', vk);
        else
            parts{k} = sprintf('%d', vk);
        end
    end
    s = sprintf('(%s, %s, %s)', parts{1}, parts{2}, parts{3});
end
