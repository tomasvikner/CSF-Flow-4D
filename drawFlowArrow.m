function arrow = drawFlowArrow(ax, x, y, direction, scale)
    dir_xy = -direction([2 3]);
    if norm(dir_xy) > 0
        dir_xy = dir_xy / norm(dir_xy) * scale;
        hold(ax, 'on');
        arrow = quiver(ax, x, y, dir_xy(1), dir_xy(2), ...
                       'Color', 'r', 'LineWidth', 2, ...
                       'MaxHeadSize', 2, 'AutoScale', 'off');
        hold(ax, 'off');
    else
        arrow = [];
    end
end