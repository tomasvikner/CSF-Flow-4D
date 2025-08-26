function m = plotMarker(ax, x, y, color)
    hold(ax, 'on');
    m(1) = plot(ax, [x-3 x+3], [y y], color, 'LineWidth', 2);
    m(2) = plot(ax, [x x], [y-3 y+3], color, 'LineWidth', 2);
    hold(ax, 'off');
end