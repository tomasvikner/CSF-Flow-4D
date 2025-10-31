function plotWaveforms(ax, vx_t, vy_t, vz_t)
    plot(ax, [vx_t vy_t vz_t], 'LineWidth', 1.5);
    legend(ax, {'vx','vy','vz'});
    title(ax, 'vx, vy, vz');
    % xlabel(ax, 'Frame');
    ylabel(ax, 'Velocity');
end