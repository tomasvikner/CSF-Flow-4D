function plotWaveforms(ax, vx_t, vy_t, vz_t)
    plot(ax, [vx_t vy_t vz_t], 'LineWidth', 1.5);
    legend(ax, {'vx (LR)','vy (AP)','vz (SI)'});
    % title(ax, 'vx, vy, vz');
    title(ax, 'Voxel waveforms');
    % xlabel(ax, 'Frame');
    ylabel(ax, 'Velocity');
end