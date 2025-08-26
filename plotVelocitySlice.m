function plotVelocitySlice(ax, vol, slice, frame, clim, titleText)
    imagesc(ax, squeeze(vol(:,:,slice,frame))');
    axis(ax, 'image');
    title(ax, titleText);
    caxis(ax, clim);
end