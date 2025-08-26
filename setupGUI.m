function [fig, g, ax, s_frame, s_slice] = setupGUI(loadFcn, updateFcn, keyFcn, ncols)
    fig = uifigure('Name','4D CSF Flow Viewer','Position',[100 100 300*ncols 900]);
    g = uigridlayout(fig,[3, ncols+1]);
    g.RowHeight = {'1x','1x','1x'};
    widths = repmat({sprintf('%gx', 1)}, 1, ncols);
    g.ColumnWidth = [widths, {80}];

    ax = gobjects(3, ncols);
    for i = 1:3
        for j = 1:ncols
            ax(i,j) = uiaxes(g);
            ax(i,j).Layout.Row = i;
            ax(i,j).Layout.Column = j;
            colormap(ax(i,j), gray);
        end
    end

    for i = 1:2
        for j = 1:ncols
            ax(i,j).XTick = [];
            ax(i,j).YTick = [];
        end
    end

    uibutton(fig, 'Text','Load 4D data', ...
        'Position',[10 870 100 30], ...
        'ButtonPushedFcn', @(btn,event) loadFcn());

    s_frame = uislider(g,'Limits',[1 10],'MajorTicks',[],'Orientation','vertical');
    s_frame.Layout.Row = 1;
    s_frame.Layout.Column = ncols+1;

    s_slice = uislider(g,'Limits',[1 10],'MajorTicks',[],'Orientation','vertical');
    s_slice.Layout.Row = 2;
    s_slice.Layout.Column = ncols+1;

    % Link slider callbacks to update function
    addlistener(s_frame,'ValueChanged',@(src,evt) updateFcn());
    addlistener(s_slice,'ValueChanged',@(src,evt) updateFcn());

    % Keyboard arrow control
    fig.WindowKeyPressFcn = @(src,event) keyFcn(src,event);

end
