function [rmax, mlag] = waveformCoupling_V3D(wf1, wf2, maxlag)

    maxframes = numel(wf1);

    % Wrap the waveform 
    wf1 = [wf1(:, end-maxlag+1:end), wf1, wf1(:, 1:maxlag)];
    wf2 = [wf2(:, end-maxlag+1:end), wf2, wf2(:, 1:maxlag)];

    wf1 = zscore(wf1);
    wf2 = zscore(wf2);

    x = 1:numel(wf1);
    xq = linspace(1, 20, 1000);
    wf1 = interp1(x,wf1,xq, "pchip");
    wf2 = interp1(x,wf2,xq, "pchip");

    maxlag = maxlag * numel(wf2) / maxframes; 

    % Cross-correlation
    [rvals, lags] = xcorr(wf2, wf1, maxlag, 'coeff');
    [~, mi] = max(abs(rvals));
    mlag = lags(mi);
    rmax = rvals(mi);

    mlag = mlag * 1000 / numel(wf1);

end