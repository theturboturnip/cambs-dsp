function qs = identify_qs(type, d)
    % This is a legacy function, it used to more than this but it just
    % calls estimate_dc_fft_qs now.
    h = histcounts(abs(d),'BinMethod','Integer');
    [~,~,qs] = estimate_dc_fft_qs(h);
end


