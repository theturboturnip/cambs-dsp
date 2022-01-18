function [ps,locs,qs_h_space_est] = estimate_dc_fft_qs(h)
    y = abs(fft(h))./length(h);
    % Recentre on zero
    y = fftshift(y);
    centre = floor(length(y)/2);

    [ps, locs] = findpeaks(y, ...
        "MinPeakHeight", max(y)/15, ...
        "MinPeakDistance", max(length(y)/100, 1));
    centred_locs = locs-centre-1;

    % Identify quantization
    if length(ps) <= 1
        % One peak = the centre
        % no other peaks => no compression
        qs_h_space_est = [];
    else
        % Assume odd amount of peaks (centre + even amounts on both sides)
        centre_peak = ceil(length(ps)/2);
        % Remove peaks to the left-of-centre
        ps_right = ps(centre_peak:end);
        centred_locs_right = centred_locs(centre_peak:end);
        
        % First scan - if we're just 1x compressed, won't have major increases in height.
        % peak_spacing(i) = location(i+1) - location(i)
        % height_diffs(i) = ps(i+1) - ps(i)
        peak_spacing = diff(centred_locs_right);
        height_diffs = diff(ps_right);
        if all(height_diffs < max(y)/5)
            % We're have no major height increases.
            % Assume the peaks are evenly spaced.
            % The FFT-space quantization interval = the spacing of peaks
            % assume the first peak_spacing is most accurate
            qs_h_space_est = [length(h)./peak_spacing(1)];
        else
            % Either not evenly spaced, or with positive height jumps.
            % We're 2x compressed, now we need to find the intervals.

            % Find the second-highest peak - this corresponds to the
            % *second* compression interval.
            % Exclude the first peak, it's the centre peak so will be the
            % highest.
            [~, i_second_highest] = max(ps_right(2:end));
            i_second_highest = i_second_highest + 1;
            q2_f_space = centred_locs_right(i_second_highest);
            q2_h_space = length(h)./q2_f_space;

            % The minor peak spacings are prop to 1/(q1_hist*q2_hist) =
            % q1_fft * q2_fft
            q1_q2_combined_q = min(peak_spacing);
            % q2_fft / (q1_fft * q2_fft) = 1/q1_fft = q1_histogram
            q1_h_space = q2_f_space / q1_q2_combined_q;

            qs_h_space_est = [q1_h_space q2_h_space];
        end
    end
    qs_h_space_est = round(qs_h_space_est);
end