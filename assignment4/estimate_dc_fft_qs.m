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
function [ps,locs,qs_est] = estimate_dc_fft_qs_old(h)
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
        qs_est = [];
    else
        % Assume odd amount of peaks (centre + even amounts on both
        % sides)
        first_peak_after_centre = ceil(length(ps)/2)+1;
        % Must at least be periodic in first peak
        qs_est = [ centred_locs(first_peak_after_centre) ];

        % See if there's a second Q-level
        % Go through subsequent peaks consecutively
        for i = first_peak_after_centre+1:length(locs)
            p = ps(i);
            % If it's significantly bigger than the first peak,
            % it's likely due to two quantizations aligning at this point.
            % Return this peak, but it's not necessarily a Q-level
            if (p - ps(first_peak_after_centre)) > max(ps)/5
                qs_est = [qs_est(1) centred_locs(i)];
                break
            end
            % Or, if the last peak location isn't a factor of this one,
            % e.g. if centred_locs(i)/qs(1) isn't close to an integer,
            % likely a second Q-level
            d = centred_locs(i)/qs_est(1);
            if abs(d - round(d)) > 0.1
                qs_est = [qs_est(1) centred_locs(i)];
                break
            end
        end
    end
end