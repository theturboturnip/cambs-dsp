function qs_est = estimate_dc_fft_qs(h)
    y = abs(fft(h))./length(h);
    % Recentre on zero
    y = fftshift(y);
    centre = floor(length(y)/2);

    % Use max(length(y)/50, 1) as minpeakdistance - for the ys of length
    % ~2000, 20 is great, but for smaller ys (and potentially AC terms)
    % don't want to be that limiting
    [ps, locs] = findpeaks(y, "MinPeakHeight", max(y)/15, "MinPeakDistance", max(length(y)/100, 1));
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
            if (p - ps(first_peak_after_centre)) > max(ps)/5%(ps(first_peak_after_centre)/5)
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
function qs_est = estimate_dc_fft_qs_new(h)
    y = abs(fft(h))./length(h);
    % Recentre on zero
    y = fftshift(y);
    centre = floor(length(y)/2);

    % Use max(length(y)/50, 1) as minpeakdistance - for the ys of length
    % ~2000, 20 is great, but for smaller ys (and potentially AC terms)
    % don't want to be that limiting
    [ps, locs] = findpeaks(y, "MinPeakHeight", max(y)/15, "MinPeakDistance", max(length(y)/100, 1));
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

        % Find the second-tallest peak after the centre
        [tallest, i_tallest] = max(ps(first_peak_after_centre:end));
        % Recentre index relative to all peaks
        i_tallest = i_tallest - 1 + first_peak_after_centre;

        % Must at least be periodic in second-tallest peak
        qs_est = [ centred_locs(i_tallest) ];

        % If that wasn't the peak directly after the centre,
        % 


        % Now, check for a lower peak for a second Q-level
        % Go through subsequent peaks consecutively
        for i = first_peak_after_centre:length(locs)
            if i == i_tallest
                continue
            end

            p = ps(i);
            % If this peak is smaller
            if (ps(i_tallest) - p) > max(ps)/5
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