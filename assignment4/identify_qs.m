function qs=identify_qs(h)
    y = abs(fft(h));
    % Recentre on zero
    y = fftshift(y);
    centre = floor(length(y)/2);

    [ps, locs] = findpeaks(y, "MinPeakHeight", max(y)/10, "MinPeakDistance", 20);
    centred_locs = locs-centre-1;

    xlim([-centre, centre]);

    % Identify quantization
    if length(ps) <= 1
        % One peak = the centre
        % no other peaks => no compression
        qs = [];
    else
        % Assume odd amount of peaks (centre + even amounts on both
        % sides)
        first_peak_after_centre = ceil(length(ps)/2)+1;
        % Must at least be periodic in first peak
        qs = [ centred_locs(first_peak_after_centre) ];

        % See if there's a second Q-level
        % Go through subsequent peaks consecutively
        for i = first_peak_after_centre+1:length(locs)
            p = ps(i);
            % If it's significantly bigger than the first peak,
            % it's likely a second Q-level
            if (p - ps(first_peak_after_centre)) > (ps(first_peak_after_centre)/5)
%                     yline(p)
%                     yline(ps(first_peak_after_centre))
                qs = [qs(1) centred_locs(i)];
                break
            end
            % Or, if the last peak location isn't a factor of this one,
            % e.g. if centred_locs(i)/qs(1) isn't close to an integer,
            % likely a second Q-level
            d = centred_locs(i)/qs(1);
            if abs(d - round(d)) > 0.1
                qs = [qs(1) centred_locs(i)];
                break
            end
        end
    end
end