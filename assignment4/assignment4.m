function assignment4(name, color_type, expected_value)
    qs = evaluate_compression(name, color_type);

    fprintf("%d,\t%s,\t%s,\t", length(qs), color_type, name)
    if length(qs) == expected_value
        fprintf("CORRECT,\t")
    else
        fprintf("INCORRECT(%d),\t", expected_value)
    end
    for q = qs
        qual_min = estimate_qual_from_q(q + 1);
        qual_max = estimate_qual_from_q(q);
        fprintf("%d,\t%3d-%-3d,\t", q, qual_min, qual_max);
    end
    fprintf("\n")
end
function qs=evaluate_compression(name, color_type)
    Di = image_DCT_params(name, color_type);
    jk = [1 1];
    d = Di{jk(1), jk(2)};
    qs = identify_qs("dc", d);
end
function qual=estimate_qual_from_q(q)
    % From cjpeg source:
    % if (quality < 50)
    %   scaling_factor = 5000 / quality;
    % else
    %   scaling_factor = 200 - quality*2;
    %
    % q = (original_q * scaling_factor + 50) / 100
    % 
    % From cjpeg source, we also know the default table's original_q for
    % (1,1) is 16
    % scaling_factor = ((q * 100) - 50)/original_q

    % cjpeg can clamp the quantization at 255 or 32767.
    % If q was clamped, we can't reverse the computation.
    % To be on the safe side, if we are given a q equal to these values
    % note they may have been clamped.
    if q == 255 || q == 32767
        warning("q-value %d may have been clamped", q)
    end

    original_q = 16;
    scaling_factor = (q * 100 - 50)/original_q;

    % for quality in [1, 50], scaling factor in [5000, 100]
    % for quality in [50, 100], scaling factor in [100, 0]
    if scaling_factor <= 100
        % Invert the bottom branch
        % scale = 200 - quality*2
        % => quality = (200 - scale)/2
        qual = (200 - scaling_factor)/2;
    else
        % Invert the top branch
        % scaling_factor = 5000 / quality
        % quality = 5000 / scaling_factor 
        qual = 5000 / scaling_factor;
    end

    qual = round(qual);
end