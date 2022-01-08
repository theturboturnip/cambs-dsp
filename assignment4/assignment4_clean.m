

imgs = [
    "test1.png" "gray";
    "test2.png" "gray";
    "test3.png" "gray";
    "test1c.png" "color";
    "test2c.png" "color";
    "test3c.png" "color";
    "./training/converted/test5.png" "color";
];
imgs_di = {};
for i_img = 1:size(imgs, 1)
    Di = image_DCT_params(imgs(i_img,1), imgs(i_img, 2));
    imgs_di{i_img} = Di;
end

%% Display data
jk_vals = [
    1 1;
%     1 2;
%     1 3;
%     2 1;
%     2 2;
%     3 1;
%     4 4;
];
figure;
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)
        histogram(d, 'BinMethod', 'Integer');
        
        h_sort = sort(unique(h));
%         ylim([0, h_sort(5)*5]);
        ylim([0, h_sort(end-5)]);
        title(imgs(i_img) + " - DCT{ " + jk(1) + "," + jk(2) + " }");
    end
end

%% FFT data
q_colors = ['m' 'r'];

figure;
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

%         h = h - mean(h);
        y = abs(fft(h));
        % Recentre on zero
        y = fftshift(y);
        centre = floor(length(y)/2);
        plot(centre-(0:(length(y)-1)), y);
        
        hold on;
        [ps, locs] = findpeaks(y, "MinPeakHeight", max(y)/10, "MinPeakDistance", 20);
        centred_locs = locs-centre-1;
        plot(centred_locs, ps, 'r*');
        hold off;
        
        xlim([-centre, centre]);
        
        hold on;
        
%         % Identify quantization
%         if length(ps) <= 1
%             % One peak = the centre
%             % no other peaks => no compression
%             qs = [];
%         else
%             % Assume odd amount of peaks (centre + even amounts on both
%             % sides)
%             first_peak_after_centre = ceil(length(ps)/2)+1;
%             % Must at least be periodic in first peak
%             qs = [ centred_locs(first_peak_after_centre) ];
%             
%             % See if there's a second Q-level
%             % Go through subsequent peaks consecutively
%             for i = first_peak_after_centre+1:length(locs)
%                 p = ps(i);
%                 % If it's significantly bigger than the first peak,
%                 % it's likely a second Q-level
%                 if (p - ps(first_peak_after_centre)) > (ps(first_peak_after_centre)/5)
% %                     yline(p)
% %                     yline(ps(first_peak_after_centre))
%                     qs = [qs(1) centred_locs(i)];
%                     break
%                 end
%                 % Or, if the last peak location isn't a factor of this one,
%                 % e.g. if centred_locs(i)/qs(1) isn't close to an integer,
%                 % likely a second Q-level
%                 d = centred_locs(i)/qs(1);
%                 if abs(d - round(d)) > 0.1
%                     qs = [qs(1) centred_locs(i)];
%                     break
%                 end
%             end
%         end
        
        qs=identify_qs(h)
        
        % Plot Q-levels
        hold on;
        for i_q = 1:length(qs)
            q = qs(i_q);
            % Find multiples of q within the window [-centre, centre]
            % Start from fix(-centre/q) (fix = round-towards-zero)
            for n = fix(-centre/q):fix(centre/q)
                xline(n*q, q_colors(i_q));
            end
        end
        hold off;
        
        t = [
            imgs(i_img) + " - FFT of DCT{ " + jk(1) + "," + jk(2) + " }";
            "qs = [" + strjoin(string(qs), ", ") + "]"
        ];
        title(t);
%         subtitle("qs = " + strjoin(string(qs), ", "));
    end
end