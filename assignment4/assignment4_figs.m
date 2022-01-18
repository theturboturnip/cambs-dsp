%% Data
imgs = [
    "test1.png" "gray";
    "test1c.png" "color";
    "test2.png" "gray";
    "test2c.png" "color";
    "test3.png" "gray";
    "test3c.png" "color";
];
img_names = [
    "Test 1",
    "Test 1C",
    "Test 2",
    "Test 2C",
    "Test 3",
    "Test 3C",
];
imgs_di = {};
for i_img = 1:size(imgs, 1)
    Di = image_DCT_params(imgs(i_img,1), imgs(i_img, 2));
    imgs_di{i_img} = Di;
end

%% Figure 1
% Histogram of DC component
figure;
jk_vals = [
    1 2;
];
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        histogram(d, 'BinMethod', 'Integer');
        
%         h_sort = sort(unique(h));
%         ylim([0, h_sort(end-5)]);
        xlim([-50 50]);
%         title(imgs(i_img));% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
sgtitle("Histogram of DC component of DCT");

%% Figure 2
% FFT of Histogram of DC component of DCT
figure;
jk_vals = [
    1 1;
];
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        y = abs(fft(h));
        % Recentre on zero
        y = fftshift(y);
        centre = floor(length(y)/2);
        plot(centre-(0:(length(y)-1)), y);
        xlim([0 1000]);
        
%         h_sort = sort(unique(h));
%         ylim([0, h_sort(end-5)]);
        title(imgs(i_img));% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
sgtitle("FFT of Histogram of DC component of DCT");

%% Figure 3
% FFT of Histogram, annotated with expected f-space quantization
figure;
jk_vals = [
    1 1;
];
qs = {
    [], [179 214], [223], [28 224], [46.5 186], []
};
q_colors = ['m' 'r'];
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        y = abs(fft(h));
        % Recentre on zero
        y = fftshift(y);
        centre = floor(length(y)/2);
        plot(centre-(0:(length(y)-1)), y);
        xlim([0 1000]);

        hold on;
        for i_q = 1:length(qs{i_img})
            q = qs{i_img}(i_q);
            % Find multiples of q within the window [-centre, centre]
            % Start from fix(-centre/q) (fix = round-towards-zero)
            for n = fix(-centre/q):fix(centre/q)
                xline(n*q, q_colors(i_q));
            end
        end
        hold off;
        
%         h_sort = sort(unique(h));
%         ylim([0, h_sort(end-5)]);
        title([imgs(i_img); ...
            "Estimated f-space Qs = [" + strjoin(string(qs{i_img}), ", ") + "]"]);% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
sgtitle("FFT of Histogram of DC component of DCT - Manually Annotated");

%% Figure 
% Peaks identified for DCT FFT
figure;
jk_vals = [
    1 1;
    1 2;
];
for i_img = 1:size(imgs, 1)
    Di = imgs_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        FigH = subplot(size(imgs, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        y = abs(fft(h))./length(h);
        % Recentre on zero
        y = fftshift(y);
        centre = floor(length(y)/2);

        hold on;
        [ps,locs,~] = estimate_dc_fft_qs(h);


        plot(centre-(0:(length(y)-1)), y,'LineWidth',2);
%         for i_p = 1:length(ps)
%             xline(locs(i_p)-centre-1, 'r');
%         end
        stem(locs-centre-1, ps, 'r');


        
%         plot(locs-centre-1, ps, 'r*');
        hold off;
        
        if jk(1) == 1 && jk(2) == 1
            xlim([-100 1000]);
        end
        
%         AxesH = axes('Parent', FigH, ...
%           'Units', 'normalized', ...
%           'Position', [0, 0, 1, 1], ...
%           'Visible', 'off', ...
%           'XLim', [0, 1], ...
%           'YLim', [0, 1], ...
%           'NextPlot', 'add');
%         TextH = text(0,1, 'Top left', ...
%           'HorizontalAlignment', 'left', ...
%           'VerticalAlignment', 'top');
%         annotation(FigH, 'textbox', [0 1 0 0], 'String', 'YourString', 'FitBoxToText', true);
        xL=xlim;
        yL=ylim;
%         text(xL(2),yL(2), ...
%             img_names(i_img) + " - DCT{ " + jk(1) + "," + jk(2) + " }", ...
%             'HorizontalAlignment','right','VerticalAlignment','top')
%         title([img_names(i_img) + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
% sgtitle("FFT Peaks");


%% Figure Extra
imgs_ex = [
    "./training/converted/test00.png" "gray";
    "./training/converted/test10.png" "gray";
    "./training/converted/test13.png" "color";
    "./training/converted/test35.png" "color";
];
imgs_ex_di = {};
for i_img = 1:size(imgs_ex, 1)
    Di = image_DCT_params(imgs_ex(i_img,1), imgs_ex(i_img, 2));
    imgs_ex_di{i_img} = Di;
end
%%
figure;
jk_vals = [
    1 1;
];
img_lines = {
    [57, 170, 595];
    [64 192];
    [322];
    [293 327];
};
for i_img = 1:size(imgs_ex, 1)
    Di = imgs_ex_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        FigH = subplot(size(imgs_ex, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        y = abs(fft(h))./length(h);
        % Recentre on zero
        y = fftshift(y);
        centre = floor(length(y)/2);

        hold on;
        [ps,locs,~] = estimate_dc_fft_qs(h);


        plot(centre-(0:(length(y)-1)), y,'LineWidth',2);
        stem(locs-centre-1, ps, 'r');

        for l = img_lines{i_img}
            xline(l, '--');
            if i_img == 4 && l < 300
                text(l, max(y), l + " ", 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
            else
                text(l, max(y), " "+l, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
            end
        end
        
        hold off;
        
        if jk(1) == 1 && jk(2) == 1
            xlim([-50 850]);
        end
        
        xL=xlim;
        yL=ylim;
    end
end