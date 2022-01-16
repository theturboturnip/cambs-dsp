%% Data
imgs = [
    "test1.png" "gray";
    "test1c.png" "color";
    "test2.png" "gray";
    "test2c.png" "color";
    "test3.png" "gray";
    "test3c.png" "color";


%     "./training/converted/test08.png" "gray"; % 0x
%     "./training/converted/test06.png" "gray"; % 1x
%     "./training/converted/test02.png" "gray"; % 2x
%     "./training/converted/test11.png" "color"; % 0x
%     "./training/converted/test17.png" "color"; % 1x
%     "./training/converted/test13.png" "color"; % 2x
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
    1 1;
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
        xlim([0 2000]);
        title(imgs(i_img));% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
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

%% Figure 4
imgs_ex = [


    "./training/converted/test08.png" "gray"; % 0x
    "./training/converted/test38.png" "gray"; % 0x
    "./training/converted/test11.png" "color"; % 0x
    "./training/converted/test29.png" "color"; % 0x

    "./training/converted/test06.png" "gray"; % 1x
    "./training/converted/test26.png" "gray"; % 1x
    "./training/converted/test33.png" "color"; % 1x
    "./training/converted/test17.png" "color"; % 1x
    
    "./training/converted/test02.png" "gray"; % 2x
    "./training/converted/test24.png" "gray"; % 2x
    "./training/converted/test35.png" "color"; % 2x
    "./training/converted/test37.png" "color"; % 2x
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
    2 2;
    1 2;
    1 3;
    2 1;
    3 1;
];
for i_img = 1:size(imgs_ex, 1)
    Di = imgs_ex_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs_ex, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        histogram(d, 'BinMethod', 'Integer');
        
%         h_sort = sort(unique(h));
%         ylim([0, h_sort(end-5)]);
%         xlim([0 2000]);
        title(imgs_ex(i_img));% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
sgtitle("Histogram of DC component of DCT");

figure;
jk_vals = [
    1 1;
    2 2;
    1 2;
    1 3;
    2 1;
    3 1;
];
for i_img = 1:size(imgs_ex, 1)
    Di = imgs_ex_di{i_img};
    for i_dct = 1:size(jk_vals,1)
        jk = jk_vals(i_dct, :);
        d = Di{jk(1), jk(2)};
        h = histcounts(abs(d),'BinMethod','Integer');
        subplot(size(imgs_ex, 1), size(jk_vals,1), (i_img-1)*size(jk_vals,1) + i_dct)

        y = abs(fft(h));
        y = fftshift(y);
        centre = floor(length(y)/2);
        plot(centre-(0:(length(y)-1)), y);
%         xlim([0 1000]);
        
%         h_sort = sort(unique(h));
%         ylim([0, h_sort(end-5)]);
        title(imgs_ex(i_img, 1));% + " - DCT{ " + jk(1) + "," + jk(2) + " }"]);
    end
end
sgtitle("FFT of Histogram of DC component of DCT");