function assignment4(name,color_type,expected_value)
    c = evaluate_compression(name,color_type);

    if c == expected_value
        fprintf("%10s:  \t CORRECT \t%d\n", name, c)
    else
        fprintf("%10s:  \tINCORRECT\t%d != actual %d\n", name, c, expected_value)
    end
end
function c=evaluate_compression(name, color_type)
    Di = image_DCT_params(name, color_type);
    jk = [1 1];
    d = Di{jk(1), jk(2)};
    h = histcounts(abs(d),'BinMethod','Integer');
    qs = identify_qs(h);
    c = length(qs);
end