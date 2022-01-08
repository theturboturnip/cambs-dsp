train_data = readtable("./training/converted/stats.csv");

for i = 1:height(train_data)
    assignment4("./training/converted/" + train_data.(1)(i), train_data.(2)(i), train_data.(3)(i))
end