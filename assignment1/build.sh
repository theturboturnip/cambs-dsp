mkdir -p ./build/
rm -r ./build/*
cp *.csv ./build/
cat ex_1a.m ex_1b.m ex_1c.m ex_1d.m > ./build/ex_complete.m
cd ./build/
echo "Building LaTeX from MATLAB..."
matlab -batch "publish('./ex_complete.m', 'latex'); exit"
echo "...done!"
echo "Rendering LaTeX..."
cd ./html/
pdflatex ./ex_complete.tex
pdflatex ./ex_complete.tex
pdflatex ./ex_complete.tex
echo "...done!"
