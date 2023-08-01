geo_path=https://zenodo.org/record/3937811
wget -O file.txt $geo_path
grep '<a href=' file.txt | sed -n 's/.*href="\([^"]*\)".*/\1/p' |sed -e '1d; $d' > test.txt
