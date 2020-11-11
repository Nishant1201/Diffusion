g++ -std=c++11 -o FTCS.out -O3 FTCS.cpp
rm FTCS_norms_Opt.out
for i in 16 32 64 128 256 512 1024 2048 4096
do ./FTCS.out $i >> FTCS_norms_Opt.out
	echo $i
done
