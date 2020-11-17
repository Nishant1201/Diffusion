g++ -std=c++11 -o BE -O3 BE.cpp
rm BE.out
for i in 64 128 256 512 1024 2048 4096
do ./BE >> BE.out $i;
done
