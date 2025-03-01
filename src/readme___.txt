
# OpenBLAS and LAPACK
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libopenblas-dev
sudo apt-get install libarpack2-dev
sudo apt-get install libsuperlu-dev

#libarmadillo
sudo apt-get install libarmadillo-dev



cd build
cmake ..
make
./PTreeCpp | tee logout.txt
