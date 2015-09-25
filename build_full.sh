
#format files
find . -name *.cpp -o -name *.h | xargs clang-format-3.4 -style=file -i

#create build folder
rm -rf build
mkdir -p build && cd build

#build
cmake ..
sudo make -j4

#covering related
#clean up
lcov --zerocounters --directory ..
#run initial
lcov --capture --initial --directory .. --output-file cover
#run tests
make test
#run final
lcov --no-checksum --directory .. --capture --output-file cover.info 
#remove external dependencies
lcov --remove cover.info 'usr/*' --output-file cover.info 
#remove test dependencies
lcov --remove cover.info '*test*' --output-file cover.info 
#generate HTML report
genhtml --highlight --legend --show-details -o html_coverage cover.info

#clean
rm -f cover cover.info
make clean
