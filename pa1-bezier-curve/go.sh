rm bezier

echo "Compiling..."
g++ bezier.cpp -o bezier -O3
echo ""
echo "Testing with test.in"
./bezier test.in output.txt
python check.py output.txt test.out
echo ""
echo "Testing with test_5x5.in"
./bezier test_5x5.in output.txt
python check.py output.txt test_5x5.out
echo ""
echo "Testing with test_7x7.in"
./bezier test_7x7.in output.txt
python check.py output.txt test_7x7.out
