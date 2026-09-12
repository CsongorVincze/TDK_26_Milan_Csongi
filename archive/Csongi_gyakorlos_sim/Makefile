CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic
TARGET = oscillon

all: $(TARGET)

$(TARGET): solver.cpp
	$(CXX) $(CXXFLAGS) solver.cpp -o $(TARGET)

run: $(TARGET)
	./$(TARGET) single

sweep: $(TARGET)
	./$(TARGET) sweep

animate: $(TARGET)
	python3 animate.py

clean:
	rm -f $(TARGET) output/snapshots.dat output/run_info.txt
	rm -rf output/sweep output/video

.PHONY: all run sweep animate clean
