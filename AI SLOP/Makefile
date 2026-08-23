CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall
SRC      = cpp/main.cpp cpp/initial_conditions.cpp cpp/solver.cpp
HDR      = cpp/params.hpp cpp/potential.hpp cpp/fields.hpp cpp/initial_conditions.hpp cpp/solver.hpp
TARGET   = bin/oscillon
DATA_DIR = output/data
IMAGE_DIR = output/images
VIDEO_DIR = output/video

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	mkdir -p $(DATA_DIR) $(IMAGE_DIR) $(VIDEO_DIR)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(DATA_DIR)/*.dat $(DATA_DIR)/run_info.txt
	rm -f $(IMAGE_DIR)/*.png $(VIDEO_DIR)/*.mp4 $(VIDEO_DIR)/*.gif

.PHONY: all run clean
