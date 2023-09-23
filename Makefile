CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++17
TARGET = rmatter
SRC = rmatter.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) -I./pcg-cpp/include $(SRC) -lpthread

clean:
	rm -f $(TARGET)
