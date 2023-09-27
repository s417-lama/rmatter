CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++20
TARGETS = rmatter converter

all: $(TARGETS)

%: %.cpp common.hpp
	$(CXX) $(CXXFLAGS) -o $@ -I./pcg-cpp/include $< -lpthread

clean:
	rm -f $(TARGET)
