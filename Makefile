CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -I src -I solvers/cadical/src -I toml/single_include
LDFLAGS = -L solvers/cadical/build -lcadical

CADICAL_LIB = solvers/cadical/build/libcadical.a

# Main binary
MAIN_SRC = src/main.cpp
MAIN_BIN = catalyst-search

# Test binaries
TEST_SRCS = $(wildcard test/test_*.cpp)
TEST_BINS = $(TEST_SRCS:.cpp=)

.PHONY: all tests clean run-tests cadical

all: $(MAIN_BIN)

cadical: $(CADICAL_LIB)

$(CADICAL_LIB):
	cd solvers/cadical && ./configure && make

$(MAIN_BIN): $(MAIN_SRC) src/*.hpp $(CADICAL_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

tests: $(TEST_BINS)

test/%: test/%.cpp src/*.hpp $(CADICAL_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run-tests: tests
	@echo "Running all tests..."
	@for test in $(TEST_BINS); do \
		echo "\n=== Running $$test ==="; \
		./$$test || exit 1; \
	done
	@echo "\n=== All tests passed ==="

clean:
	rm -f $(MAIN_BIN) $(TEST_BINS)
