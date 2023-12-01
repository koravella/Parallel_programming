BINS := $(patsubst %.cpp, %, $(wildcard ./programs/*.cpp)) check
LINEAR := linear.cpp
DEPS := linear.h
LIBS := -fopenmp


CPPFLAGS = -I./

.PHONY: all clean

all: $(addsuffix .out, $(BINS))

%.out: %.cpp $(LINEAR) $(DEPS)
	$(CXX) $(LINEAR) $< -o $@ $(LIBS) $(CPPFLAGS)

clean:
	rm $(addsuffix .out, $(BINS))
