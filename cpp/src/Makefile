TARGET_EXEC = FEM
CXX = g++

#------------------------------------------------
EIGEN_INCLUDE_PATH = /opt/homebrew/include/eigen3
# EIGEN_INCLUDE_PATH = /usr/local/include/eigen3
# EIGEN_INCLUDE_PATH = /usr/include/eigen3

#------------------------------------------------
CXXFLAGS = -std=c++17 -Wall --pedantic-errors -O2
CPPFLAGS = -I$(EIGEN_INCLUDE_PATH)

#------------------------------------------------
BUILD_DIR = obj

SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(SRCS:%=$(BUILD_DIR)/%.d)

$(TARGET_EXEC): $(OBJS) all.h.gch
	$(CXX) $(CXXFLAGS) -include all.h $(OBJS) -o $@

-include $(DEPS)

all.h.gch: all.h
	$(CXX) $(CXXFLAGS) -x c++-header -o $@ $<

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -MMD -MP $< -o $@ 

clean:
	rm -f $(TARGET_EXEC)
	rm -f ./all.h.gch
	rm -rf $(BUILD_DIR)

.PHONY: clean

MKDIR_P = mkdir -pv
