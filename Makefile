CC = g++
CXX = g++

TARGET_EXEC ?= build_geometry.out

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name "*.s" -or -name "*.cc" ! -name "*test*")
TEST_SRCS := $(shell find $(SRC_DIRS) -name "*.s" -or -name "*.cc" ! -name "main.cc")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
TEST_OBJS := $(TEST_SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

TESTLDFLAGS ?= -l gtest -l pthread

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))


CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -std=c++11

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cc.o: %.cc
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

all: $(BUILD_DIR)/$(TARGET_EXEC)

test: $(TEST_OBJS)
	$(CC) $(TEST_OBJS) -o $(BUILD_DIR)/$@ $(TESTLDFLAGS)
	./$(BUILD_DIR)/$@

-include $(DEPS)

MKDIR_P ?= mkdir -p
