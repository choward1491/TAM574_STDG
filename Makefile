# rule to print variables using command line
print-% : ; @echo $* = $($*)

# Main macros for compilation related objects
COMPILER:= g++
D_CFLAGS:= -O0 -g -Wall -std=c++11 -DUNIT_TEST
R_CFLAGS:= -O2 -std=c++11 -DUNIT_TEST
LFLAGS  := 
MEXEC   := sim_exec
BIN     := bin
DOBJ    := objs_d
ROBJ    := objs_r
SRC     := src
INCLUDES:= $(addprefix -I,$(SRC)) -I/usr/include/ -I/usr/local/include/

# setup source files, object files
src    := $(foreach SRC, $(SRC), $(wildcard $(SRC)/*.cpp))
objs   := $(notdir $(patsubst %.cpp, %.o, $(src)))
d_objs := $(patsubst %.o, $(DOBJ)/%.o, $(objs))
r_objs := $(patsubst %.o, $(ROBJ)/%.o, $(objs))

# setup vpath
VPATH = $(SRC)

# compile shit
all: build

rebuild: clean build

build: debug release

debug: $(d_objs)
	mkdir -p $(BIN)
	$(COMPILER) $^ -o $(BIN)/$(MEXEC)_d $(LFLAGS)

release: $(r_objs)
	mkdir -p $(BIN)
	$(COMPILER) $^ -o $(BIN)/$(MEXEC) $(LFLAGS)

$(DOBJ)/%.o: %.cpp
	mkdir -p $(DOBJ)
	$(COMPILER) $(D_CFLAGS) $(INCLUDES) -c $< -o $@

$(ROBJ)/%.o: %.cpp
	mkdir -p $(ROBJ)
	$(COMPILER) $(R_CFLAGS) $(INCLUDES) -c $< -o $@

clean: 
	rm -rf $(DOBJ) $(ROBJ) $(BIN)/$(MEXEC)*
