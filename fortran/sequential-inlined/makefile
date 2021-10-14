SRC_DIR = src
SRC = spacetime.f90 boundary.f90 utils.f90 main.f90
MOD_DIR = modules
MOD = $(patsubst %.f90, $(MOD_DIR)/%.mod, $(SRC))
OBJ_DIR = objects
OBJ = $(patsubst %.f90, $(OBJ_DIR)/%.o, $(SRC))
MAIN = main

FF = ftn
FFLAGS = -O3 -finline-functions -J $(MOD_DIR)
FLINK =


all: main


main: $(OBJ)
	$(FF) $(FFLAGS)	-o $(MAIN) 	$(OBJ) $(FLINK)


$(OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(MOD_DIR)
	@mkdir -p $(OBJ_DIR)
	$(FF) $(FFLAGS)	-o $@ 	-c $< 	$(FLINK)


.PHONY: clean distclean


clean:
	$(RM) -r $(MOD_DIR)


distclean: clean
	$(RM) -r $(OBJ_DIR) $(MAIN)