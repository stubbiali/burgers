FF = gfortran
FFLAGS = -O4
FLINK =

SRC_DIR = src
SRC = spacetime.f90 boundary.f90 utils.f90 main.f90
MOD_DIR = modules
MOD = $(patsubst %.f90, $(MOD_DIR)/%.mod, $(SRC))
OBJ_DIR = objects
OBJ = $(patsubst %.f90, $(OBJ_DIR)/%.o, $(SRC))
MAIN = main


all: main


main: $(OBJ)
	$(FF) $(FFLAGS) -o $(MAIN) $(OBJ) $(FLINK)


$(OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(MOD_DIR)
	@mkdir -p $(OBJ_DIR)
	$(FF) $(FFLAGS)	-c $< 	-o $@ 	-J $(MOD_DIR)


.PHONY: clean distclean


clean:
	rm -rf $(MOD_DIR)


distclean: clean
	rm -rf $(OBJ_DIR) $(MAIN)