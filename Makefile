CC = g++ -fdiagnostics-color=always
FLAGS = -std=c++17 -g -pg -O3 -Wall -I.
LIBS = -pthread

BUILD = build
OBJ = build/obj

GRAPH = graph
MATCH = matching
UTILS = utils
MODEL = model
DMS = DecisionMakingSystem
BUILD_TOOLS = build/tools

all : dir $(BUILD)/csm

dir: $(OBJ)

$(OBJ) :
	mkdir -p $(OBJ)

#################### start ####################

$(BUILD)/csm: $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/csmpp.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o $(OBJ)/logistic_regression.o \
		$(OBJ)/common.o $(OBJ)/DecisionMakingSystem.o
	$(CC) $(FLAGS) $(OBJ)/main.o \
		$(OBJ)/matching.o \
		$(OBJ)/csmpp.o \
		$(OBJ)/graph.o $(OBJ)/induced_graph.o \
		$(OBJ)/globals.o $(OBJ)/logistic_regression.o \
		$(OBJ)/common.o $(OBJ)/DecisionMakingSystem.o \
		-o $(BUILD)/csm $(LIBS)

$(OBJ)/main.o: $(MATCH)/main.cpp \
		$(UTILS)/CLI11.hpp \
		$(UTILS)/globals.h $(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/csmpp.h
	$(CC) -c $(FLAGS) $(MATCH)/main.cpp -o $(OBJ)/main.o

#################### matching ####################
$(OBJ)/csmpp.o: $(MATCH)/csmpp.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(UTILS)/globals.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h \
		$(MATCH)/csmpp.h \
		$(DMS)/DecisionMakingSystem.h
	$(CC) -c $(FLAGS) $(MATCH)/csmpp.cpp \
		-o $(OBJ)/csmpp.o

$(OBJ)/matching.o: $(MATCH)/matching.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/graph.h \
		$(MATCH)/matching.h
	$(CC) -c $(FLAGS) $(MATCH)/matching.cpp \
		-o $(OBJ)/matching.o

#################### graph ####################

$(OBJ)/graph.o: $(GRAPH)/graph.cpp \
		$(UTILS)/types.h $(UTILS)/utils.h \
		$(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(GRAPH)/graph.cpp \
		-o $(OBJ)/graph.o

$(OBJ)/induced_graph.o: $(GRAPH)/induced_graph.cpp \
		$(UTILS)/types.h \
		$(GRAPH)/induced_graph.h $(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(GRAPH)/induced_graph.cpp \
		-o $(OBJ)/induced_graph.o

#################### utils ####################

$(OBJ)/globals.o: $(UTILS)/globals.cpp $(UTILS)/globals.h
	$(CC) -c $(FLAGS) $(UTILS)/globals.cpp \
		-o $(OBJ)/globals.o

#################### DMS ####################

$(OBJ)/DecisionMakingSystem.o: $(DMS)/DecisionMakingSystem.cpp $(DMS)/DecisionMakingSystem.h $(MODEL)/logistic_regression.h $(UTILS)/types.h $(UTILS)/globals.h $(GRAPH)/graph.h
	$(CC) -c $(FLAGS) $(DMS)/DecisionMakingSystem.cpp \
		-o $(OBJ)/DecisionMakingSystem.o

#################### model ####################

$(OBJ)/logistic_regression.o: $(MODEL)/logistic_regression.cpp $(MODEL)/logistic_regression.h $(MODEL)/common.h
	$(CC) -c $(FLAGS) $(MODEL)/logistic_regression.cpp \
		-o $(OBJ)/logistic_regression.o

$(OBJ)/common.o: $(MODEL)/common.cpp $(MODEL)/common.h
	$(CC) -c $(FLAGS) $(MODEL)/common.cpp \
		-o $(OBJ)/common.o

#################### end ####################

.PHONY: clean

clean: 
	rm -r ${BUILD}
