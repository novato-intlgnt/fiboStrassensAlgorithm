CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic -march=native
TARGET = main
DATA = results.dat
SCRIPT = plot.gp
GRAPH = ./proofs/graph.png

all: $(GRAPH)

$(GRAPH): $(DATA) $(SCRIPT)
	gnuplot $(SCRIPT)
	@echo "Gráfico generado: $(GRAPH)"

$(DATA): $(TARGET)
	./$(TARGET) > $(DATA)

$(TARGET): main.cc
	$(CXX) $(CXXFLAGS) main.cc -o $(TARGET)

clean:
	rm -f $(TARGET) $(DATA) $(GRAPH)
