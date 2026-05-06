CC = gcc
#CFLAGS = -O2 -std=c11 -Wall -Wextra -pedantic -march=native -lm
CFLAGS = -O2 -std=c11 -Wall -Wextra -pedantic -march=native
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

$(TARGET): main.c
	$(CC) $(CFLAGS) main.c -o $(TARGET) -lm

clean:
	rm -f $(TARGET) $(DATA) $(GRAPH)
