#
# Makefile for non-Microsoft compilers 
# c++11
#

CC=g++
CFLAGS=-std=c++11

all: es es_4 ews

es: es.cpp graph_es.h
	${CC} ${CFLAGS} es.cpp -o es

es_4: es_circle4.cpp
	${CC} ${CFLAGS} es_circle4.cpp -o es_4

ews: ews.cpp graph_ews.h
	${CC} ${CFLAGS} ews.cpp -o ews

clean:
	\rm *.o  es es_4 ews
