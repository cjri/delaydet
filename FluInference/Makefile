CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I  /usr/local/include/gsl/ 
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
ONE		= multi_haps_3501.o utilities_haps.o get_files.o
ONEDD		= multi_haps_3501_dd.o utilities_haps_dd.o get_files.o

one: $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_haps_3501 $(LD_FLAGS)
multi_haps_3501.o: multi_haps_3501.cpp
	$(CC) $(CC_FLAGS) -c multi_haps_3501.cpp

onedd: $(ONEDD)
	$(CC) $(CC_FLAGS) $(ONEDD) -o run_haps_3501_dd $(LD_FLAGS)
multi_haps_3501_dd.o: multi_haps_3501_dd.cpp
	$(CC) $(CC_FLAGS) -c multi_haps_3501_dd.cpp
