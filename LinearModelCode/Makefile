CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/
LD_FLAGS        = -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
SE_OBJECTS	= simulation_run_exp.o utilities_comp.o
ID_OBJECTS	= inference_det.o utilities_comp.o
IDD_OBJECTS	= inference_delaydet.o utilities_comp.o
IS_OBJECTS	= inference_stoch.o utilities_comp.o


sim_e:  $(SE_OBJECTS)
	$(CC) $(CC_FLAGS) $(SE_OBJECTS) -o comp_sim_exp $(LD_FLAGS)
simulation_run_exp.o: simulation_run_exp.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp.cpp
utilities_comp.o: utilities_comp.cpp
	$(CC) $(CC_FLAGS) -c utilities_comp.cpp

infd:	$(ID_OBJECTS)
	$(CC) $(CC_FLAGS) $(ID_OBJECTS) -o comp_inf_det $(LD_FLAGS)
inference_det.o: inference_det.cpp
	$(CC) $(CC_FLAGS) -c inference_det.cpp

infdd:	$(IDD_OBJECTS)
	$(CC) $(CC_FLAGS) $(IDD_OBJECTS) -o comp_inf_ddd $(LD_FLAGS)
inference_delaydet.o: inference_delaydet.cpp
	$(CC) $(CC_FLAGS) -c inference_delaydet.cpp

infs:	$(IS_OBJECTS)
	$(CC) $(CC_FLAGS) $(IS_OBJECTS) -o comp_inf_stoch $(LD_FLAGS)
inference_stoch.o: inference_stoch.cpp
	$(CC) $(CC_FLAGS) -c inference_stoch.cpp

