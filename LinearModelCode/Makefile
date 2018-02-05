CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/
LD_FLAGS        = -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
SE_OBJECTS	= simulation_run_exp.o utilities_comp.o
ST_OBJECTS	= simulation_run_exp_time.o utilities_comp.o
ID_OBJECTS	= inference_det.o utilities_comp.o
IDD_OBJECTS	= inference_delaydet.o utilities_comp.o
IDR_OBJECTS	= inference_det_rapid.o utilities_comp.o
IDDR_OBJECTS	= inference_delaydet_rapid.o utilities_comp.o
IS_OBJECTS	= inference_stoch.o utilities_comp.o
GR_OBJECTS	= simulation_run_exp_grow.o utilities_comp.o
DO_OBJECTS	= simulation_run_exp_down.o utilities_comp.o
LD_OBJECTS	= inference_det_land.o utilities_land.o
LDD_OBJECTS	= inference_ddet_land.o utilities_land.o

sim_e:  $(SE_OBJECTS)
	$(CC) $(CC_FLAGS) $(SE_OBJECTS) -o comp_sim_exp $(LD_FLAGS)
simulation_run_exp.o: simulation_run_exp.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp.cpp
utilities_comp.o: utilities_comp.cpp
	$(CC) $(CC_FLAGS) -c utilities_comp.cpp

grow: 	$(GR_OBJECTS)
	$(CC) $(CC_FLAGS) $(GR_OBJECTS) -o comp_sim_grow $(LD_FLAGS)
simulation_run_exp_grow.o: simulation_run_exp_grow.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_grow.cpp

down:   $(DO_OBJECTS)
	$(CC) $(CC_FLAGS) $(DO_OBJECTS) -o comp_sim_down $(LD_FLAGS)
simulation_run_exp_down.o: simulation_run_exp_down.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_down.cpp


infd:	$(ID_OBJECTS)
	$(CC) $(CC_FLAGS) $(ID_OBJECTS) -o comp_inf_det $(LD_FLAGS)
inference_det.o: inference_det.cpp
	$(CC) $(CC_FLAGS) -c inference_det.cpp

infdr: 	$(IDR_OBJECTS)
	$(CC) $(CC_FLAGS) $(IDR_OBJECTS) -o comp_inf_det_rap $(LD_FLAGS)
inference_det_rapid.o: inference_det_rapid.cpp
	$(CC) $(CC_FLAGS) -c inference_det_rapid.cpp

infdd:	$(IDD_OBJECTS)
	$(CC) $(CC_FLAGS) $(IDD_OBJECTS) -o comp_inf_ddd $(LD_FLAGS)
inference_delaydet.o: inference_delaydet.cpp
	$(CC) $(CC_FLAGS) -c inference_delaydet.cpp

infddr: $(IDDR_OBJECTS)
	$(CC) $(CC_FLAGS) $(IDDR_OBJECTS) -o comp_inf_ddr_rap $(LD_FLAGS)
inference_delaydet_rapid.o: inference_delaydet_rapid.cpp
	$(CC) $(CC_FLAGS) -c inference_delaydet_rapid.cpp

infdl:  $(LD_OBJECTS)
	$(CC) $(CC_FLAGS) $(LD_OBJECTS) -o comp_inf_det_land $(LD_FLAGS)
inference_det_land.o: inference_det_land.cpp
	$(CC) $(CC_FLAGS) -c inference_det_land.cpp
utilities_land.o: utilities_land.cpp
	$(CC) $(CC_FLAGS) -c utilities_land.cpp

infddl: $(LDD_OBJECTS)
	$(CC) $(CC_FLAGS) $(LDD_OBJECTS) -o comp_inf_ddet_land $(LD_FLAGS)
inference_ddet_land.o: inference_ddet_land.cpp
	$(CC) $(CC_FLAGS) -c inference_ddet_land.cpp

infs:	$(IS_OBJECTS)
	$(CC) $(CC_FLAGS) $(IS_OBJECTS) -o comp_inf_stoch $(LD_FLAGS)
inference_stoch.o: inference_stoch.cpp
	$(CC) $(CC_FLAGS) -c inference_stoch.cpp

times:  $(ST_OBJECTS)
	$(CC) $(CC_FLAGS) $(ST_OBJECTS) -o comp_sim_exp_time $(LD_FLAGS)
simulation_run_exp_time.o: simulation_run_exp_time.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_exp_time.cpp

