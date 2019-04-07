# Makefile for parallel simulated annealer

PREFIX=${PARSECDIR}/pkgs/kernels/canneal/inst/${PARSECPLAT}

TARGET=canneal
LIBS:=$(LIBS) -lm 
CXXFLAGS = -O0 
NTHREADS = 1
COMPILER_LIST = g++ "/opt/intel/bin/icc" clang++-5.0
SSEFLAGS = -mssse3

ifdef OPT
	echo $(NTHREADS)\
    CXXFLAGS := -O3
endif

ifdef THREAD
     ifeq ($(shell test $(THREAD) -ge 1; echo $$? ), 0) 
        CXXFLAGS+= -DENABLE_THREADS -pthread -std=gnu++11
        NTHREADS:= $(THREAD)
     endif
endif

ifdef COMPILER
    CXX := $(COMPILER)
endif

.PHONY: all test clean
all:
	$(CXX) $(CXXFLAGS) annealer_thread.cpp -c -o annealer_thread.o
	$(CXX) $(CXXFLAGS) rng.cpp -c -o rng.o
	$(CXX) $(CXXFLAGS) netlist.cpp -c -o netlist.o
	$(CXX) $(CXXFLAGS) main.cpp -c -o main.o
	$(CXX) $(CXXFLAGS) netlist_elem.cpp -c -o netlist_elem.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) *.o $(LIBS) -o $(TARGET)_orig
	$(CXX) $(CXXFLAGS) annealer_thread.cpp -c -o annealer_thread.o
	$(CXX) $(CXXFLAGS) rng.cpp -c -o rng.o
	$(CXX) $(CXXFLAGS) netlist.cpp -c -o netlist.o
	$(CXX) $(CXXFLAGS) main.cpp -c -o main.o
	$(CXX) $(CXXFLAGS) $(SSEFLAGS) netlist_elem_sse.cpp -c -o netlist_elem.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) *.o $(LIBS) -o $(TARGET)_ssse
	$(CXX) $(CXXFLAGS) annealer_thread.cpp -c -o annealer_thread.o
	$(CXX) $(CXXFLAGS) rng.cpp -c -o rng.o
	$(CXX) $(CXXFLAGS) netlist.cpp -c -o netlist.o
	$(CXX) $(CXXFLAGS) main.cpp -c -o main.o
	$(CXX) $(CXXFLAGS) $(SSEFLAGS) netlist_elem_opt.cpp -c -o netlist_elem.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) *.o $(LIBS) -o $(TARGET)_opt

clean:
	rm -f *.o $(TARGET)_orig $(TARGET)_ssse $(TARGET)_opt

install:
	mkdir -p $(PREFIX)/bin
	cp -f $(TARGET) $(PREFIX)/bin/$(TARGET)

test:
	@echo "ORIG Optimal 0"
	@for compiler in $(COMPILER_LIST);do \
		echo madir -p result_$${compiler}_orig; \
		mkdir -p result_$${compiler}_orig; \
		echo make COMPILER=$$compiler THREAD=1; \
		make --no-print-directory THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_orig $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_orig/result_O0_$$nthread.txt"; \
			(time ./$(TARGET)_orig $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_orig/result_O0_$$nthread.txt; \
		done \
	done
	
	@echo "ORIG Optimal 3"
	@for compiler in $(COMPILER_LIST);do \
		echo make OPT=1 COMPILER=$$compiler THREAD=1; \
		make --no-print-directory OPT=1 THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_orig $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_orig/result_O3_$$nthread.txt"; \
			(time ./$(TARGET)_orig $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_orig/result_O3_$$nthread.txt; \
		done \
	done
	
	@echo "SSSE Optimal 0"
	@for compiler in $(COMPILER_LIST);do \
		echo madir -p result_$${compiler}_ssse; \
		mkdir -p result_$${compiler}_ssse; \
		echo make COMPILER=$$compiler THREAD=1; \
		make --no-print-directory THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_ssse $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_ssse/result_O0_$$nthread.txt"; \
			(time ./$(TARGET)_ssse $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_ssse/result_O0_$$nthread.txt; \
		done \
	done
	
	@echo "SSSE Optimal 3"
	@for compiler in $(COMPILER_LIST);do \
		echo make OPT=1 COMPILER=$$compiler THREAD=1; \
		make --no-print-directory OPT=1 THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_ssse $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_ssse/result_O3_$$nthread.txt"; \
			(time ./$(TARGET)_ssse $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_ssse/result_O3_$$nthread.txt; \
		done \
	done
	
	@echo "OPT Optimal 0"
	@for compiler in $(COMPILER_LIST);do \
		echo madir -p result_$${compiler}_opt; \
		mkdir -p result_$${compiler}_opt; \
		echo make COMPILER=$$compiler THREAD=1; \
		make --no-print-directory THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_opt $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_opt/result_O0_$$nthread.txt"; \
			(time ./$(TARGET)_opt $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_opt/result_O0_$$nthread.txt; \
		done \
	done
	
	@echo "OPT Optimal 3"
	@for compiler in $(COMPILER_LIST);do \
		echo make OPT=1 COMPILER=$$compiler THREAD=1; \
		make --no-print-directory OPT=1 THREAD=1 COMPILER=$$compiler; \
		for nthread in 1 32;do \
			echo time ./$(TARGET)_opt $$nthread 15000 2000 ../inputs/2500000.nets 6000 ">> ./result_$${compiler}_opt/result_O3_$$nthread.txt"; \
			(time ./$(TARGET)_opt $$nthread 15000 2000 ../inputs/2500000.nets 6000 ;) 2>> ./result_$${compiler}_opt/result_O3_$$nthread.txt; \
		done \
	done
