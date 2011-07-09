SYSTEM     = x86-64_RHEL3.0_3.2

COPT =  -Wall -Wno-sign-compare -O2 -fPIC -fexceptions 
#COPT =  -g -Wall -Wno-sign-compare -fPIC -fexceptions 

# Create a list of the object files using macro substitutions
OBJECTS = $(SOURCES:.cpp=.o)

# redefine inbuilt macros and define new ones
CC    = g++
CFLAGS = $(COPT) 


#.SUFFIXES:       # remove inbuilt definitions
#.SUFFIXES: .cpp .o # define default compilation procedure
#.cpp.o:            # .o files depend on .c files
#	$(CC) $(CFLAGS) $*.c # start with a tab character!!

# Default target
all: multiply splinter preproc getmers checkss center cluster countCenters filtern
#	@echo "compilation done"

filtern: filtern.o
	$(CC) filtern.o -o filtern
filtern.o : filtern.cpp
	$(CC) -c $(COPT) filtern.cpp -o filtern.o
multiply: multiply.o
	$(CC) multiply.o -o multiply
multiply.o : multiply.cpp
	$(CC) -c $(COPT) multiply.cpp -o multiply.o
splinter: splinter.o
	$(CC) splinter.o -o splinter
splinter.o : splinter.cpp
	$(CC) -c $(COPT) splinter.cpp -o splinter.o
preproc: preproc.o
	$(CC) preproc.o -o preproc
preproc.o : preproc.cpp
	$(CC) -c $(COPT) preproc.cpp -o preproc.o
getmers: getmers.o
	$(CC) getmers.o -o getmers
getmers.o : getmers.cpp
	$(CC) -c $(COPT) getmers.cpp -o getmers.o
checkss: checkss.o
	$(CC) checkss.o -o checkss
checkss.o : checkss.cpp
	$(CC) -c $(COPT) checkss.cpp -o checkss.o
cluster: cluster.o
	$(CC) cluster.o -o cluster  -lpthread
cluster.o : cluster.cpp
	$(CC) -c $(COPT) cluster.cpp -o cluster.o
center: center.o
	$(CC) center.o -o center -lpthread
center.o : center.cpp
	$(CC) -c $(COPT) center.cpp -o center.o
countCenters: countCenters.o
	$(CC) countCenters.o -o countCenters -lpthread
countCenters.o : countCenters.cpp
	$(CC) -c $(COPT) countCenters.cpp -o countCenters.o

# Target deleting unwanted files
clean:
	rm -f *.o *~ core mppcore
