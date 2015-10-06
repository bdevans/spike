program_NAME := Spike
program_C_SRCS := $(filter-out ._*, $(wildcard *.c))	 #$(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_SRCS := $(program_C_SRCS) $(program_CXX_SRCS)
program_HDRS := $(wildcard *.h)
program_ANALYSIS := $(wildcard *.m)
program_DOCUMENTOR := doxygen
program_DOC_CONFIG := Doxyfile_$(program_NAME)
program_DOC_GEN := $(program_DOCUMENTOR) $(program_DOC_CONFIG)
program_ADDITIONAL = $(MAKEFILE_LIST) $(program_DOC_CONFIG)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o} #No C++ files
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := /opt/local/include 	# $(shell gsl-config --prefix)/include
program_LIBRARY_DIRS := /opt/local/lib		# $(shell gsl-config --prefix)/lib
program_LIBRARIES := m gomp gsl
program_ARCHIVE := $(program_NAME).tbz
program_REPO_UPDATE := svn update
program_REPO_CHECKOUT := svn checkout
program_REPO_COMMIT := svn commit
program_REPO_SRCS := https://evans@mac0.cns.ox.ac.uk:443/svn/SpikeNet
program_TEST_ARGS := -f defaults.m
# program_LIBRARIES := -lm -lgomp -lgsl

# target: dependencies
# [tab] system command
# *** Line endings must be UNIX (LF) ***

# $@	The file name of the target.
# $<	The name of the first dependency.
# $*	The part of a filename which matched a suffix rule.
# $?	The names of all the dependencies newer than the target separated by spaces.
# $^    The names of all the dependencies separated by spaces, but with duplicate names removed.
# $+    The names of all the dependencies separated by spaces with duplicate names included and in the same order as in the rule.

# \ 	Continuation character so commands run in same shell
# =		Recursively expanded
# :=	Simply expanded

# Boilerplate in case the standard $(MAKE) variable is not defined
#ifeq ($(MAKE),)
#	MAKE := make
#endif
MAKE ?= make # N.B. Does nothing if MAKE is defined but empty

OPT_FLAGS := -O3 -xT
DBG_FLAGS := -Wall -D DEBUG=3
PRO_FLAGS := -pg # Profiling flags

CFLAGS += -fopenmp

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

# These targets should execute their build rules, even if a newer file with a name matching the target exists.
.PHONY: all clean distclean newdoc doc tar project untar hash inspect

all: $(program_NAME) # This should be the first target so that 'make' is equivalent to 'make all'

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) $(OPT_FLAGS) -o $(program_NAME)

debug: $(program_NAME)	# private?
	$(OPT_FLAGS) := $(DBG_FLAGS)
	
profile: $(program_NAME)
	$(OPT_FLAGS) := $(PRO_FLAGS)

install:
	# Copy all relevant binaries and parameter files to a suitable working directory. Default: '/usr/local/' if DESTDIR is not provided

list:
	cd .; \
	ls; \
	echo $(program_C_SRCS); \
	echo $(program_CXX_SRCS); \
	echo $(program_HDRS); \
	echo $(program_ANALYSIS);

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	find . -name ._\* -exec rm -f {} \;

# Using the pipe `|' means : normal-prerequisites | order-only-prerequisites 
fresh : | clean clearscr all

clearscr:
	clear

distclean: clean	# Get source from svn e.g. sccs get $@
	$(program_REPO_CHECKOUT) $(program_REPO_SRCS)

update: clean
	$(program_REPO_UPDATE) $(program_REPO_SRCS)

commit: clean
	# Add svn commit changes here / change to git?

newdoc: $(program_SRCS)	# Generate a new config file with doxygen -g <Doxyfile>
	$(program_DOCUMENTOR) -g $(program_DOC_CONFIG)
	# Replace tags in the config file 
	# string='PATTERN_TO_MATCH         = REPLACE_ME '; printf '<%s>\n' "${string%=*}= REPLACE_ME"
	# PROJECT_NAME
	# OPTIMIZE_OUTPUT_FOR_C
	make doc # This uses a second instance of make http://stackoverflow.com/questions/3267145/makefile-execute-another-target

doc: $(program_SRCS)
	$(program_DOC_GEN)

tar: clean $(program_ARCHIVE) hash

project : program_ADDITIONAL += $(program_NAME).xcodeproj

project : tar
#	rm -f $(program_ARCHIVE)

$(program_ARCHIVE) : $(program_C_SRCS) $(program_HDRS) $(program_ANALYSIS) $(program_ADDITIONAL)
	tar -cjf $(program_ARCHIVE) $(program_C_SRCS) $(program_HDRS) $(program_ANALYSIS) $(program_ADDITIONAL)

### http://stackoverflow.com/questions/2148892/conditionally-appending-to-a-variable-inside-a-makefile-target
LIST = item1

targetMain: 
# DO all the work in here
	echo $(LIST)

targetA: LIST+=itemA
targetB: LIST+=itemB

targetA targetB: targetMain
###

#TAR := tar -cjf $(program_NAME).tbz $(program_C_SRCS) $(program_HDRS) $(program_ANALYSIS) $(program_ADDITIONAL) #$(MAKEFILE_LIST)

#tar: $(program_NAME).tbz
#	@- $(RM) $(program_NAME).tbz
#	@- $(RM) $(program_OBJS)
#	find . -name ._\* -exec rm -f {} \;
#	$(TAR)

#project: $(program_ADDITIONAL) += $(program_NAME).xcodeproj
#project: tar

#project : clean tar # $(program_NAME).tbz
#	find . -name ._\* -exec rm -f {} \;
#	$(TAR) $(program_NAME).xcodeproj
# 	$(program_ADDITIONAL) += $(program_NAME).xcodeproj
#	tar -cjf $(program_NAME).tbz $(program_C_SRCS) $(program_HDRS) $(program_ANALYSIS) $(MAKEFILE_LIST) $(program_NAME).xcodeproj
# Could just build original archive then include the xcodeproj directory

untar:     # $(program_NAME).tbz
	tar -xvf $(program_ARCHIVE)

hash: # SHA algorithms: 1 (default), 224, 256, 384, 512
	md5 $(program_ARCHIVE)
	shasum -a 512 $(program_ARCHIVE)

inspect:	# Inspect program for dynamically linked libraries
	otool -L ./$(program_NAME) # ldd ./$(program_NAME) # for linux

checklib:	# Test for the presence of GSL libraries on the system
	gsl-config --prefix --version

test:	# Run test(s) on compiled binary
	./$(program_NAME) $(program_TEST_ARGS)

gridtest:
	# Submit job(s) to Xgrid wrapped in a bash script to check it works with the sandbox

# www.cs.duke.edu/~ola/courses/programming/Makefiles/node6.html
PRINTCMD = enscript -2rG -pprintout.ps

print: $(program_SRCS) $(program_HDRS)
	echo Spooling $? to printer using $(PRINTCMD)
	$(PRINTCMD) $?
	print printout.ps
	rm printout.ps

# To see the rules and macros being used type: make -p
help:
	$(MAKE) -v
	$(MAKE) -p

# The following lines generate a list of dependencies for header files but only works with C not C++
depend:
	makedepend -- $(CPPFLAGS) $(CFLAGS) -- $(program_C_SRCS)

# Don't place anything below this line, since
# the make depend program will overwrite it
# DO NOT DELETE THIS LINE -- make depend depends on it.
