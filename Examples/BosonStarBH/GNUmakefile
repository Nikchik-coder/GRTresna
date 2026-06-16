# -*- Mode: Makefile -*-

ifndef CHOMBO_HOME
    $(error Please define CHOMBO_HOME - see installation instructions.)
endif

makefiles += releasedExamples_AMRPoisson_execVariableCoefficient

cxxcppflags = -DUSE_CTTKHybrid -DUSE_COMPLEX_SCALAR_MATTER

ebase = Main_BosonStarBH

LibNames = AMRElliptic AMRTools BoxTools

INPUT = params_boson_star_test.txt

src_dirs := ../../Source \
            ../../Source/Core \
            ../../Source/Matter \
            ../../Source/Methods \
            ../../Source/Tools \
            ../../Source/Variables \
            ../../Source/TaggingCriteria \
            ../../Source/Operator \
            ../../Source/Operator/SolverOperator

include $(CHOMBO_HOME)/mk/Make.example
