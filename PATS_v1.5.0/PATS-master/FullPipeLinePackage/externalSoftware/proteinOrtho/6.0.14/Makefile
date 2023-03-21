# This makefile is part of Proteinortho.
##########################################

# Run 'make' for compiling everything in the current directory (using the installed version of lapack in e.g. /usr/lib/, you can install lapack with e.g. apt-get install libatlas3-base or liblapack3)
# Run 'make STATIC=TRUE' for a static version
# Run 'make USELAPACK=FALSE' for a version without(!) LAPACK (only power iteration is used)
# Run 'make USEPRECOMPILEDLAPACK=FALSE' for directly recompiling the provided lapack version 3.8.0 and linking dynamically

# Run 'make CXX=g++-7' for using the g++-7 compiler. See Flags below for more informations
# Run 'make CXX=clang++' for using the clang compiler

# Run 'make install' for installing the compiled files to /usr/local/bin
# Run 'make install PREFIX=/home/paul/bin/' for local installation

############ OPTIONS: ##########################################
## STATIC=TRUE : enable static compiling (default:FALSE) 
## CXX=g++ : the g++ compiler
## CXXFLAGS = compiler flags passed to g++
## CXXLIBRARY = the path to the libs like lapack,... (dont forget the -L)
## CXXINCLUDE = include path (^) (dont forget the -I)
## PREFIX = the installation prefix (only for make install)
################################################################

##########################
## enviroment variables ##
##########################

# ALIAS for PREFIX
INSTALLDIR=/usr/local/bin
# ALIAS for LAPACK
USELAPACK=TRUE
# compile statically
STATIC=FALSE

ifdef static
STATIC=$(static)
endif
ifdef PREFIX
INSTALLDIR=$(PREFIX)
endif
ifdef prefix
INSTALLDIR=$(prefix)
endif
ifdef installdir
INSTALLDIR=$(installdir)
endif
ifdef LAPACK
USELAPACK=$(LAPACK)
endif
ifeq ($(STATIC),true)
STATIC=TRUE
endif

USEPRECOMPILEDLAPACK=TRUE

UNAME_S=$(shell uname -s)_$(shell uname -m)
# output dir of make (make install moves these to PREFIX)
BUILDDIR=src/BUILD/$(UNAME_S)

CC=cc
CXX=g++

CXXFLAGS_PO=-Wall -O3 -std=c++11 -Wno-unused-result

IS_COLOR_COMPATIBLE:=$(shell tput color 2>/dev/null)
ifdef IS_COLOR_COMPATIBLE
RED=\033[1;31m
GREEN=\033[1;32m
ORANGE=\033[1;33m
NC=\033[0m
endif

##############
# MAKEFILE : #
##############

dir_guard=@if [ ! -d $(BUILDDIR) ]; then echo "Creating build directory ..."; mkdir -p $(BUILDDIR); fi

.PHONY: all
all:$(BUILDDIR)/proteinortho_extract_from_graph.pl $(BUILDDIR)/proteinortho_compareProteinorthoGraphs.pl $(BUILDDIR)/proteinortho_grab_proteins.pl $(BUILDDIR)/proteinortho_formatUsearch.pl $(BUILDDIR)/proteinortho_do_mcl.pl $(BUILDDIR)/proteinortho2tree.pl $(BUILDDIR)/proteinortho2html.pl $(BUILDDIR)/proteinortho2xml.pl $(BUILDDIR)/proteinortho_singletons.pl $(BUILDDIR)/proteinortho_summary.pl $(BUILDDIR)/proteinortho_ffadj_mcs.py $(BUILDDIR)/proteinortho_clustering $(BUILDDIR)/proteinortho_history.pl $(BUILDDIR)/proteinortho_graphMinusRemovegraph $(BUILDDIR)/proteinortho_cleanupblastgraph $(BUILDDIR)/proteinortho_treeBuilderCore
	@echo "[100%] $(GREEN)Everything is compiled with no errors.$(NC)"

$(BUILDDIR)/proteinortho_extract_from_graph.pl: src/proteinortho_extract_from_graph.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_compareProteinorthoGraphs.pl: src/proteinortho_compareProteinorthoGraphs.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_grab_proteins.pl: src/proteinortho_grab_proteins.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_do_mcl.pl: src/proteinortho_do_mcl.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_formatUsearch.pl: src/proteinortho_formatUsearch.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho2html.pl: src/proteinortho2html.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho2xml.pl: src/proteinortho2xml.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_singletons.pl: src/proteinortho_singletons.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho2tree.pl: src/proteinortho2tree.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_ffadj_mcs.py: src/proteinortho_ffadj_mcs.py
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_history.pl: src/proteinortho_history.pl
	$(dir_guard)
	@cp $< $@

$(BUILDDIR)/proteinortho_summary.pl: src/proteinortho_summary.pl
	$(dir_guard)
	@cp $< $@

echoENV:
	@echo -n "CC = "
	@echo $(CC)
	@echo -n "CCFLAGS = "
	@echo $(CCFLAGS)
	@echo -n "CXX = "
	@echo $(CXX)
	@echo -n "CXXFLAGS = "
	@echo $(CXXFLAGS)
	@echo -n "LDFLAGS = "
	@echo $(LDFLAGS)
	@echo -n "LDLIBS = "
	@echo $(LDLIBS)

# 1. Try to compile statically with LAPACK
# 2. try to compile dynamically with the given lapack lib in src/

$(BUILDDIR)/proteinortho_clustering: src/proteinortho_clustering.cpp
	$(dir_guard)
	@echo "[ 10%] Prepare proteinortho_clustering ..."
ifeq ($(USELAPACK),TRUE)
ifeq ($(USEPRECOMPILEDLAPACK),TRUE)
ifeq ($(STATIC),TRUE)
	@echo "[ 20%] Building **proteinortho_clustering** with LAPACK (static linking)";
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -static -Wl,--allow-multiple-definition -llapack -lblas -lgfortran -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lquadmath && ([ $$? -eq 0 ] ) || ( \
		echo "......$(ORANGE)static linking failed, now I try dynamic linking.$(NC)"; \
		$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -llapack -lblas -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
			echo "......$(ORANGE)dynamic linking failed too, now I try dynamic linking without -WL,-whole-archive (this should now work for OSX).$(NC)"; \
			$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -llapack -lblas -pthread -lpthread && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
				echo "......$(ORANGE)dynamic linking failed (without -WL,-whole-archive) too too, now I try to openblas.$(NC)"; \
				$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -lopenblas -pthread -lpthread && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
					echo "......$(ORANGE)last linking failed too too too, now I try to recompile lapack (v.3.8.0) and then compile proteinortho_clustering with dynamic linking.$(NC)"; \
					echo "......[ 33%] Extracting the LAPACK library"; \
					if [ ! -d src/lapack-3.8.0 ]; then cd src; tar -xzvf lapack-3.8.0.tar.gz > /dev/null 2>&1; cd ..; fi; \
					echo "......[ 66%] Compiling the LAPACK library (using cmake + make)"; \
					if [ ! -f src/lapack-3.8.0/build/lib/liblapack.a ]; then mkdir src/lapack-3.8.0/build 2> /dev/null; cd src/lapack-3.8.0/build; cmake .. -DCMAKE_CXX_COMPILER=$(CXX) > /dev/null; make all -j4 > /dev/null 2>&1; cd ../../.. ; fi; \
					echo "......[ 99%] Building **proteinortho_clustering** with LAPACK (dynamic linking)"; \
					$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp -o $@ $< -Isrc/lapack-3.8.0/build/include/ -Lsrc/lapack-3.8.0/build/lib/ -llapack -lblas $(LDFLAGS) $(LDLIBS) -lgfortran && echo "......OK dynamic linking was successful for proteinortho_clustering!" || ( echo "" ) ; ) ) ) )
else
	@echo "[ 20%] Building **proteinortho_clustering** with LAPACK (dynamic linking)";
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -llapack -lblas -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
		echo "......$(ORANGE)dynamic linking failed too, now I try dynamic linking without -WL,-whole-archive (this should now work for OSX).$(NC)"; \
		$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -llapack -lblas -pthread -lpthread && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
			echo "......$(ORANGE)dynamic linking failed (without -WL,-whole-archive) too too, now I try to openblas.$(NC)"; \
			$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp  -o $@ $< $(LDFLAGS) $(LDLIBS) -lopenblas -pthread -lpthread && ([ $$? -eq 0 ] && echo "......OK dynamic linking was successful for proteinortho_clustering!";) || ( \
				echo "......$(ORANGE)last linking failed too too too, now I try to recompile lapack (v.3.8.0) and then compile proteinortho_clustering with dynamic linking.$(NC)"; \
				echo "......[ 33%] Extracting the LAPACK library"; \
				if [ ! -d src/lapack-3.8.0 ]; then cd src; tar -xzvf lapack-3.8.0.tar.gz > /dev/null 2>&1; cd ..; fi; \
				echo "......[ 66%] Compiling the LAPACK library (using cmake + make)"; \
				if [ ! -f src/lapack-3.8.0/build/lib/liblapack.a ]; then mkdir src/lapack-3.8.0/build 2> /dev/null; cd src/lapack-3.8.0/build; cmake .. -DCMAKE_CXX_COMPILER=$(CXX) > /dev/null; make all -j4 > /dev/null 2>&1; cd ../../.. ; fi; \
				echo "......[ 99%] Building **proteinortho_clustering** with LAPACK (dynamic linking)"; \
				$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp -o $@ $< -Isrc/lapack-3.8.0/build/include/ -Lsrc/lapack-3.8.0/build/lib/ -llapack -lblas $(LDFLAGS) $(LDLIBS) -lgfortran && echo "......OK dynamic linking was successful for proteinortho_clustering!" || ( echo "" ) ; ) ) )
endif
	@if [ ! -e $(BUILDDIR)/proteinortho_clustering ]; then echo "proteinortho_clustering compilation failed. Please visit https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes"; false; fi
else
	@echo "[ 15%] Extracting the LAPACK library";
	@if [ ! -d src/lapack-3.8.0 ]; then cd src; tar -xzvf lapack-3.8.0.tar.gz > /dev/null 2>&1; cd ..; fi;
	@echo "[ 18%] Compiling the LAPACK library (using cmake + make)";
	@if [ ! -f src/lapack-3.8.0/build/lib/liblapack.a ]; then mkdir src/lapack-3.8.0/build; cd src/lapack-3.8.0/build; cmake .. -DCMAKE_CXX_COMPILER=$(CXX) > /dev/null 2>&1; make all -j4 > /dev/null 2>&1; cd ../../.. ; fi;
	@echo "[ 20%] Building **proteinortho_clustering** with LAPACK (dynamic linking)";
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp -o $@ $< -Isrc/lapack-3.8.0/build/include/ -Lsrc/lapack-3.8.0/build/lib/ -llapack -lblas $(LDFLAGS) $(LDLIBS) -lgfortran;
endif
endif
ifeq ($(USELAPACK),FALSE)
	@echo "[ 20%] Building **proteinortho_clustering** WITHOUT(!) LAPACK";
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp -o $@ src/proteinortho_clustering_nolapack.cpp $(LDFLAGS) $(LDLIBS) -static && ([ $$? -eq 0 ] ) || ( \
		echo "......$(ORANGE)static linking failed of proteinortho_clustering_nolapack, now i switch to dynamic linking.$(NC)"; \
		$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -fopenmp -o $@ src/proteinortho_clustering_nolapack.cpp $(LDFLAGS) $(LDLIBS) && echo "......OK dynamic linking was successful for proteinortho_clustering_nolapack!"; )
endif

$(BUILDDIR)/proteinortho_cleanupblastgraph: src/cleanupblastgraph.cpp
	$(dir_guard)
	@echo "[ 50%] Building **cleanupblastgraph**"
ifeq ($(STATIC),TRUE)
	@$(CXX) -std=c++11 $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< $(LDFLAGS) $(LDLIBS) -static && ([ $$? -eq 0 ] ) || ( \
		echo "......$(ORANGE)static linking failed of cleanupblastgraph, now i switch to dynamic linking.$(NC)"; \
		$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< && echo "......OK dynamic linking was successful of cleanupblastgraph!"; )
else
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< $(LDFLAGS) $(LDLIBS)
endif

$(BUILDDIR)/proteinortho_graphMinusRemovegraph: src/graphMinusRemovegraph.cpp
	$(dir_guard)
	@echo "[ 25%] Building **graphMinusRemovegraph**"
ifeq ($(STATIC),TRUE)
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< $(LDFLAGS) $(LDLIBS) -static && ([ $$? -eq 0 ] ) || ( \
		echo "......$(ORANGE)static linking failed of graphMinusRemovegraph, now i switch to dynamic linking.$(NC)"; \
		$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< && echo "......OK dynamic linking was successful of graphMinusRemovegraph!"; )
else
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_PO) -o $@ $< $(LDFLAGS) $(LDLIBS)
endif

$(BUILDDIR)/proteinortho_treeBuilderCore: src/po_tree.c
	$(dir_guard)
	@echo "[ 75%] Building **po_tree**"
ifeq ($(STATIC),TRUE)
	@$(CC) $(CCFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS) -static && ([ $$? -eq 0 ] ) || ( \
		echo "......$(ORANGE)static linking failed of po_tree, now i switch to dynamic linking.$(NC)"; \
		$(CC) $(CCFLAGS) -o $@ $< && echo "......OK dynamic linking was successful of po_tree!"; )
else
	@$(CC) $(CCFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
endif

.PHONY: install
install: proteinortho6.pl proteinortho $(BUILDDIR)/proteinortho_extract_from_graph.pl $(BUILDDIR)/proteinortho_formatUsearch.pl $(BUILDDIR)/proteinortho_compareProteinorthoGraphs.pl $(BUILDDIR)/proteinortho_do_mcl.pl $(BUILDDIR)/proteinortho2html.pl $(BUILDDIR)/proteinortho2xml.pl $(BUILDDIR)/proteinortho_clustering $(BUILDDIR)/proteinortho_singletons.pl $(BUILDDIR)/proteinortho_ffadj_mcs.py $(BUILDDIR)/proteinortho2tree.pl $(BUILDDIR)/proteinortho_history.pl $(BUILDDIR)/proteinortho_cleanupblastgraph $(BUILDDIR)/proteinortho_graphMinusRemovegraph $(BUILDDIR)/proteinortho_treeBuilderCore $(BUILDDIR)/proteinortho_grab_proteins.pl $(BUILDDIR)/proteinortho_summary.pl $(BUILDDIR)/proteinortho_history.pl 
	@echo "INSTALLING everything to $(INSTALLDIR)"
	@install -v $^ $(INSTALLDIR);
	@echo "$(GREEN)Everything installed successfully to $(INSTALLDIR).$(NC)"
	@echo "If needed you can add $(INSTALLDIR) to \$$PATH with 'export PATH=\$$PATH:$(INSTALLDIR)'."

.PHONY: test
test: proteinortho6.pl test_clean test_step2 test_step3 test_clean2
	@echo "[TEST] All tests $(GREEN)passed$(NC)"

.PHONY: test_step2
test_step2: proteinortho6.pl
	@echo "[TEST] 1. basic proteinortho6.pl -step=2 test. (algorithms that are not present are skipped)"
	@echo -n " [1/12] -p=blastp+ test: "
	@if [ "$(shell which blastp)" = "" ]; then\
		echo "$(ORANGE)blastp missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_blastp -p=blastp+ test/*.faa; \
		set -e ; ./src/chk_test.pl test_blastp.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [2/12] -p=blastp+ synteny (PoFF) test: "
	@if [ "$(shell which blastp)" = "" ]; then\
		echo "$(ORANGE)blastp missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_synteny -synteny -singles -p=blastp+ test/*.faa; \
		set -e ; ./src/chk_test.pl test_synteny.proteinortho.tsv; \
		set -e ; ./src/chk_test.pl test_synteny.poff.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [3/12] -p=diamond test: "
	@if [ "$(shell which diamond)" = "" ]; then\
		echo "$(ORANGE)diamond missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_diamond -p=diamond test/*.faa; \
		set -e ; ./src/chk_test.pl test_diamond.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [4/12] -p=diamond (--moresensitive) test (subparaBlast): "
	@if [ "$(shell which diamond)" = "" ]; then\
		echo "$(ORANGE)diamond missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_diamondmoresensitive -p=diamond -subparaBlast="--more-sensitive" test/*.faa; \
		set -e ; ./src/chk_test.pl test_diamondmoresensitive.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [5/12] -p=lastp (lastal) test: "
	@if [ "$(shell which lastal)" = "" ]; then\
		echo "$(ORANGE)lastal missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_lastp -p=lastp test/*.faa; \
		set -e ; ./src/chk_test.pl test_lastp.proteinortho.tsv ; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [6/12] -p=topaz test: "
	@if [ "$(shell which topaz)" = "" ]; then\
		echo "$(ORANGE)topaz missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_topaz -p=topaz test/*.faa; \
		set -e ; ./src/chk_test.pl test_topaz.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [7/12] -p=usearch test: "
	@if [ "$(shell which usearch)" = "" ]; then\
		echo "$(ORANGE)usearch missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_usearch -p=usearch test/*.faa; \
		set -e ; ./src/chk_test.pl test_usearch.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [8/12] -p=ublast test: "
	@if [ "$(shell which usearch)" = "" ]; then\
		echo "$(ORANGE)usearch missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_ublast -p=ublast test/*.faa; \
		set -e ; ./src/chk_test.pl test_ublast.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [9/12] -p=rapsearch test: "
	@if [ "$(shell which rapsearch)" = "" ]; then\
		echo "$(ORANGE)rapsearch missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_rapsearch -p=rapsearch test/*.faa; \
		set -e ; ./src/chk_test.pl test_rapsearch.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [10/12] -p=blatp (blat) test: "
	@if [ "$(shell which blat)" = "" ]; then\
		echo "$(ORANGE)blat missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_blatp -p=blatp test/*.faa; \
		set -e ; ./src/chk_test.pl test_blatp.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [11/12] -p=mmseqsp (mmseqs) test: "
	@if [ "$(shell which mmseqs)" = "" ]; then\
		echo "$(ORANGE)mmseqs missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_mmseqsp -p=mmseqsp test/*.faa; \
		set -e ; ./src/chk_test.pl test_blatp.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

	@echo -n " [12/12] -p=blastp (legacy blast) test: "
	@if [ "$(shell which blastall)" = "" ]; then\
		echo "$(ORANGE)blastall missing, skipping...$(NC)"; \
	else \
		./proteinortho6.pl -silent -force -project=test_blastall -p=blastp test/*.faa; \
		set -e ; ./src/chk_test.pl test_blastall.proteinortho.tsv; \
		echo "$(GREEN)passed$(NC)"; \
	fi

.PHONY: test_step3
test_step3: proteinortho6.pl test_step2
	@echo "[TEST] 2. -step=3 tests (proteinortho_clustering) "
	@echo -n " [1/2] various test functions of proteinortho_clustering (if this fails, try make clean first): "; \
	$(BUILDDIR)/proteinortho_clustering -test> /dev/null 2>&1
	@echo "$(GREEN)passed$(NC)"
	@echo -n " [2/2] Compare results of 'with lapack' and 'without lapack': "; \
	$(BUILDDIR)/proteinortho_clustering -epsilon 0 test_blastp.blast-graph> /dev/null 2>&1; \
	sort remove.graph -o test.A> /dev/null 2>&1; \
	$(BUILDDIR)/proteinortho_clustering -epsilon 0 -lapack 0 test_blastp.blast-graph> /dev/null 2>&1; \
	sort remove.graph -o test.B> /dev/null 2>&1; \
	set -e ; diff test.A test.B;
	@echo "$(GREEN)passed$(NC)"

.PHONY: test_clean
test_clean:
	@echo "[TEST] Clean up all test files..."; \
	rm -rf proteinortho_cache_test_* test.* test_* test/C.faa.* test/E.faa.* test/C2.faa.* test/L.faa.* test/M.faa.*> /dev/null 2>&1;

.PHONY: test_clean2
test_clean2:
	@echo "[TEST] Clean up all test files..."; \
	rm -rf proteinortho_cache_test_* test.* test_* test/C.faa.* test/E.faa.* test/C2.faa.* test/L.faa.* test/M.faa.*> /dev/null 2>&1;

.PHONY: clean
clean:
	rm -rf src/BUILD test/C.faa.* test/E.faa.* test/C2.faa.* test/L.faa.* test/M.faa.*
	rm -rf src/lapack-3.8.0/
