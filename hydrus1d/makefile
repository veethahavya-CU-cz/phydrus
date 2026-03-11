# Makefile for HYDRUS-1D — build from repository root

SRCDIR    := src
MCMCDIR   := src/mcmc
BINDIR    := bin
CC        := gfortran
FFLAGS    := -g -ffpe-summary=none -std=legacy
LDFLAGS   ?=
# macOS: invoke `make STRIP_RPATH=/path/to/rpath` to strip a specific rpath
STRIP_RPATH ?=
UNAME_S   := $(shell uname -s)

FILES := HYDRUS.FOR INPUT.FOR HYSTER.FOR MATERIAL.FOR OUTPUT.FOR SINK.FOR \
         SOLUTE.FOR TEMPER.FOR TIME.FOR WATFLOW.FOR

SRCS      := $(addprefix $(SRCDIR)/,$(FILES))
OBJS      := $(SRCS:.FOR=.o)

MCMC_SRCS := $(addprefix $(MCMCDIR)/,$(FILES))
MCMC_OBJS := $(MCMC_SRCS:.FOR=.o)

.PHONY: default all mcmc clean f2py

default: $(BINDIR)/hydrus

all: $(BINDIR)/hydrus $(BINDIR)/hydrus-mcmc

mcmc: $(BINDIR)/hydrus-mcmc

$(BINDIR)/hydrus: $(OBJS)
	@mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)
	@if [ "$(UNAME_S)" = "Darwin" ] && [ -n "$(STRIP_RPATH)" ]; then \
		echo "Removing rpath $(STRIP_RPATH) from $@"; \
		install_name_tool -delete_rpath $(STRIP_RPATH) $@ || true; \
	fi
	@rm -f $(OBJS)

$(BINDIR)/hydrus-mcmc: $(MCMC_OBJS)
	@mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) -o $@ $(MCMC_OBJS)
	@if [ "$(UNAME_S)" = "Darwin" ] && [ -n "$(STRIP_RPATH)" ]; then \
		echo "Removing rpath $(STRIP_RPATH) from $@"; \
		install_name_tool -delete_rpath $(STRIP_RPATH) $@ || true; \
	fi
	@rm -f $(MCMC_OBJS)

$(SRCDIR)/%.o: $(SRCDIR)/%.FOR
	$(CC) $(FFLAGS) -c -o $@ $<

$(MCMCDIR)/%.o: $(MCMCDIR)/%.FOR
	$(CC) $(FFLAGS) -c -o $@ $<

f2py:
	f2py -c $(SRCS) -m hydrus

clean:
	rm -f $(OBJS) $(MCMC_OBJS) $(BINDIR)/hydrus $(BINDIR)/hydrus-mcmc
