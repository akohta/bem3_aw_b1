# makefile for GNU C compiler
CC     =gcc
CFLAGS =-O3 -Wall
LDFLAGS=-lpng -lz -lgsl -lgslcblas
MPFLAGS=-fopenmp
MP_LD  =-lgomp
# link mkl lapack library, static linking 
LDF_MKL =-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl#mkl lapack link library (static linking) 
# mkl compile option
COP_MKL =-DMKL_ILP64 -m64 -I${MKLROOT}/include
SRCDIR1=d3b1_src
SRCDIR2=maw_src
SRCDIR3=com_src
OBJDIR =obj
TARGET1=aw_d3b1_bv_solver
TRGSRC1=aw_d3b1_bv_solver.c
TARGET2=example1.out
TRGSRC2=example1.c
TARGET3=example2.out
TRGSRC3=example2.c
TARGET4 =example3.out
TRGSRC4 =example3.c

SRCS1=$(wildcard $(SRCDIR1)/*.c)
OBJS1=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS1)) ))
HEAD1=$(wildcard $(SRCDIR1)/*.h)

SRCS2=$(wildcard $(SRCDIR2)/*.c)
OBJS2=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS2)) ))
HEAD2=$(wildcard $(SRCDIR2)/*.h)

SRCS3=$(wildcard $(SRCDIR3)/*.c)
OBJS3=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS3)) ))
HEAD3=$(wildcard $(SRCDIR3)/*.h)

TRGOBJ1=$(OBJS1) $(OBJS2) $(OBJS3)
TRGOBJ2=$(filter-out $(OBJDIR)/$(TARGET1).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC2))
TRGOBJ3=$(filter-out $(OBJDIR)/$(TARGET1).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC3))
TRGOBJ4=$(filter-out $(OBJDIR)/$(TARGET1).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC4))

all : directories $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

directories:
	@mkdir -p $(OBJDIR)

$(TARGET1) : $(TRGOBJ1) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LDF_MKL) $(MP_LD)

$(TARGET2) : $(TRGOBJ2) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LDF_MKL) $(MP_LD)

$(TARGET3) : $(TRGOBJ3) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LDF_MKL) $(MP_LD)

$(TARGET4) : $(TRGOBJ4) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LDF_MKL) $(MP_LD)

$(OBJDIR)/%.o : $(SRCDIR1)/%.c
	$(CC) $(CFLAGS) $(COP_MKL) $(MPFLAGS) -I$(SRCDIR2) -I$(SRCDIR3) -c $< -o $@

$(OBJDIR)/%.o : $(SRCDIR2)/%.c
	$(CC) $(CFLAGS) -I$(SRCDIR3) -c $< -o $@

$(OBJDIR)/%.o : $(SRCDIR3)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) $(COP_MKL) $(MPFLAGS) -I$(SRCDIR1) -I$(SRCDIR2) -I$(SRCDIR3) -c $<

clean:
	@rm -rf $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(OBJDIR) ./*.o

$(OBJS1) : $(HEAD1) $(HEAD2) $(HEAD3)
$(OBJS2) : $(HEAD2) $(HEAD3)
$(OBJS3) : $(HEAD3)
