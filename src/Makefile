# Directories.
TOPDIR = ..
BINDIR = $(TOPDIR)/bin
DATADIR = $(TOPDIR)/data
OBJDIR = ./obj
DIRS = $(BINDIR) $(OBJDIR)

vpath %.dat $(DATADIR)

CFLAGS += -c
CFLAGS += -Wall
CFLAGS += -O2
CFLAGS += -march=native -mtune=native
CFLAGS += -fomit-frame-pointer

LDFLAGS = -lm

SOURCES = crypto_hash_sha512.c \
					crypto_stream.c \
					randombytes.c \
					fastrandombytes.c \
					shred.c \
					convert.c \
					pack.c \
					pol.c \
					params.c \
					pqntrusign.c

BENCHOBJ = bench.o
TESTOBJ = sanity.o

OBJECTS= $(addprefix $(OBJDIR)/, $(SOURCES:.c=.o))
BENCH=$(BINDIR)/bench
TEST=$(BINDIR)/test


.PHONY : all test dirs

all: dirs $(BENCH)

test: $(TEST)

dirs : $(DIRS)
$(DIRS):
	-mkdir -p $@

$(BENCH): $(OBJECTS) $(BENCHOBJ)
	$(CC) $(LDFLAGS) $(OBJECTS) $(BENCHOBJ) -o $@

$(TEST): $(OBJECTS) $(TESTOBJ)
	$(CC) $(LDFLAGS) $(OBJECTS) $(TESTOBJ) -o $@

$(OBJDIR)/%.o : %.c | $(OBJDIR)
	$(CC) $(CFLAGS) $< -o $@

.PHONY : clean
clean:
	-rm -rf $(OBJECTS) $(TESTOBJ) $(BENCHOBJ)
	-rmdir $(OBJDIR)
	-rm -f $(BENCH) $(TEST)



