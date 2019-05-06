
LIBS=seq_file string_buffer gsl cJSON

all: $(LIBS)

clean:
	cd seq_file && make clean
	cd string_buffer && make clean

seq_file: seq_file/Makefile
	cd seq_file && git pull && make

string_buffer: string_buffer/Makefile
	cd string_buffer && git pull && make

gsl : gsl-1.16/configure.ac

#
# Clone
#

seq_file/Makefile:
	git clone https://github.com/noporpoise/seq_file.git

string_buffer/Makefile:
	git clone https://github.com/noporpoise/string_buffer.git

gsl-1.16/configure.ac:
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	tar xfz gsl-1.16.tar.gz
	rm gsl-1.16.tar.gz
	cd gsl-1.16 && ./configure && make

.PHONY: all clean $(LIBS)
