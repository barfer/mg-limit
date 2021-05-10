PROGRAM = mackeyGlass
CC = g++
CFLAGS = -O2

USERDIR = /mnt/c/Users/ferin

CAPDBINLIB = ${USERDIR}/lib/capd/bin
CAPDLIBS = `${CAPDBINLIB}/capd-config --cflags --libs`

INCLUDE = -I.

HEADERFILES = main.h tools.h branch.h segment.h

${PROGRAM}: ${PROGRAM}.cpp segment.o branch.o tools.o pugixml.o ${HEADERFILES}
	${CC} ${CFLAGS} ${INCLUDE} ${PROGRAM}.cpp -o ${PROGRAM} segment.o branch.o tools.o pugixml.o ${CAPDLIBS}

pugixml.o: pugixml.hpp
	${CC} ${CFLAGS} ${INCLUDE} -o pugixml.o -c pugixml.cpp

segment.o: segment.cpp ${HEADERFILES}
	${CC} ${CFLAGS} ${INCLUDE} -o segment.o -c segment.cpp ${CAPDLIBS}

branch.o: branch.cpp ${HEADERFILES}
	${CC} ${CFLAGS} ${INCLUDE} -o branch.o -c branch.cpp ${CAPDLIBS}

tools.o: tools.cpp ${HEADERFILES}
	${CC} ${CFLAGS} ${INCLUDE} -o tools.o -c tools.cpp ${CAPDLIBS}

clean:
	rm *~ ${PROGRAM}

tar:
	tar cvfz ${PROGRAM}.tgz ${PROGRAM}.cpp ${HEADERFILES} Makefile

