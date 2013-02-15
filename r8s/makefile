# Make file for r8s
# Updated April 2003 for Linux and Mac OS 10.2. No other builds currently supported
#

FC = gfortran  # Fortran compiler; standard on Linux, may have to go fetch this on OS X! (get it from Fink)
CC = gcc  # GNU C compiler
LIBS =  -lgfortran -lm  # FORTRAN to C library, and standard C library 
#LPATH = -L/usr/local/gfortran/lib  
#LPATH = -L/usr/local/lib/gfortran  
LPATH = -L/usr/local/gfortran/lib # correct location as of August 2011. 
CFLAGS = -g -std=c99 -pedantic  # for debugging, etc. 
#CFLAGS =  # usual case 


OBJS =  DrawTree.o TreeSim.o moment.o powell.o ConstrOpt.o \
MyUtilities.o WuLi.o DistrFuncs.o NRCvectorUtils.o distance.o penalty.o \
ObjFunc.o main.o relativeRates.o ReadNexusFile2.o  memory.o root3taxa.o \
TimeAlgorithms.o nextToken2.o storeNexusFile.o MinND.o TreeUtils.o structures.o \
blas.o tn.o TNwrapper.o continuousML.o ancestral.o covarion.o


#r8s: ${OBJS}
#	${CC} -o r8s  ${OBJS} ${LIBS}
# Use the following for MAC OS X
r8s: ${OBJS}
	${CC} -o r8s  ${LPATH} ${OBJS} ${LIBS}

# DO NOT DELETE


covarion.o:covarion.h
ancestral.o:ancestral.h
continuousML.o:continuousML.h
TNwrapper.o:TNwrapper.h
DrawTree.o:DrawTree.h TreeUtils.h /usr/include/string.h /usr/include/ctype.h /usr/include/stdio.h /usr/include/stdlib.h
TreeSim.o:TreeSim.h TreeUtils.h /usr/include/stdio.h /usr/include/stdlib.h /usr/include/math.h 
moment.o:/usr/include/math.h moment.h MyUtilities.h
ConstrOpt.o: /usr/include/stdio.h /usr/include/stdlib.h
ConstrOpt.o: /usr/include/math.h Maximize.h
ConstrOpt.o: NRCvectorUtils.h ConstrOpt.h MyUtilities.h /usr/include/string.h
ConstrOpt.o: /usr/include/ctype.h
DistrFuncs.o: DistrFuncs.h /usr/include/math.h
MinND.o: Maximize.h NRCvectorUtils.h
MyUtilities.o: MyUtilities.h /usr/include/stdio.h /usr/include/string.h
MyUtilities.o: /usr/include/stdlib.h
MyUtilities.o: /usr/include/ctype.h memory.h
NRCvectorUtils.o: NRCvectorUtils.h /usr/include/stdio.h /usr/include/stdlib.h
NRCvectorUtils.o: memory.h
ObjFunc.o: /usr/include/stdio.h /usr/include/stdlib.h 
ObjFunc.o: ObjFunc.h TreeUtils.h /usr/include/string.h /usr/include/ctype.h
ObjFunc.o: structures.h TimeAlgorithms.h ConstrOpt.h objective.h ObjFuncHeader.h
ObjFunc.o: /usr/include/math.h MyUtilities.h
ReadNexusFile2.o: WuLi.h nexus.h /usr/include/stdio.h
ReadNexusFile2.o: /usr/include/string.h /usr/include/stdlib.h
ReadNexusFile2.o:  /usr/include/ctype.h TreeUtils.h
ReadNexusFile2.o: structures.h MyUtilities.h memory.h ObjFunc.h
ReadNexusFile2.o: TimeAlgorithms.h TreeSim.h 
TimeAlgorithms.o: TreeUtils.h /usr/include/string.h /usr/include/stdio.h
TimeAlgorithms.o: /usr/include/ctype.h /usr/include/stdlib.h
TimeAlgorithms.o: structures.h TimeAlgorithms.h
TimeAlgorithms.o: memory.h /usr/include/math.h penalty.h
TreeUtils.o: TreeUtils.h /usr/include/string.h /usr/include/stdio.h
TreeUtils.o: /usr/include/ctype.h /usr/include/stdlib.h
TreeUtils.o: structures.h memory.h
WuLi.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
WuLi.o: nexus.h /usr/include/string.h
WuLi.o: /usr/include/ctype.h TreeUtils.h structures.h WuLi.h
WuLi.o: MyUtilities.h memory.h DistrFuncs.h distance.h
distance.o: distance.h nexus.h /usr/include/stdio.h /usr/include/string.h
distance.o: /usr/include/stdlib.h /usr/include/ctype.h
distance.o: TreeUtils.h structures.h /usr/include/math.h
main.o: /usr/include/stdio.h /usr/include/stdlib.h
main.o: /usr/include/string.h storeNexusFile.h nexus.h /usr/include/ctype.h
main.o: TreeUtils.h structures.h
memory.o: /usr/include/stdlib.h memory.h
memory.o: /usr/include/errno.h 
nextToken2.o: nexus.h /usr/include/stdio.h /usr/include/string.h
nextToken2.o: /usr/include/stdlib.h
nextToken2.o: /usr/include/ctype.h TreeUtils.h structures.h MyUtilities.h
nextToken2.o: memory.h
#optimizeNRCLR.o: NRCvectorUtils.h optimizeNRCLR.h /usr/include/math.h
#optimizeNRCLR.o: /usr/include/stdio.h /usr/include/stdlib.h
penalty.o: TreeUtils.h /usr/include/string.h /usr/include/stdio.h
penalty.o: /usr/include/ctype.h /usr/include/stdlib.h
penalty.o: structures.h penalty.h /usr/include/math.h
relativeRates.o: nexus.h /usr/include/stdio.h /usr/include/string.h
relativeRates.o: /usr/include/stdlib.h
relativeRates.o: /usr/include/ctype.h TreeUtils.h structures.h
relativeRates.o: relativeRates.h WuLi.h MyUtilities.h root3taxa.h
relativeRates.o: distance.h
root3taxa.o: TreeUtils.h /usr/include/string.h /usr/include/stdio.h
root3taxa.o: /usr/include/ctype.h /usr/include/stdlib.h
root3taxa.o: structures.h memory.h root3taxa.h
storeNexusFile.o: /usr/include/stdio.h /usr/include/stdlib.h
storeNexusFile.o: storeNexusFile.h MyUtilities.h
storeNexusFile.o: /usr/include/string.h /usr/include/ctype.h
structures.o: structures.h /usr/include/stdlib.h
structures.o: /usr/include/stdio.h /usr/include/string.h MyUtilities.h
structures.o: /usr/include/ctype.h memory.h
#BFGSwrapper.o: BFGSwrapper.h
