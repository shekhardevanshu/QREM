
default: app

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

app: main.o buildH.o
		$(LINK.C) -o $@ $^ ${SLEPC_EPS_LIB}
		${RM} app.o

main.o: main.c buildH.h
		$(LINK.C) -c main.c ${SLEPC_EPS_LIB}

buildH.o: buildH.c buildH.h
		$(LINK.C) -c buildH.c ${SLEPC_EPS_LIB}


