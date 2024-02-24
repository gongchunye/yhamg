DIRS = src test

all:
	@ \
	for i in ${DIRS}; \
	do \
	  if [ -d $$i ]; \
	  then (cd $$i; make -j); \
	  fi; \
	done

clean:
	@ \
	for i in ${DIRS}; \
	do \
	  if [ -d $$i ]; \
	  then (cd $$i; make clean); \
	  fi; \
	done
