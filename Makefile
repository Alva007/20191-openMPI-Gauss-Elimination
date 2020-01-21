compile-cycle: cycle-inverse.c
	mpicc -o target/cycle-inverse cycle-inverse.c

run-cycle: 
	mpirun --oversubscribe -np 8 --hostfile mpi_hosts target/cycle-inverse --test 1200

compile-inverse: inverse.c
	mpicc -o target/inverse inverse.c

run-inverse: 
	mpirun --oversubscribe -np 8 target/inverse --test 0

create-test: maketest.c
	./test/makeTest

compile-test: ./test/maketest.c
	gcc -o ./test/makeTest ./test/maketest.c
