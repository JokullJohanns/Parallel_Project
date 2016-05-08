default: program


program:
	mkdir bin
	mpicc kMeans_MPI.c -lm -o bin/kMeans_MPI.out
	gcc kMeans_serial.c -lm -o bin/kMeans_serial.out

clean:
	-rm -f bin/kMeans_serial.out
	-rm -f bin/kMeans_MPI.out
	rmdir bin