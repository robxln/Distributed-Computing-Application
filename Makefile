build:
	mpicc -o distributed_computing_app distributed_computing_app.c -Wall

clean:
	rm -rf distributed_computing_app
