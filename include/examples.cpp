

Mesh my_mesh(128);
int tid, i, i_n = 0;
printf("Length = %d\n", my_mesh.length);

#pragma omp parallel private(tid, i, i_n)
{	
	#pragma omp for nowait
	for (i = 0; i < my_mesh.length; i++)
	{
		my_mesh[i]=0.;
		i_n = i;
	}
	/* Obtain and print thread id */
	tid = omp_get_thread_num();
	printf("Hello World from thread = %d\n", tid);
	printf("Thread %d ended up with i = %i.\n", tid, i_n);
}

