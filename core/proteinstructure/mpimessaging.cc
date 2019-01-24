
#ifdef WITH_MPI

using namespace std;

#include "mpimessaging.h"

/*
 * Returns a MPI_Datatype used to communicate the edge
 * connections that need to be evaluated by secondary processes
 */
MPI_Datatype mpi_type_update_score_send()
{
	MPI_Datatype sendtype;
	MPI_Type_contiguous(4, MPI_INT,&sendtype);
    MPI_Type_commit(&sendtype);

	return sendtype;
}

/*
 * Returns a data type that can be used to send the
 * scores calculated by secondary processes back to the
 * main process
 */
MPI_Datatype mpi_type_update_score_receive()
{
	MPI_Datatype receivetype;
        int blocklengths[5]={1,1,1,1,12};
	MPI_Datatype types[5]={MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint displacements[5];
	MPI_Aint intex, doubleex;

	MPI_Type_extent(MPI_INT, &intex);
	MPI_Type_extent(MPI_DOUBLE, &doubleex);
	
	displacements[0] = (MPI_Aint) 0;
	displacements[1] = intex;
	displacements[2] = intex + doubleex;
	displacements[3] = intex + (2 * doubleex);
	displacements[4] = (2 * intex) + (2 * doubleex);

	MPI_Type_struct(5, blocklengths, displacements, types, &receivetype);
	MPI_Type_commit(&receivetype);

	return receivetype;
}

/*
 * Returns a data type used to exchange the values of spanning
 * tree edges between the main and auxiliary ranks
 */
MPI_Datatype mpi_type_edge()
{
	MPI_Datatype edgetype;
	MPI_Type_contiguous(3, MPI_INT,&edgetype);
    MPI_Type_commit(&edgetype);

	return edgetype;
}

#endif
