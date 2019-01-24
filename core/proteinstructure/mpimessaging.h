#ifndef _MPIMESSAGING_H_
#define _MPIMESSAGING_H_

using namespace std;

#ifdef WITH_MPI

#include "mpi.h"
#include "soroban_score.h"

/* Give an integer value to the tags used for mpi messaging (one for each message type sent) */
#define MPI_MESSAGE_UPDATE_SCORE_SEND_TAG 50
#define MPI_MESSAGE_UPDATE_SCORE_RECEIVE_TAG 51
#define MPI_MESSAGE_FINISH_TAG 52

// The first tag MPI_MESSAGE_UPDATE_PRE_POPULATION_TAG is received by the
// main Recv instruction. It just notifies that it should ignore the body
// of the current message and then wait for a MPI_MESSAGE_UPDATE_POPULATION_SIZE_TAG
// followed by a MPI_MESSAGE_UPDATE_POPULATION_TAG
// that will contain the edge information that constitutes a population update
#define MPI_MESSAGE_UPDATE_PRE_POPULATION_TAG 53
#define MPI_MESSAGE_UPDATE_POPULATION_SIZE_TAG 54
#define MPI_MESSAGE_UPDATE_POPULATION_TAG 55

// Message sent from the main to aux ranks to indicate them that they should do clustering calculations
#define MPI_MESSAGE_CLUSTER_SEND_TAG 56
// Message from aux ranks to main notifying the clustering results
#define MPI_MESSAGE_CLUSTER_RECEIVE_TAG 57

/*
 * Mirrors the contents of md_edge in mdgraph.h
 * When the message is sent this will actually be sent
 * as a series of these elements (as many as we have edges
 * for the spanning tree)
 */
typedef struct mpi_message_update_score_send_t
{
	/* 
	 * identifies the specific prediction
	 * that this edge belongs to. This is added for convenience
	 * so that when the score is returned to the rank
	 * that coordinates it knows what prediction to update the
	 * score to
	 */
	int population_index;
	int receptor, ligand, prediction_number;
} mpi_message_update_score_send_t;

/*
 * Contains the scoring elements that have been calculated
 * by a secondary process, and now they need to be integrated
 * to the rest of the scores in the main rank process
 */
typedef struct mpi_message_update_score_receive_t
{
	/*
	 * This index indicates that this score corresponds
	 * to the element at population[population_index]
	 */
	int population_index;
	int clashes;
	double score;
	double pairwise_score;
	soroban_t soro_score;
} mpi_message_update_score_receive_t;

/*
 * When a population has been altered it is necessary to
 * notify the auxiliary ranks.
 * It is necessary to notify, for each population element,
 * all the edges that form the element's spanning tree.
 *
 * This type defines the information that describes one edge.
 * A total of populationSize * (numEdges * thisTypeSize) bytes
 * will be sent in this type of message
 */
typedef struct mpi_message_edge_t
{
	int receptor, ligand, prediction_number;
} mpi_message_edge_t;

/*
 * Returns a MPI_Datatype used to communicate the edge
 * connections that need to be evaluated by secondary processes
 */
MPI_Datatype mpi_type_update_score_send();

/*
 * Returns a data type that can be used to send the
 * scores calculated by secondary processes back to the
 * main process
 */
MPI_Datatype mpi_type_update_score_receive();

/*
 * Returns a data type used to exchange the values of spanning
 * tree edges between the main and auxiliary ranks
 */
MPI_Datatype mpi_type_edge();

#endif

#endif
