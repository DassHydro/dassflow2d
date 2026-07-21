#ifndef _ADMPI_ADMPI_H_
#define _ADMPI_ADMPI_H_

#include <mpi.h>

/**
 * does the request originate with a  send or a receive 
 */
enum adMPI_Request_origin_E { 
  WAIT_FOR_SEND,
  WAIT_FOR_RECV
};

#ifdef ADMPI_FORTRANCOMPATIBLE
typedef int adMPI_Request_origin;
#else 
typedef enum adMPI_Request_origin_E adMPI_Request_origin;
#endif 

/**
 * MPI_Request augmented with extra information 
 */ 
struct adMPI_Request_S {
  // MPI_Isend / MPI_Irecv  buf  parameter 
  void *buf;
  // The corresponding adjoint buffer
  void *adjointBuf ;
  // MPI_Isend / MPI_Irecv  count  parameter 
  int count;
  // MPI_Isend / MPI_Irecv  datatype  parameter 
  MPI_Datatype datatype;
  // MPI_Isend / MPI_Irecv  dst or src  parameter 
  int endPoint;
  // MPI_Isend / MPI_Irecv  tag parameter 
  int tag;
  // MPI_Isend / MPI_Irecv  comm  parameter 
  MPI_Comm comm;
  // temporary adjoint buffer
  void *adjointTempBuf;
  // the count of the adjoint buffer size in terms of the original data type
  int adjointCount;
  // the "plain" request returned by MPI_Isend / MPI_Irecv
  MPI_Request plainRequest;
  // the "plain" request returned in the bwd sweep by the MPI_Isend / MPI_Irecv in the adMPI_Wait_bwd
  MPI_Request bwRequest;
  // the "plain" request returned in tangent mode by the shadow MPI_Isend / MPI_Irecv
  MPI_Request shadowRequest;
  // MPI_Isend  / MPI_Irecv sets this to tell the Wait what it is waiting for
  enum adMPI_Request_origin_E origin;
};

#ifdef ADMPI_FORTRANCOMPATIBLE
typedef MPI_Request adMPI_Request;
#else 
typedef struct adMPI_Request_S adMPI_Request;
#endif

/** The type required for TANGENT user-given reduction functions,
 * that are passed e.g. to adMPI_Reduce_d in Tapenade tangent diff MPI code. */
typedef void (TLM_userFunctionF) (void*, void*, void*, void*, int*, MPI_Datatype*, MPI_Datatype*) ;

/** The type required for ADJOINT user-given reduction functions,
 * that are passed e.g. to adMPI_Reduce_b in Tapenade adjoint diff MPI code. */
typedef void (ADJ_userFunctionF) (void*, void*, void*, void*, int*, MPI_Datatype*, MPI_Datatype*) ;

// =Turn=

/** Special primitive that must be made visible to (the adjoint diff of) the application.
 * Specifies that the primal "buf" variable has "adjointBuf" as corresponding adjoint.
 */
void adMPI_Turn(void* buf, void* adjointBuf) ;

// =Init=

/**
 * this wrapper variant of MPI_Init has no adjoint transformation / trace functionality; to be used outside of the transformed/traced code section
 */
int adMPI_Init(int* argc, char*** argv) ;


// =Finalize=

/**
 * this wrapper variant of MPI_Finalize has no adjoint transformation / trace functionality; to be used outside of the transformed/traced code section
 */
int adMPI_Finalize(void) ;

// =Buffer_attach=

// =Buffer_detach=

// =Type_contiguous=

int adMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype* newtype) ;

// =Type_create_struct=

int adMPI_Type_create_struct(int count,
                               int array_of_blocklengths[],
                               MPI_Aint array_of_displacements[],
                               MPI_Datatype array_of_types[],
                               MPI_Datatype *newtype) ;

// =Type_commit=

int adMPI_Type_commit(MPI_Datatype *datatype) ;

// =Type_create_resized=

int adMPI_Type_create_resized(MPI_Datatype oldtype,
                                MPI_Aint lb,
                                MPI_Aint extent,
                                MPI_Datatype *newtype) ;

// =Type_free=

int adMPI_Type_free(MPI_Datatype *datatype) ;

// =Op_create=

int adMPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op) ;

// =Op_free=

int adMPI_Op_free(MPI_Op *op) ;

// =Send=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param dest
 * @param tag
 * @param comm
 * @return
 */

int MPI_Send_d(void* buf,  void* shadowbuf,
               int count,
               MPI_Datatype datatype, MPI_Datatype shadowdatatype,
               int dest,
               int tag,
               MPI_Comm comm) ;

int MPI_Send_fwd(void* buf, 
                 int count, 
                 MPI_Datatype datatype, 
                 int dest, 
                 int tag,
                 MPI_Comm comm) ;

int MPI_Send_bwd(void* buf,
                 int count, 
                 MPI_Datatype datatype, 
                 int dest, 
                 int tag,
                 MPI_Comm comm) ;

// =Isend=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param dest
 * @param tag
 * @param comm
 * @param request see \ref requests
 * @return
 */
int adMPI_Isend(void* buf, 
		int count, 
		MPI_Datatype datatype, 
		int dest, 
		int tag,
		MPI_Comm comm, 
		adMPI_Request* request) ;

int adMPI_Isend_d(void* buf, void* shadowbuf,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                    int dest,
                    int tag,
                    MPI_Comm comm,
                    adMPI_Request* request) ;

int adMPI_Isend_fwd(void* buf,
		   int count, 
		   MPI_Datatype datatype, 
		   int dest, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) ;

int adMPI_Isend_bwd(void* buf,
		   int count, 
		   MPI_Datatype datatype, 
		   int dest, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) ;

// =Bsend=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param dest
 * @param tag
 * @param comm
 * @return
 */

// =Rsend=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param dest
 * @param tag
 * @param comm
 * @return
 */

// =Recv=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param src
 * @param tag
 * @param comm
 * @param status
 * @return
 */

int MPI_Recv_d(void* buf, void* shadowbuf,
                  int count,
                  MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                  int src,
                  int tag,
                  MPI_Comm comm,
                  MPI_Status* status) ;

int MPI_Recv_fwd(void* buf, 
		 int count,
		 MPI_Datatype datatype, 
		 int src, 
		 int tag,
		 MPI_Comm comm,
		 MPI_Status* status) ;  

int MPI_Recv_bwd(void* buf, 
		 int count,
		 MPI_Datatype datatype, 
		 int src, 
		 int tag,
		 MPI_Comm comm,
		 MPI_Status* status) ;  
                  
// =Irecv=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param src
 * @param tag
 * @param comm
 * @param request see \ref requests
 * @return
 */
int adMPI_Irecv(void* buf, 
		int count, 
		MPI_Datatype datatype, 
		int src, 
		int tag,
		MPI_Comm comm, 
		adMPI_Request* request) ;

int adMPI_Irecv_d(void* buf, void* shadowbuf,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                    int source,
                    int tag,
                    MPI_Comm comm,
                    adMPI_Request* request) ;

int adMPI_Irecv_fwd(void* buf,
		   int count,
		   MPI_Datatype datatype,
		   int source,
		   int tag,
		   MPI_Comm comm,
		   adMPI_Request* request) ;

int adMPI_Irecv_bwd(void* buf, 
		   int count, 
		   MPI_Datatype datatype, 
		   int source, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) ;

// =Wait=

/**
 * before we start reverse we need to make sure there are no pending requests in our userIF bookkeeping 
 */
int adMPI_Wait(adMPI_Request *request, MPI_Status *status) ;

int adMPI_Wait_d(adMPI_Request *request, MPI_Status *status) ;

int adMPI_Wait_fwd(adMPI_Request *request, MPI_Status *status) ;

int adMPI_Wait_bwd(adMPI_Request *request, MPI_Status *status) ;

// =Waitall=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param count
 * @param requests see \ref requests
 * @param statuses
 * @return
 */
int adMPI_Waitall(int count, adMPI_Request requests[], MPI_Status statuses[]) ;

// =Barrier=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param comm
 * @return
 */

int MPI_Barrier_d(MPI_Comm comm) ;

int MPI_Barrier_fwd(MPI_Comm comm) ;

int MPI_Barrier_bwd(MPI_Comm comm) ;

// =Gather=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sendbuf
 * @param sendcnt
 * @param sendtype see \ref datatypes
 * @param recvbuf
 * @param recvcnt
 * @param recvtype see \ref datatypes
 * @param root
 * @param comm
 * @return
 */

int MPI_Gather_d(void *sendbuf, void *shadowsendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                    void *recvbuf, void *shadowrecvbuf,
                    int recvcnt,
                    MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                    int root,
                    MPI_Comm comm) ;

int MPI_Gather_fwd(void *sendbuf,
		   int sendcnt,
		   MPI_Datatype sendtype,
		   void *recvbuf,
		   int recvcnt,
		   MPI_Datatype recvtype,
		   int root,
		   MPI_Comm comm) ;

int MPI_Gather_bwd(void *sendbuf,
		   int sendcnt,
		   MPI_Datatype sendtype,
		   void *recvbuf,
		   int recvcnt,
		   MPI_Datatype recvtype,
		   int root,
		   MPI_Comm comm) ;

// =Gatherv=


int MPI_Gatherv_d(void *sendbuf, void *shadowsendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                     void *recvbuf, void *shadowrecvbuf,
                     int *recvcnts,
                     int *displs,
                     MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                     int root,
                     MPI_Comm comm) ;

int MPI_Gatherv_fwd(void *sendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype,
                    void *recvbuf,
                    int *recvcnts,
                    int *displs,
                    MPI_Datatype recvtype,
                    int root,
                    MPI_Comm comm) ;

int MPI_Gatherv_bwd(void *sendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype,
                    void *recvbuf,
                    int *recvcnts,
                    int *displs,
                    MPI_Datatype recvtype,
                    int root,
                    MPI_Comm comm) ;

// =Allgather=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sendbuf
 * @param sendcount
 * @param sendtype see \ref datatypes
 * @param recvbuf
 * @param recvcount
 * @param recvtype see \ref datatypes
 * @param comm
 * @return
 */

int MPI_Allgather_d(void *sendbuf, void *shadowsendbuf,
                       int sendcount,
                       MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                       void *recvbuf, void *shadowrecvbuf,
                       int recvcount,
                       MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                       MPI_Comm comm) ;

int MPI_Allgather_fwd(void *sendbuf,
                      int sendcount,
                      MPI_Datatype sendtype,
                      void *recvbuf,
                      int recvcount,
                      MPI_Datatype recvtype,
                      MPI_Comm comm) ;

int MPI_Allgather_bwd(void *sendbuf,
                      int sendcount,
                      MPI_Datatype sendtype,
                      void *recvbuf,
                      int recvcount,
                      MPI_Datatype recvtype,
                      MPI_Comm comm) ;

// =Allgatherv=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sendbuf
 * @param sendcnt
 * @param sendtype see \ref datatypes
 * @param recvbuf
 * @param recvcnts
 * @param displs
 * @param recvtype see \ref datatypes
 * @param comm
 * @return
 */

int MPI_Allgatherv_d(void *sendbuf, void *shadowsendbuf,
                        int sendcnt,
                        MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                        void *recvbuf, void *shadowrecvbuf,
                        int *recvcnts,
                        int *displs,
                        MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                        MPI_Comm comm) ;

int MPI_Allgatherv_fwd(void *sendbuf,
                       int sendcnt,
                       MPI_Datatype sendtype,
                       void *recvbuf,
                       int *recvcnts,
                       int *displs,
                       MPI_Datatype recvtype,
                       MPI_Comm comm) ;

int MPI_Allgatherv_bwd(void *sendbuf,
                       int sendcnt,
                       MPI_Datatype sendtype,
                       void *recvbuf,
                       int *recvcnts,
                       int *displs,
                       MPI_Datatype recvtype,
                       MPI_Comm comm) ;

// =Scatter=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sendbuf
 * @param sendcnt
 * @param sendtype see \ref datatypes
 * @param recvbuf
 * @param recvcnt
 * @param recvtype see \ref datatypes
 * @param root
 * @param comm
 * @return
 */

int MPI_Scatter_d(void *sendbuf, void *shadowsendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                     void *recvbuf, void *shadowrecvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                     int root,
                     MPI_Comm comm) ;

int MPI_Scatter_fwd(void *sendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                    MPI_Comm comm) ;

int MPI_Scatter_bwd(void *sendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                    MPI_Comm comm) ;

// =Scatterv=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sendbuf
 * @param sendcnts
 * @param displs
 * @param sendtype see \ref datatypes
 * @param recvbuf
 * @param recvcnt
 * @param recvtype see \ref datatypes
 * @param root
 * @param comm
 * @return
 */

int MPI_Scatterv_d(void *sendbuf, void *shadowsendbuf,
                      int *sendcnts,
                      int *displs,
                      MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                      void *recvbuf, void *shadowrecvbuf,
                      int recvcnt,
                      MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                      int root, MPI_Comm comm) ;

int MPI_Scatterv_fwd(void *sendbuf,
                     int *sendcnts,
                     int *displs,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) ;

int MPI_Scatterv_bwd(void *sendbuf,
                     int *sendcnts,
                     int *displs,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) ;

// =Bcast=


/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param buf
 * @param count
 * @param datatype see \ref datatypes
 * @param root
 * @param comm
 * @return
 */

int MPI_Bcast_d(void* buf, void* shadowbuf,
                   int count,
                   MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                   int root,
                   MPI_Comm comm) ;

int MPI_Bcast_fwd(void* buf,
                   int count,
                   MPI_Datatype datatype,
                   int root,
                   MPI_Comm comm) ;

int MPI_Bcast_bwd(void* buf,
                   int count,
                   MPI_Datatype datatype,
                   int root,
                   MPI_Comm comm) ;

// =Reduce=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sbuf
 * @param rbuf
 * @param count
 * @param datatype see \ref datatypes
 * @param op
 * @param root
 * @param comm
 * @return
 */

/** Tangent diff of MPI_Reduce. */
int MPI_Reduce_d(void* sbuf, void* sbufd,
                    void* rbuf, void* rbufd,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype datatyped,
                    MPI_Op op, TLM_userFunctionF* uopd,
                    int root,
                    MPI_Comm comm) ;

/** Adjoint diff of MPI_Reduce, forward sweep. */
int MPI_Reduce_fwd(void* sbuf,
                   void* rbuf,
                   int count,
                   MPI_Datatype datatype,
                   MPI_Op op,
                   int root,
                   MPI_Comm comm) ;

/** Adjoint diff of MPI_Reduce, backward sweep. */
int MPI_Reduce_bwd(void* sbuf, void* sbufb,
		   void* rbuf, void* rbufb,
		   int count,
		   MPI_Datatype datatype, MPI_Datatype datatypeb,
		   MPI_Op op, TLM_userFunctionF* uopb,
                   int root,
                   MPI_Comm comm) ;

// =Allreduce=

/**
 * all parameters as in the MPI standard with exceptions as listed
 * @param sbuf
 * @param rbuf
 * @param count
 * @param datatype see \ref datatypes
 * @param op
 * @param comm
 * @return
 */

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_d(void* sbuf, void* sbufd,
                       void* rbuf, void* rbufd,
                       int count,
                       MPI_Datatype datatype, MPI_Datatype datatyped,
                       MPI_Op op, TLM_userFunctionF* uopd,
                       MPI_Comm comm) ;

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_fwd(void* sbuf,
                       void* rbuf,
                       int count,
                       MPI_Datatype datatype,
                       MPI_Op op,
                       MPI_Comm comm) ;

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_bwd(void* sbuf, void* sbufb,
                       void* rbuf, void* rbufb,
                       int count,
                       MPI_Datatype datatype, MPI_Datatype datatypeb,
                       MPI_Op op, TLM_userFunctionF* uopb,
                       MPI_Comm comm) ;

// =Comm_size=

// =Comm_rank=

// =Comm_dup=

/**
 * In addition to MPI_Comm_dup(), creates and registers a shadow comm
 * Same as MPI_Comm_dup but manages the shadow communicators if code is differentiated in tangent mode with ST-AD with shadow variables (e.g. Tapenade)
 */
int MPI_Comm_dup_d(MPI_Comm comm, MPI_Comm *dupComm) ;

// =Comm_split=

/**
 * In addition to MPI_Comm_split(), creates and registers a shadow comm
 * Same as MPI_Comm_split but manages the shadow communicators if code is differentiated in tangent mode with ST-AD with shadow variables (e.g. Tapenade)
 */
int MPI_Comm_split_d(MPI_Comm comm, int color, int key, MPI_Comm *dupComm) ;

// =Comm_create=

/**
 * In addition to MPI_Comm_create(), creates and registers a shadow comm
 * Same as MPI_Comm_create but manages the shadow communicators if code is differentiated in tangent mode with ST-AD with shadow variables (e.g. Tapenade)
 */
int MPI_Comm_create_d(MPI_Comm comm, MPI_Group group, MPI_Comm *dupComm) ;

// =Comm_free=

/**
 * In addition to MPI_Comm_free(), frees the duplicate shadow comm
 * Same as MPI_Comm_free but manages the shadow communicators if code is differentiated in tangent mode with ST-AD with shadow variables (e.g. Tapenade)
 */
int MPI_Comm_free_d(MPI_Comm *comm) ;

#endif
