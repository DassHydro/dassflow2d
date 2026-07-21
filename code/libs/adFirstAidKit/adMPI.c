/*
#######################################################################
#  This code is largely derived from the AdjoinableMPI library        #
#  Copyright (C) 2012-2014 Laurent Hascoet, Michel Schanen, Jean Utke #
#  AdjoinableMPI was released under the MIT License                   #
#######################################################################
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <malloc.h>
#include <assert.h>
#include <string.h>
#include <stddef.h>

#include "adMPI.h"

#include "adStack.h"

// =utilities=

/** The global MPI_COMM_WORLD_D */
extern MPI_Comm adMPI_COMM_WORLD_SHADOW;

#ifdef ADMPI_FORTRANCOMPATIBLE
/** Must be defined in the fortranSupport.F of the particular AD tool */
void admpi_fortransetupbindings_() ;
#endif

int adMPI_isActiveType(MPI_Datatype datatype) {
  return (datatype==MPI_DOUBLE
          ||
          datatype==MPI_FLOAT
#ifdef ADMPI_FORTRANCOMPATIBLE
          ||
          datatype==MPI_DOUBLE_PRECISION
          ||
          datatype==MPI_REAL
          ||
          datatype==MPI_REAL8
          ||
          datatype==MPI_REAL16
#endif
          ) ;
}

/**
 * user-defined type data
 * only one instance of derivedTypeData exists at once
 * get pointer from getDTypeData, add new stuff with addDTypeData
 */
typedef struct {
  int size;
  int preAlloc;
  int* num_actives;
  /* displacements of first/last active blocks */
  MPI_Aint* first_active_blocks;
  MPI_Aint* last_active_blocks;
  /* need to know last active block length to find last active element */
  int* last_active_block_lengths;
  MPI_Datatype* derived_types;
  int* counts;
  int** arrays_of_blocklengths;
  MPI_Aint** arrays_of_displacements;
  MPI_Datatype** arrays_of_types;
  MPI_Aint* lbs;
  MPI_Aint* extents;
  /* corresponding typemaps packed for sending */
  MPI_Datatype* packed_types;
  int** arrays_of_p_blocklengths;
  MPI_Aint** arrays_of_p_displacements;
  MPI_Datatype** arrays_of_p_types;
  MPI_Aint* p_extents;
} derivedTypeData;

derivedTypeData* getDTypeData() {
  static derivedTypeData* dat = NULL;
  if (dat==NULL) {
    derivedTypeData* newdat = (derivedTypeData*)malloc(sizeof(derivedTypeData));
    newdat->size = 0;
    newdat->preAlloc = 0;
    newdat->num_actives = NULL;
    newdat->first_active_blocks = NULL;
    newdat->last_active_blocks = NULL;
    newdat->last_active_block_lengths = NULL;
    newdat->derived_types = NULL;
    newdat->counts = NULL;
    newdat->arrays_of_blocklengths = NULL;
    newdat->arrays_of_displacements = NULL;
    newdat->arrays_of_types = NULL;
    newdat->lbs = NULL;
    newdat->extents = NULL;
    newdat->packed_types = NULL;
    newdat->arrays_of_p_blocklengths = NULL;
    newdat->arrays_of_p_displacements = NULL;
    newdat->arrays_of_p_types = NULL;
    newdat->p_extents = NULL;
    dat = newdat;
  }
  return dat;
}

void releaseDTypeData() {
  int i;
  derivedTypeData* dat = getDTypeData();
  for (i=0;i<dat->size;i++) {
    free(dat->arrays_of_blocklengths[i]);
    free(dat->arrays_of_displacements[i]);
    free(dat->arrays_of_types[i]);
    free(dat->arrays_of_p_blocklengths[i]);
    free(dat->arrays_of_p_displacements[i]);
    free(dat->arrays_of_p_types[i]);
    if (dat->packed_types[i]!=MPI_DATATYPE_NULL) MPI_Type_free(dat->packed_types+i);
  }
  free(dat->num_actives);
  free(dat->first_active_blocks);
  free(dat->last_active_blocks);
  free(dat->last_active_block_lengths);
  free(dat->derived_types);
  free(dat->counts);
  free(dat->arrays_of_blocklengths);
  free(dat->arrays_of_displacements);
  free(dat->arrays_of_types);
  free(dat->lbs);
  free(dat->extents);
  free(dat->packed_types);
  free(dat->arrays_of_p_blocklengths);
  free(dat->arrays_of_p_displacements);
  free(dat->arrays_of_p_types);
  free(dat->p_extents);
  free(dat);
}

int derivedTypeIdx(MPI_Datatype datatype) {
  int i;
  derivedTypeData* dtdata = getDTypeData();
  for (i=0;i<dtdata->size;i++) {
    if (dtdata->derived_types[i]==datatype) return i;
  }
  return -1;
}

int isDerivedType(int dt_idx) {
  return dt_idx!=-1;
}

/**
 * @param dat
 * @param count
 * @param array_of_blocklengths
 * @param array_of_displacements
 * @param array_of_types
 * @param lower_bound
 * @param extent
 * @param array_of_p_blocklengths
 * @param array_of_p_displacements
 * @param array_of_p_types
 * @param p_extent
 * @param newtype
 * @param packed_type
 * addDTypeData takes derived type data and adds a new entry; returns
 * position of new type in data struct; returns -1 if struct contains no active types;
 * doubles data struct size every time there's overflow
 */
void addDTypeData(derivedTypeData* dat,
		  int count,
		  int array_of_blocklengths[],
		  MPI_Aint array_of_displacements[],
		  MPI_Datatype array_of_types[],
		  MPI_Aint lower_bound,
		  MPI_Aint extent,
		  int array_of_p_blocklengths[],
		  MPI_Aint array_of_p_displacements[],
		  MPI_Datatype array_of_p_types[],
		  MPI_Aint p_extent,
		  MPI_Datatype* newtype,
		  MPI_Datatype* packed_type) {
  assert(dat);
  int i, dt_idx;
  int num_actives=0, fst_ablk_set=0;
  MPI_Aint fst_active_blk=0, lst_active_blk=0, lst_active_blk_len=0;
  for (i=0;i<count;i++) {
    if (adMPI_isActiveType(array_of_types[i])) {
      num_actives += array_of_blocklengths[i];
      if (!fst_ablk_set) {
	fst_active_blk = array_of_displacements[i];
	fst_ablk_set = 1;
      }
      lst_active_blk = array_of_displacements[i];
      lst_active_blk_len = array_of_blocklengths[i];
      continue;
    }
    dt_idx = derivedTypeIdx(array_of_types[i]);
    if (isDerivedType(dt_idx)) {
      num_actives += dat->num_actives[dt_idx]*array_of_blocklengths[i];
      if (!fst_ablk_set) {
	fst_active_blk = array_of_displacements[i] + dat->first_active_blocks[dt_idx];
	fst_ablk_set = 1;
      }
      lst_active_blk = array_of_displacements[i] + (array_of_blocklengths[i]-1)*dat->extents[dt_idx] + dat->last_active_blocks[dt_idx];
      lst_active_blk_len = dat->last_active_block_lengths[dt_idx];
    }
  }
  if (dat->preAlloc == dat->size) {
    dat->preAlloc += 16;
    dat->num_actives = realloc(dat->num_actives, (dat->preAlloc)*sizeof(int));
    dat->first_active_blocks = realloc(dat->first_active_blocks, (dat->preAlloc)*sizeof(MPI_Aint));
    dat->last_active_blocks = realloc(dat->last_active_blocks, (dat->preAlloc)*sizeof(MPI_Aint));
    dat->last_active_block_lengths = realloc(dat->last_active_block_lengths, (dat->preAlloc)*sizeof(int));
    dat->derived_types = realloc(dat->derived_types,
				 (dat->preAlloc)*sizeof(MPI_Datatype));
    dat->counts = realloc(dat->counts, (dat->preAlloc)*sizeof(int));
    dat->arrays_of_blocklengths = realloc(dat->arrays_of_blocklengths,
					  (dat->preAlloc)*sizeof(int*));
    dat->arrays_of_displacements = realloc(dat->arrays_of_displacements,
					   (dat->preAlloc)*sizeof(MPI_Aint*));
    dat->arrays_of_types = realloc(dat->arrays_of_types,
				   (dat->preAlloc)*sizeof(MPI_Datatype*));
    dat->lbs = realloc(dat->lbs, (dat->preAlloc)*sizeof(MPI_Aint));
    dat->extents = realloc(dat->extents, (dat->preAlloc)*sizeof(MPI_Aint));
    dat->packed_types = realloc(dat->packed_types,
				(dat->preAlloc)*sizeof(MPI_Datatype));
    dat->arrays_of_p_blocklengths = realloc(dat->arrays_of_p_blocklengths,
					    (dat->preAlloc)*sizeof(int*));
    dat->arrays_of_p_displacements = realloc(dat->arrays_of_p_displacements,
					     (dat->preAlloc)*sizeof(MPI_Aint*));
    dat->arrays_of_p_types = realloc(dat->arrays_of_p_types,
				     (dat->preAlloc)*sizeof(MPI_Datatype*));
    dat->p_extents = realloc(dat->p_extents, (dat->preAlloc)*sizeof(MPI_Aint));
  }
  dat->num_actives[dat->size] = num_actives;
  dat->first_active_blocks[dat->size] = fst_active_blk;
  dat->last_active_blocks[dat->size] = lst_active_blk;
  dat->last_active_block_lengths[dat->size] = lst_active_blk_len;
  dat->derived_types[dat->size] = *newtype;
  dat->counts[dat->size] = count;
  dat->arrays_of_blocklengths[dat->size] = malloc(count*sizeof(int));
  memcpy(dat->arrays_of_blocklengths[dat->size], array_of_blocklengths, count*sizeof(int));
  dat->arrays_of_displacements[dat->size] = malloc(count*sizeof(MPI_Aint));
  memcpy(dat->arrays_of_displacements[dat->size], array_of_displacements, count*sizeof(MPI_Aint));
  dat->arrays_of_types[dat->size] = malloc(count*sizeof(MPI_Datatype));
  memcpy(dat->arrays_of_types[dat->size], array_of_types, count*sizeof(MPI_Datatype));
  dat->lbs[dat->size] = lower_bound;
  dat->extents[dat->size] = extent;
  dat->packed_types[dat->size] = *packed_type;
  dat->arrays_of_p_blocklengths[dat->size] = malloc(count*sizeof(int));
  memcpy(dat->arrays_of_p_blocklengths[dat->size], array_of_p_blocklengths, count*sizeof(int));
  dat->arrays_of_p_displacements[dat->size] = malloc(count*sizeof(MPI_Aint));
  memcpy(dat->arrays_of_p_displacements[dat->size], array_of_p_displacements, count*sizeof(MPI_Aint));
  dat->arrays_of_p_types[dat->size] = malloc(count*sizeof(MPI_Datatype));
  memcpy(dat->arrays_of_p_types[dat->size], array_of_p_types, count*sizeof(MPI_Datatype));
  dat->p_extents[dat->size] = p_extent;
  dat->size += 1;
}

/**
 * user-defined reduction op data
 * only one instance of userDefinedOpData exists at once
 * get pointer from getUOpData, add new stuff with addUOpData
 */
typedef struct {
  int size;
  int preAlloc;
  MPI_Op* ops;
  MPI_User_function** functions;
  int* commutes;
} userDefinedOpData;

userDefinedOpData* getUOpData() {
  static userDefinedOpData* dat = NULL;
  if (dat==NULL) {
    userDefinedOpData* newdat = (userDefinedOpData*)malloc(sizeof(userDefinedOpData));
    newdat->size = 0;
    newdat->preAlloc = 0;
    newdat->ops = NULL;
    newdat->functions = NULL;
    newdat->commutes = NULL;
    dat = newdat;
  }
  return dat;
}

void releaseUOpData() {
  userDefinedOpData* dat = getUOpData();
  free(dat->ops);
  free(dat->functions);
  free(dat->commutes);
  free(dat);
}

int userDefinedOpIdx(MPI_Op op) {
  int i;
  userDefinedOpData* uopdata = getUOpData();
  for (i=0;i<uopdata->size;i++) {
    if (uopdata->ops[i]==op) return i;
  }
  return -1;
}

int isUserDefinedOp(int uop_idx) {
  return uop_idx!=-1;
}

/**
 * @param dat
 * @param op a user-defined operation
 * @param function
 * @param commute
 * takes user-defined op  and adds a new entry;
 * doubles data struct size every
 * time there's overflow
 */
void addUOpData(userDefinedOpData* dat,
		MPI_Op* op,
		MPI_User_function* function,
		int commute) {
  assert(dat);
  if (dat->preAlloc == dat->size) {
    dat->preAlloc += 16;
    dat->ops = realloc(dat->ops,(dat->preAlloc)*sizeof(MPI_Op));
    dat->functions = realloc(dat->functions,(dat->preAlloc)*sizeof(MPI_User_function*));
    dat->commutes = realloc(dat->commutes,(dat->preAlloc)*sizeof(int));
  }
  dat->ops[dat->size] = *op;
  dat->functions[dat->size] = function;
  dat->commutes[dat->size] = commute;
  dat->size += 1;
}

void* ADMPI_BOTTOM_F = NULL ;
void* ADMPI_IN_PLACE_F = NULL ;
void* ADMPI_STATUS_IGNORE_F = NULL ;
void* ADMPI_STATUSES_IGNORE_F = NULL ;
void* ADMPI_ERRCODES_IGNORE_F = NULL ;
void* ADMPI_ARGV_NULL_F = NULL ;
void* ADMPI_ARGVS_NULL_F = NULL ;

#ifdef ADMPI_FORTRANCOMPATIBLE
void admpi_sendfortranbindings_(void *bt_f, void *ip_f, void *si_f, void *ssi_f,
                                      void *ei_f, void *av_f, void *avs_f) {
  ADMPI_BOTTOM_F = bt_f ;
  ADMPI_IN_PLACE_F = ip_f ;
  ADMPI_STATUS_IGNORE_F = si_f ;
  ADMPI_STATUSES_IGNORE_F = ssi_f ;
  ADMPI_ERRCODES_IGNORE_F = ei_f ;
  ADMPI_ARGV_NULL_F = av_f ;
  ADMPI_ARGVS_NULL_F = avs_f ;
}
#endif

/* Variables that hold, on the C side, the value of binding constants from the Fortran side. */
extern void* ADMPI_BOTTOM_F ;
extern void* ADMPI_IN_PLACE_F ;
extern void* ADMPI_STATUS_IGNORE_F ;
extern void* ADMPI_STATUSES_IGNORE_F ;
extern void* ADMPI_ERRCODES_IGNORE_F ;
extern void* ADMPI_ARGV_NULL_F ;
extern void* ADMPI_ARGVS_NULL_F ;

struct adMPI_Request_stack {
  struct adMPI_Request_stack *next_p;
  void *buf ;
  void *adjointBuf ;
  int count ;
  MPI_Datatype datatype ;
  int endPoint ;
  int tag ;
  MPI_Comm comm;
  enum adMPI_Request_origin_E origin;
} ;

static struct adMPI_Request_stack* requestStackTop = NULL ;

struct RequestListItem { 
  struct adMPI_Request_S admpiRequest; /*[llh] I'd rather put *admpiRequest to not copy */
  struct RequestListItem *next_p;
  struct RequestListItem *prev_p;
};

static struct RequestListItem* requestListHead=0;
static struct RequestListItem* requestListTail=0;
static struct RequestListItem* unusedRequestStack=0;

struct ADMPI_ShadowComm_list {
  struct ADMPI_ShadowComm_list *next_p;
  MPI_Comm comm ;
  MPI_Comm shadowComm ;
} ;

struct ADMPI_ShadowComm_list * adMPI_SHADOWCOMMLIST = NULL ;

/**
 * Register the info that the shadow communicator "dupComm"
 * has been created for the new communicator "comm"
 */
void adMPI_addShadowComm(MPI_Comm comm, MPI_Comm dupComm) {
  struct ADMPI_ShadowComm_list *newCell =
    (struct ADMPI_ShadowComm_list *)malloc(sizeof(struct ADMPI_ShadowComm_list)) ;
  newCell->next_p = adMPI_SHADOWCOMMLIST ;
  newCell->comm = comm;
  newCell->shadowComm = dupComm;
  adMPI_SHADOWCOMMLIST = newCell;
}

/**
 * Get the shadow communicator used to separate the communication graph of
 * (tangent-)diff variables from the communication graph of original variables
 */
MPI_Comm adMPI_getShadowComm(MPI_Comm comm) {
  struct ADMPI_ShadowComm_list * inShadowCommList = adMPI_SHADOWCOMMLIST ;
  while (inShadowCommList!=NULL && inShadowCommList->comm!=comm) {
    inShadowCommList = inShadowCommList->next_p ;
  }
  if (inShadowCommList) {
    return inShadowCommList->shadowComm ;
  } else {
    /* Nothing found about "comm": this is wrong!! fallback return comm */
    return comm ;
  }
}

/**
 * Removes the info about the shadow communicator associated to "comm".
 */
void adMPI_delShadowComm(MPI_Comm comm) {
  struct ADMPI_ShadowComm_list ** toinShadowCommList = &adMPI_SHADOWCOMMLIST ;
  while (*toinShadowCommList!=NULL && (*toinShadowCommList)->comm!=comm) {
    toinShadowCommList = &((*toinShadowCommList)->next_p) ;
  }
  if (*toinShadowCommList!=NULL) {
    struct ADMPI_ShadowComm_list *cell = *toinShadowCommList ;
    *toinShadowCommList = cell->next_p ;
    free(cell);
  }
}

void* adMPI_allocateTempBuf(int adjointCount, MPI_Datatype datatype, MPI_Comm comm) {
  size_t s=0;
  int dt_idx = derivedTypeIdx(datatype);
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION)
    s=sizeof(double);
  else if (datatype==MPI_FLOAT || datatype==MPI_REAL)
    s=sizeof(float);
  else if (datatype==MPI_REAL8)
    s=8;
  else if (datatype==MPI_REAL16)
    s=16;
  else if (isDerivedType(dt_idx))
    s = getDTypeData()->p_extents[dt_idx];
  else
    MPI_Abort(comm, MPI_ERR_TYPE);
  return (void*)malloc(adjointCount*s);
}

void adMPI_releaseAdjointTempBuf(void *tempBuf) { 
  free(tempBuf) ;
}

void* adMPI_allocateTempActiveBuf(int count,
                                        MPI_Datatype datatype,
                                        MPI_Comm comm) {
  MPI_Aint lb, extent ;
  int rc = MPI_Type_get_extent(datatype, &lb, &extent) ;
  assert(rc==MPI_SUCCESS);
  void* ptr = NULL ;
  int size = ((char*)extent)-((char*)lb) ;
  ptr = malloc(count*size) ;
  assert(ptr);
  return ptr;
}

void *adMPI_copyActiveBuf(void* source,
                                void* target,
                                int count,
                                MPI_Datatype datatype,
                                MPI_Comm comm) {
  MPI_Aint lb,extent ;
  int rc = MPI_Type_get_extent(datatype, &lb, &extent) ;
  assert(rc==MPI_SUCCESS);
  int size = extent - lb ;
  memcpy(target, source, count*size) ;
  return source;
}

/**
 * Push the contents of buffer somewhere
 */
void adMPI_pushBuffer(int count, MPI_Datatype datatype, MPI_Comm comm,
                            void* buffer) {
    MPI_Aint lb, extent ;
    MPI_Type_get_extent(datatype,&lb,&extent);
    int length = count*(int)extent ;
    pushCharacterArray((char*)buffer, length) ;
}

/**
 * Pop the contents of buffer from somewhere
 */
void adMPI_popBuffer(int count, MPI_Datatype datatype, MPI_Comm comm,
                           void* buffer) {
    MPI_Aint lb, extent ;
    MPI_Type_get_extent(datatype,&lb,&extent);
    int length = count*(int)extent ;
    popCharacterArray((char*)buffer, length) ;
}

void adMPI_pushRequest(struct adMPI_Request_S  *admpiRequest) { 
  struct adMPI_Request_stack* newTop =
    (struct adMPI_Request_stack*)malloc(sizeof(struct adMPI_Request_stack)) ;
  newTop->next_p = requestStackTop ;
  newTop->buf = admpiRequest->buf ;
  newTop->adjointBuf = admpiRequest->adjointBuf ;
  newTop->count = admpiRequest->count ;
  newTop->datatype = admpiRequest->datatype ;
  newTop->endPoint = admpiRequest->endPoint ;
  newTop->tag = admpiRequest->tag ;
  newTop->comm = admpiRequest->comm ;
  newTop->origin = admpiRequest->origin ;
  requestStackTop = newTop ;
}

void adMPI_popRequest(struct adMPI_Request_S  *admpiRequest) { 
  struct adMPI_Request_stack* oldTop = requestStackTop ;
  admpiRequest->buf = oldTop->buf ;
  admpiRequest->adjointBuf = oldTop->adjointBuf ;
  admpiRequest->count = oldTop->count ;
  admpiRequest->datatype = oldTop->datatype ;
  admpiRequest->endPoint = oldTop->endPoint ;
  admpiRequest->tag = oldTop->tag ;
  admpiRequest->comm = oldTop->comm ;
  admpiRequest->origin = oldTop->origin ;
  requestStackTop = oldTop->next_p ;
  free(oldTop) ;
}

void adMPI_tangentMultiply(int count, MPI_Datatype datatype, MPI_Comm comm,
                                 void *source, void *tangentSource,
                                 void* target, void* tangentTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtd = (double*)tangentTarget ;
    double* src = (double*)source ;
    double* srcd = (double*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgtd) tgtd[i] = tgtd[i]*src[i] + tgt[i]*srcd[i] ;
      tgt[i] = tgt[i]*src[i] ;
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtd = (float*)tangentTarget ;
    float* src = (float*)source ;
    float* srcd = (float*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgtd) tgtd[i] = tgtd[i]*src[i] + tgt[i]*srcd[i] ;
      tgt[i] = tgt[i]*src[i] ;
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/** This is the tangent of assignment target=MIN(source,target).
 * If targetd is NULL, does only the original assignment */
void adMPI_tangentMin(int count, MPI_Datatype datatype, MPI_Comm comm,
                            void *source, void *tangentSource,
                            void* target, void* tangentTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtd = (double*)tangentTarget ;
    double* src = (double*)source ;
    double* srcd = (double*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgt[i] > src[i]) {
        if (tgtd) tgtd[i] = srcd[i] ;
        tgt[i] = src[i] ;
      }
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtd = (float*)tangentTarget ;
    float* src = (float*)source ;
    float* srcd = (float*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgt[i] > src[i]) {
        if (tgtd) tgtd[i] = srcd[i] ;
        tgt[i] = src[i] ;
      }
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/** This is the tangent of assignment target=MAX(source,target).
 * If targetd is NULL, does only the original assignment */
void adMPI_tangentMax(int count, MPI_Datatype datatype, MPI_Comm comm,
                            void *source, void *tangentSource,
                            void* target, void* tangentTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtd = (double*)tangentTarget ;
    double* src = (double*)source ;
    double* srcd = (double*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgt[i] < src[i]) {
        if (tgtd) tgtd[i] = srcd[i] ;
        tgt[i] = src[i] ;
      }
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtd = (float*)tangentTarget ;
    float* src = (float*)source ;
    float* srcd = (float*)tangentSource ;
    for (i=0 ; i<count ; ++i) {
      if (tgt[i] < src[i]) {
        if (tgtd) tgtd[i] = srcd[i] ;
        tgt[i] = src[i] ;
      }
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/**
 * This is the adjoint of assignment target=source*target
 */
void adMPI_adjointMultiply(int count, MPI_Datatype datatype, MPI_Comm comm,
                                 void *source, void *adjointSource,
                                 void* target, void* adjointTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtb = (double*)adjointTarget ;
    double* src = (double*)source ;
    double* srcb = (double*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      srcb[i] += tgt[i]*tgtb[i] ;
      tgtb[i] *= src[i] ;
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtb = (float*)adjointTarget ;
    float* src = (float*)source ;
    float* srcb = (float*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      srcb[i] += tgt[i]*tgtb[i] ;
      tgtb[i] *= src[i] ;
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/**
 * This is the adjoint of assignment target=MIN(source,target)
 */
void adMPI_adjointMin(int count, MPI_Datatype datatype, MPI_Comm comm,
                                 void *source, void *adjointSource,
                                 void* target, void* adjointTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtb = (double*)adjointTarget ;
    double* src = (double*)source ;
    double* srcb = (double*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      if (src[i]<tgt[i]) {
        srcb[i] += tgtb[i] ;
        tgtb[i] = 0.0 ;
      }
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtb = (float*)adjointTarget ;
    float* src = (float*)source ;
    float* srcb = (float*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      if (src[i]<tgt[i]) {
        srcb[i] += tgtb[i] ;
        tgtb[i] = 0.0 ;
      }
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/**
 * This is the adjoint of assignment target=MAX(source,target)
 */
void adMPI_adjointMax(int count, MPI_Datatype datatype, MPI_Comm comm,
                                 void *source, void *adjointSource,
                                 void* target, void* adjointTarget) {
  int i ;
  if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double* tgt = (double*)target ;
    double* tgtb = (double*)adjointTarget ;
    double* src = (double*)source ;
    double* srcb = (double*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      if (src[i]>tgt[i]) {
        srcb[i] += tgtb[i] ;
        tgtb[i] = 0.0 ;
      }
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float* tgt = (float*)target ;
    float* tgtb = (float*)adjointTarget ;
    float* src = (float*)source ;
    float* srcb = (float*)adjointSource ;
    for (i=0 ; i<count ; ++i) {
      if (src[i]>tgt[i]) {
        srcb[i] += tgtb[i] ;
        tgtb[i] = 0.0 ;
      }
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/**
 * Increment the given buffer "target", which holds an adjoint variable,
 * with the given additional adjoint value found in "source".
 */
void adMPI_incrementAdjoint(int adjointCount, MPI_Datatype datatype, MPI_Comm comm, void* target, void *source) { 
  int dt_idx = derivedTypeIdx(datatype);
  if (isUserDefinedOp(dt_idx)) {
    derivedTypeData* dat = getDTypeData();
    MPI_Aint lb, extent ;
    MPI_Type_get_extent(datatype,&lb,&extent);
    MPI_Aint*   fieldOffsets = dat->arrays_of_displacements[dt_idx] ;
    int*   fieldBlocklengths = dat->arrays_of_blocklengths[dt_idx] ;
    MPI_Datatype* fieldTypes = dat->arrays_of_types[dt_idx] ;
    int nbfields = dat->counts[dt_idx] ;
    int i,j ;
    for (i=0 ; i<adjointCount ; ++i) {
      for (j=0 ; j<nbfields ; ++j) {
        adMPI_incrementAdjoint(fieldBlocklengths[j], fieldTypes[j], comm,
                                     target+fieldOffsets[j], source+fieldOffsets[j]) ;
      }
      target += extent ;
      source += extent ;
    }
  } else if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double *vb = (double *)target ;
    double *nb = (double *)source ;
    int i ;
    for (i=0 ; i<adjointCount ; ++i) {
      *vb = *vb + *nb ;
      ++vb ;
      ++nb ;
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float *vb = (float *)target ;
    float *nb = (float *)source ;
    int i ;
    for (i=0 ; i<adjointCount ; ++i) {
      *vb = *vb + *nb ;
      ++vb ;
      ++nb ;
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

/**
 * Reset to zero the given buffer "target", which holds an adjoint variable.
 */
void adMPI_nullifyAdjoint(int adjointCount, MPI_Datatype datatype, MPI_Comm comm,
                                void* target) {
  int dt_idx = derivedTypeIdx(datatype);
  if (isUserDefinedOp(dt_idx)) {
    derivedTypeData* dat = getDTypeData();
    MPI_Aint lb, extent ;
    MPI_Type_get_extent(datatype,&lb,&extent);
    MPI_Aint*   fieldOffsets = dat->arrays_of_displacements[dt_idx] ;
    int*   fieldBlocklengths = dat->arrays_of_blocklengths[dt_idx] ;
    MPI_Datatype* fieldTypes = dat->arrays_of_types[dt_idx] ;
    int nbfields = dat->counts[dt_idx] ;
    int i,j ;
    for (i=0 ; i<adjointCount ; ++i) {
      for (j=0 ; j<nbfields ; ++j) {
        adMPI_nullifyAdjoint(fieldBlocklengths[j], fieldTypes[j], comm,
                                   target+fieldOffsets[j]) ;
      }
      target += extent ;
    }
  } else if (datatype==MPI_DOUBLE || datatype==MPI_DOUBLE_PRECISION || datatype==MPI_REAL8) {
    double *vb = (double *)target ;
    int i ;
    for (i=0 ; i<adjointCount ; ++i) {
      *vb = 0.0 ;
      ++vb ;
    }
  } else if (datatype==MPI_FLOAT || datatype==MPI_REAL) {
    float *vb = (float *)target ;
    int i ;
    for (i=0 ; i<adjointCount ; ++i) {
      *vb = 0.0 ;
      ++vb ;
    }
  } else {
    printf("Unknown size of MPI Type %i\n", datatype) ;
    MPI_Abort(comm, MPI_ERR_TYPE);
  }
}

MPI_Datatype adMPI_FW_rawType(MPI_Datatype datatype) {
  int dt_idx = derivedTypeIdx(datatype);
  if (datatype==MPI_DOUBLE) return MPI_DOUBLE;
  else if (datatype==MPI_FLOAT) return MPI_FLOAT;
  else if (isDerivedType(dt_idx)) return getDTypeData()->packed_types[dt_idx];
  else return datatype;
}

MPI_Datatype adMPI_BW_rawType(MPI_Datatype datatype) {
  int dt_idx = derivedTypeIdx(datatype);
  if (datatype==MPI_DOUBLE) return MPI_DOUBLE;
  else if (datatype==MPI_FLOAT) return MPI_FLOAT;
  else if (isDerivedType(dt_idx)) return MPI_DOUBLE;
  else return datatype;
}

/**
 * \param admpiRequest is added (by deep copy) to the internal bookkeeping
 * using the already set value of member plainRequest as key
 */
void adMPI_putRequest(struct adMPI_Request_S *admpiRequest) { 
  struct RequestListItem* inList_p=0;
  if (! unusedRequestStack) { 
    unusedRequestStack=(struct RequestListItem*)malloc(sizeof(struct RequestListItem));
    assert(unusedRequestStack);
    unusedRequestStack->next_p=0; 
    unusedRequestStack->prev_p=0; 
  }
  /* get it from the unused stack */
  inList_p=unusedRequestStack;
  unusedRequestStack=inList_p->prev_p;
  inList_p->prev_p=0;
  /* add it to the list */
  if (!requestListHead) requestListHead=inList_p;
  if (requestListTail) { 
    requestListTail->next_p=inList_p;
    inList_p->prev_p=requestListTail;
  }
  requestListTail=inList_p;
  inList_p->admpiRequest=*admpiRequest;
}

/**
 * \param request is used as key to look up the associated adMPI_Request_S instance which then is deep copied 
 * \param admpiRequest pointer to the structure into which the values are copied 
 * the information is removed from the internal bookkeeping data
 */
void adMPI_getRequest(MPI_Request *request, struct adMPI_Request_S *admpiRequest) { 
  struct RequestListItem* inList_p=requestListHead;
  while(inList_p) { 
    if (inList_p->admpiRequest.plainRequest==*request) break;
    inList_p=inList_p->next_p;
  }
  assert(inList_p);
  *admpiRequest = inList_p->admpiRequest;
  /* remove inList_p from the list */
  if (requestListHead==inList_p) { 
    requestListHead=inList_p->next_p;
    if (requestListHead) requestListHead->prev_p=0;
    inList_p->next_p=0;
  }
  if (requestListTail==inList_p) { 
    requestListTail=inList_p->prev_p;
    if (requestListTail) requestListTail->next_p=0;
    inList_p->prev_p=0;
  }
  if (inList_p->next_p && inList_p->prev_p) {
    inList_p->prev_p->next_p=inList_p->next_p;
    inList_p->next_p->prev_p=inList_p->prev_p;
    inList_p->next_p=0; 
    inList_p->prev_p=0;
  }
  /* add it to the unused stack */
  if (unusedRequestStack) { 
    inList_p->prev_p=unusedRequestStack;
  }
  unusedRequestStack=inList_p;
}

/** Common implementation for tangent and adjoint of reductions (Reduce and Allreduce).
 *  Manages arbitrary reduction operation, with optimized implementation for MPI_SUM.
 *  Used so far only for the shadowed case (i.e. Association-by-Name, i.e. Tapenade)
 *  but could also be used for Association-by-Address in principle.
 *  When no adjoint is required (sbufb null), pass split_mode 0, otherwise:
 *   -pass split_mode 0 to obtain a joint adjoint reduction-driver and a joint adjoint reduction-op.
 *   -pass split_mode 1 to obtain a split adjoint reduction-driver and a joint adjoint reduction-op.
 *   -pass split_mode 2 to obtain a split adjoint reduction-driver and a split adjoint reduction-op.
 *  Pass all_mode 1 to perform an "Allreduce" ; in that case we suggest that caller passes root=0.
 */
int PEDESTRIAN_ADMPI_Reduce(void* sbuf, void* sbufd, void* sbufb,
                    void* rbuf, void* rbufd, void* rbufb,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype datatyped, MPI_Datatype datatypeb,
                    MPI_Op op, TLM_userFunctionF* uopd, ADJ_userFunctionF* uopb,
                    int split_mode, int all_mode,
                    int root,
                    MPI_Comm comm) {
  if (count == 0) return MPI_SUCCESS;
  int rc, rank ;
  MPI_Comm_rank(comm,&rank) ;
  int reduceTgt = (rbufd!=NULL) ;
  int reduceAdj = (rbufb!=NULL) ;
  MPI_Comm shadowcomm = comm ;
  if (reduceTgt)
    shadowcomm = adMPI_getShadowComm(comm) ;

  if (uopd || uopb || op!=MPI_SUM) {

    int uop_idx = userDefinedOpIdx(op);
    userDefinedOpData* uopdata = (isUserDefinedOp(uop_idx)?getUOpData():NULL) ;
    int is_commutative = (uopdata?uopdata->commutes[uop_idx]:1) ;
    if (uopdata) {
      if (reduceTgt) assert(uopd) ;
      if (reduceAdj) assert(uopb) ;
    }

    void *exch_buf=NULL ;
    int switched = 0 ;

    void *obuf=NULL, *obufd=NULL ;

    int dt_idx = derivedTypeIdx(datatype);
    MPI_Aint lb = (isDerivedType(dt_idx)?getDTypeData()->lbs[dt_idx]:0) ;
    int dt_idxd = 0 ;
    MPI_Aint lbd = 0 ;
    obuf =
      adMPI_allocateTempActiveBuf(count,datatype,comm);
    obuf = (void*)((char*)obuf - lb);
    if (reduceTgt) {
      dt_idxd = derivedTypeIdx(datatyped);
      lbd = (isDerivedType(dt_idxd)?getDTypeData()->lbs[dt_idxd]:0) ;
      obufd =
        adMPI_allocateTempActiveBuf(count,datatyped,shadowcomm);
      obufd = (void*)((char*)obufd - lbd);
    }

    void *orig_rbuf=NULL, *orig_rbufd=NULL ;
    if (all_mode) {
      orig_rbuf = rbuf ;
      if (reduceTgt) {orig_rbufd = rbufd ;}
    }

    if (sbuf==MPI_IN_PLACE) {
      if (rank != root) {
        exch_buf = rbuf ;
        rbuf = adMPI_allocateTempActiveBuf(count,datatype,comm);
        rbuf = (void*)((char*)rbuf - lb);
        adMPI_copyActiveBuf(exch_buf, rbuf, count, datatype, comm);
        if (reduceTgt) {
          exch_buf = rbufd ;
          rbufd = adMPI_allocateTempActiveBuf(count,datatyped,shadowcomm);
          rbufd = (void*)((char*)rbufd - lbd);
          adMPI_copyActiveBuf(exch_buf, rbufd, count, datatyped, shadowcomm);
        }
      }
    } else {  /* Standard case: sbuf != MPI_IN_PLACE */
      if (rank != root) {
        rbuf = adMPI_allocateTempActiveBuf(count,datatype,comm);
        rbuf = (void*)((char*)rbuf - lb);
        if (reduceTgt) {
          rbufd = adMPI_allocateTempActiveBuf(count,datatyped,shadowcomm);
          rbufd = (void*)((char*)rbufd - lbd);
        }
      }
      adMPI_copyActiveBuf(sbuf, rbuf, count, datatype, comm);
      if (reduceTgt)
        adMPI_copyActiveBuf(sbufd, rbufd, count, datatyped, shadowcomm);
    }

    MPI_Status status;
    int comm_size ;
    MPI_Comm_size(comm,&comm_size);
    int other, action ;
    int maskup = 0xffffffff ;
    int mask   = 0x1;

    if (!reduceAdj || split_mode==0) {
     while (mask < comm_size) {
      if ((rank&mask) == 0) { /* Typical action is RECV */
        other = (rank==root?root&maskup:rank) | mask ;
        if (other >= comm_size)
          action = 0/*NOACTION*/ ;
        else if ((other&maskup) == (root&maskup)) {
          other = root ;
          action = 1/*SEND*/ ;
        } else {
          action = -1/*RECV*/ ;
        }
      } else { /* mask&rank == 1: Typical action is SEND */
        other = (rank==root?root&maskup:rank) & ~mask ;
        if ((other&maskup) == (root&maskup)) other = root ;
        if (rank==root)
          action = -1/*RECV*/ ;
        else
          action = 1/*SEND*/ ;
      }
      maskup = maskup & ~mask ;
      mask<<=1;

      if (action==1/*SEND*/) {
        /* TODO Not sure this "(..., 11, comm)" is correct. Would better use shadowcomm ? */
        rc = MPI_Send(rbuf, count, datatype, other, 11, comm) ;
        assert(rc==MPI_SUCCESS);
        if (reduceTgt) {
          rc = MPI_Send(rbufd, count, datatyped, other, 11, shadowcomm) ;
          assert(rc==MPI_SUCCESS);
        }
	break;
      } else if (action==-1/*RECV*/) {
        rc = MPI_Recv(obuf, count, datatype, other, 11, comm, &status);
        assert(rc==MPI_SUCCESS);
        if (reduceTgt) {
          rc = MPI_Recv(obufd, count, datatyped, other, 11, shadowcomm, &status);
          assert(rc==MPI_SUCCESS);
        }
        if (is_commutative || (other<rank)) {
          /* Save obuf and rbuf for future use in the adjoint sweep */
          adMPI_pushBuffer(count,datatype,comm,obuf) ;
          if (split_mode!=2) {
            adMPI_pushBuffer(count,datatype,comm,rbuf) ;
          }
          if (isUserDefinedOp(uop_idx)) {
            if (reduceTgt)
              (*uopd)(obuf, obufd, rbuf, rbufd, &count, &datatype, &datatyped);
            else
              (*(uopdata->functions[uop_idx]))(obuf, rbuf, &count, &datatype) ;
          } else {
            if (op==MPI_PROD) {
              adMPI_tangentMultiply
                (count, datatype, comm, obuf, obufd, rbuf, rbufd) ;
            } else if (op==MPI_MIN) {
              adMPI_tangentMin
                (count, datatype, comm, obuf, obufd, rbuf, rbufd) ;
            } else if (op==MPI_MAX) {
              adMPI_tangentMax
                (count, datatype, comm, obuf, obufd, rbuf, rbufd) ;
            } else {
              printf(__FILE__ ": tangent adMPI reduction not yet implemented for std op==%i\n",uop_idx) ;
            }
          }
        } else {
          /* Save obuf and rbuf for future use in the adjoint sweep */
          adMPI_pushBuffer(count,datatype,comm,rbuf) ;
          if (split_mode!=2) {
            adMPI_pushBuffer(count,datatype,comm,obuf) ;
          }
          if (reduceTgt)
            (*uopd)(rbuf, rbufd, obuf, obufd, &count, &datatype, &datatyped);
          else
            (*(uopdata->functions[uop_idx]))(rbuf, obuf, &count, &datatype) ;
          exch_buf = obuf ; obuf = rbuf ; rbuf = exch_buf ;
          if (reduceTgt) {
            exch_buf = obufd ; obufd = rbufd ; rbufd = exch_buf ;
          }
          switched = ~switched ;
        }
      }
     }

     if (switched) {
      if (!reduceAdj) { /* Adjoint joint mode does not need to return a correct rbuf */
        exch_buf = obuf ; obuf = rbuf ; rbuf = exch_buf ;
        if (rank==root)
          adMPI_copyActiveBuf(obuf, rbuf, count, datatype, comm) ;
      }
      if (reduceTgt) {
        exch_buf = obufd ; obufd = rbufd ; rbufd = exch_buf ;
        if (rank==root)
          adMPI_copyActiveBuf(obufd, rbufd, count, datatyped, shadowcomm) ;
      }
     }

    }

    if (all_mode) {
      if (!reduceAdj) { /* Adjoint joint mode does not need to return a correct rbuf */
        /* ?? adMPI_pushBuffer(count,datatype,comm,orig_rbuf) ; */
        rc=MPI_Bcast(orig_rbuf, count, datatype, root, comm) ;
        assert(rc==MPI_SUCCESS) ;
      }
      if (reduceTgt) {
        rc=MPI_Bcast(orig_rbufd, count, datatyped, root, comm) ;
        assert(rc==MPI_SUCCESS) ;
      }
    }

    if (split_mode!=0) {
      if (!reduceAdj) {
        adMPI_pushBuffer(1,MPI_INT,comm,&switched) ;
        adMPI_pushBuffer(1,MPI_INT,comm,&maskup) ;
        adMPI_pushBuffer(1,MPI_INT,comm,&mask) ;
      } else {
        adMPI_popBuffer(1,MPI_INT,comm,&mask) ;
        adMPI_popBuffer(1,MPI_INT,comm,&maskup) ;
        adMPI_popBuffer(1,MPI_INT,comm,&switched) ;
      }
    }

    if (all_mode) {
      if (reduceAdj) {
        if (rank != root)
          rc=MPI_Reduce(rbufb, rbufb, count, datatypeb, MPI_SUM, root, comm) ;
        else
          rc=MPI_Reduce(MPI_IN_PLACE, rbufb, count, datatypeb, MPI_SUM, root, comm) ;
        assert(rc==MPI_SUCCESS) ;
        if (rank != root)
          adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
      }
    }

    if (reduceAdj) {
      void *rbufb_initial=NULL ;
      int dt_idxb = derivedTypeIdx(datatypeb);
      MPI_Aint lbb = (isDerivedType(dt_idxb)?getDTypeData()->lbs[dt_idxb]:0) ;
      void *obufb =
        adMPI_allocateTempActiveBuf(count,datatypeb,comm);
      obufb = (void*)((char*)obufb - lbb);
      adMPI_nullifyAdjoint(count,datatypeb,comm,obufb);
      if (rank != root) {
        rbufb_initial = rbufb ; 
        rbufb = adMPI_allocateTempActiveBuf(count,datatypeb,comm);
        rbufb = (void*)((char*)rbufb - lbb);
        adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
      }
      if (switched && rank==root) {
        adMPI_copyActiveBuf(rbufb, obufb, count, datatypeb, comm) ;
        adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
        exch_buf = obufb ; obufb = rbufb ; rbufb = exch_buf ;
      }
      while (mask!=0x1) {
        mask>>=1 ;
        maskup = maskup | mask ;
        if ((rank&mask) == 0) { /* Typical action fw is RECV */
          other = (rank==root?root&maskup:rank) | mask ;
          if (other >= comm_size)
            action = 0/*NOACTION*/ ;
          else if ((other&maskup) == (root&maskup)) {
            other = root ;
            action = 1/* fw SEND*/ ;
          } else {
            action = -1/* fw RECV*/ ;
          }
        } else { /* mask&rank == 1: Typical action is SEND */
          other = (rank==root?root&maskup:rank) & ~mask ;
          if ((other&maskup) == (root&maskup)) other = root ;
          if (rank==root)
            action = -1/* fw RECV*/ ;
          else
            action = 1/* fw SEND*/ ;
        }

        if (action==1/* fw SEND*/) {
          rc = MPI_Recv(obufb, count, datatypeb, other, 11, comm, &status);
          assert(rc==MPI_SUCCESS);
          adMPI_incrementAdjoint(count,datatypeb,comm,rbufb,obufb) ;
        } else if (action==-1/* fw RECV*/) {
          if (is_commutative || (other<rank)) {
            /* Retrieve obuf and rbuf for the adjoint call */
            if (split_mode!=2) {
              adMPI_popBuffer(count,datatype,comm,rbuf) ;
            }
            adMPI_popBuffer(count,datatype,comm,obuf) ;
            adMPI_nullifyAdjoint(count,datatypeb,comm,obufb);
            if (isUserDefinedOp(uop_idx)) {
              (*uopb)(obuf, obufb, rbuf, rbufb, &count, &datatype, &datatypeb) ;
            } else {
              if (op==MPI_PROD) {
                adMPI_adjointMultiply
                  (count, datatype, comm, obuf, obufb, rbuf, rbufb) ;
              } else if (op==MPI_MIN) {
                adMPI_adjointMin
                  (count, datatype, comm, obuf, obufb, rbuf, rbufb) ;
              } else if (op==MPI_MAX) {
                adMPI_adjointMax
                  (count, datatype, comm, obuf, obufb, rbuf, rbufb) ;
              } else {
                printf(__FILE__ ": adjoint adMPI reduction not yet implemented for std op==%i\n",uop_idx) ;
              }
            }
          } else {
            exch_buf = obuf ; obuf = rbuf ; rbuf = exch_buf ;
            exch_buf = obufb ; obufb = rbufb ; rbufb = exch_buf ;
            /* Retrieve obuf and rbuf for the adjoint call */
            if (split_mode!=2) {
              adMPI_popBuffer(count,datatype,comm,obuf) ;
            }
            adMPI_popBuffer(count,datatype,comm,rbuf) ;
            (*uopb)(rbuf, rbufb, obuf, obufb, &count, &datatype, &datatypeb) ;
          }
          rc = MPI_Send(obufb, count, datatypeb, other, 11, comm);
          assert(rc==MPI_SUCCESS);
        }
        adMPI_nullifyAdjoint(count,datatypeb,comm,obufb);
      }
      if (sbuf==MPI_IN_PLACE) {
        if (rank != root) {
          adMPI_incrementAdjoint(count,datatypeb,comm,rbufb_initial,rbufb) ;
        }
      } else {
        adMPI_incrementAdjoint(count,datatypeb,comm,sbufb,rbufb) ;
        adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
      }
      if (rank!=root) free(rbufb);
      free(obufb);
    }

    free(obuf);
    if (reduceTgt) free(obufd);
    if (rank!=root) {
      free(rbuf);
      if (reduceTgt) free(rbufd);
    }
    rc = MPI_SUCCESS;

  } else {  /* i.e. op==MPI_SUM and no user-given derivative */

    if (!reduceAdj) {
      if (all_mode)
        rc=MPI_Allreduce(sbuf,
                    rbuf,
                    count,
                    datatype,
                    op,
                    comm);
      else
        rc=MPI_Reduce(sbuf,
                    rbuf,
                    count,
                    datatype,
                    op,
                    root,
                    comm);
      assert(rc==MPI_SUCCESS);
    }
    if (reduceTgt) {
      if (all_mode)
        rc=MPI_Allreduce(sbufd,
                    rbufd,
                    count,
                    datatyped,
                    op,
                    shadowcomm);
      else
        rc=MPI_Reduce(sbufd,
                    rbufd,
                    count,
                    datatyped,
                    op,
                    root,
                    shadowcomm);
      assert(rc==MPI_SUCCESS);
    }
    if (reduceAdj) {
      if (all_mode) {
        rc=MPI_Allreduce(MPI_IN_PLACE,
                         rbufb,
                         count,
                         datatypeb,
                         op,
                         shadowcomm);
        assert(rc==MPI_SUCCESS) ;
        if (sbuf!=MPI_IN_PLACE) {
          adMPI_incrementAdjoint(count,datatypeb,comm,sbufb,rbufb) ;
          adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
        }
      } else {
        int dt_idxb = derivedTypeIdx(datatypeb);
        MPI_Aint lbb = (isDerivedType(dt_idxb)?getDTypeData()->lbs[dt_idxb]:0) ;
        void *tmp_bufb =
          adMPI_allocateTempActiveBuf(count,datatypeb,comm);
        tmp_bufb = (void*)((char*)tmp_bufb - lbb);
        if (rank==root) {
          adMPI_copyActiveBuf(rbufb,tmp_bufb, count, datatypeb, comm);
        }
        rc=MPI_Bcast(tmp_bufb, count, datatypeb, root, comm) ;
        assert(rc==MPI_SUCCESS) ;
        if (sbuf==MPI_IN_PLACE) {
          if (rank != root)
            adMPI_incrementAdjoint(count,datatypeb,comm,rbufb,tmp_bufb) ;
        } else {
          if (rank == root)
            adMPI_nullifyAdjoint(count,datatypeb,comm,rbufb);
          adMPI_incrementAdjoint(count,datatypeb,comm,sbufb,tmp_bufb) ;
        }
        free(tmp_bufb);
      }
      rc = MPI_SUCCESS ;
    }
  }
  return rc;
}

// =Turn=

/** Adds into the request-to-buffer association list the associated
 * adjoint buffer <tt>adjointBuf</tt> of non-diff buffer <tt>buf</tt>
 * This should be done upon turn from FW sweep to BW sweep. */
void adMPI_Turn(void* buf, void* adjointBuf) {
  struct adMPI_Request_stack* inStack = requestStackTop ;
  while (inStack!=NULL) {
    if (inStack->buf==buf) {
      inStack->adjointBuf = adjointBuf ;
    }
    inStack = inStack->next_p ;
  }
}

// =Init=

int adMPI_Init(int* argc, 
              char*** argv) {
  int rc = MPI_Init(argc, argv);
  MPI_Comm worldDup ;
  int rc2 = MPI_Comm_dup(MPI_COMM_WORLD, &worldDup) ;
  assert(rc2==MPI_SUCCESS);
  adMPI_SHADOWCOMMLIST = NULL ;
  adMPI_addShadowComm(MPI_COMM_WORLD, worldDup) ;
  return rc ;
}

// =Finalize=

int adMPI_Finalize(void) {
  releaseDTypeData();
  releaseUOpData();
  return MPI_Finalize();
}

// =Buffer_attach=

// =Buffer_detach=

// =Type_contiguous=

int adMPI_Type_contiguous(int count,
                            MPI_Datatype oldtype,
                            MPI_Datatype* newtype) {
  int rc;
  rc = MPI_Type_contiguous(count,
                            oldtype,
                            newtype);
  assert(rc==MPI_SUCCESS);
  MPI_Datatype type, temp_packed_type, packed_type;
  MPI_Aint array_of_displacements[1] = {(MPI_Aint)0};
  int s=0, is_active, dt_idx;
  MPI_Aint p_mapsize, extent, lb;
  is_active = adMPI_isActiveType(oldtype);
  dt_idx = derivedTypeIdx(oldtype);
  derivedTypeData* dtd = getDTypeData();
  if (is_active) {
    type = MPI_DOUBLE;
    s = sizeof(double);
  }
  else if (isDerivedType(dt_idx)) {
    type = dtd->packed_types[dt_idx];
    s = dtd->p_extents[dt_idx];
  }
  else {
    type = oldtype;
    if (oldtype==MPI_DOUBLE) s = sizeof(double);
    else if (oldtype==MPI_INT) s = sizeof(int);
    else if (oldtype==MPI_FLOAT) s = sizeof(float);
    else if (oldtype==MPI_CHAR) s = sizeof(char);
    else assert(0);
  }
  p_mapsize = count*s;
  MPI_Type_get_extent(*newtype,&lb,&extent);
  rc = MPI_Type_contiguous(count,
                            type,
                            &temp_packed_type);
  assert(rc==MPI_SUCCESS);
  rc = MPI_Type_create_resized(temp_packed_type,
                                0,
                                (MPI_Aint)p_mapsize,
                                &packed_type);
  addDTypeData(dtd,
               1,
               &count,
               array_of_displacements,
               &oldtype,
               lb,
               extent,
               &count,
               array_of_displacements,
               &type,
               p_mapsize,
               newtype,
               &packed_type);
  return rc;
}

// =Type_create_struct=

int adMPI_Type_create_struct(int count,
                               int array_of_blocklengths[],
                               MPI_Aint array_of_displacements[],
                               MPI_Datatype array_of_types[],
                               MPI_Datatype *newtype) {
  int i, rc;
  rc = MPI_Type_create_struct(count,
                               array_of_blocklengths,
                               array_of_displacements,
                               array_of_types,
                               newtype);
  assert(rc==MPI_SUCCESS);
  MPI_Datatype temp_packed_type, packed_type;
  int array_of_p_blocklengths[count];
  MPI_Aint array_of_p_displacements[count];
  MPI_Datatype array_of_p_types[count], datatype;
  int s=0, is_active, is_derived, dt_idx;
  MPI_Aint p_mapsize=0, extent, lb;
  derivedTypeData* dat = getDTypeData();
  for (i=0;i<count;i++) {
    datatype = array_of_types[i];
    is_active = adMPI_isActiveType(datatype);
    dt_idx = derivedTypeIdx(datatype);
    is_derived = isDerivedType(dt_idx);
    array_of_p_blocklengths[i] = array_of_blocklengths[i];
    array_of_p_displacements[i] = p_mapsize;
    if (is_active) {
      array_of_p_types[i] = MPI_DOUBLE;
      s = sizeof(double);
    }
    else if (is_derived) {
      array_of_p_types[i] = dat->packed_types[dt_idx];
      s = dat->p_extents[dt_idx];
    }
    else {
      array_of_p_types[i] = array_of_types[i];
      if (array_of_types[i]==MPI_DOUBLE) s = sizeof(double);
      else if (array_of_types[i]==MPI_INT) s = sizeof(int);
      else if (array_of_types[i]==MPI_FLOAT) s = sizeof(float);
      else if (array_of_types[i]==MPI_CHAR) s = sizeof(char);
      else assert(0);
    }
    p_mapsize += array_of_blocklengths[i]*s;
  }
  MPI_Type_get_extent(*newtype,&lb,&extent);
  rc = MPI_Type_create_struct(count,
                               array_of_p_blocklengths,
                               array_of_p_displacements,
                               array_of_p_types,
                               &temp_packed_type);
  assert(rc==MPI_SUCCESS);
  rc = MPI_Type_create_resized(temp_packed_type,
                                0,
                                (MPI_Aint)p_mapsize,
                                &packed_type);
  addDTypeData(dat,
               count,
               array_of_blocklengths,
               array_of_displacements,
               array_of_types,
               lb,
               extent,
               array_of_p_blocklengths,
               array_of_p_displacements,
               array_of_p_types,
               p_mapsize,
               newtype,
               &packed_type);
  MPI_Type_free(&temp_packed_type);
  return rc;
}

// =Type_commit=

int adMPI_Type_commit(MPI_Datatype *datatype) {
  int dt_idx = derivedTypeIdx(*datatype);
  if (isDerivedType(dt_idx)) MPI_Type_commit(&(getDTypeData()->packed_types[dt_idx]));
  return MPI_Type_commit(datatype);
}

// =Type_create_resized=

int adMPI_Type_create_resized(MPI_Datatype oldtype,
                                MPI_Aint lb,
                                MPI_Aint extent,
                                MPI_Datatype *newtype) {
  int rc;
  rc = MPI_Type_create_resized(oldtype,
                               lb,
                               extent,
                               newtype);
  int dt_idx = derivedTypeIdx(oldtype);
  if (isDerivedType(dt_idx)) {
    derivedTypeData* dtd = getDTypeData();
    dtd->lbs[dt_idx] = lb;
    dtd->extents[dt_idx] = extent;
    dtd->derived_types[dt_idx] = *newtype;
  }
  return rc;
}

// =Type_free=

int adMPI_Type_free(MPI_Datatype *datatype) {
  int dt_idx = derivedTypeIdx(*datatype);
  if (isDerivedType(dt_idx)) MPI_Type_free(&(getDTypeData()->packed_types[dt_idx]));
  return MPI_Type_free(datatype);
}

// =Op_create=

int adMPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op) {
  int rc;
  rc = MPI_Op_create(function, commute, op);
  if (!(rc==MPI_SUCCESS)) assert(0);
  userDefinedOpData* dat = getUOpData();
  addUOpData(dat, op, function, commute);
  return rc;
}

// =Op_free=

int adMPI_Op_free(MPI_Op *op) {
  return MPI_Op_free(op);
}

// =Send=

/**
 * Tangent Send, with separate shadow (i.e. tangent) buffer.
 */
int MPI_Send_d(void* buf,  void* shadowbuf,
                   int count,
                   MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                   int dest,
                   int tag,
                   MPI_Comm comm) {
  int rc = MPI_Send(buf, count, datatype, dest, tag, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Send(shadowbuf, count, shadowdatatype, dest, tag, shadowcomm) ;
  assert(rc==MPI_SUCCESS);
  return rc ;
}

int MPI_Send_fwd(void* buf, 
                  int count, 
                  MPI_Datatype datatype, 
                  int dest, 
                  int tag,
                  MPI_Comm comm) {
  int rc=0;
  double* mappedbuf=NULL;
  int dt_idx = derivedTypeIdx(datatype);
  int is_derived = isDerivedType(dt_idx);
  if(adMPI_isActiveType(datatype)) {
      mappedbuf=buf;
  } else if(is_derived) {
      mappedbuf=adMPI_allocateTempBuf(count,datatype,comm);
  } else {
      mappedbuf=buf;
  }
  rc=MPI_Send(mappedbuf,
		count,
		adMPI_FW_rawType(datatype),
		/* if derived then need to replace typemap */
		dest,
		tag,
		comm);
  if (is_derived) free(mappedbuf);
  if (rc==MPI_SUCCESS && (adMPI_isActiveType(datatype) || is_derived)) {
      pushInteger4(dest);
      pushInteger4(tag) ;
  }
  return rc;
}

int MPI_Send_bwd(void* buf,
                  int count, 
                  MPI_Datatype datatype, 
                  int dest, 
                  int tag,
                  MPI_Comm comm) {
  int rc=0;
  popInteger4(&tag) ;
  popInteger4(&dest) ;
  MPI_Datatype mappedtype = adMPI_BW_rawType(datatype);
  void *tempBuf = adMPI_allocateTempBuf(count,mappedtype,comm) ;
  rc=MPI_Recv(tempBuf,
                  count,
                  mappedtype,
                  dest,
                  tag,
                  comm,
                  MPI_STATUS_IGNORE) ;
  adMPI_incrementAdjoint(count, mappedtype, comm, buf, tempBuf);
  free(tempBuf);
  return rc;
}

// =Isend=

int adMPI_Isend(void* buf, 
		int count, 
		MPI_Datatype datatype, 
		int dest, 
		int tag,
		MPI_Comm comm, 
		adMPI_Request* request) { 
  return MPI_Isend(buf,
		   count,
		   datatype,
		   dest,
		   tag,
		   comm,
#ifdef ADMPI_FORTRANCOMPATIBLE
		   request
#else 
		   &(request->plainRequest)
#endif 
		   );
}

/**
 * Tangent Isend, with separate shadow (i.e. tangent) buffer.
 */
int adMPI_Isend_d(void* buf, void* shadowbuf,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                    int dest,
                    int tag,
                    MPI_Comm comm,
                    adMPI_Request* request) {
  int rc = 0 ;
  MPI_Comm shadowcomm ;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  admpiRequest->plainRequest=*request;
#else 
  admpiRequest=request;
#endif
  /* adMPI_Isend of non-differentiable values must use a dummy shadowbuf */
  int sendShadow = 1 ;
  if (!adMPI_isActiveType(datatype)) {
    shadowbuf = NULL ;
    sendShadow = 0 ;
  }
  /* fill in the other info. */
  admpiRequest->endPoint=dest;
  admpiRequest->tag=tag;
  admpiRequest->count=count;
  admpiRequest->datatype=datatype;
  admpiRequest->comm=comm;
  admpiRequest->origin=WAIT_FOR_SEND;
  admpiRequest->buf = shadowbuf ;
  admpiRequest->adjointBuf=shadowbuf ;
  rc = MPI_Isend(buf, count, datatype, dest, tag, comm,
                     &(admpiRequest->plainRequest)) ;
  if (sendShadow) {
    assert(rc==MPI_SUCCESS);
    shadowcomm = adMPI_getShadowComm(comm) ;
    rc = MPI_Isend(shadowbuf, count, shadowdatatype, dest, tag, shadowcomm,
                   &(admpiRequest->shadowRequest)) ;
  }
#ifdef ADMPI_FORTRANCOMPATIBLE
  *request = admpiRequest->plainRequest ;
  adMPI_putRequest(admpiRequest);
#endif
  return rc ;
}

int adMPI_Isend_fwd(void* buf,
		   int count, 
		   MPI_Datatype datatype, 
		   int dest, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) { 
  int rc=0;
  double* mappedbuf=NULL;
  if(adMPI_isActiveType(datatype)) {
      mappedbuf=buf;
  } else {
      mappedbuf=buf;
  }
  rc= MPI_Isend(mappedbuf,
		  count,
		  datatype,
		  dest,
		  tag,
		  comm,
#ifdef ADMPI_FORTRANCOMPATIBLE
		  request
#else 
		  &(request->plainRequest)
#endif 
		  );
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  admpiRequest->plainRequest=*request;
#else 
  admpiRequest=request;
#endif
  /* fill in the other info */
  admpiRequest->endPoint=dest;
  admpiRequest->tag=tag;
  admpiRequest->count=count;
  admpiRequest->datatype=datatype;
  admpiRequest->comm=comm;
  admpiRequest->origin=WAIT_FOR_SEND;
  admpiRequest->buf = buf ;
#ifdef ADMPI_FORTRANCOMPATIBLE
  adMPI_putRequest(admpiRequest);
#endif
  return rc;
}

int adMPI_Isend_bwd(void* buf,
		   int count, 
		   MPI_Datatype datatype, 
		   int dest, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) { 
  int rc=0;
  MPI_Request *plainRequest;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  plainRequest=request;
  adMPI_getRequest(plainRequest,admpiRequest);
  plainRequest=&(admpiRequest->bwRequest);
#else 
  admpiRequest=request;
  plainRequest=&(admpiRequest->plainRequest);
#endif
  assert(admpiRequest->origin==WAIT_FOR_SEND) ;
  if (!adMPI_isActiveType(admpiRequest->datatype)) {
    /* If the type was passive, there is no BW communication: */
    rc = MPI_SUCCESS ;
  } else { 
    rc=MPI_Wait(plainRequest, MPI_STATUS_IGNORE);
    void* adjointTarget = (admpiRequest->adjointBuf?admpiRequest->adjointBuf:buf) ;
    adMPI_incrementAdjoint(admpiRequest->adjointCount,
                                 admpiRequest->datatype,
                                 admpiRequest->comm,
                                 adjointTarget,
                                 admpiRequest->adjointTempBuf);
    free(admpiRequest->adjointTempBuf);
  }
  return rc;
}

// =Bsend=

// =Rsend=

// =Recv=

/**
 * Tangent Recv, with separate shadow (i.e. tangent) buffer.
 */
int MPI_Recv_d(void* buf, void* shadowbuf,
                  int count,
                  MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                  int src,
                  int tag,
                  MPI_Comm comm,
                  MPI_Status* status) {
  MPI_Status status1 ;
  int rc = MPI_Recv(buf, count, datatype, src, tag, comm, &status1) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Recv(shadowbuf, count, shadowdatatype,
                    (src==MPI_ANY_SOURCE?status1.MPI_SOURCE:src),
                    (tag==MPI_ANY_TAG?status1.MPI_TAG:tag),
                    shadowcomm, status) ;
  assert(rc==MPI_SUCCESS);
  return rc ;
}

int MPI_Recv_fwd(void* buf, 
		 int count,
		 MPI_Datatype datatype, 
		 int src, 
		 int tag,
		 MPI_Comm comm,
		 MPI_Status* status) { 
  int rc=0;
  MPI_Status myStatus;
  double* mappedbuf=NULL;
  int dt_idx = derivedTypeIdx(datatype);
  int is_derived = isDerivedType(dt_idx);
  if(adMPI_isActiveType(datatype)) {
    mappedbuf=buf;
  } else if(is_derived) {
    mappedbuf=adMPI_allocateTempBuf(count,datatype,comm);
  } else {
    mappedbuf=buf;
  }
  rc=MPI_Recv(mappedbuf,
		count,
		adMPI_FW_rawType(datatype),
		src,
		tag,
		comm,
		&myStatus); /* because status as passed in may be MPI_STATUS_IGNORE */
  if (rc==MPI_SUCCESS && (adMPI_isActiveType(datatype) || is_derived)) {
      if (is_derived) {
	free(mappedbuf);
      }
      if(tag==MPI_ANY_TAG) tag=myStatus.MPI_TAG;
      if(src==MPI_ANY_SOURCE) src=myStatus.MPI_SOURCE;
      pushInteger4(src) ;
      pushInteger4(tag) ;
  }
  if (status!=MPI_STATUS_IGNORE) *status=myStatus;
  return rc;
}  

int MPI_Recv_bwd(void* buf, 
		 int count,
		 MPI_Datatype datatype, 
		 int src, 
		 int tag,
		 MPI_Comm comm,
		 MPI_Status* status) {
  int rc=0;
  popInteger4(&tag) ;
  popInteger4(&src) ;
  MPI_Datatype mappedtype = adMPI_BW_rawType(datatype);
  rc=MPI_Send(buf,
              count,
              mappedtype,
              src,
              tag,
              comm);
  adMPI_nullifyAdjoint(count,mappedtype,comm, buf);
  return rc;
}
                  
// =Irecv=

int adMPI_Irecv(void* buf, 
		int count, 
		MPI_Datatype datatype, 
		int src, 
		int tag,
		MPI_Comm comm, 
		adMPI_Request* request) { 
  return MPI_Irecv(buf,
		   count,
		   datatype,
		   src,
		   tag,
		   comm,
#ifdef ADMPI_FORTRANCOMPATIBLE
		   request
#else 
		   &(request->plainRequest)
#endif 
		   );
}

/**
 * Tangent Irecv, with separate shadow (i.e. tangent) buffer.
 */
int adMPI_Irecv_d(void* buf, void* shadowbuf,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                    int source,
                    int tag,
                    MPI_Comm comm,
                    adMPI_Request* request) {

  int rc=0;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  admpiRequest->plainRequest=*request;
#else 
  admpiRequest=request;
#endif
  /* adMPI_Irecv of non-differentiable values must use a dummy shadowbuf */
  if (!adMPI_isActiveType(datatype)) {
    shadowbuf = NULL ;
  }
  /* fill in the info needed to Recv the shadowbuf later.*/
  admpiRequest->endPoint=source;
  admpiRequest->tag=tag;
  admpiRequest->count=count;
  admpiRequest->datatype=shadowdatatype;
  admpiRequest->comm=comm;
  admpiRequest->origin=WAIT_FOR_RECV;
  admpiRequest->adjointBuf=shadowbuf ;
  rc= MPI_Irecv(buf,
                count,
                datatype,
                source,
                tag,
                comm,
                &(admpiRequest->plainRequest));
#ifdef ADMPI_FORTRANCOMPATIBLE
  *request = admpiRequest->plainRequest ;
  adMPI_putRequest(admpiRequest);
#endif
  return rc;
}

int adMPI_Irecv_fwd(void* buf,
		   int count,
		   MPI_Datatype datatype,
		   int source,
		   int tag,
		   MPI_Comm comm,
		   adMPI_Request* request) {
  int rc=0;
  double* mappedbuf=NULL;
  if(adMPI_isActiveType(datatype)) {
    mappedbuf=buf;
  } else {
    mappedbuf=buf;
  }
  rc= MPI_Irecv(mappedbuf,
		  count,
		  datatype,
		  source,
		  tag,
		  comm,
#ifdef ADMPI_FORTRANCOMPATIBLE
                  request
#else
                  &(request->plainRequest)
#endif
		  );
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  admpiRequest->plainRequest=*request;
#else 
  admpiRequest=request;
#endif
  /* fill in the other info */
  admpiRequest->endPoint=source;
  admpiRequest->tag=tag;
  admpiRequest->count=count;
  admpiRequest->datatype=datatype;
  admpiRequest->comm=comm;
  admpiRequest->origin=WAIT_FOR_RECV;
  admpiRequest->buf = buf ;
#ifdef ADMPI_FORTRANCOMPATIBLE
  adMPI_putRequest(admpiRequest);
#endif
  return rc;
}

int adMPI_Irecv_bwd(void* buf, 
		   int count, 
		   MPI_Datatype datatype, 
		   int source, 
		   int tag,
		   MPI_Comm comm, 
		   adMPI_Request* request) {
  int rc=0;
  MPI_Request *plainRequest;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  plainRequest=request;
  adMPI_getRequest(plainRequest,admpiRequest);
  plainRequest =&(admpiRequest->bwRequest);
#else
  plainRequest=&(request->plainRequest) ;
  admpiRequest=request;
#endif
  assert(admpiRequest->origin==WAIT_FOR_RECV) ;
  if (!adMPI_isActiveType(admpiRequest->datatype)) {
    /* If the type was passive, there is no BW communication: */
    rc = MPI_SUCCESS ;
  } else { 
    rc=MPI_Wait(plainRequest, MPI_STATUS_IGNORE);
    void* adjointTarget = (admpiRequest->adjointBuf?admpiRequest->adjointBuf:buf) ;
    adMPI_nullifyAdjoint(admpiRequest->adjointCount,
                               admpiRequest->datatype,
                               admpiRequest->comm,
                               adjointTarget);
  }
  return rc;
}

// =Wait=

int adMPI_Wait(adMPI_Request *request, MPI_Status *status) { 
  return MPI_Wait(
#ifdef ADMPI_FORTRANCOMPATIBLE
		   request
#else 
		   &(request->plainRequest)
#endif 
		   ,status);
}

/**
 * Tangent Wait, with separate shadow (i.e. tangent) buffer.
 */
int adMPI_Wait_d(adMPI_Request *request,
                  MPI_Status *status) {
  int rc=0;
  MPI_Status status1 ;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  adMPI_getRequest(request,admpiRequest);
#else 
  admpiRequest=request;
#endif 
  rc=MPI_Wait(&(admpiRequest->plainRequest), &status1);
  switch(admpiRequest->origin) { 
  case WAIT_FOR_SEND: {
    if (admpiRequest->adjointBuf) {
      assert(rc==MPI_SUCCESS);
      rc=MPI_Wait(&(admpiRequest->shadowRequest), status);
    }
    break ;
  }
  case WAIT_FOR_RECV: { 
    assert(rc==MPI_SUCCESS);
    if (admpiRequest->adjointBuf) {
      MPI_Comm shadowcomm = adMPI_getShadowComm(admpiRequest->comm) ;
      rc = MPI_Recv(admpiRequest->adjointBuf, admpiRequest->count, admpiRequest->datatype,
                    (admpiRequest->endPoint==MPI_ANY_SOURCE?status1.MPI_SOURCE:admpiRequest->endPoint),
                    (admpiRequest->tag==MPI_ANY_TAG?status1.MPI_TAG:admpiRequest->tag),
                    shadowcomm, status) ;
    }
    break ;
  }
  default:
    rc=MPI_Abort(admpiRequest->comm, MPI_ERR_ARG);
    break ;
  }
  return rc;
}

int adMPI_Wait_fwd(adMPI_Request *request,
		 MPI_Status *status) { 
  // [llh]: there used to be a version with an extra argument "buf", maybe
  //    in case admpiRequest's shadowBuf was incorrect or not set
  int rc=0;
  MPI_Request *plainRequest;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
  plainRequest=request;
  adMPI_getRequest(plainRequest,admpiRequest);
#else 
  plainRequest=&(request->plainRequest);
  admpiRequest=request;
#endif 
  rc=MPI_Wait(plainRequest,
	      status);
  if (rc==MPI_SUCCESS) {
    if (adMPI_isActiveType(admpiRequest->datatype)) {
      if(admpiRequest->tag==MPI_ANY_TAG) admpiRequest->tag=status->MPI_TAG;
      if(admpiRequest->endPoint==MPI_ANY_SOURCE) admpiRequest->endPoint=status->MPI_SOURCE;
    }
    adMPI_pushRequest(admpiRequest);
  }
  return rc;
}

int adMPI_Wait_bwd(adMPI_Request *request,
		 MPI_Status *status) {
  // [llh]: there used to be a version with an extra argument "buf", maybe
  //    in case admpiRequest.adjointBuf was incorrect or not set, e.g. by adMPI_Turn
  int rc=0;
  struct adMPI_Request_S *admpiRequest;
#ifdef ADMPI_FORTRANCOMPATIBLE
  struct adMPI_Request_S admpiRequestInst;
  admpiRequest=&admpiRequestInst;
#else 
  admpiRequest=request;
#endif 
  /* pop request  */
  adMPI_popRequest(admpiRequest);
  MPI_Request bwRequest ;
  if (adMPI_isActiveType(admpiRequest->datatype)) {
    switch(admpiRequest->origin) { 
    case WAIT_FOR_SEND: {
      admpiRequest->adjointCount=admpiRequest->count;
      admpiRequest->adjointTempBuf =
        adMPI_allocateTempBuf(admpiRequest->adjointCount,
                              admpiRequest->datatype,
                              admpiRequest->comm) ;
      rc=MPI_Irecv(admpiRequest->adjointTempBuf,
                   admpiRequest->adjointCount,
                   admpiRequest->datatype,
                   admpiRequest->endPoint,
                   admpiRequest->tag,
                   admpiRequest->comm,
                   &bwRequest/*(admpiRequest->plainRequest)*/);
      break;
    }
    case WAIT_FOR_RECV: {
      admpiRequest->adjointCount=admpiRequest->count;
      /* If this assert triggers, it means that association buf/adjointBuf was not done.
       * This can result from a missing adMPI_Turn() in Tapenade? */
      if (admpiRequest->adjointBuf==NULL) printf(__FILE__ ": missing adMPI_Turn()\n") ;
      assert(admpiRequest->adjointBuf!=NULL) ;
      rc=MPI_Isend(admpiRequest->adjointBuf,
                   admpiRequest->adjointCount,
                   admpiRequest->datatype,
                   admpiRequest->endPoint,
                   admpiRequest->tag,
                   admpiRequest->comm,
                   &bwRequest/*(admpiRequest->plainRequest)*/);
      break;
    }
    default:  
      rc=MPI_Abort(admpiRequest->comm, MPI_ERR_TYPE);
      break;
    }
#ifdef ADMPI_FORTRANCOMPATIBLE 
    admpiRequest->plainRequest=*request;
#else
    admpiRequest->plainRequest=bwRequest;
#endif
    admpiRequest->bwRequest=bwRequest;
  } else {
#ifdef ADMPI_FORTRANCOMPATIBLE 
    /* Even for a passive buffer, we must adMPI_putRequest a dummy admpiRequest with the key "request" */
    admpiRequest->plainRequest=*request;
#endif
  }
#ifdef ADMPI_FORTRANCOMPATIBLE
    adMPI_putRequest(admpiRequest);
#endif
  return rc;
}

// =Waitall=

int adMPI_Waitall(int count, 
		  adMPI_Request requests[], 
		  MPI_Status statuses[]) { 
#ifndef ADMPI_FORTRANCOMPATIBLE
  int i; 
  /* extract original requests */
  MPI_Request * origRequests=(MPI_Request*)malloc(count*sizeof(MPI_Request));
  assert(origRequests);
  for (i=0;i<count;++i) origRequests[i]=requests[i].plainRequest; 
#endif 
  return MPI_Waitall(count,
#ifdef ADMPI_FORTRANCOMPATIBLE
		     requests,
#else
		     origRequests,
#endif
		     statuses);
}

// =Barrier=

int MPI_Barrier_d(MPI_Comm comm){
  int rc=0;
  rc=MPI_Barrier(comm);
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc=MPI_Barrier(shadowcomm);
  return rc;
}

int MPI_Barrier_fwd(MPI_Comm comm){
  int rc=0;
  rc=MPI_Barrier(comm);
  return rc;
}

int MPI_Barrier_bwd(MPI_Comm comm){
  int rc;
  comm=0;
  rc=MPI_Barrier(comm);
  return rc;
}

// =Gather=

int MPI_Gather_d(void *sendbuf, void *shadowsendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                    void *recvbuf, void *shadowrecvbuf,
                    int recvcnt,
                    MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                    int root,
                    MPI_Comm comm) {
  int rc = MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Gather(shadowsendbuf, sendcnt, shadowsendtype, shadowrecvbuf, recvcnt, shadowrecvtype, root, shadowcomm) ;
  return rc;
}

int MPI_Gather_fwd(void *sendbuf,
		   int sendcnt,
		   MPI_Datatype sendtype,
		   void *recvbuf,
		   int recvcnt,
		   MPI_Datatype recvtype,
		   int root,
		   MPI_Comm comm) {
  void *rawSendBuf=sendbuf, *rawRecvBuf=recvbuf;
  int rc=MPI_SUCCESS;
  int isInPlace=(sendbuf==MPI_IN_PLACE);
  int myRank, myCommSize;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (!isInPlace && adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (!isInPlace && adMPI_isActiveType(sendtype))
      rawSendBuf=sendbuf;
    if (myRank==root) {
      if (adMPI_isActiveType(recvtype))
        rawRecvBuf=recvbuf;
    }
    rc=MPI_Gather(rawSendBuf,
		  sendcnt,
		  sendtype,
		  rawRecvBuf,
		  recvcnt,
		  recvtype,
		  root,
		  comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(recvtype)) {
      pushInteger4(myRank==root ? myCommSize : 0);
    }
  }
  return rc;
}

int MPI_Gather_bwd(void *sendbuf,
		   int sendcnt,
		   MPI_Datatype sendtype,
		   void *recvbuf,
		   int recvcnt,
		   MPI_Datatype recvtype,
		   int root,
		   MPI_Comm comm) {
  int rc=MPI_SUCCESS;
  int commSizeForRootOrNull, rTypeSize,i;
  popInteger4(&commSizeForRootOrNull) ;
  void *tempBuf = 0;
  if (sendcnt>0) tempBuf = adMPI_allocateTempBuf(sendcnt,sendtype,comm) ;
  else {
    if (commSizeForRootOrNull) 
      tempBuf=MPI_IN_PLACE;
    else 
      tempBuf=0;
  }
  rc=MPI_Scatter(recvbuf,
		 recvcnt,
		 recvtype,
		 tempBuf,
		 sendcnt,
		 sendtype,
		 root,
		 comm);
  adMPI_incrementAdjoint(sendcnt, sendtype, comm, sendbuf, tempBuf);
  if (commSizeForRootOrNull) {
    MPI_Type_size(recvtype,&rTypeSize);
    for (i=0;i<commSizeForRootOrNull;++i) { 
      if (! (i==root && sendcnt==0)) { /* don't nullify the segment if "in place" on root */
	void *recvbufSegment=(char*)recvbuf+(i*recvcnt*rTypeSize);
	adMPI_nullifyAdjoint(recvcnt,recvtype,comm,
							 recvbufSegment);
      }
    }
  }
  if (tempBuf!=MPI_IN_PLACE && tempBuf!=0) free(tempBuf);
  return rc;
}

// =Gatherv=

int MPI_Gatherv_d(void *sendbuf, void *shadowsendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                     void *recvbuf, void *shadowrecvbuf,
                     int *recvcnts,
                     int *displs,
                     MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                     int root,
                     MPI_Comm comm) {
  int rc = MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Gatherv(shadowsendbuf, sendcnt, shadowsendtype, shadowrecvbuf, recvcnts, displs, shadowrecvtype, root, shadowcomm) ;
  return rc;
}

int MPI_Gatherv_fwd(void *sendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype,
                    void *recvbuf,
                    int *recvcnts,
                    int *displs,
                    MPI_Datatype recvtype,
                    int root,
                    MPI_Comm comm) {
  void *rawSendBuf=sendbuf, *rawRecvBuf=recvbuf;
  int rc=MPI_SUCCESS;
  int isInPlace=(sendbuf==MPI_IN_PLACE);
  int myRank, myCommSize;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (!isInPlace && adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (!isInPlace && adMPI_isActiveType(sendtype))
      rawSendBuf=sendbuf;
    if (myRank==root) {
      if (adMPI_isActiveType(recvtype))
        rawRecvBuf=recvbuf;
    }
    rc=MPI_Gatherv(rawSendBuf,
                   sendcnt,
                   sendtype,
                   rawRecvBuf,
                   recvcnts,
                   displs,
                   recvtype,
                   root,
                   comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(recvtype)) {
      pushInteger4(myRank==root ? myCommSize : 0) ;
    }
  }
  return rc;
}

int MPI_Gatherv_bwd(void *sendbuf,
                    int sendcnt,
                    MPI_Datatype sendtype,
                    void *recvbuf,
                    int *recvcnts,
                    int *displs,
                    MPI_Datatype recvtype,
                    int root,
                    MPI_Comm comm) {
  int i;
  int rc=MPI_SUCCESS;
  int myRank, commSizeForRootOrNull, rTypeSize;
  popInteger4(&commSizeForRootOrNull) ;
  MPI_Comm_rank(comm, &myRank);
  void *tempBuf = 0;
  if (sendcnt>0) tempBuf = adMPI_allocateTempBuf(sendcnt,sendtype,comm) ;
  else {
    if (commSizeForRootOrNull) 
      tempBuf=MPI_IN_PLACE;
    else 
      tempBuf=0;
  }
  rc=MPI_Scatterv(recvbuf,
                  recvcnts,
                  displs,
                  recvtype,
                  tempBuf,
                  sendcnt,
                  sendtype,
                  root,
                  comm);
  adMPI_incrementAdjoint(sendcnt, sendtype, comm, sendbuf, tempBuf);
  if (commSizeForRootOrNull) {
    MPI_Type_size(recvtype,&rTypeSize);
    for (i=0;i<commSizeForRootOrNull;++i) {
      if (! (i==root && sendcnt==0)) { /* don't nullify the segment if "in place" on root */
	void* recvbufSegment=(char*)recvbuf+(rTypeSize*displs[i]); /* <----------  very iffy! */
	adMPI_nullifyAdjoint(recvcnts[i],recvtype,comm,
							 recvbufSegment);
      }
    }
  }
  if (tempBuf!=MPI_IN_PLACE && tempBuf!=0) free(tempBuf);
  return rc;
}

// =Allgather=

int MPI_Allgather_d(void *sendbuf, void *shadowsendbuf,
                       int sendcount,
                       MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                       void *recvbuf, void *shadowrecvbuf,
                       int recvcount,
                       MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                       MPI_Comm comm) {
  int rc = MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Allgather(shadowsendbuf, sendcount, shadowsendtype, shadowrecvbuf, recvcount, shadowrecvtype, shadowcomm) ;
  return rc;
}

int MPI_Allgather_fwd(void *sendbuf,
                      int sendcount,
                      MPI_Datatype sendtype,
                      void *recvbuf,
                      int recvcount,
                      MPI_Datatype recvtype,
                      MPI_Comm comm) {
  void *rawSendBuf=NULL, *rawRecvBuf=NULL;
  int rc=MPI_SUCCESS;
  int myRank, myCommSize;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (adMPI_isActiveType(sendtype))
      rawSendBuf=sendbuf;
    else rawSendBuf=sendbuf;
    if (adMPI_isActiveType(recvtype))
      rawRecvBuf=recvbuf;
    else rawRecvBuf=recvbuf;
    rc=MPI_Allgather(rawSendBuf,
                     sendcount,
                     sendtype,
                     rawRecvBuf,
                     recvcount,
                     recvtype,
                     comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(recvtype)) {
      pushInteger4(myCommSize);
    }
  }
  return rc;
}

int MPI_Allgather_bwd(void *sendbuf,
                      int sendcount,
                      MPI_Datatype sendtype,
                      void *recvbuf,
                      int recvcount,
                      MPI_Datatype recvtype,
                      MPI_Comm comm) {
  int rc=MPI_SUCCESS, rootPlaceholder;
  int commSizeForRootOrNull, rTypeSize, *recvcounts,i;
  popInteger4(&commSizeForRootOrNull) ;
  recvcounts=(int*)malloc(sizeof(int)*commSizeForRootOrNull);
  for (i=0;i<commSizeForRootOrNull;++i) recvcounts[i]=sendcount;
  void *tempBuf = adMPI_allocateTempBuf(sendcount,sendtype,comm);
  rc=MPI_Reduce_scatter(recvbuf,
                        tempBuf,
                        recvcounts,
                        MPI_DOUBLE,
                        MPI_SUM,
                        comm);
  adMPI_incrementAdjoint(sendcount, sendtype, comm, sendbuf, tempBuf);
  if (commSizeForRootOrNull) {
    MPI_Type_size(recvtype,&rTypeSize);
    adMPI_nullifyAdjoint(recvcount*commSizeForRootOrNull,
                               recvtype,comm,recvbuf);
  }
  free(tempBuf);
  if (recvcounts) free((void*)recvcounts);
  return rc;
}

// =Allgatherv=

int MPI_Allgatherv_d(void *sendbuf, void *shadowsendbuf,
                        int sendcnt,
                        MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                        void *recvbuf, void *shadowrecvbuf,
                        int *recvcnts,
                        int *displs,
                        MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                        MPI_Comm comm) {
  int rc = MPI_Allgatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Allgatherv(shadowsendbuf, sendcnt, shadowsendtype, shadowrecvbuf, recvcnts, displs, shadowrecvtype, shadowcomm) ;
  return rc;
}

int MPI_Allgatherv_fwd(void *sendbuf,
                       int sendcnt,
                       MPI_Datatype sendtype,
                       void *recvbuf,
                       int *recvcnts,
                       int *displs,
                       MPI_Datatype recvtype,
                       MPI_Comm comm) {
  void *rawSendBuf=NULL, *rawRecvBuf=NULL;
  int rc=MPI_SUCCESS;
  int myRank, myCommSize;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (adMPI_isActiveType(sendtype))
      rawSendBuf=sendbuf;
    else rawSendBuf=sendbuf;
    if (adMPI_isActiveType(recvtype))
      rawRecvBuf=recvbuf;
    else rawRecvBuf=recvbuf;
    rc=MPI_Allgatherv(rawSendBuf,
                      sendcnt,
                      sendtype,
                      rawRecvBuf,
                      recvcnts,
                      displs,
                      recvtype,
                      comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(recvtype)) {
      pushInteger4(myCommSize) ;
    }
  }
  return rc;
}

int MPI_Allgatherv_bwd(void *sendbuf,
                       int sendcnt,
                       MPI_Datatype sendtype,
                       void *recvbuf,
                       int *recvcnts,
                       int *displs,
                       MPI_Datatype recvtype,
                       MPI_Comm comm) {
  int i;
  int rc=MPI_SUCCESS;
  int myRank, commSizeForRootOrNull, rTypeSize,rootPlaceholder;
  popInteger4(&commSizeForRootOrNull) ;
  MPI_Comm_rank(comm, &myRank);
  void *tempBuf = adMPI_allocateTempBuf(recvcnts[myRank],sendtype,comm) ;
  rc=MPI_Reduce_scatter(recvbuf,
                        tempBuf,
                        recvcnts,
                        MPI_DOUBLE,
                        MPI_SUM,
                        comm);
  adMPI_incrementAdjoint(sendcnt, sendtype, comm, sendbuf, tempBuf);
  MPI_Type_size(recvtype,&rTypeSize);
  for (i=0;i<commSizeForRootOrNull;++i) {
    void* buf=(char*)recvbuf+(rTypeSize*displs[i]); /* <----------  very iffy! */
    adMPI_nullifyAdjoint(recvcnts[i],recvtype,comm,buf);
  }
  free(tempBuf);
  return rc;
}

// =Scatter=

int MPI_Scatter_d(void *sendbuf, void *shadowsendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                     void *recvbuf, void *shadowrecvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                     int root,
                     MPI_Comm comm){
  int rc = MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Scatter(shadowsendbuf, sendcnt, shadowsendtype, shadowrecvbuf, recvcnt, shadowrecvtype, root, shadowcomm) ;
  return rc;
}

int MPI_Scatter_fwd(void *sendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) {
  int rc=MPI_SUCCESS;
  int myRank, myCommSize;
  int isInPlace=(recvbuf==MPI_IN_PLACE);
  void *rawSendBuf=sendbuf, *rawRecvBuf=recvbuf;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (!isInPlace && adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (myRank==root) {
      if (adMPI_isActiveType(sendtype))
        rawSendBuf=sendbuf;
    }
    if (!isInPlace && adMPI_isActiveType(recvtype))
      rawRecvBuf=recvbuf;
    rc=MPI_Scatter(rawSendBuf,
                   sendcnt,
                   sendtype,
                   rawRecvBuf,
                   recvcnt,
                   recvtype,
                   root,
                   comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(sendtype)) {
      pushInteger4(myRank==root ? myCommSize : 0);
    }
  }
  return rc;
}

int MPI_Scatter_bwd(void *sendbuf,
                     int sendcnt,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) {
  int rc=MPI_SUCCESS;
  int commSizeForRootOrNull,i,rTypeSize;
  popInteger4(&commSizeForRootOrNull) ;
  void *tempBuf = NULL;
  if (commSizeForRootOrNull>0) tempBuf=adMPI_allocateTempBuf(sendcnt*commSizeForRootOrNull,sendtype,comm);
  rc=MPI_Gather(recvbuf,
		recvcnt,
		recvtype,
		tempBuf,
                sendcnt,
		sendtype,
		root,
		comm);
  adMPI_nullifyAdjoint(recvcnt,recvtype,comm, recvbuf);
  if (commSizeForRootOrNull>0) MPI_Type_size(recvtype,&rTypeSize);
  for (i=0;i<commSizeForRootOrNull;++i) {
    if (! (i==root && recvcnt==0)) { /* don't increment the segment if "in place" on root */
      void *tempBufSegment=(char*)tempBuf+i*sendcnt*rTypeSize;
      void *sendBufSegment=(char*)sendbuf+i*sendcnt*rTypeSize;
      adMPI_incrementAdjoint(sendcnt,
                             sendtype,
                             comm,
                             sendBufSegment,
                             tempBufSegment);
    }
  }
  if (commSizeForRootOrNull>0 && tempBuf)free(tempBuf);
  return rc;
}

// =Scatterv=

int MPI_Scatterv_d(void *sendbuf, void *shadowsendbuf,
                      int *sendcnts,
                      int *displs,
                      MPI_Datatype sendtype, MPI_Datatype shadowsendtype,
                      void *recvbuf, void *shadowrecvbuf,
                      int recvcnt,
                      MPI_Datatype recvtype, MPI_Datatype shadowrecvtype,
                      int root, MPI_Comm comm){
  int rc = MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Scatterv(shadowsendbuf, sendcnts, displs, shadowsendtype, shadowrecvbuf, recvcnt, shadowrecvtype, root, shadowcomm) ;
  return rc;
}

int MPI_Scatterv_fwd(void *sendbuf,
                     int *sendcnts,
                     int *displs,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) {
  int rc=MPI_SUCCESS;
  int myRank, myCommSize;
  int isInPlace=(recvbuf==MPI_IN_PLACE);
  void *rawSendBuf=sendbuf, *rawRecvBuf=recvbuf;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &myCommSize);
  if (!isInPlace && adMPI_isActiveType(sendtype)!=adMPI_isActiveType(recvtype)) {
    rc=MPI_Abort(comm, MPI_ERR_ARG);
  } else {
    if (myRank==root) {
      if (adMPI_isActiveType(sendtype))
        rawSendBuf=sendbuf;
    }
    if (!isInPlace && adMPI_isActiveType(recvtype))
      rawRecvBuf=recvbuf;
    rc=MPI_Scatterv(rawSendBuf,
                    sendcnts,
                    displs,
                    sendtype,
                    rawRecvBuf,
                    recvcnt,
                    recvtype,
                    root,
                    comm);
    if (rc==MPI_SUCCESS && adMPI_isActiveType(sendtype)) {
      pushInteger4(myRank==root ? myCommSize : 0) ;
    }
  }
  return rc;
}

int MPI_Scatterv_bwd(void *sendbuf,
                     int *sendcnts,
                     int *displs,
                     MPI_Datatype sendtype,
                     void *recvbuf,
                     int recvcnt,
                     MPI_Datatype recvtype,
                     int root,
                     MPI_Comm comm) {
  int rc=MPI_SUCCESS;
  int sendSize=0,i, typeSize;
  int myRank, commSizeForRootOrNull, *tempDispls;
  popInteger4(&commSizeForRootOrNull) ;
  MPI_Comm_rank(comm, &myRank);
  tempDispls=(int*)malloc(sizeof(int)*commSizeForRootOrNull);
  for (i=0;i<commSizeForRootOrNull;++i) {
    tempDispls[i]=sendSize;
    sendSize+=sendcnts[i];
  }
  void *tempBuf = NULL;
  if (commSizeForRootOrNull>0) tempBuf=adMPI_allocateTempBuf(sendSize,sendtype,comm);
  rc=MPI_Gatherv(recvbuf,
                 recvcnt,
                 recvtype,
                 tempBuf,
                 sendcnts,
                 tempDispls,
                 sendtype,
                 root,
                 comm);
  adMPI_nullifyAdjoint(recvcnt,recvtype,comm,recvbuf);
  if (commSizeForRootOrNull>0) {
    MPI_Type_size(sendtype,&typeSize);
    for (i=0;i<commSizeForRootOrNull;++i) {
      if (! (i==root && recvcnt==0)) { /* don't increment the segment if "in place" on root */
        void* buf=(char*)sendbuf+(typeSize*displs[i]); /* <----------  very iffy! */
        void* sourceBuf=(char*)tempBuf+(typeSize*tempDispls[i]);
        adMPI_incrementAdjoint(sendcnts[i], sendtype, comm, buf, sourceBuf);
      }
    }
    free(tempBuf);
  }
  if (tempDispls) free((void*)tempDispls);
  return rc;
}

// =Bcast=

int MPI_Bcast_d(void* buf, void* shadowbuf,
                   int count,
                   MPI_Datatype datatype, MPI_Datatype shadowdatatype,
                   int root,
                   MPI_Comm comm){
  int rc = MPI_Bcast(buf, count, datatype, root, comm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowcomm = adMPI_getShadowComm(comm) ;
  rc = MPI_Bcast(shadowbuf, count, shadowdatatype, root, shadowcomm) ;
  return rc;
}

int MPI_Bcast_fwd(void* buf,
                   int count,
                   MPI_Datatype datatype,
                   int root,
                   MPI_Comm comm) {
  int rc=0;
  double* mappedbuf=NULL;
  int dt_idx = derivedTypeIdx(datatype);
  int is_derived = isDerivedType(dt_idx);
  if(adMPI_isActiveType(datatype)) {
    mappedbuf=buf;
  }
  else if(is_derived) {
    mappedbuf=adMPI_allocateTempBuf(count,datatype,comm);
  }
  else {
    mappedbuf=buf;
  }
  rc=MPI_Bcast(mappedbuf,
               count,
               adMPI_FW_rawType(datatype),
               root,
               comm);
  if (rc==MPI_SUCCESS && (adMPI_isActiveType(datatype) || is_derived )) {
    if (is_derived) {
      free(mappedbuf);
    }
  }
  return rc;
}

int MPI_Bcast_bwd(void* buf,
                   int count,
                   MPI_Datatype datatype,
                   int root,
                   MPI_Comm comm) {
  int rc,rank;
  MPI_Comm_rank(comm,&rank);
  MPI_Datatype mappedtype = adMPI_BW_rawType(datatype);
  void *tempBuf = adMPI_allocateTempBuf(count,datatype,comm);
  rc=MPI_Reduce(buf,
                tempBuf,
                count,
                mappedtype,
                MPI_SUM,
                root,
                comm);
  adMPI_nullifyAdjoint(count, mappedtype, comm, buf);
  if (rank==root) {
    adMPI_incrementAdjoint(count, mappedtype, comm, buf, tempBuf);
  }
  free(tempBuf);
  return rc;
}

// =Reduce=

/** Tangent diff of MPI_Reduce.
 This version for shadowed (i.e. Association-by-Name) :
 */
int MPI_Reduce_d(void* sbuf, void* sbufd,
                    void* rbuf, void* rbufd,
                    int count,
                    MPI_Datatype datatype, MPI_Datatype datatyped,
                    MPI_Op op, TLM_userFunctionF* uopd,
                    int root,
                    MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, sbufd, NULL,
                         rbuf, rbufd, NULL,
                         count,
                         datatype, datatyped, datatype,
                         op, uopd, NULL,
                         0, 0,
                         root,
                         comm) ;
}

/** Adjoint diff of MPI_Reduce, forward sweep. */
int MPI_Reduce_fwd(void* sbuf,
                   void* rbuf,
                   int count,
                   MPI_Datatype datatype,
                   MPI_Op op,
                   int root,
                   MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, NULL, NULL,
                         rbuf, NULL, NULL,
                         count,
                         datatype, datatype, datatype, 
                         op, NULL, NULL,
                         1, 0,
                         root,
                         comm) ;
}

/** Adjoint diff of MPI_Reduce, backward sweep.
 [llh 16/10/2013] This version for shadowed (i.e. Association-by-Name) : */
int MPI_Reduce_bwd(void* sbuf, void* sbufb,
		   void* rbuf, void* rbufb,
		   int count,
		   MPI_Datatype datatype, MPI_Datatype datatypeb,
		   MPI_Op op, TLM_userFunctionF* uopb,
                   int root,
                   MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, NULL, sbufb,
                         rbuf, NULL, rbufb,
                         count,
                         datatype, datatype, datatypeb,
                         op, NULL, uopb,
                         1, 0,
                         root,
                         comm) ;
}

// =Allreduce=

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_d(void* sbuf, void* sbufd,
                       void* rbuf, void* rbufd,
                       int count,
                       MPI_Datatype datatype, MPI_Datatype datatyped,
                       MPI_Op op, TLM_userFunctionF* uopd,
                       MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, sbufd, NULL,
                         rbuf, rbufd, NULL,
                         count,
                         datatype, datatyped, datatype,
                         op, uopd, NULL,
                         0, 1,
                         0,
                         comm) ;
}

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_fwd(void* sbuf,
                       void* rbuf,
                       int count,
                       MPI_Datatype datatype,
                       MPI_Op op,
                       MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, NULL, NULL,
                         rbuf, NULL, NULL,
                         count,
                         datatype, datatype, datatype, 
                         op, NULL, NULL,
                         1, 1,
                         0,
                         comm) ;
}

/** Adjoint forward sweep of MPI_Allreduce */
int MPI_Allreduce_bwd(void* sbuf, void* sbufb,
                       void* rbuf, void* rbufb,
                       int count,
                       MPI_Datatype datatype, MPI_Datatype datatypeb,
                       MPI_Op op, TLM_userFunctionF* uopb,
                       MPI_Comm comm) {
  return PEDESTRIAN_ADMPI_Reduce(sbuf, NULL, sbufb,
                         rbuf, NULL, rbufb,
                         count,
                         datatype, datatype, datatypeb,
                         op, NULL, uopb,
                         1, 1,
                         0,
                         comm) ;
}

// =Comm_size=

// =Comm_rank=

// =Comm_dup=

/**
 * In addition to MPI_Comm_dup(), creates and registers a shadow comm
 */
int MPI_Comm_dup_d(MPI_Comm comm, MPI_Comm *dupComm) {
  int rc = MPI_Comm_dup(comm, dupComm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowDupComm ;
  rc = MPI_Comm_dup(*dupComm, &shadowDupComm) ;
  adMPI_addShadowComm(*dupComm, shadowDupComm) ;
  return rc;
}

// =Comm_split=

/**
 * In addition to MPI_Comm_split(), creates and registers a shadow comm
 */
int MPI_Comm_split_d(MPI_Comm comm, int color, int key, MPI_Comm *dupComm) {
  int rc = MPI_Comm_split(comm, color, key, dupComm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowDupComm ;
  rc = MPI_Comm_dup(*dupComm, &shadowDupComm) ;
  adMPI_addShadowComm(*dupComm, shadowDupComm) ;
  return rc;
}

// =Comm_create=

/**
 * In addition to MPI_Comm_create(), creates and registers a shadow comm
 */
int MPI_Comm_create_d(MPI_Comm comm, MPI_Group group, MPI_Comm *dupComm) {
  int rc = MPI_Comm_create(comm, group, dupComm) ;
  assert(rc==MPI_SUCCESS);
  MPI_Comm shadowDupComm ;
  rc = MPI_Comm_dup(*dupComm, &shadowDupComm) ;
  adMPI_addShadowComm(*dupComm, shadowDupComm) ;
  return rc;
}

// =Comm_free=

/**
 * In addition to MPI_Comm_free(), frees the duplicate shadow comm
 */
int MPI_Comm_free_d(MPI_Comm *comm) {
  if (comm) {
    MPI_Comm shadowComm = adMPI_getShadowComm(*comm) ;
    adMPI_delShadowComm(*comm) ;
    if (shadowComm!=*comm) MPI_Comm_free(&shadowComm) ;
    return MPI_Comm_free(comm) ;
  } else
    return 0;
}


// =FORTRAN INTERFACES=

//// =Turn=

void admpi_turn_(double *v, double *vb) {
  adMPI_Turn(v, vb) ;
}

//// =Init=

void admpi_init_(int* err_code) {
  *err_code = adMPI_Init(0, 0);
#ifdef ADMPI_FORTRANCOMPATIBLE
  admpi_fortransetupbindings_() ;
#endif
}

//// =Finalize=

void admpi_finalize_(int* err_code) {
  *err_code = adMPI_Finalize();
}

//// =Send=

void mpi_send_d_(void* buf, void* shadowbuf,
                    int *count,
                    MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                    int *dest, 
                    int *tag,
                    int *commF,
                    int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Send_d(buf, shadowbuf, *count, datatype, shadowdatatype,
                            *dest, *tag,  commC);
}

void mpi_send_fwd_(void* buf, 
                   int *count, 
                   MPI_Fint *datatypeF, 
                   int *dest, 
                   int *tag,
                   int *commF,
                   int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Send_fwd(buf, *count, datatype,
                           *dest, *tag,  commC);
}

void mpi_send_bwd_(void* buf,
                   int *count,
                   MPI_Fint *datatypeF,
                   int *dest, 
                   int *tag,
                   int *commF,
                   int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Send_bwd(buf, *count, datatype,
                           *dest, *tag,  commC);
}

//// =Isend=

void admpi_isend_(void* buf,
                 int *count,
                 MPI_Fint *datatypeF,
                 int *dest,
                 int *tag,
                 int *commF,
                 MPI_Fint *requestF,
                 int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Isend(buf, *count, datatype,
                         *dest, *tag,  commC,
                         &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_isend_d_(void* buf, void* shadowbuf,
                     int *count,
                     MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                     int *dest,
                     int *tag,
                     int *commF,
                     MPI_Fint *requestF,
                     int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Isend_d(buf, shadowbuf, *count, datatype, shadowdatatype,
                             *dest, *tag,  commC, &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_isend_fwd_(void* buf,
                    int *count,
                    MPI_Fint *datatypeF,
                    int *dest,
                    int *tag,
                    int *commF,
                    MPI_Fint *requestF,
                    int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Isend_fwd(buf, *count, datatype,
                            *dest, *tag,  commC, &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_isend_bwd_(void* buf,
                    int *count,
                    MPI_Fint *datatypeF,
                    int *dest,
                    int *tag,
                    int *commF,
                    MPI_Fint *requestF,
                    int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Isend_bwd(buf, *count, datatype,
                            *dest, *tag,  commC, &request);
}

//// =Bsend=

//// =Rsend=

//// =Recv=

void mpi_recv_d_(void* buf, void* shadowbuf,
                    int *count,
                    MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                    int* src,
                    int* tag,
                    int* commF,
                    int* status,
                    int* err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Recv_d(buf, shadowbuf, *count, datatype, shadowdatatype,
                            *src, *tag,  commC,
                           (MPI_Status*)status);
}

void mpi_recv_fwd_(void* buf,
                   int *count,
                   MPI_Fint *datatypeF,
                   int *src,
                   int *tag,
                   int *commF,
                   int *status,
                   int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Recv_fwd(buf, *count, datatype,
                           *src, *tag,  commC,
                           (MPI_Status*)status);
}

void mpi_recv_bwd_(void* buf,
                   int *count,
                   MPI_Fint *datatypeF,
                   int* src,
                   int* tag,
                   int* commF,
                   int* status,
                   int* err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Recv_bwd(buf, *count, datatype,
                           *src, *tag,  commC,
                           (MPI_Status*)status);
}

//// =Irecv=

void admpi_irecv_(void* buf,
                 int *count,
                 MPI_Fint *datatypeF,
                 int *source,
                 int *tag,
                 int *commF,
                 MPI_Fint *requestF,
                 int *err_code){
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Irecv(buf, *count, datatype,
                         *source, *tag,  commC,
                         &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_irecv_d_(void* buf, void* shadowbuf,
                     int *count,
                     MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                     int *source,
                     int *tag,
                     int *commF,
                     MPI_Fint *requestF,
                     int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Irecv_d(buf, shadowbuf, *count, datatype, shadowdatatype,
                             *source, *tag,  commC, &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_irecv_fwd_(void* buf,
                    int *count,
                    MPI_Fint *datatypeF,
                    int *source,
                    int *tag,
                    int *commF,
                    MPI_Fint *requestF,
                    int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Irecv_fwd(buf, *count, datatype,
                            *source, *tag,  commC, &request);
  *requestF = MPI_Request_c2f(request);
}

void admpi_irecv_bwd_(void* buf,
                    int *count,
                    MPI_Fint *datatypeF,
                    int *source,
                    int *tag,
                    int *commF,
                    MPI_Fint *requestF,
                    int *err_code) {
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Request request   = MPI_Request_f2c(*requestF);
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = adMPI_Irecv_bwd(buf, *count, datatype,
                            *source, *tag,  commC, &request);
}

//// =Wait=

void admpi_wait_( MPI_Fint *requestF, MPI_Fint *statusF, int* err_code) {
  MPI_Request request;
  request = MPI_Request_f2c( *requestF );
  if (statusF == ADMPI_STATUS_IGNORE_F ) {
    *err_code = adMPI_Wait( &request,  MPI_STATUS_IGNORE );
  } else if (statusF == ADMPI_STATUSES_IGNORE_F ) {
    *err_code = adMPI_Wait( &request,  MPI_STATUSES_IGNORE );
  } else {
    MPI_Status status;
    MPI_Status_f2c( statusF, &status );
    *err_code = adMPI_Wait( &request,  &status );
    MPI_Status_c2f( &status, statusF ) ;
  }
}

void admpi_wait_d_(MPI_Fint *requestF, MPI_Fint *statusF, int* err_code) {
  MPI_Request request = MPI_Request_f2c( *requestF );
  if (statusF == ADMPI_STATUS_IGNORE_F ) {
    *err_code = adMPI_Wait_d( &request,  MPI_STATUS_IGNORE );
  } else if (statusF == ADMPI_STATUSES_IGNORE_F ) {
    *err_code = adMPI_Wait_d( &request,  MPI_STATUSES_IGNORE );
  } else {
    MPI_Status status;
    MPI_Status_f2c( statusF, &status );
    *err_code = adMPI_Wait_d( &request,  &status );
    MPI_Status_c2f( &status, statusF ) ;
  }
}

void admpi_wait_fwd_(MPI_Fint *requestF, MPI_Fint *statusF, int* err_code) {
  MPI_Request request = MPI_Request_f2c( *requestF );
  if (statusF == ADMPI_STATUS_IGNORE_F ) {
    *err_code = adMPI_Wait_fwd( &request,  MPI_STATUS_IGNORE );
  } else if (statusF == ADMPI_STATUSES_IGNORE_F ) {
    *err_code = adMPI_Wait_fwd( &request,  MPI_STATUSES_IGNORE );
  } else {
    MPI_Status status;
    MPI_Status_f2c( statusF, &status );
    *err_code = adMPI_Wait_fwd( &request,  &status );
    MPI_Status_c2f( &status, statusF ) ;
  }
}

void admpi_wait_bwd_(MPI_Fint *requestF, MPI_Fint *statusF, int* err_code) {
  MPI_Request request = MPI_Request_f2c( *requestF );
  *err_code = adMPI_Wait_bwd( &request,  MPI_STATUS_IGNORE );
  *requestF = MPI_Request_c2f(request);
}

//// =Barrier=

void mpi_barrier_d_(int *commF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Barrier_d(commC) ;
}

void mpi_barrier_fwd_(int *commF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Barrier_fwd(commC) ;
}

void mpi_barrier_bwd_(int *commF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Barrier_bwd(commC) ;
}

//// =Gather=

void mpi_gather_d_(void* sbuf, void* shadowsbuf, int *scount, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcount, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gather_d(sbuf, shadowsbuf, *scount, stype, shadowstype, rbuf, shadowrbuf, *rcount, rtype, shadowrtype, *root, commC) ;
}

void mpi_gather_fwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gather_fwd(sbuf, *scount, stype, rbuf, *rcount, rtype, *root, commC) ;
}

void mpi_gather_bwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gather_bwd(sbuf, *scount, stype, rbuf, *rcount, rtype, *root, commC) ;
}

//// =Gatherv=

void mpi_gatherv_d_(void* sbuf, void* shadowsbuf, int *scount, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcounts, int *displs, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gatherv_d(sbuf, shadowsbuf, *scount, stype, shadowstype, rbuf, shadowrbuf, rcounts, displs, rtype, shadowrtype, *root, commC) ;
}

void mpi_gatherv_fwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int* rcounts, int *displs, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gatherv_fwd(sbuf, *scount, stype, rbuf, rcounts, displs, rtype, *root, commC) ;
}

void mpi_gatherv_bwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcounts, int *displs, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Gatherv_bwd(sbuf, *scount, stype, rbuf, rcounts, displs, rtype, *root, commC) ;
}

//// =Allgather=

void mpi_allgather_d_(void* sbuf, void* shadowsbuf, int *scount, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcount, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgather_d(sbuf, shadowsbuf, *scount, stype, shadowstype, rbuf, shadowrbuf, *rcount, rtype, shadowrtype, commC) ;
}

void mpi_allgather_fwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgather_fwd(sbuf, *scount, stype, rbuf, *rcount, rtype, commC) ;
}

void mpi_allgather_bwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgather_bwd(sbuf, *scount, stype, rbuf, *rcount, rtype, commC) ;
}

//// =Allgatherv=

void mpi_allgatherv_d_(void* sbuf, void* shadowsbuf, int *scount, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcounts, int *displs, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgatherv_d(sbuf, shadowsbuf, *scount, stype, shadowstype, rbuf, shadowrbuf, rcounts, displs, rtype, shadowrtype, commC) ;
}

void mpi_allgatherv_fwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcounts, int *displs, MPI_Fint *rtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgatherv_fwd(sbuf, *scount, stype, rbuf, rcounts, displs, rtype, commC) ;
}

void mpi_allgatherv_bwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcounts, int *displs, MPI_Fint *rtypeF, int *commF, int* err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allgatherv_bwd(sbuf, *scount, stype, rbuf, rcounts, displs, rtype, commC) ;
}

//// =Scatter=

void mpi_scatter_d_(void* sbuf, void* shadowsbuf, int *scount, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcount, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  if (shadowrbuf==ADMPI_IN_PLACE_F) shadowrbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatter_d(sbuf, shadowsbuf, *scount, stype, shadowstype, rbuf, shadowrbuf, *rcount, rtype, shadowrtype, *root, commC) ;
}

void mpi_scatter_fwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatter_fwd(sbuf, *scount, stype, rbuf, *rcount, rtype, *root, commC) ;
}

void mpi_scatter_bwd_(void* sbuf, int *scount, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatter_bwd(sbuf, *scount, stype, rbuf, *rcount, rtype, *root, commC) ;
}

//// =Scatterv=

void mpi_scatterv_d_(void* sbuf, void* shadowsbuf, int *scounts, int *displs, MPI_Fint *stypeF, MPI_Fint *shadowstypeF, void* rbuf, void* shadowrbuf, int *rcount, MPI_Fint *rtypeF, MPI_Fint *shadowrtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  if (shadowrbuf==ADMPI_IN_PLACE_F) shadowrbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype shadowstype = MPI_Type_f2c(*shadowstypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Datatype shadowrtype = MPI_Type_f2c(*shadowrtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatterv_d(sbuf, shadowsbuf, scounts, displs, stype, shadowstype, rbuf, shadowrbuf, *rcount, rtype, shadowrtype, *root, commC) ;
}

void mpi_scatterv_fwd_(void* sbuf, int *scounts, int *displs, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatterv_fwd(sbuf, scounts, displs, stype, rbuf, *rcount, rtype, *root, commC) ;
}

void mpi_scatterv_bwd_(void* sbuf, int *scounts, int *displs, MPI_Fint *stypeF, void* rbuf, int *rcount, MPI_Fint *rtypeF, int *root, int *commF, int* err_code) {
  if (rbuf==ADMPI_IN_PLACE_F) rbuf = MPI_IN_PLACE;
  MPI_Datatype stype = MPI_Type_f2c(*stypeF) ;
  MPI_Datatype rtype = MPI_Type_f2c(*rtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Scatterv_bwd(sbuf, scounts, displs, stype, rbuf, *rcount, rtype, *root, commC) ;
}

//// =Bcast=

void mpi_bcast_d_(void* buf, void* shadowbuf, int *count, MPI_Fint *typeF, MPI_Fint *shadowtypeF, int *root, int *commF, int* err_code) {
  MPI_Datatype type = MPI_Type_f2c(*typeF) ;
  MPI_Datatype shadowtype = MPI_Type_f2c(*shadowtypeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Bcast_d(buf, shadowbuf, *count, type, shadowtype, *root, commC) ;
}

void mpi_bcast_fwd_(void* buf, int *count, MPI_Fint *typeF, int *root, int *commF, int* err_code) {
  MPI_Datatype type = MPI_Type_f2c(*typeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Bcast_fwd(buf, *count, type, *root, commC) ;
}

void mpi_bcast_bwd_(void* buf, int *count, MPI_Fint *typeF, int *root, int *commF, int* err_code) {
  MPI_Datatype type = MPI_Type_f2c(*typeF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Bcast_bwd(buf, *count, type, *root, commC) ;
}

//// =Reduce=

void mpi_reduce_d_(void* sbuf, void* shadowsbuf,
                      void* rbuf, void* shadowrbuf,
                      int *count,
                      MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                      MPI_Fint *opF, void* uopdF,
                      int *root,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  TLM_userFunctionF* uopd = 0 /*???(uopdF)*/ ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Reduce_d(sbuf, shadowsbuf,
                              rbuf, shadowrbuf,
                              *count,
                              datatype, shadowdatatype,
                              op, uopd,
                              *root, commC) ;
}

void mpi_reduce_fwd_(void* sbuf, void* rbuf,
                      int *count,
                      MPI_Fint *datatypeF,
                      MPI_Fint *opF,
                      int *root,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Reduce_fwd(sbuf, rbuf, *count, datatype,
                              op, *root, commC) ;
}

void mpi_reduce_bwd_(void* sbuf, void* sbufb,
                      void* rbuf, void* rbufb,
                      int *count,
                      MPI_Fint *datatypeF, MPI_Fint *datatypebF,
                      MPI_Fint *opF, void* uopbF,
                      int *root,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (sbufb==ADMPI_IN_PLACE_F) sbufb = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype datatypeb = MPI_Type_f2c(*datatypebF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  TLM_userFunctionF* uopb = 0 /*???(uopbF)*/ ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Reduce_bwd(sbuf, sbufb,
                              rbuf, rbufb,
                              *count,
                              datatype, datatypeb,
                              op, uopb,
                              *root, commC) ;
}

//// =Allreduce=

void mpi_allreduce_d_(void* sbuf, void* shadowsbuf,
                      void* rbuf, void* shadowrbuf,
                      int *count,
                      MPI_Fint *datatypeF, MPI_Fint *shadowdatatypeF,
                      MPI_Fint *opF, void* uopdF,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (shadowsbuf==ADMPI_IN_PLACE_F) shadowsbuf = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype shadowdatatype = MPI_Type_f2c(*shadowdatatypeF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  TLM_userFunctionF* uopd = 0 /*???(uopdF)*/ ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allreduce_d(sbuf, shadowsbuf,
                              rbuf, shadowrbuf,
                              *count,
                              datatype, shadowdatatype,
                              op, uopd,
                              commC) ;
}

void mpi_allreduce_fwd_(void* sbuf, void* rbuf,
                      int *count,
                      MPI_Fint *datatypeF,
                      MPI_Fint *opF,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allreduce_fwd(sbuf, rbuf, *count, datatype,
                              op, commC) ;
}

void mpi_allreduce_bwd_(void* sbuf, void* sbufb,
                      void* rbuf, void* rbufb,
                      int *count,
                      MPI_Fint *datatypeF, MPI_Fint *datatypebF,
                      MPI_Fint *opF, void* uopbF,
                      int *commF,
                      int *err_code) {
  if (sbuf==ADMPI_IN_PLACE_F) sbuf = MPI_IN_PLACE;
  if (sbufb==ADMPI_IN_PLACE_F) sbufb = MPI_IN_PLACE;
  MPI_Datatype datatype = MPI_Type_f2c(*datatypeF) ;
  MPI_Datatype datatypeb = MPI_Type_f2c(*datatypebF) ;
  MPI_Op op = MPI_Op_f2c(*opF) ;
  TLM_userFunctionF* uopb = 0 /*???(uopbF)*/ ;
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Allreduce_bwd(sbuf, sbufb,
                              rbuf, rbufb,
                              *count,
                              datatype, datatypeb,
                              op, uopb,
                              commC) ;
}

//// =Comm_size=

//// =Comm_rank=

//// =Comm_dup=

void mpi_comm_dup_d_(MPI_Fint *commF, MPI_Fint *dupCommF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  MPI_Comm dupCommC ;
  *err_code = MPI_Comm_dup_d(commC, &dupCommC) ;
  *dupCommF = MPI_Comm_c2f(dupCommC) ;
}

//// =Comm_split=

void mpi_comm_split_d_(MPI_Fint *commF, int *color, int *key, MPI_Fint *dupCommF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  MPI_Comm dupCommC ;
  *err_code = MPI_Comm_split_d(commC, *color, *key, &dupCommC) ;
  *dupCommF = MPI_Comm_c2f(dupCommC) ;
}

//// =Comm_create=

void mpi_comm_create_d_(MPI_Fint *commF, MPI_Fint *groupF, MPI_Fint *dupCommF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  MPI_Comm dupCommC ;
  MPI_Group group = MPI_Group_f2c(*groupF) ;
  *err_code = MPI_Comm_create_d(commC, group, &dupCommC) ;
  *dupCommF = MPI_Comm_c2f(dupCommC) ;
}

//// =Comm_free=

void mpi_comm_free_d_(MPI_Fint *commF, int* err_code) {
  MPI_Comm commC = MPI_Comm_f2c( *commF ) ;
  *err_code = MPI_Comm_free_d(&commC) ;
}
