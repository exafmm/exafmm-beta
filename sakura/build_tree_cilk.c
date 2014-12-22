#include <cilk/cilk.h>
#include <math.h> 
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#define MAXBINS 8
#define MAXBINS64 64
#define NP 512
#define NPL 64
#define NP3 64
#define NPS 16
#define th 16
#define th2 64
#define DIM 3
#define NPD 32
#define LDIM 12
#define MASK64 0x3F
#define MASK32 0x07

void decode_morton_code(int *x, int *y, int *z, uint64_t mcode);

void relocate_data_radix6(uint32_t (*restrict pointIds), uint32_t (*restrict index), 
			  uint16_t (*restrict zcodes), uint16_t (*restrict codes), 
			  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = index[j];
      jj++;
      str[ii]=jj;
    }
  }
}

void relocate_data_radix6(uint32_t (*restrict pointIds), uint32_t (*restrict index), 
			  uint32_t (*restrict zcodes), uint32_t (*restrict codes), 
			  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = index[j];
      jj++;
      str[ii]=jj;
    }
  }
}

void relocate_data_radix6(uint32_t (*restrict pointIds), uint32_t (*restrict index), 
			  uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
			  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = index[j];
      jj++;
      str[ii]=jj;
    }
  }
}

void relocate_data_radix6_noindex(uint32_t (*restrict pointIds), 
				  uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
				  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = P + j;
      jj++;
      str[ii]=jj;
    }
  }
}

void relocate_data_radix6_noindex(uint32_t (*restrict pointIds),
				  uint32_t (*restrict zcodes), uint32_t (*restrict codes),
				  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = P + j;
      jj++;
      str[ii]=jj;
    }
  }
}

void relocate_data_radix6_noindex(uint32_t (*restrict pointIds),
				  uint16_t (*restrict zcodes), uint16_t (*restrict codes),
				  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = P + j;
      jj++;
      str[ii]=jj;
    }
  }
}



void parallel_cpy(uint16_t (*restrict codes), uint16_t (*restrict zcodes), 
		  uint32_t (*restrict pointIds), uint32_t (*restrict index), int N){

  codes[0:N] = zcodes[0:N];
  pointIds[0:N] = index[0:N];

}

void parallel_cpy(uint32_t (*restrict codes), uint32_t (*restrict zcodes), 
		  uint32_t (*restrict pointIds), uint32_t (*restrict index), int N){

  codes[0:N] = zcodes[0:N];
  pointIds[0:N] = index[0:N];

}

void parallel_cpy(uint64_t (*restrict codes), uint64_t (*restrict zcodes), 
		  uint32_t (*restrict pointIds), uint32_t (*restrict index), int N){

  codes[0:N] = zcodes[0:N];
  pointIds[0:N] = index[0:N];

}

void find_bin_sizes(uint16_t (*restrict morton_codes), 
		    int (*restrict BinSizes), int sft, int M){

#pragma ivdep
  for(int i=0; i<M; i++){
    uint32_t ii = (morton_codes[i] >> sft) & 0x3F;
    BinSizes[ii]++;
  }

}

void find_bin_sizes(uint32_t (*restrict morton_codes), 
		    int (*restrict BinSizes), int sft, int M){

#pragma ivdep
  for(int i=0; i<M; i++){
    uint32_t ii = (morton_codes[i] >> sft) & 0x3F;
    BinSizes[ii]++;
  }

}

void find_bin_sizes(uint64_t (*restrict morton_codes), 
		    int (*restrict BinSizes), int sft, int M){

#pragma ivdep
  for(int i=0; i<M; i++){
    uint32_t ii = (morton_codes[i] >> sft) & 0x3F;
    BinSizes[ii]++;
  }

}


void binning(uint64_t (*restrict morton_codes), 
	     int (*restrict BinSizes), int M){

#pragma ivdep
  for(int i=0; i<M; i++){
    BinSizes[morton_codes[i]]++;
  }

}

void bin_sort_serial_radix6_bitmap_old(uint64_t *zcodes, uint64_t *codes, 
				       uint32_t *pointIds, uint32_t *index, 
				       uint32_t *bit_map,
				       int N, 
				       int sft, int lv, int stop, 
				       int population_threshold){

  int BinSizes[MAXBINS64];
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;

  bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level

  if(N<=population_threshold || sft < stop){ // Base case. The node is a leaf

    pointIds[0:N] = index[0:N]; // Exchange the pernutation vector
    codes[0:N] = zcodes[0:N]; 
    return;
  }
  else{

    BinSizes[:] = 0;

    // Find which child each point belongs to 
#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      BinSizes[ii]++;
    }

    // scan prefix (must change this code)  

    int offset = 0;
#pragma ivdep
    for(int i=0; i<MAXBINS64; i++){
      int ss = BinSizes[i];
      BinSizes[i] = offset;
      offset += ss;
    }

#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      pointIds[BinSizes[ii]] = index[j];
      codes[BinSizes[ii]] = zcodes[j];
      BinSizes[ii]++;
    }

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers 
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    /* Call the function to split the lower levels */
    int cursor = 0;
    int super_box_start = 0;
    int super_box_size = 0;
    offset = 0; 
    //#pragma cilk grainsize = 4
    int i = 0;

    //////////////
    
    for(int i=0; i<MAXBINS; i++){
      
      int super_box_size = BinSizes[(i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(i+1)*MAXBINS-1];
      
      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;

      if(super_box_size > th){
	//#pragma cilk grainsize = 4
	for(int j=0; j<MAXBINS; j++){
	  int size = BinSizes[i*MAXBINS+j] - offset;

	  if( size>0 ){
	      
	    cilk_spawn bin_sort_serial_radix6_bitmap_old(&zcodes[offset], 
							 &codes[offset], 
							 &pointIds[offset], 
							 &index[offset], 
							 &bit_map[offset], 
							 size, sft-6, 
							 lv+1, stop,
							 population_threshold);
	      
	  }
	  offset += size;
	}
      }
      else{ // must terminate
	  
	cilk_spawn parallel_cpy(&codes[offset], 
				&zcodes[offset], 
				&pointIds[offset], 
				&index[offset], super_box_size);	  
	  
	offset += super_box_size;
      }

    }
    cilk_sync;    
      
  }
  
}

/* The function performing the first recursion for the radix 6 sorting mathod 
   In this version the scan prefix is parallel*/
void bin_sort_radix6_bitmap_old(uint64_t *zcodes, uint64_t *codes, 
				uint32_t *pointIds, uint32_t *index,
				uint32_t *bit_map,
				int N, int sft, int lv, int stop, 
				int population_threshold){

  int BinSizes[NP*MAXBINS64];
  int str[NP*MAXBINS64];
  uint32_t Sizes[MAXBINS64];
  uint32_t acm_sizes[MAXBINS64];
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;  

  BinSizes[:] = 0;
  str[:] = 0;
  Sizes[:] = 0;
  acm_sizes[:] = 0;

  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
    //Bsizes[0] = N; // record the bin size
  }


  if(N<population_threshold || sft<stop){

    pointIds[0:N] = index[0:N];
    //level[0] = lv-1;
    return;

  }
  else{

    int M = (int)ceil((float)N / (float)NP);

    // Find which child each point belongs to
    //#pragma cilk grainsize = 64
    cilk_for(int i=0; i<NP; i++){
#pragma ivdep
      for(int j=0; j<M; j++){
	if(i*M+j<N){
	  uint32_t ii = (zcodes[i*M + j]>>sft) & 0x3F;
	  BinSizes[i*MAXBINS64 + ii]++;
	}
      }
    }

    // Merge the results from different threads
    int dd = 0;
    for(int i=0; i<MAXBINS64; i++){
      str[i] = dd;
      acm_sizes[i] = dd;
#pragma ivdep
      for(int j=1; j<NP; j++){
	str[j*MAXBINS64+i] = str[(j-1)*MAXBINS64+i] + BinSizes[(j-1)*MAXBINS64+i];
	Sizes[i] += BinSizes[(j-1)*MAXBINS64+i];
      }
      dd = str[(NP-1)*MAXBINS64+i] + BinSizes[(NP-1)*MAXBINS64 + i];
      Sizes[i] += BinSizes[(NP-1)*MAXBINS64 + i];
    }

    // Relocate the indices
    for(int i=0; i<NP; i++){
      cilk_spawn relocate_data_radix6(pointIds, &index[i*M], 
				      &zcodes[i*M], codes, 
				      &str[i*MAXBINS64], i*M, 
				      M, N, sft);
    }
    cilk_sync;

    //swap the index pointers  

    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    // Spawn the children

    int cursor = 0;
    int super_box_start = 0;

    for(int i=0; i<MAXBINS; i++){
      int super_box_size = 0;
      for(int j=0; j<MAXBINS; j++){
	if( Sizes[i*MAXBINS + j]>0 ){
	   
	  cilk_spawn bin_sort_serial_radix6_bitmap_old(&zcodes[acm_sizes[i*MAXBINS + j]], 
						       &codes[acm_sizes[i*MAXBINS + j]], 
						       &pointIds[acm_sizes[i*MAXBINS + j]],
						       &index[acm_sizes[i*MAXBINS + j]], 
						       &bit_map[acm_sizes[i*MAXBINS + j]], 
						       Sizes[i*MAXBINS + j], sft-6, 
						       lv+1, stop,
						       population_threshold);
	   
	   


	}	
	super_box_size += Sizes[i*MAXBINS + j]; 
      }
      bit_map[cursor] |= 1;
      cursor += super_box_size;
    }

  }
  cilk_sync;
   

}


/*Function that performs the recursion after the first. 
  In this version the scan prefix is serial*/
void bin_sort_serial_radix6_bitmap(uint64_t (*restrict zcodes), 
				   uint64_t (*restrict codes), 
				   uint32_t (*restrict pointIds), 
				   uint32_t (*restrict index), 
				   uint32_t (*restrict bit_map),
				   int N, 
				   int sft, int lv, int stop, 
				   int population_threshold){

  int BinSizes[MAXBINS64] = {0};
  int BinCursor[MAXBINS64] = {0};
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;


  bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level

  if(N<population_threshold || sft < stop){ // Base case. The node is a leaf

    pointIds[0:N] = index[0:N]; // Exchange the pernutation vector
    codes[0:N] = zcodes[0:N]; 

    return;
  }
  else{

    // Find which child each point belongs to 
#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      BinSizes[ii]++;
    }

    // scan prefix (must change this code)  
    int offset = 0;
#pragma ivdep
    for(int i=0; i<MAXBINS64; i++){
      int ss = BinSizes[i];
      BinCursor[i] = offset;
      offset += ss;
      BinSizes[i] = offset;
    }

#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      pointIds[BinCursor[ii]] = index[j];
      codes[BinCursor[ii]] = zcodes[j];
      BinCursor[ii]++;
    }

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers 
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    /* Call the function to split the lower levels */
    int cursor = 0;
    int super_box_start = 0;
    int super_box_size = 0;
    offset = 0; 
    int i = 0;

    for(int i=0; i<MAXBINS; i++){
      
      int super_box_size = BinSizes[(i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(i+1)*MAXBINS-1];
      
      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;

      if(super_box_size > th){
	//#pragma cilk grainsize = 4
	for(int j=0; j<MAXBINS; j++){
	    
	  int size = BinSizes[i*MAXBINS+j] - offset;

	  if( size>0 ){
	      
	    cilk_spawn bin_sort_serial_radix6_bitmap(&zcodes[offset], 
						     &codes[offset], 
						     &pointIds[offset], 
						     &index[offset], 
						     &bit_map[offset], 
						     size, sft-6, 
						     lv+1, stop,
						     population_threshold);
	      
	  }
	  offset += size;
	}
      }
      else{ // must terminate
	  
	cilk_spawn parallel_cpy(&codes[offset], 
				&zcodes[offset], 
				&pointIds[offset], 
				&index[offset], super_box_size);	  
	  
	  
	offset += super_box_size;
      }

    }
    cilk_sync;    
      
  }
  
}

/* The function performing the first recursion for the radix 6 sorting mathod 
   In this version the scan prefix is parallel*/
void bin_sort_radix6_bitmap(uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
			    uint32_t (*restrict pointIds), uint32_t (*restrict index),
			    uint32_t (*restrict bit_map),
			    int N, int sft, int lv, int stop, 
			    int population_threshold){

  int BinSizes[NP*MAXBINS64];
  int BinCursor[NP*MAXBINS64];
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;  

  BinSizes[:] = 0;

  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
  }

  if(N<population_threshold || sft<stop){
    pointIds[0:N] = index[0:N];
    return;
  }
  else{

    int M = (int)ceil((float)N / (float)NP);
    // Find which child each point belongs to
    for(int i=0; i<NP; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      cilk_spawn find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size); 
    }
    cilk_sync;

    // Merge the results from different threads
    int dd = 0;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
#pragma ivdep
      for(int j=0; j<NP; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinCursor[j*MAXBINS64 + i] = offset;
	offset += size; 
	BinSizes[j*MAXBINS64 + i] = offset;
      }
    }

    // Relocate the indices
    for(int i=0; i<NP; i++){
      /*
	relocate_data_radix6(pointIds, &index[i*M], 
	&zcodes[i*M], codes, 
	&BinCursor[i*MAXBINS64], i*M, 
	M, N, sft);
      */
       
      cilk_spawn relocate_data_radix6_noindex(pointIds, 
					      &zcodes[i*M], codes, 
					      &BinCursor[i*MAXBINS64], i*M, 
					      M, N, sft);
       
    }
    cilk_sync;

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    // Spawn the children
    int cursor = 0;
    int super_box_start = 0;
    //int start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NP-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NP-1)*MAXBINS64 + (i+1)*MAXBINS-1];

      for(int j=0; j<MAXBINS; j++){
	int start = (i==0 && j==0) ? 0 : BinSizes[(NP-1)*MAXBINS64 + i*MAXBINS + j-1];
	int size = BinSizes[(NP-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	if( size > 0 ){
	   
	  cilk_spawn bin_sort_serial_radix6_bitmap(&zcodes[start], 
						   &codes[start], 
						   &pointIds[start],
						   &index[start], 
						   &bit_map[start], 
						   size, sft-6, 
						   lv+1, stop,
						   population_threshold);
	   
	   


	}	
	//start += size;
      }
      bit_map[cursor] |= 1;
      cursor += super_box_size;
    }
     
  }
  cilk_sync;
   
}

// Dense version for the qube only
void bin_sort_serial_radix6_dense(uint32_t (*restrict zcodes), 
				  uint32_t (*restrict codes), 
				  uint32_t (*restrict pointIds), 
				  uint32_t (*restrict index), 
				  uint32_t (*restrict pointer2data),
				  int N, int sft, int lv, int stop, 
 				  int binId,
				  int dataoffset){
  int BinSizes[MAXBINS64] = {0};
  int BinCursor[MAXBINS64] = {0};
  uint32_t *tmp_ptr;
  uint32_t *tmp_code;

  if(sft < stop){ // Base case. The node is a leaf
    pointer2data[binId] = dataoffset;
    codes[0:N] = zcodes[0:N];
    pointIds[0:N] = index[0:N];
    return;
  }
  else{

    int sftr = MAX(sft, 0);
    int mask = (sft < 0) ? 0x07 : MASK64;
#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sftr) & mask;
      BinSizes[ii]++;
    }
    int offset = 0;
#pragma ivdep
    for(int i=0; i<MAXBINS64; i++){
      int ss = BinSizes[i];
      BinCursor[i] = offset;
      offset += ss;
      BinSizes[i] = offset;
    }

#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sftr) & mask;
      pointIds[BinCursor[ii]] = index[j];
      codes[BinCursor[ii]] = zcodes[j];
      BinCursor[ii]++;
    }

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers 
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    /* Call the function to split the lower levels */
    int cursor = 0;
    int super_box_start = 0;
    int super_box_size = 0;
    offset = 0; 


    int BINS = (sft>0) ? MAXBINS64 : 8;

    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(i+1)*MAXBINS-1];
      for(int j=0; j<MAXBINS; j++){
	int size = BinCursor[i*MAXBINS+j] - offset;
	if( size>0 ){
	  bin_sort_serial_radix6_dense(&zcodes[offset], 
				       &codes[offset], 
				       &pointIds[offset], 
				       &index[offset], 
				       pointer2data, 
				       size, sft-6, 
				       lv+1, stop,
				       binId*MAXBINS64 + i*MAXBINS + j,
				       dataoffset + offset);
	      
	}
	offset += size;
      }
    }
  }
}

void bin_sort_radix6_dense(uint32_t (*restrict zcodes), uint32_t (*restrict codes), 
			   uint32_t (*restrict pointIds), uint32_t (*restrict index),
			   uint32_t (*restrict pointer2data),
			   int N, int sft, int lv, int stop, int binId, int dataoffset){

  int BinSizes[NPD*MAXBINS64];
  int BinCursor[NPD*MAXBINS64];
  uint32_t *tmp_ptr;
  uint32_t *tmp_code;  

  BinSizes[:] = 0;

  //printf("stop: %d\n", stop);

  if(sft<stop){
    pointIds[0:N] = index[0:N];
    return;
  }
  else{

    int M = (int)ceil((float)N / (float)NPD);
    // Find which child each point belongs to
    for(int i=0; i<NPD; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      cilk_spawn find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size); 

    }
    cilk_sync;

    // Merge the results from different threads
    int dd = 0;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
#pragma ivdep
      for(int j=0; j<NPD; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinCursor[j*MAXBINS64 + i] = offset;
	offset += size; 
	BinSizes[j*MAXBINS64 + i] = offset;
      }
    }
    for(int i=0; i<NPD; i++){
      cilk_spawn relocate_data_radix6_noindex(pointIds, 
					      &zcodes[i*M], codes, 
					      &BinCursor[i*MAXBINS64], i*M, 
					      M, N, sft);	

    }
    cilk_sync;
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;
    int cursor = 0;
    int super_box_start = 0;
    int start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NPD-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NPD-1)*MAXBINS64 + (i+1)*MAXBINS-1];
      for(int j=0; j<MAXBINS; j++){
	int start = (i==0 && j==0) ? 0 : BinSizes[(NPD-1)*MAXBINS64 + i*MAXBINS + j-1];
	int size = BinSizes[(NPD-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	if( size > 0 ){
	  cilk_spawn bin_sort_serial_radix6_dense(&zcodes[start], 
						  &codes[start], 
						  &pointIds[start],
						  &index[start], 
						  pointer2data, 
						  size, sft-6, 
						  lv+1, stop,
						  binId*MAXBINS64 + i*MAXBINS + j,
						  dataoffset + start);
	}	
	start += size;
      }
    }
    cilk_sync;
  }
}


void bin_sort_dense_singlepass(float *Y, float *X, uint32_t *keys, 
			       uint32_t *permutation_vector,
			       int N, int maxLevel){

  int  nbins = 1 << (3*maxLevel);
  int *BinSizes = (int*)calloc(nbins,sizeof(int));

  int Imin = __sec_reduce_min(keys[0:N]);

  for(int i=0; i<N; i++){
    BinSizes[keys[i]-Imin]++;
  }
  int offset = 0;
  for(int i=0; i<nbins; i++){
    int ss = BinSizes[i];
    BinSizes[i] = offset;
    offset += ss;
  }

  for(int i=0; i<N; i++){
    int ii = BinSizes[keys[i]-Imin];
    Y[ii*LDIM:LDIM] = X[i*LDIM:LDIM];
    permutation_vector[ii] = i;
    BinSizes[keys[i]-Imin]++;       
  }

  free(BinSizes);

}

void bin_sort_serial_radix6_bitmap_small(float (*restrict Y), float (*restrict X), 
					 uint64_t (*restrict zcodes), 
					 uint64_t (*restrict codes), 
					 uint32_t (*restrict pointIds), 
					 uint32_t (*restrict index), 
					 uint32_t (*restrict bit_map),
					 int N, 
					 int sft, int lv, int stop, 
					 int population_threshold){

  int BinSizes[MAXBINS64] = {0};
  int BinCursor[MAXBINS64] = {0};
#ifndef MSTEPS
  uint32_t *tmp_ptr;
#else
  float *X_tmp;
#endif
  uint64_t *tmp_code;

  bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
  if(N<population_threshold || sft < stop){ // Base case. The node is a leaf
    codes[0:N] = zcodes[0:N]; 
    
#ifndef MSTEPS
    for(uint32_t i=0; i<N; i++){
      Y[(uint32_t)i*LDIM:LDIM] = X[((uint32_t)index[i]*LDIM):LDIM];
    }
#else
    Y[0:N*LDIM] = X[0:N*LDIM];
#endif

    return;
  }
  else{

    // Find which child each point belongs to 
    //#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      BinSizes[ii]++;
    }

    // scan prefix (must change this code)  
    int offset = 0;
    //#pragma ivdep
    for(int i=0; i<MAXBINS64; i++){
      int ss = BinSizes[i];
      BinCursor[i] = offset;
      offset += ss;
      BinSizes[i] = offset;
    }

    //#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
#ifdef MSTEPS
      Y[(uint32_t)BinCursor[ii]*LDIM:LDIM] = X[(uint32_t)j*LDIM:LDIM]; // data relocation
#else
      pointIds[BinCursor[ii]] = index[j];
#endif
      codes[BinCursor[ii]] = zcodes[j];
      BinCursor[ii]++;
    }

    //swap the index pointers  
#ifndef MSTEPS
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;
#else
    X_tmp = Y;
    Y = X;
    X = X_tmp;
#endif

    //swap the code pointers 
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    /* Call the function to split the lower levels */
    int cursor = 0;
    int super_box_start = 0;
    int super_box_size = 0;
    offset = 0; 
    int i = 0;

    for(int i=0; i<MAXBINS; i++){
      
      int super_box_size = BinSizes[(i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(i+1)*MAXBINS-1];
      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;

      if(super_box_size > th){
	//#pragma cilk grainsize = 4
	for(int j=0; j<MAXBINS; j++){
	    
	  int size = BinSizes[i*MAXBINS+j] - offset;

	  if( size>0 ){
	      
#ifndef MSTEPS
	    bin_sort_serial_radix6_bitmap_small(&Y[(uint32_t)offset*LDIM], 
						X, &zcodes[offset], 
						&codes[offset], 
						&pointIds[offset], 
						&index[offset], 
						&bit_map[offset], 
						size, sft-6, 
						lv+1, stop,
						population_threshold);
#else
	    bin_sort_serial_radix6_bitmap_small(&Y[(uint32_t)offset*LDIM], 
						&X[(uint32_t)offset*LDIM], 
						&zcodes[offset], 
						&codes[offset], 
						&pointIds[offset], 
						&index[offset], 
						&bit_map[offset], 
						size, sft-6, 
						lv+1, stop,
						population_threshold);
#endif
	      
	  }
	  offset += size;
	}
      }
      else{ // must terminate
	  
#ifndef MSTEPS
	for(int i=0; i<super_box_size; i++){
	  Y[(uint32_t)(offset+i)*LDIM:LDIM] = X[(uint32_t)index[offset+i]*LDIM:LDIM];
	}
#else
	Y[0:(uint32_t)super_box_size*LDIM] = X[0:(uint32_t)super_box_size*LDIM];
#endif
	  
	/*
	  parallel_cpy(&codes[offset], 
	  &zcodes[offset], 
	  &pointIds[offset], 
	  &index[offset], super_box_size);	  
	*/
	  
	offset += super_box_size;
      }

    }
      
  }
  
}


void local_data_perm(float *Y, float *X, uint32_t *index, int offset, int M, int N){

  for(uint32_t i=0; i<M; i++){
    if(offset+i<N){
      Y[i*LDIM:LDIM] = X[(uint32_t)index[i]*LDIM:LDIM];
    }
			  
  }

}

void local_data_perm(float *Y, float *X, 
		     uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
		     int* str, int P, int M, int N, int sft){
   
#pragma ivdep
  for(uint32_t j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      uint32_t jj = str[ii];
      codes[jj] = zcodes[j];
      //pointIds[jj] = index[j];
      Y[jj*LDIM:LDIM] = X[j*LDIM:LDIM];
      jj++;
      str[ii]=jj;
    }
  }

}

void bin_sort_radix6_bitmap_small(float (*restrict Y), float (*restrict X), 
				  uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
				  uint32_t (*restrict pointIds), 
				  uint32_t (*restrict index),
				  uint32_t (*restrict bit_map),
				  int N, int sft, int lv, int stop, 
				  int population_threshold){
  int BinSizes[NPS*MAXBINS64];
  int BinCursor[NPS*MAXBINS64];
#ifndef MSTEPS
  uint32_t *tmp_ptr;
#else
  float *X_tmp;
#endif
  uint64_t *tmp_code;
  BinSizes[:] = 0;
  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
  }
  if(N<population_threshold || sft<stop){
    pointIds[0:N] = index[0:N];
    return;
  }
  else{

    int M = (int)ceil((float)N / (float)NPS);
    // Find which child each point belongs to
    for(int i=0; i<NPS; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      cilk_spawn find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size); 

    }
    cilk_sync;

    // Merge the results from different threads
    int dd = 0;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
      //#pragma ivdep
      for(int j=0; j<NPS; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinCursor[j*MAXBINS64 + i] = offset;
	offset += size; 
	BinSizes[j*MAXBINS64 + i] = offset;
      }
    }

    // Relocate the indices
    for(int i=0; i<NPS; i++){

#ifndef MSTEPS
      cilk_spawn relocate_data_radix6_noindex(pointIds, 
					      &zcodes[i*M], codes, 
					      &BinCursor[i*MAXBINS64], i*M, 
					      M, N, sft);
#else

      cilk_spawn local_data_perm(Y, X,  
				 &zcodes[i*M], codes, 
				 &BinCursor[i*MAXBINS64], i*M, 
				 M, N, sft);

#endif

    }
    cilk_sync;

    //swap the index pointers  
#ifndef MSTEPS
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;
#else
    X_tmp = Y;
    Y = X;
    X = X_tmp;
#endif

    //swap the code pointers
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    // Spawn the children
    int cursor = 0;
    int super_box_start = 0;
    //int start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NPS-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NPS-1)*MAXBINS64 + (i+1)*MAXBINS-1];
      //int cursor = (i>0) * BinSizes[(NP-1)*MAXBINS64 + i*MAXBINS-1];

      for(int j=0; j<MAXBINS; j++){
	int start = (i==0 && j==0) ? 0 : BinSizes[(NPS-1)*MAXBINS64 + i*MAXBINS + j-1];
	int size = BinSizes[(NPS-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	if( size > 0 ){
	   
#ifndef MSTEPS
	  cilk_spawn bin_sort_serial_radix6_bitmap_small(&Y[(uint32_t)start*LDIM], X, 
							 &zcodes[start], 
							 &codes[start], 
							 &pointIds[start],
							 &index[start], 
							 &bit_map[start], 
							 size, sft-6, 
							 lv+1, stop,
							 population_threshold);
#else
	  cilk_spawn bin_sort_serial_radix6_bitmap_small(&Y[(uint32_t)start*LDIM], 
							 &X[(uint32_t)start*LDIM], 
							 &zcodes[start], 
							 &codes[start], 
							 &pointIds[start],
							 &index[start], 
							 &bit_map[start], 
							 size, sft-6, 
							 lv+1, stop,
							 population_threshold);
#endif

	}	
	//start += size;
      }
      bit_map[cursor] |= 1;
      cursor += super_box_size;
    }
     
  }
  cilk_sync;
   
}



// Plummer distribution
void bin_sort_radix6_bitmap_plummer_recursion(uint64_t (*restrict zcodes), 
					      uint64_t (*restrict codes), 
					      uint32_t (*restrict pointIds), 
					      uint32_t (*restrict index),
					      uint32_t (*restrict bit_map),
					      int N, int sft, int lv, int stop, 
					      int population_threshold){



  int BinSizes[NPL*MAXBINS64] = {0};
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;  

  //BinSizes[:] = 0;

  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
  }

  if(N<=population_threshold || sft<stop){

    pointIds[0:N] = index[0:N];
    codes[0:N] = zcodes[0:N]; 
    return;

  }
  else{

    int M = (int)ceil((float)N / (float)NPL);

    // Find which child each point belongs to
    cilk_for(int i=0; i<NPL; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size);  
    }

    // Merge the results from different threads
    int dd = 0;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
#pragma ivdep
      for(int j=0; j<NPL; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinSizes[j*MAXBINS64 + i] = offset;
	offset += size; 
      }
    }

    // Relocate the indices
    cilk_for(int i=0; i<NPL; i++){
       
      relocate_data_radix6(pointIds, &index[i*M], 
			   &zcodes[i*M], codes, 
			   &BinSizes[i*MAXBINS64], i*M, 
			   M, N, sft);
       
      /*
	relocate_data_radix6_noindex(pointIds, 
	&zcodes[i*M], codes, 
	&BinSizes[i*MAXBINS64], i*M, 
	M, N, sft);
      */
    }
    //cilk_sync;

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    // Spawn the children
    int cursor = 0;
    int super_box_start = 0;
    int start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NPL-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NPL-1)*MAXBINS64 + (i+1)*MAXBINS-1];

      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;

      if(super_box_size > th){
	for(int j=0; j<MAXBINS; j++){
	  int size = BinSizes[(NPL-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	  if( size > 0 ){	     
	    cilk_spawn bin_sort_radix6_bitmap_plummer_recursion(&zcodes[start], 
								&codes[start], 
								&pointIds[start],
								&index[start], 
								&bit_map[start], 
								size, sft-6, 
								lv+1, stop,
								population_threshold);	   
	  }	
	  start += size;
	}
      }
      else{
	cilk_spawn parallel_cpy(&codes[start], 
				&zcodes[start], 
				&pointIds[start], 
				&index[start], super_box_size);	  
	 
	start += super_box_size;
      }
    }
     
  }
  cilk_sync;
}

// caller function
void bin_sort_radix6_bitmap_plummer(uint64_t (*restrict zcodes), 
				    uint64_t (*restrict codes), 
				    uint32_t (*restrict pointIds), 
				    uint32_t (*restrict index),
				    uint32_t (*restrict bit_map),
				    int N, int sft, int lv, int stop, 
				    int population_threshold){



  int BinSizes[NPL*MAXBINS64] = {0};
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;  

  //BinSizes[:] = 0;

  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1)); // Set the bit of the corresponding level
  }

  if(N<=population_threshold || sft<stop){

    pointIds[0:N] = index[0:N];
    codes[0:N] = zcodes[0:N]; 
    return;

  }
  else{

    int M = (int)ceil((float)N / (float)NPL);

    // Find which child each point belongs to
    cilk_for(int i=0; i<NPL; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size);  
    }

    // Merge the results from different threads
    int dd = 0;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
#pragma ivdep
      for(int j=0; j<NPL; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinSizes[j*MAXBINS64 + i] = offset;
	offset += size; 
      }
    }

    // Relocate the indices
    cilk_for(int i=0; i<NPL; i++){
      /*
	relocate_data_radix6(pointIds, &index[i*M], 
	&zcodes[i*M], codes, 
	&BinSizes[i*MAXBINS64], i*M, 
	M, N, sft);
      */
      relocate_data_radix6_noindex(pointIds, 
				   &zcodes[i*M], codes, 
				   &BinSizes[i*MAXBINS64], i*M, 
				   M, N, sft);
    }
    //cilk_sync;

    //swap the index pointers  
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    //swap the code pointers
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    // Spawn the children
    int cursor = 0;
    int super_box_start = 0;
    int start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NPL-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NPL-1)*MAXBINS64 + (i+1)*MAXBINS-1];

      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;

      if(super_box_size > th){
	for(int j=0; j<MAXBINS; j++){
	  int size = BinSizes[(NPL-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	  if( size > 0 ){	     
	    cilk_spawn bin_sort_radix6_bitmap_plummer_recursion(&zcodes[start], 
								&codes[start], 
								&pointIds[start],
								&index[start], 
								&bit_map[start], 
								size, sft-6, 
								lv+1, stop,
								population_threshold);	   
	  }	
	  start += size;
	}
      }
      else{
	cilk_spawn parallel_cpy(&codes[start], 
				&zcodes[start], 
				&pointIds[start], 
				&index[start], super_box_size);	  
	 
	start += super_box_size;
      }
    }
     
  }
  cilk_sync;
}

void count_bins_per_level_bitmap(int (*restrict nodes_per_level), 
				 uint32_t (*restrict bit_map), uint32_t (*restrict masks), 
				 int N, int L){
  nodes_per_level[0:L] = 0;
  for(int i=0; i<N; i++){
    for(int j=0; j<L; j++){
      if((bit_map[i] & masks[j]) > 0){
	nodes_per_level[j]++;
      }
    }
  }  
}

int count_bins_bitmap_wrapper(int (*restrict nodes_per_level), int (*restrict nodes_block_first),
			      uint32_t (*restrict bit_map), int N, int L){
  uint32_t *masks = (uint32_t *)malloc(L*sizeof(uint32_t));
  cilk_for(int i=0; i<L; i++){
    masks[i] = 1 << i;
  }
  int M = (int)ceil((float)N / (float)NP3); 
  int *nodes_per_level_buff = (int*)malloc(NP3*L*sizeof(int));
  for(int i=0; i<NP3; i++){
    int size = ((i+1)*M<N) ? M : N - i*M;
    cilk_spawn count_bins_per_level_bitmap(&nodes_per_level_buff[i*L], 
					   &bit_map[i*M], masks, 
					   size, L);
  }
  cilk_sync;
  nodes_per_level[0:L] = 0;
  nodes_block_first[0:L] = 0;
  for(int i=1; i<NP3; i++){
    nodes_block_first[i*L:L] = nodes_block_first[(i-1)*L:L] + 
      nodes_per_level_buff[(i-1)*L:L];
  }
  nodes_per_level[0:L] = nodes_block_first[(NP3-1)*L:L] + 
    nodes_per_level_buff[(NP3-1)*L:L];
  int height = L;
  for(int i=0; i<L; i++){
    if(nodes_per_level[i] == 0){
      height = i;
      break;
    }
  } 
  
  free(masks);
  free(nodes_per_level_buff);
  
  return height;
}

void parent_children_connection_singlepass(int (**restrict node_pointers),
					   int (**restrict num_children),
					   int (**restrict node_codes), 
					   int (*restrict child_counter),
					   int (*restrict remaining_child),
					   int (*restrict block_first), 
					   uint32_t (*restrict masks),
					   uint32_t (*restrict bit_map), 
					   uint64_t (*restrict leaf_morton_codes),
					   int N,
					   int L, int Base, int maxL, int maxlev){
  int cursor[20];
  cursor[0:L] = block_first[0:L];
  child_counter[0:L] = 0;
  remaining_child[0:L] = 0;  
  for(int i=0; i<N; i++){
    for(int j=0; j<L-1; j++){
      if( (bit_map[i] & masks[j]) > 0 ){
	node_pointers[j][cursor[j]] = Base + i;
	decode_morton_code(&node_codes[j][3*cursor[j]],
			   &node_codes[j][3*cursor[j]+1],
			   &node_codes[j][3*cursor[j]+2],
			   leaf_morton_codes[i] >> (3*(maxlev - j -1)) );
	if(cursor[j]>block_first[j]){
	  num_children[j][cursor[j]-1] = child_counter[j+1];
	}else{
	  remaining_child[j] = child_counter[j+1];
	}
	child_counter[j+1] = 0;
	cursor[j]++;
	child_counter[j]++;
      }
    }
    if( (bit_map[i] & masks[L-1]) > 0 ){       
      node_pointers[L-1][cursor[L-1]] = Base + i;
      decode_morton_code(&node_codes[L-1][3*cursor[L-1]],
			 &node_codes[L-1][3*cursor[L-1]+1],
			 &node_codes[L-1][3*cursor[L-1]+2],
			 leaf_morton_codes[i] >> (3*(maxlev - L)));
      cursor[L-1]++;
      child_counter[L-1]++;
    }
  }
  for(int j=0; j<L-1; j++){
    if(cursor[j]>block_first[j]){
      num_children[j][cursor[j]-1] = child_counter[j+1];
    }else{
      remaining_child[j] = child_counter[j+1];
    }
  }
}

void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					int (**restrict node_codes), 
					int (*restrict nodes_block_first),
					uint32_t (*restrict bit_map), 
					uint64_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev){
  int *child_counter = (int *)malloc(L*NP3*sizeof(int));
  int *remaining_child = (int *)malloc(L*NP3*sizeof(int));
  uint32_t *masks = (uint32_t *)malloc(L*sizeof(uint32_t));  
  int M = (int)ceil((float)N / (float)NP3);
  cilk_for(int i=0; i<L; i++){
    masks[i] = 1 << i;
  }
  for(int i=0; i<NP3; i++){  
    int size = ((i+1)*M < N) ? M : N - i*M;
    cilk_spawn parent_children_connection_singlepass(node_pointers,
						     num_children,
						     node_codes, 
						     &child_counter[i*L],
						     &remaining_child[i*L],
						     &nodes_block_first[i*maxL],
						     masks,
						     &bit_map[i*M], 
						     &leaf_morton_codes[i*M],
						     size,
						     L, i*M, maxL, maxlev);
  }
  cilk_sync;
  for(int i=0; i<NP3; i++){
    for(int j=0; j<L; j++){ 
      if(nodes_block_first[i*maxL+j]-1 >= 0){
	num_children[j][nodes_block_first[i*maxL+j]-1] += remaining_child[i*L + j];
      }
    }
  }
  free(child_counter);
  free(remaining_child);
  free(masks);
}

void first_child_position(int *children_first, int *num_children, int N){
  children_first[0] = num_children[0];
  for(int i=1; i<N; i++){ 
    children_first[i] = children_first[i-1] + num_children[i];
  }
}

void first_child_position_wrapper(int **children_first, 
				  int **num_children, 
				  int *nodes_per_level, 
				  int L){
  cilk_for(int i=0; i<L; i++){
    if(nodes_per_level[i]>0){
      first_child_position(children_first[i], num_children[i], nodes_per_level[i]);
    }
  }
}

int cmpfunc_tree (const void * a, const void * b)
{
  return ( *(uint64_t*)a - *(uint64_t*)b );
}

void build_tree(float *Y, float *X, uint64_t (*restrict mcodes), 
		uint64_t (*restrict scodes), 
		uint32_t (*restrict permutation_vector), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist){
    bin_sort_radix6_bitmap(mcodes, scodes, permutation_vector, 
			   index, bit_map, 
			   N, 3*(maxlev-2), 
			   0, 3*(maxlev-maxheight), population_threshold);
}

