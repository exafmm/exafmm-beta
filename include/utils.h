/* datasets.c file */
void cube(float *X, int N);
void cubeT(float *X, int N);
void cubeTL(float *X, int N);
void octant(float *X, int N);
void octantT(float *X, int N);
void octantTL(float *X, int N);
void create_dataset(float *X, int N, int dist);
void create_dataset_T(float *X, int N, int dist);
void create_dataset_TL(float *X, int N, int dist);

/* quantization.c file */
void quantize(uint *codes, float *X, float low, float step, int N);
void quantizeT(uint *codes, float *X, float* low, float* step, int N);
void quantizeTL(uint *codes, float *X, float* low, float* step, int N);
void compute_quantization_codes(uint* codes, float *X, int N, int nbins);
void compute_quantization_codes_T(uint* codes, float *X, int N, int nbins);
void compute_quantization_codes_TL(uint* codes, float *X, int N, int nbins);

/* build_tree_cilk.c file */
unsigned long int mortonEncode_LUT(unsigned int x, unsigned int y, unsigned int z);
void morton_encoding(unsigned long int* mcode, uint *codes, int N);
void morton_encoding_T(unsigned long int* mcodes, uint *codes, int N);
void morton_encoding_TL(unsigned long int* mcodes, uint *codes, int N);
void bin_sort2(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv);
void bin_sort_radix6(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv, int stop);
void bin_sort_radix6_hr(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv, int stop, uint* permut);
void bin_sort_serial_radix6(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv, int stop);


/* data_rearrange.c file*/
void relocate_par(float *Y, float *X, uint *index, int N);
void rearrange_dataT(float *Y, float *X, uint *Ids, int N);
void permute(float *Y, float *X, uint *Ids, int N);
void rearrange_dataT_nlev(float *Y, float *X, uint *Ids, int N, int nbins, uint *bsizes, uint *baccum_sizes);
void rearrange_data(float *Y, float *X, uint *Ids, int N);



