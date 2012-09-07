#ifndef _MACROS_H_
#define _MACROS_H_

#define NLEAF 16
#define NCRIT 32
#define NTHREAD 128
#define NBLOCK 512
#define WARP_SIZE 32
#define MAXLEVELS 30
#define LMEM_STACK_SIZE 2048
#define NWARP (NTHREAD / WARP_SIZE)

#if NLEAF == 8
#define NLEAF2 3
#define LEAFBIT 29
#define BODYMASK 0x1FFFFFFF
#define INVBMASK 0xE0000000
#elif NLEAF == 16
#define NLEAF2 4
#define LEAFBIT 28
#define BODYMASK 0x0FFFFFFF
#define INVBMASK 0xF0000000
#elif NLEAF == 32
#define NLEAF2 5
#define LEAFBIT 27
#define BODYMASK 0x07FFFFFF
#define INVBMASK 0xF8000000
#elif NLEAF == 64
#define NLEAF2 6
#define LEAFBIT 26
#define BODYMASK 0x03FFFFFF
#define INVBMASK 0xFC000000
#elif NLEAF == 128
#define NLEAF2 7
#define LEAFBIT 25
#define BODYMASK 0x01FFFFFF
#define INVBMASK 0xFE000000
#else
#error "Please choose correct NLEAF available in node_specs.h"
#endif

#if NCRIT == 8
#define NCRIT2 3
#define CRITBIT 29
#define CRITMASK 0x1FFFFFFF
#define INVCMASK 0xE0000000
#elif NCRIT == 16
#define NCRIT2 4
#define CRITBIT 28
#define CRITMASK 0x0FFFFFFF
#define INVCMASK 0xF0000000
#elif NCRIT == 32
#define NCRIT2 5
#define CRITBIT 27
#define CRITMASK 0x07FFFFFF
#define INVCMASK 0xF8000000
#elif NCRIT == 64
#define NCRIT2 6
#define CRITBIT 26
#define CRITMASK 0x03FFFFFF
#define INVCMASK 0xFC000000
#elif NCRIT == 128
#define NCRIT2 7
#define CRITBIT 25
#define CRITMASK 0x01FFFFFF
#define INVCMASK 0xFE000000
#else
#error "Please choose correct NCRIT available in node_specs.h"
#endif

#if NTHREAD == 8
#define NTHREAD2 3
#elif NTHREAD == 16
#define NTHREAD2 4
#elif NTHREAD == 32
#define NTHREAD2 5
#elif NTHREAD == 64
#define NTHREAD2 6
#elif NTHREAD == 96
#define NTHREAD2 7
#elif NTHREAD == 128
#define NTHREAD2 7
#elif NTHREAD == 256
#define NTHREAD2 8
#else
#error "Please choose correct NTHREAD available in node_specs.h"
#endif

#if WARP_SIZE == 16
#define WARP_SIZE2 4
#elif WARP_SIZE == 32
#define WARP_SIZE2 5
#else
#error "Please choose correct WARP_SIZE available in node_specs.h"
#endif

#if NCRIT > 2*WARP_SIZE
#error "NCRIT in include/node_specs.h must be <= WARP_SIZE"
#endif

#if NCRIT < NLEAF
#error "Fatal, NCRIT < NLEAF. Please check that NCRIT >= NLEAF"
#endif

#define ALIGN(a, b) ((a - 1) / b + 1)

#endif
