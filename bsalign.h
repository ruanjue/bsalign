#ifndef BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H
#define BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H

/**
 *
 * bsalign.h
 *
 * Jue Ruan <ruanjue@gmail.com>
 *
 * References:
 * Myers,G. 1999. "A fast bit-vector algorithm for approximate string matching based on dynamic programming." J. ACM, 46, 395¨C41.
 * Farrar, Michael. 2007. "Striped Smith-Waterman Speeds Database Searches Six Times over Other SIMD Implementations." Bioinformatics 23 (2): 156¨C61.
 * Suzuki, Hajime, and Masahiro Kasahara. 2017. "Acceleration of Nucleotide Semi-Global Alignment with Adaptive Banded Dynamic Programming" BioRxiv, September, 130633.
 * Suzuki, Hajime, and Masahiro Kasahara. 2018. "Introducing Difference Recurrence Relations for Faster Semi-Global Alignment of Long Sequences." BMC Bioinformatics 19 (Suppl 1).
 * Li, Heng. 2018. "Minimap2: Pairwise Alignment for Nucleotide Sequences." Bioinformatics 34 (18): 3094¨C3100.
 *
 * My Algorithm = global/overlap alignment + striped vectorization + difference recurrence relation + adaptive banding + active F-loop
 * 
 * I compute the score matrix in the way of one row by one row, that is y always increases 1 in next call.
 * x is resorted into striped vector based on W and B.
 * B: number of values in a SIMD word, W = bandwidth / B.
 * When B = 4 and W = 2, there have
 * normal array     [0, 1, 2, 3, 4, 5, 6, 7]
 * striped blocks   [0, 2, 4, 6; 1, 3, 5, 7]
 * running blocks [0, 1; 2, 3; 4, 5; 6, 7]
 * In implementation, I use xint to store 16 int8_t, B = 16.
 * H, E, F, Q, G are the absolute scores
 * u, e, f, q, g are the relative scores
 * S is score for bases matching
 *
 * Important formulas:
 * H(x, y) = max(H(x - 1, y - 1) + S(x, y), E(x, y), Q(x, y), F(x, y), G(x, y))
 * u(x, y) = H(x, y) - H(x - 1, y)
 * h(x, y) = H(x, y) - H(x - 1, y - 1)
 * # vertical cross two rows
 * E(x, y) = max(E(x, y - 1) + gap_e, H(x, y - 1) + gap_o + gap_e)
 * e(x, y) = E(x, y) - H(x, y)
 * Q(x, y) = max(Q(x, y - 1) + gap_e2, H(x, y - 1) + gap_o2 + gap_e2)
 * q(x, y) = Q(x, y) - H(x, y)
 * # horizontal along a row
 * F(x, y) = max(F(x - 1, y) + gap_e, H(x - 1, y) + gap_o + gap_e)
 * f(x, y) = F(x, y) - H(x - 1, y - 1)
 * G(x, y) = max(G(x - 1, y) + gap_e2, H(x - 1, y) + gap_o2 + gap_e2)
 * g(x, y) = G(x, y) - H(x - 1, y - 1)
 *
 * Main steps
 * 1, shifting
 *  shuffling the previous row into be well aligned with current row, all in striped. See banded_striped_epi8_seqalign_piece2_row_mov
 *
 * 2, first pass of current row to obtain f(x, y) and g(x, y) of the last striped block
 *  S(x, y) can be load from prepared striped query profile
 *  load e and q from previous row
 *  f(0, y) = g(0, y) = MIN_INF
 *  h(x, y) = max(S(x, y), e(x, y), q(x, y), f(x, y), g(x, y))
 *  f(x + 1, y) = max(f(x, y) + gap_e, h(x, y) + gap_o + gap_e) + u(x, y - 1) // please note u(x, y - 1) = H(x, y - 1) - H(x - 1, y - 1)
 *  g(x + 1, y) = max(g(x, y) + gap_e2, h(x, y) + gap_o2 + gap_e2) + u(x, y - 1)
 *  sum_u[] = sum(u([], y - 1)) // summing u of each running block in previous row
 *
 * 3, active F-loop
 *  F was calculated within each running block, need to check whether F can update Hs/Fs in next one or more running blocks, call this F-penetration
 *  Imagining a bigger F that continous updates all following Hs/Fs to the very ending, traditional lazy F-loop algorithm will perform badly
 *  When doing global alignment, such like cases often happen. Active F-loop first updates all Fs at the very beggining of each running block,
 *  then there must be at most W loops instead of 16 * W.
 *  f[0] = MIN_INF // f[0] = {f(0, y), f(W, y), f(2 * W, y), ..., f((B - 1) * W, y)} in striped coordinate, f(0 ... B - 1, y) in normal coordinate
 *  f[i] = max(f[i], f[i - 1] + sum_u[i] + 16 * gap_e)
 *  g[0] = MIN_INF
 *  g[i] = max(g[i], g[i - 1] + sum_u[i] + 16 * gap_e2)
 *  Now, there will be no F-penetration
 *
 * 4, second pass of row to calculate scores
 *  ... // calculate h(x, y), and update e(x, y + 1), q(x, y + 1), f(x + 1, y), and g(x + 1, y)
 *  u(x, y) = h(x, y) - (h(x - 1, y) - u(x - 1, y)) // update u in the loop of striped blocks
 *
 * 5, find max score within a row
 *  H(x, y) = H(x - 1, y) + u(x, y)
 *
 * 6, adaptive band
 * if sum(H[0 .. W - 1]) > sum(H[15 * W .. 16 * W - 1])  row_offset = row_offset + 0 // in normal coordinate
 * if sum(H[0 .. W - 1]) == sum(H[15 * W .. 16 * W - 1]) row_offset = row_offset + 1
 * if sum(H[0 .. W - 1]) < sum(H[15 * W .. 16 * W - 1])  row_offset = row_offset + 2
 * keep the bandwidth, but shift the offset of band in this row for next call, see 1) shifting
 *
 */

/*
 * To use bsalign.h in your program, please copy bsalign.h, list.h, sort.h and mem_share.h together
 */

#include "list.h"
#include "sort.h"
#ifdef __AVX2__
#include <immintrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#endif

#define SEQALIGN_MODE_GLOBAL	0
#define SEQALIGN_MODE_OVERLAP	1
#define SEQALIGN_MODE_EXTEND	2
#define SEQALIGN_MODE_EDIT		3

#define SEQALIGN_BT_M	0
#define SEQALIGN_BT_D	1
#define SEQALIGN_BT_I	2

#define SEQALIGN_BT1_DE	4
#define SEQALIGN_BT1_IE	8

#define SEQALIGN_BT2_D1	1
#define SEQALIGN_BT2_I1	2
#define SEQALIGN_BT2_D2	3
#define SEQALIGN_BT2_I2	4
#define SEQALIGN_BT2_DE1	8
#define SEQALIGN_BT2_IE1	16
#define SEQALIGN_BT2_DE2	32
#define SEQALIGN_BT2_IE2	64

#define SEQALIGN_SCORE_EPI8_MIN	(-(MAX_B1 >> 1))
#define SEQALIGN_SCORE_EPI8_MAX	(MAX_B1 >> 1)
#define SEQALIGN_SCORE_MIN	(-(MAX_B4 >> 1))
#define SEQALIGN_SCORE_MAX	(MAX_B4 >> 1)

#ifdef __AVX2__

//#pragma message("Choose AVX2 in " __FILE__ ". Just a message, ignore it")
#define WORDSIZE	32
#define WORDSHIFT	5
typedef __m256i	xint;
#define mm_load	_mm256_load_si256
#define mm_loadu	_mm256_loadu_si256
#define mm_store	_mm256_store_si256
#define mm_or	_mm256_or_si256
#define mm_xor	_mm256_xor_si256
#define mm_and	_mm256_and_si256
#define mm_andnot	_mm256_andnot_si256
#define mm_set1_epi8	_mm256_set1_epi8
#define mm_set1_epi16	_mm256_set1_epi16
#define mm_set1_epi32	_mm256_set1_epi32
#define mm_set1_epi64x	_mm256_set1_epi64x
#define mm_srli	_mm256_srli_si256
#define mm_slli	_mm256_slli_si256
#define mm_srli_epi64	_mm256_srli_epi64
#define mm_slli_epi64	_mm256_slli_epi64
#define mm_insert_epi8	_mm256_insert_epi8
#define mm_extract_epi16	_mm256_extract_epi16
#define mm_extract_epi32	_mm256_extract_epi32
#define mm_adds_epi8	_mm256_adds_epi8
#define mm_adds_epi16	_mm256_adds_epi16
#define mm_add_epi32	_mm256_add_epi32
#define mm_subs_epi8	_mm256_subs_epi8
#define mm_subs_epi16	_mm256_subs_epi16
#define mm_sub_epi32	_mm256_sub_epi32
#define mm_cmpeq_epi8	_mm256_cmpeq_epi8
#define mm_cmpgt_epi8	_mm256_cmpgt_epi8
#define mm_cmpeq_epi8	_mm256_cmpeq_epi8
#define mm_cmpgt_epi32	_mm256_cmpgt_epi32
#define mm_movemask_epi8	_mm256_movemask_epi8
#define mm_max_epi8	_mm256_max_epi8
#define mm_max_epi16	_mm256_max_epi16
#define mm_max_epi32	_mm256_max_epi32
#define mm_blendv	_mm256_blendv_epi8
#define mm_cvtepi8lo_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(a))
#define mm_cvtepi8hi_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
#define mm_cvtepi16lo_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(a))
#define mm_cvtepi16hi_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
#define mm_cvtepi8x0_epi32(a)	_mm256_cvtepi8_epi32(_mm256_castsi256_si128(a))
#define mm_cvtepi8x1_epi32(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 8)))
#define mm_cvtepi8x2_epi32(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
#define mm_cvtepi8x3_epi32(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 24)))
#define mm_packs_epi32	_mm256_packs_epi32
#define mm_packs_epi16	_mm256_packs_epi16

#else

//#pragma message("Choose SSE4.2 in " __FILE__ ". Just a message, ignore it")
#define WORDSIZE	16
#define WORDSHIFT	4
typedef __m128i	xint;
#define mm_load	_mm_load_si128
#define mm_loadu	_mm_loadu_si128
#define mm_store	_mm_store_si128
#define mm_or	_mm_or_si128
#define mm_xor	_mm_xor_si128
#define mm_and	_mm_and_si128
#define mm_andnot	_mm_andnot_si128
#define mm_set1_epi8	_mm_set1_epi8
#define mm_set1_epi16	_mm_set1_epi16
#define mm_set1_epi32	_mm_set1_epi32
#define mm_set1_epi64x	_mm_set1_epi64x
#define mm_srli	_mm_srli_si128
#define mm_slli	_mm_slli_si128
#define mm_srli_epi64	_mm_srli_epi64
#define mm_slli_epi64	_mm_slli_epi64
#define mm_insert_epi8	_mm_insert_epi8
#define mm_extract_epi16	_mm_extract_epi16
#define mm_extract_epi32	_mm_extract_epi32
#define mm_adds_epi8	_mm_adds_epi8
#define mm_adds_epi16	_mm_adds_epi16
#define mm_add_epi32	_mm_add_epi32
#define mm_subs_epi8	_mm_subs_epi8
#define mm_subs_epi16	_mm_subs_epi16
#define mm_sub_epi32	_mm_sub_epi32
#define mm_cmpeq_epi8	_mm_cmpeq_epi8
#define mm_cmpgt_epi8	_mm_cmpgt_epi8
#define mm_cmpgt_epi32	_mm_cmpgt_epi32
#define mm_movemask_epi8	_mm_movemask_epi8
#define mm_max_epi8	_mm_max_epi8
#define mm_max_epi16	_mm_max_epi16
#define mm_max_epi32	_mm_max_epi32
#define mm_blendv	_mm_blendv_epi8
#define mm_cvtepi8lo_epi16(a)	_mm_cvtepi8_epi16(a)
#define mm_cvtepi8hi_epi16(a)	_mm_cvtepi8_epi16(_mm_srli_si128(a, 8))
#define mm_cvtepi16lo_epi32(a)	_mm_cvtepi16_epi32(a)
#define mm_cvtepi16hi_epi32(a)	_mm_cvtepi16_epi32(_mm_srli_si128(a, 8))
#define mm_cvtepi8x0_epi32(a)	_mm_cvtepi8_epi32(a)
#define mm_cvtepi8x1_epi32(a)	_mm_cvtepi8_epi32(_mm_srli_si128(a, 4))
#define mm_cvtepi8x2_epi32(a)	_mm_cvtepi8_epi32(_mm_srli_si128(a, 8))
#define mm_cvtepi8x3_epi32(a)	_mm_cvtepi8_epi32(_mm_srli_si128(a, 12))
#define mm_packs_epi32	_mm_packs_epi32
#define mm_packs_epi16	_mm_packs_epi16

#endif

#define mm_sr1bit(x) mm_or(mm_slli_epi64(mm_srli(x, 8), 63), mm_srli_epi64(x, 1))
#define mm_sl1bit(x) mm_or(mm_srli_epi64(mm_slli(x, 8), 63), mm_slli_epi64(x, 1))


typedef struct {
	int score;
	int qb, qe;
	int tb, te;
	int mat, mis, ins, del, aln;
} seqalign_result_t;

/**
 * Basic function referings for DNA alignment using edit distance, there are implementations in this file
 */

//#define DEBUG_ED
//#define DEBUG_ED_FULL

#define striped_epi2_seqedit_getval(xs, W, pos) (((xs)[((((pos) % (W)) * WORDSIZE) + (((pos) / (W)) >> 3))] >> (((((pos) / ((W)))) & 0x7))) & 0x1)
#define striped_epi2_seqedit_qprof_size(qlen) (roundup_times(qlen, WORDSIZE * 8) / 2)
typedef void (*striped_epi2_seqedit_set_query_prof_func)(u1i *qseq, u4i qlen, b1i *qprof);
typedef void (*striped_epi2_seqedit_row_init_func)(b1i *us[2], u4i W);
typedef void (*striped_epi2_seqedit_row_cal_func)(b1i *us[2][2], b1i *hs, b1i *qprof, u4i W, u1i base);
typedef seqalign_result_t (*striped_epi2_seqedit_backtrace_func)(b1i *uts[2], u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, u4v *cigars);
static inline seqalign_result_t striped_epi2_seqedit_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int verbose);

/**
 * Basic function referings for global/extend/overlap DNA sequence alignment, there are implementations in this file
 */

#define banded_striped_epi8_pos2idx(bandwidth, pos) ((((pos) % (bandwidth >> WORDSHIFT)) << WORDSHIFT) + ((pos) / ((bandwidth) >> WORDSHIFT)))

static inline void banded_striped_epi8_seqalign_set_score_matrix(b1i matrix[16], b1i mat, b1i mis){ u4i i; for(i=0;i<16;i++) matrix[i] = ((i ^ (i >> 2)) & 0x3)? mis : mat; }

#define banded_striped_epi8_seqalign_qprof_size(qlen, bandwidth) (((num_max(qlen, bandwidth) + 1) * 4) * WORDSIZE)
#define banded_striped_epi8_seqalign_get_qprof_value(qprof, pos, base) (qprof)[((pos) * 4 + (base)) * WORDSIZE]

// prepare query profile
typedef void (*banded_striped_epi8_seqalign_set_query_prof_func)(u1i *qseq, u4i qlen, b1i *qprof, u4i bandwidth, b1i score_matrix[16]);

// prepare the rows[-1]
// mode = SEQALIGN_MODE_GLOBAL/SEQALIGN_MODE_EXTEND/SEQALIGN_MODE_OVERLAP
// length of ubegs = WORDSIZE + 1
typedef void (*banded_striped_epi8_seqalign_piecex_row_init_func)(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs, int *begs, int mode, u4i bandwidth, b1i max_nt_score, b1i min_nt_score, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// W = bandwidth / WORDSIZE, WORDSIZE = 16 (SSEx 128) or 32 (AVX2 512)
// mov <= W
// converting from us[1] to us[0]
// make the two row aligned
// ubegs is the absolute scores for the first striped block
// this function will update ubegs
typedef void (*banded_striped_epi8_seqalign_piecex_row_mov_func)(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs, u4i W, u4i mov, int piecewise);

// core func to update row scores, from us[0] to us[1]
// ph is the score of H(-1, y - 1)
// rh is the revised score of H(-1, y - 1)
// ubegs[i] - ubegs[i-1] is used in F-penetration
// ubegs will be updated
// @return: score of H(-1, y), note that H(-1, y) is useful to restore all scores of row in row_max
typedef int (*banded_striped_epi8_seqalign_piecex_row_cal_func)(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh);

// us[0] and us[1] is two row-offset-aligned progenitors
// cmp [0] and [1], merges into [0], records original in os[u, e, q], {us, es, qs, ubegs}
// os stores array of the origins, 0 or 1, respectively for u, e, q
// this function is used to merge multiple progenitors into current node in graph based alignment, program should call row_mov(i)+row_cal(i) before row_merge(1 .. n in pairwise)
typedef void (*banded_striped_epi8_seqalign_piecex_row_merge_func)(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *os[3], b1i osmarks[2], u4i W, int piecewise);

typedef int (*banded_striped_epi8_seqalign_piecex_row_verify_func)(int rbeg, int W, int ph, int rh, int hh, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// find max score and return the normal position in band
typedef u4i (*banded_striped_epi8_seqalign_row_max_func)(b1i *us, int *ubegs, u4i W, int *max_score);

// @return suggestion for band moving
//  0: move band toward left
//  1: keep band in diagonal
//  2: move band toward right
typedef int (*banded_striped_epi8_seqalign_band_mov_func)(b1i *us, int *ubegs, u4i W, u4i tidx, u4i qbeg, u4i qlen);

// combine multiple input bands into a max weight output band, used in graph alignment
// every band contributes its weight in the manner of max(v1, .., vn), and find a max weight region(band)
// the weights in a running block are all set to (ubegs[i] + ubegs[i+1]) / 2 to fast calculate
typedef u4i (*banded_striped_epi8_seqalign_band_comb_func)(u4i cnt, u4i *qoffs, int **ubegs, b1v *mempool, u4i W, u4i qlen);

// get the absolute score of a position
typedef int (*banded_striped_epi8_seqalign_getscore_func)(b1i *us, int *ubegs, u4i W, u4i pos);

// backtrace
// rs->qe, rs->te and rs->score MUST be set before call this function
// begs provides the band's offset of rows, it is continously suming of mov in banded_striped_epi8_seqalign_piecex_row_mov_func
//typedef u4i (*banded_striped_epi8_seqalign_piecex_trace_func)(u1i *qseq, u1i *tseq, b1i **ups, b1i **vps, int *begs, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, seqalign_result_t *rs, u4v *cigars);
typedef u4i (*banded_striped_epi8_seqalign_piecex_backtrace_func)(u1i *qseq, u1i *tseq, b1i *bs, int *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, u4v *cigars);


typedef u4i (*seqalign_cigar2alnstr_func)(u1i *qseq, u1i *tseq, seqalign_result_t *rs, u4v *cigars, char *alnstr[3], u4i length);

// implementation of overlap alignment for two sequences
// bandwidth should be times of WORDSIZE
static inline seqalign_result_t banded_striped_epi8_seqalign_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int mode, u4i bandwidth, b1i matrix[16], b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, int verbose);

static inline void striped_epi2_seqedit_set_query_prof(u1i *qseq, u4i qlen, b1i *qprof){
	b1i *qp;
	u4i xlen, x, y, pos, i, j, k, W;
	xlen = roundup_times(qlen, WORDSIZE * 8); // 4 base into 1 byte
	W = xlen / (WORDSIZE * 8);
	memset(qprof, 0, 4 * W * WORDSIZE);
	qp = qprof;
	for(i=0;i<4;i++){
		for(j=0;j<W;j++){
			for(k=0;k<WORDSIZE;k++){
				y = 8 * k * W + j;
				for(x=0;x<8;x++){
					pos = y + x * W;
					if(pos < qlen && qseq[pos] == i){
						qp[0] |= 0x1 << ((x & 0x7));
#ifdef DEBUG_ED
						if(striped_epi2_seqedit_getval(qprof + i * W * WORDSIZE, W, pos) != 0x1){
							fprintf(stderr, " -- something wrong i=%d j=%d k=%d x=%d pos=%d in %s -- %s:%d --\n", i, j, k, x, pos, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
							abort();
						}
#endif
					}
				}
				qp ++;
			}
		}
	}
}

static inline void striped_epi2_seqedit_row_init(b1i *us[2], u4i W){
	xint BIT0, BIT1;
	u4i i;
	// 01
	BIT0 = mm_set1_epi8(0x00);
	BIT1 = mm_set1_epi8(0xFF);
	for(i=0;i<W;i++){
		mm_store(((xint*)us[0]) + i, BIT0);
		mm_store(((xint*)us[1]) + i, BIT1);
	}
}

static inline void fprint_striped_epi1_word(FILE *out, xint v){
	u4i i, j, x;
	u1i xs[WORDSIZE];
	mm_store((xint*)xs, v);
	for(i=0;i<WORDSIZE;i++){
		for(j=0;j<8;j++){
			x = (xs[i] >> (((j)&0x7))) & 0x1;
			fprintf(out, "%d", x);
		}
		fprintf(out, " ");
	}
}

static u1i *QSEQ = NULL;
static u4i  QLEN = 0;
static int SCORE = 0;

static inline void print_epi2_row(u1i tbase, b1i *us[2], u4i W){
	u4i j, b1, b2;
	int score;
	score = SCORE;
	fprintf(stdout, "[%c]\t", "ACGTN"[tbase]);
	for(j=0;j<W*WORDSIZE*8;j++){
		b1 = striped_epi2_seqedit_getval(us[0], W, j);
		b2 = striped_epi2_seqedit_getval(us[1], W, j);
		if(b1 == 0 && b2 == 1) score ++;
		else if(b1 == 1 && b2 == 0) score --;
		fprintf(stdout, "%3d%c%03d:%d%d ", j, j < QLEN? "ACGTN"[QSEQ[j]] : '*', score, b1, b2);
	}
	fprintf(stdout, "\n");
}

// global DNA edit distance
static inline void striped_epi2_seqedit_row_cal(b1i *us[2][2], b1i *hs, b1i *qprof, u4i W, u1i base){
	xint s, u1, u2, u3, u4, v1, v2, h, h2, BIT0, BIT1, BMSK;
	u4i i, running;
	BIT0 = mm_set1_epi8(0x00);
	BIT1 = mm_set1_epi8(0xFF);
	BMSK = mm_srli(mm_set1_epi8(0x1), WORDSIZE - 1);
	v1 = BIT0;
	v2 = BIT1;
	for(i=0;i<W;i++){
	inline void row_cal_sub1(){
		// set s = 0/1 = 10/00
		s  = mm_load(((xint*)qprof) + base * W + i);
		// u = H(x, y - 1) - H(x - 1, y - 1) = -1/0/1 = 10/00/01
		// v = H(x - 1, y) - H(x - 1, y - 1) = 10/00/01
		u1  = mm_load(((xint*)us[0][0]) + i);
		u2  = mm_load(((xint*)us[0][1]) + i);
	}
		row_cal_sub1();
	inline void row_cal_sub2(){
		// H(x, y) = H(x - 1, y - 1) + min(s, u + 1, v + 1)
		// h = H(x, y) - H(x - 1, y - 1) = min(s, u + 1, v + 1) = 0/1
		//    s       u       v   =    h
		//  0(01)  -1(10)  -1(10) =  0(00)
		//  0(01)  -1(10)   0(00) =  0(00)
		//  0(01)  -1(10)   1(01) =  0(00)
		//  0(01)   0(00)  -1(10) =  0(00)
		//  0(01)   0(00)   0(00) =  0(00)
		//  0(01)   0(00)   1(01) =  0(00)
		//  0(01)   1(01)  -1(10) =  0(00)
		//  0(01)   1(01)   0(00) =  0(00)
		//  0(01)   1(01)   1(01) =  0(00)
		//  1(00)  -1(10)   0(00) =  0(00)
		//  1(00)  -1(10)   1(01) =  0(00)
		//  1(00)   0(00)  -1(10) =  0(00)
		//  1(00)   0(00)   0(00) =  1(01)
		//  1(00)   0(00)   1(01) =  1(01)
		//  1(00)   1(01)  -1(10) =  0(00)
		//  1(00)   1(01)   0(00) =  1(01)
		//  1(00)   1(01)   1(01) =  1(01)
		//    ^^      ^^      ^^       ^^
		//    ab      cd      ef       xy
		// b <- s
		// c <- u1
		// e <- v1
		// x = 0
		// y = ~(b | c | e)
		h  = mm_andnot(mm_or(s, mm_or(u1, v1)), BIT1);
		// manipulate bit by bit for u' = h - v
		//    h       v   =    u'
		//  0(00)   0(00) =  0(00)
		//  0(00)   1(01) = -1(10)
		//  0(00)  -1(10) =  1(01)
		//  1(01)   0(00) =  1(01)
		//  1(01)   1(01) =  0(00)
		//    ^^      ^^       ^^
		//    ab      cd       xy
		// b <- h
		// c <- v1
		// d <- v2
		// x <- u3
		// y <- u4
		// x = d & (~(b))
		// y = d ^ (b | c | d)
		// the same will be with v' = h - u
		u3 = mm_andnot(h, v2);
		u4 = mm_xor(v2, mm_or(h, mm_or(v1, v2)));
#ifdef DEBUG_ED_FULL
		fprintf(stdout, "i  = %d\n", i);
		fprintf(stdout, "s  = "); fprint_striped_epi1_word(stdout, s); fprintf(stdout, "\n");
		fprintf(stdout, "u1 = "); fprint_striped_epi1_word(stdout, u1); fprintf(stdout, "\n");
		fprintf(stdout, "u2 = "); fprint_striped_epi1_word(stdout, u2); fprintf(stdout, "\n");
		fprintf(stdout, "v1 = "); fprint_striped_epi1_word(stdout, v1); fprintf(stdout, "\n");
		fprintf(stdout, "v2 = "); fprint_striped_epi1_word(stdout, v2); fprintf(stdout, "\n");
		fprintf(stdout, "h  = "); fprint_striped_epi1_word(stdout,  h); fprintf(stdout, "\n");
		fprintf(stdout, "u3 = "); fprint_striped_epi1_word(stdout, u3); fprintf(stdout, "\n");
		fprintf(stdout, "u4 = "); fprint_striped_epi1_word(stdout, u4); fprintf(stdout, "\n");
#endif
		v1 = mm_andnot(h, u2);
		v2 = mm_xor(u2, mm_or(h, mm_or(u1, u2)));
	}
		row_cal_sub2();
	inline void row_cal_sub3(){
		mm_store(((xint*)hs) + i, h);
		mm_store(((xint*)us[1][0]) + i, u3);
		mm_store(((xint*)us[1][1]) + i, u4);
	}
		row_cal_sub3();
		//print_epi2_row(base, us[1], W);
	}
	running = 1;
	while(running){
		//fprintf(stdout, "\nv1 = "); fprint_striped_epi1_word(stdout, v1); fprintf(stdout, "\n");
		//fprintf(stdout, "v2 = "); fprint_striped_epi1_word(stdout, v2); fprintf(stdout, "\n");
		v1 = mm_sl1bit(v1);
		v2 = mm_or(mm_sl1bit(v2), BMSK);
		//fprintf(stdout, "v1 = "); fprint_striped_epi1_word(stdout, v1); fprintf(stdout, "\n");
		//fprintf(stdout, "v2 = "); fprint_striped_epi1_word(stdout, v2); fprintf(stdout, "\n");
		//fprintf(stdout, "#\n");
		for(i=0;i<W;i++){
			s   = mm_load(((xint*)qprof) + base * W + i);
			h2  = mm_load(((xint*)hs) + i);
			u1  = mm_load(((xint*)us[0][0]) + i);
			u2  = mm_load(((xint*)us[0][1]) + i);
			h   = mm_andnot(mm_or(s, mm_or(u1, v1)), BIT1);
			//fprintf(stdout, "hold = "); fprint_striped_epi1_word(stdout, h2); fprintf(stdout, "\n");
			//fprintf(stdout, "hnew = "); fprint_striped_epi1_word(stdout,  h); fprintf(stdout, "\n");
			u3 = mm_andnot(h, v2);
			u4 = mm_xor(v2, mm_or(h, mm_or(v1, v2)));
			v1 = mm_andnot(h, u2);
			v2 = mm_xor(u2, mm_or(h, mm_or(u1, u2)));
			mm_store(((xint*)hs) + i, h);
			mm_store(((xint*)us[1][0]) + i, u3);
			mm_store(((xint*)us[1][1]) + i, u4);
#ifdef DEBUG_ED_FULL
			fprintf(stdout, "i  = %d\n", i);
			fprintf(stdout, "s  = "); fprint_striped_epi1_word(stdout, s); fprintf(stdout, "\n");
			fprintf(stdout, "u1 = "); fprint_striped_epi1_word(stdout, u1); fprintf(stdout, "\n");
			fprintf(stdout, "u2 = "); fprint_striped_epi1_word(stdout, u2); fprintf(stdout, "\n");
			fprintf(stdout, "v1 = "); fprint_striped_epi1_word(stdout, v1); fprintf(stdout, "\n");
			fprintf(stdout, "v2 = "); fprint_striped_epi1_word(stdout, v2); fprintf(stdout, "\n");
			fprintf(stdout, "h  = "); fprint_striped_epi1_word(stdout,  h); fprintf(stdout, "\n");
			fprintf(stdout, "u3 = "); fprint_striped_epi1_word(stdout, u3); fprintf(stdout, "\n");
			fprintf(stdout, "u4 = "); fprint_striped_epi1_word(stdout, u4); fprintf(stdout, "\n");
#endif
			//print_epi2_row(base, us[1], W);
			//when there is nothing to update, break
			if(!mm_movemask_epi8(mm_andnot(mm_cmpeq_epi8(h2, h), BIT1))){
				running = 0;
				break;
			}
		}
	}
}

static inline seqalign_result_t striped_epi2_seqedit_backtrace(b1i *uts[2], u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, u4v *cigars){
	seqalign_result_t rs;
	u4i W, cg, op;
	int x, y, u1, u2, u3, u4;
	ZEROS(&rs);
	rs.qe = qlen;
	rs.te = tlen;
	W = roundup_times(qlen, WORDSIZE * 8) / (WORDSIZE * 8);
	x = qlen - 1;
	y = tlen - 1;
	if(cigars) clear_u4v(cigars);
	cg = op = 0;
	while(x >= 0 && y >= 0){
		if(qseq[x] == tseq[y]){
			rs.mat ++;
			op = 0;
			x --;
			y --;
		} else {
			u3 = striped_epi2_seqedit_getval(uts[0] + (y + 1) * W * WORDSIZE, W, x);
			u4 = striped_epi2_seqedit_getval(uts[1] + (y + 1) * W * WORDSIZE, W, x);
			if(u3 == 0 && u4 == 1){ // H(x - 1, y) - H(x - 1, y - 1) + 1 == 0
				rs.ins ++;
				op = 1; // I
				x --;
			} else {
				u1 = striped_epi2_seqedit_getval(uts[0] + (y + 0) * W * WORDSIZE, W, x);
				u2 = striped_epi2_seqedit_getval(uts[1] + (y + 0) * W * WORDSIZE, W, x);
				if(u1 == 1 && u2 == 0){
					rs.del ++;
					op = 2; // D
					y --;
				} else {
					rs.mis ++;
					op = 0;
					x --;
					y --;
				}
			}
		}
		/*
		if(op == 0){
			fprintf(stdout, "%c %c\n", "ACGT*"[qseq[x + 1]], "ACGT*"[tseq[y + 1]]);
		} else if(op == 1){
			fprintf(stdout, "%c %c\n", "ACGT*"[qseq[x + 1]], '-');
		} else {
			fprintf(stdout, "%c %c\n", '-', "ACGT*"[tseq[y + 1]]);
		}
		*/
		if(op == (cg & 0xf)){
			cg += 0x10;
		} else {
			if(cg && cigars) push_u4v(cigars, cg);
			cg = 0x10 | op;
		}
	}
	rs.qb = x + 1;
	rs.tb = y + 1;
	rs.aln = rs.mat + rs.mis + rs.ins + rs.del;
	if(rs.qb){
		op = 1;
		if(op == (cg & 0xf)){
			cg += 0x10 * rs.qb;
		} else {
			if(cg && cigars) push_u4v(cigars, cg);
			cg = (0x10 * rs.qb) | op;
		}
		rs.qb = 0;
	}
	if(rs.tb){
		op = 2;
		if(op == (cg & 0xf)){
			cg += 0x10 * rs.tb;
		} else {
			if(cg && cigars) push_u4v(cigars, cg);
			cg = (0x10 * rs.tb) | op;
		}
		rs.tb = 0;
	}
	if(cg && cigars) push_u4v(cigars, cg);
	if(cigars) reverse_u4v(cigars);
	return rs;
}

static inline seqalign_result_t striped_epi2_seqedit_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int verbose){
	striped_epi2_seqedit_set_query_prof_func set_qprof;
	striped_epi2_seqedit_row_init_func       row_init;
	striped_epi2_seqedit_row_cal_func        row_cal;
	striped_epi2_seqedit_backtrace_func      backtrace;
	seqalign_result_t rs;
	b1i *memp, *mempb, *qprof, *uts[2], *us[2][2], *hs;
	u8i mpsize;
	u4i i, W;
	set_qprof = striped_epi2_seqedit_set_query_prof;
	row_cal   = striped_epi2_seqedit_row_cal;
	row_init  = striped_epi2_seqedit_row_init;
	backtrace = striped_epi2_seqedit_backtrace;
	W = roundup_times(qlen, WORDSIZE * 8) / (WORDSIZE * 8);
	//fprintf(stdout, "qlen=%d\tW=%d\n", qlen, W);
	mpsize = 0;
	mpsize += striped_epi2_seqedit_qprof_size(qlen); // qprof[]
	mpsize += 2 * W * WORDSIZE * ((tlen + 1)); // uts[]
	mpsize += W * WORDSIZE; // hs[]
	if(mempool){
		clear_and_encap_b1v(mempool, mpsize);
		memp = mempool->buffer + WORDSIZE;
		mempb = NULL;
	} else {
		mempb = malloc(mpsize);
		memp = mempb + WORDSIZE;
	}
	qprof  = memp; memp += striped_epi2_seqedit_qprof_size(qlen);
	uts[0] = memp; memp += W * WORDSIZE * (tlen + 1);
	uts[1] = memp; memp += W * WORDSIZE * (tlen + 1);
	hs     = memp; memp += W * WORDSIZE;
	set_qprof(qseq, qlen, qprof);
	us[0][0] = uts[0];
	us[0][1] = uts[1];
	QSEQ = qseq;
	QLEN = qlen;
	row_init(us[0], W);
	for(i=0;i<tlen;i++){
		SCORE = i + 1;
		us[0][0] = uts[0] + (i + 0) * W * WORDSIZE;
		us[0][1] = uts[1] + (i + 0) * W * WORDSIZE;
		us[1][0] = uts[0] + (i + 1) * W * WORDSIZE;
		us[1][1] = uts[1] + (i + 1) * W * WORDSIZE;
		row_cal(us, hs, qprof, W, tseq[i]);
		if(verbose){
			int vals[2][2] = {{0, 1}, {-1, 2}};
			int j, b1, b2, u, u2, v, v2, score;
			score = i + 1;
			v2 = 1;
			fprintf(stdout, "[%04d:%c]\t", i, "ACGTN"[tseq[i]]);
			for(j=0;j<Int(qlen);j++){
				b1 = striped_epi2_seqedit_getval(us[0][0], W, j);
				b2 = striped_epi2_seqedit_getval(us[0][1], W, j);
				u = vals[b1][b2];
				v = v2;
				if(qseq[j] == tseq[i] || u == -1 || v == -1){
					u2 = 0 - v;
					v2 =  0 - u;
				} else {
					u2 = 1 - v;
					v2  = 1 - u;
				}
				b1 = striped_epi2_seqedit_getval(us[1][0], W, j);
				b2 = striped_epi2_seqedit_getval(us[1][1], W, j);
				if(u2 != vals[b1][b2]){
					fflush(stdout);
					fprintf(stderr, " -- j=%d s=%c%c u=%d v=%d u2=%d should be %d in %s -- %s:%d --\n", j, "ACGTN"[qseq[j]], "ACGTN"[tseq[i]], u, v, vals[b1][b2], u2, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				if(b1 == 0 && b2 == 1) score ++;
				else if(b1 == 1 && b2 == 0) score --;
				fprintf(stdout, "%c%03d:%c:%c ", "ACGTN"[qseq[j]], score, "-0+"[vals[b1][b2] + 1], "-0+"[v2 + 1]);
			}
			fprintf(stdout, "\n");
		}
	}
	rs = backtrace(uts, qseq, qlen, tseq, tlen, cigars);
	if(mempb) free(mempb);
	return rs;
}

static void banded_striped_epi8_seqalign_piecex_row_init(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs, int *begs, int mode, u4i bandwidth, b1i max_nt, b1i min_nt, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	xint ZERO, MIN, GAP;
	u4i k, W, xp;
	int s, t;
	W = bandwidth / WORDSIZE;
	ZERO = mm_set1_epi8(0);
	MIN = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	if(mode == SEQALIGN_MODE_GLOBAL || mode == SEQALIGN_MODE_EXTEND){
		if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1 && (gapo1 - gapo2) / (gape1 - gape2) < Int(bandwidth)){
			xp = (gapo2 - gapo1) / (gape1 - gape2);
			GAP = mm_set1_epi8(gape2);
			for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, GAP);
			for(k=0;k<WORDSIZE;k++) ubegs[k] = gape2 * W;
			us[1][0] = gapo1 + gape1 + min_nt - max_nt;
			ubegs[0] += us[1][0] - gape2;
			for(k=1;k<=xp;k++){
				us[1][banded_striped_epi8_pos2idx(bandwidth, k)] = gape1;
				ubegs[k / W] += gape1 - gape2;
			}
		} else {
			GAP = mm_set1_epi8(gape1);
			for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, GAP);
			us[1][0] = gapo1 + gape1 + min_nt - max_nt;
			for(k=0;k<WORDSIZE;k++) ubegs[k] = gape1 * W;
			ubegs[0] += us[1][0] - gape1;
		}
		s = 0;
		for(k=0;k<WORDSIZE;k++){
			t = ubegs[k];
			ubegs[k] = s;
			s += t;
		}
		ubegs[k] = s;
		begs[-1] = - Int(bandwidth);
	} else {
		for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, ZERO);
		begs[-1] = 0;
		memset(ubegs, 0, (WORDSIZE + 1) * sizeof(int));
	}
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		for(k=0;k<W;k++) mm_store(((xint*)es[1]) + k, MIN);
		for(k=0;k<W;k++) mm_store(((xint*)qs[1]) + k, MIN);
	} else if(gapo1){
		for(k=0;k<W;k++) mm_store(((xint*)es[1]) + k, MIN);
	} else {
	}
}

static inline void banded_striped_epi8_seqalign_set_query_prof(u1i *qseq, u4i qlen, b1i *qprof, u4i bandwidth, b1i mtx[16]){
	b1i *qp;
	u4i xlen, x, pos, i, j, W, b;
	W = bandwidth / WORDSIZE;
	xlen = num_max(qlen, bandwidth); // In case of bandwidth > qlen
	for(x=0;x<=xlen;x++){ // leave the last block all SEQALIGN_SCORE_EPI8_MIN for accessing [W]
		qp = qprof + x * 4 * WORDSIZE;
		for(j=0;j<WORDSIZE;j++){
			pos = x + j * W;
			if(pos < qlen){
				b = qseq[pos];
				b *= 4;
				for(i=0;i<4;i++){
					qp[i * WORDSIZE + j] = mtx[b + i];
				}
			} else {
				for(i=0;i<4;i++){
					qp[i * WORDSIZE + j] = SEQALIGN_SCORE_EPI8_MIN;
				}
			}
		}
	}
}

static inline void banded_striped_epi8_seqalign_piecex_row_mov(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs, u4i W, u4i mov, int piecewise){
	xint u, x, UBS[4];
	u4i i, div;
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = mm_load(((xint*)us[1]) + i + mov);
		mm_store(((xint*)us[0]) + i, u);
	}
	if(piecewise){
		for(i=0;i<div;i++){
			x = mm_load(((xint*)es[1]) + i + mov);
			mm_store(((xint*)es[0]) + i, x);
		}
	}
	if(piecewise == 2){
		for(i=0;i<div;i++){
			x = mm_load(((xint*)qs[1]) + i + mov);
			mm_store(((xint*)qs[0]) + i, x);
		}
	}
	if(!mov) return;
	{
		u = mm_load(((xint*)us[1]) + i + mov - W);
		u = mm_srli(u, 1);
		u = mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, WORDSIZE - 1);
	}
	for(i=div;i<W;i++){
		mm_store(((xint*)us[0]) + i, u);
		u = mm_load(((xint*)us[1]) + (i + mov - W + 1) % W);
		u = mm_srli(u, 1);
	}
	{
		UBS[0] = mm_load(((xint*)ubegs) + 0);
		UBS[1] = mm_load(((xint*)ubegs) + 1);
		UBS[2] = mm_load(((xint*)ubegs) + 2);
		UBS[3] = mm_load(((xint*)ubegs) + 3);
	}
	for(i=0;i<mov;i++){
		u = mm_load(((xint*)us[1]) + i);
		UBS[0] = mm_add_epi32(UBS[0], mm_cvtepi8x0_epi32(u));
		UBS[1] = mm_add_epi32(UBS[1], mm_cvtepi8x1_epi32(u));
		UBS[2] = mm_add_epi32(UBS[2], mm_cvtepi8x2_epi32(u));
		UBS[3] = mm_add_epi32(UBS[3], mm_cvtepi8x3_epi32(u));
	}
	{
		mm_store(((xint*)ubegs) + 0, UBS[0]);
		mm_store(((xint*)ubegs) + 1, UBS[1]);
		mm_store(((xint*)ubegs) + 2, UBS[2]);
		mm_store(((xint*)ubegs) + 3, UBS[3]);
		ubegs[WORDSIZE] += SEQALIGN_SCORE_EPI8_MIN;
	}
	if(piecewise){
		for(;i<W;i++){
			x = mm_load(((xint*)es[1]) + i + mov - W);
			x = mm_srli(x, 1);
			mm_store(((xint*)es[0]) + i, x);
		}
	}
	if(piecewise == 2){
		for(;i<W;i++){
			x = mm_load(((xint*)qs[1]) + i + mov - W);
			x = mm_srli(x, 1);
			mm_store(((xint*)qs[0]) + i, x);
		}
	}
}

static inline void banded_striped_epi8_seqalign_piecex_row_merge(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *os[3], b1i osmarks[2], u4i W, int piecewise){
	xint s[2][4], m[2][4], u[2], x[2][4], e[2], c[4], mrk[2];
	u4i i, j, k, p;
	mrk[0] = mm_set1_epi8(osmarks[0]);
	mrk[1] = mm_set1_epi8(osmarks[1]);
	for(k=0;k<2;k++){
		for(j=0;j<4;j++){
			s[k][j] = mm_load(((xint*)ubegs[k]) + j);
		}
	}
	p = 0;
	for(j=0;j<4;j++){
		m[p][j] = mm_max_epi32(s[0][j], s[1][j]);
		mm_store(((xint*)ubegs[0]) + j, m[p][j]);
	}
	ubegs[0][WORDSIZE] = num_max(ubegs[0][WORDSIZE], ubegs[1][WORDSIZE]);
	for(i=0;i<W;i++){
		p = !p;
		for(k=0;k<2;k++){
			u[k] = mm_load(((xint*)us[k]) + i);
			s[k][0] = mm_add_epi32(s[k][0], mm_cvtepi8x0_epi32(u[k]));
			s[k][1] = mm_add_epi32(s[k][1], mm_cvtepi8x1_epi32(u[k]));
			s[k][2] = mm_add_epi32(s[k][2], mm_cvtepi8x2_epi32(u[k]));
			s[k][3] = mm_add_epi32(s[k][3], mm_cvtepi8x3_epi32(u[k]));
		}
		for(j=0;j<4;j++){
			c[j] = mm_cmpgt_epi32(s[1][j], s[0][j]);
			m[p][j] = mm_max_epi32(s[0][j], s[1][j]);
			m[!p][j] = mm_sub_epi32(m[p][j], m[!p][j]);
		}
		u[0] = mm_packs_epi32(m[!p][0], m[!p][1]);
		u[1] = mm_packs_epi32(m[!p][2], m[!p][3]);
		u[0] = mm_packs_epi16(u[0], u[1]);
		mm_store(((xint*)us[0]) + i, u[0]);
		c[0] = mm_packs_epi32(c[0], c[1]);
		c[1] = mm_packs_epi32(c[2], c[3]);
		c[0] = mm_packs_epi16(c[0], c[1]);
		c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
		mm_store(((xint*)os[0]) + i, c[0]);
		if(__builtin_constant_p(piecewise == 0)) continue;
		for(k=0;k<2;k++){
			e[k] = mm_load(((xint*)es[k]) + i);
			x[k][0] = mm_add_epi32(s[k][0], mm_cvtepi8x0_epi32(e[k]));
			x[k][1] = mm_add_epi32(s[k][1], mm_cvtepi8x1_epi32(e[k]));
			x[k][2] = mm_add_epi32(s[k][2], mm_cvtepi8x2_epi32(e[k]));
			x[k][3] = mm_add_epi32(s[k][3], mm_cvtepi8x3_epi32(e[k]));
		}
		for(j=0;j<4;j++){
			c[j] = mm_cmpgt_epi32(x[1][j], x[0][j]);
			m[!p][j] = mm_sub_epi32(m[p][j], m[!p][j]);
		}
		e[0] = mm_packs_epi32(m[!p][0], m[!p][1]);
		e[1] = mm_packs_epi32(m[!p][2], m[!p][3]);
		e[0] = mm_packs_epi16(e[0], e[1]);
		mm_store(((xint*)es[0]) + i, e[0]);
		c[0] = mm_packs_epi32(c[0], c[1]);
		c[1] = mm_packs_epi32(c[2], c[3]);
		c[0] = mm_packs_epi16(c[0], c[1]);
		c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
		mm_store(((xint*)os[1]) + i, c[0]);
		if(__builtin_constant_p(piecewise == 1)) continue;
		for(k=0;k<2;k++){
			e[k] = mm_load(((xint*)qs[k]) + i);
			x[k][0] = mm_add_epi32(s[k][0], mm_cvtepi8x0_epi32(e[k]));
			x[k][1] = mm_add_epi32(s[k][1], mm_cvtepi8x1_epi32(e[k]));
			x[k][2] = mm_add_epi32(s[k][2], mm_cvtepi8x2_epi32(e[k]));
			x[k][3] = mm_add_epi32(s[k][3], mm_cvtepi8x3_epi32(e[k]));
		}
		for(j=0;j<4;j++){
			c[j] = mm_cmpgt_epi32(x[1][j], x[0][j]);
			m[!p][j] = mm_sub_epi32(m[p][j], m[!p][j]);
		}
		e[0] = mm_packs_epi32(m[!p][0], m[!p][1]);
		e[1] = mm_packs_epi32(m[!p][2], m[!p][3]);
		e[0] = mm_packs_epi16(e[0], e[1]);
		mm_store(((xint*)qs[0]) + i, e[0]);
		c[0] = mm_packs_epi32(c[0], c[1]);
		c[1] = mm_packs_epi32(c[2], c[3]);
		c[0] = mm_packs_epi16(c[0], c[1]);
		c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
		mm_store(((xint*)os[2]) + i, c[0]);
	}
}

static inline __attribute__((always_inline)) int banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(xint h, xint u, xint v, b1i *us[2], b1i *fs, int *ubegs, u4i W, int ph, int rh, int bt0){
	u4i i;
	UNUSED(W);
	// revise the first striped block of u
	v = mm_subs_epi8(h, u); // v(x - 1, y) = h(x - 1, y) - u(x - 1, y - 1)
	mm_store((xint*)fs, v);
	for(i=1;i<=WORDSIZE;i++){
		ubegs[i] += fs[i - 1];
	}
	v = mm_slli(v, 1);
	u = mm_load(((xint*)us[1]) + 0);
	u = mm_subs_epi8(u, v); // because I previously set v to zero, now update it, u(x, y) = h(x, y) - v(x - 1, y)
	mm_store(((xint*)us[1]) + 0, u);
	// shift score to fit EPI8
	if(bt0 == SEQALIGN_BT_M){
		if(us[1][0] == 0){
			rh = ph - 1;
			us[1][0] = 1;
		} else {
			rh = ph;
		}
	} else {
		rh = ph + us[1][0];
		us[1][0] = 0;
	}
	ubegs[0] = rh;
	return rh;
}

static inline __attribute__((always_inline)) xint banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(u4i W, xint f, b1i *fs, int *ubegs, b1i gape){
	int i, s, t;
	f = mm_slli(f, 1);
	mm_store((xint*)fs, f);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape;
	s = t + fs[0] - (ubegs[1] - ubegs[0]);
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s) fs[i] = s;
		s = t + fs[i] - (ubegs[i + 1] - ubegs[i]);
	}
	f = mm_load((xint*)fs);
	return f;
}

static inline int banded_striped_epi8_seqalign_piece0_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, f, u, v, c, d, z;
	xint I, D, GapE;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(es);
	UNUSED(qs);
	UNUSED(gapo1);
	UNUSED(gapo2);
	UNUSED(gape2);
	UNUSED(mov);
	I     = mm_set1_epi8(SEQALIGN_BT_I);
	D     = mm_set1_epi8(SEQALIGN_BT_D);
	GapE  = mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + gape1; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
		bt0 = SEQALIGN_BT_M;
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
		bt0 = SEQALIGN_BT_D;
	}
	// preparing initial f for each running block
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	h = mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;i++){
		u = mm_load(((xint*)us[0]) + i);
		// max h, e, f
		e = mm_adds_epi8(u, GapE);
		h = mm_max_epi8(e, h);
		h = mm_max_epi8(f, h);
		// preparing next f
		f = mm_adds_epi8(h, GapE);
		f = mm_subs_epi8(f, u);
		// preparing next h
		h = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	f = banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(W, f, fs, ubegs, gape1);
	// main loop
	// don't use the h from last W - 1 to calculate v = h - u, because this h maybe updated when F-penetration
	// will revise v after this loop
	v = mm_set1_epi8(0);
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	z = mm_insert_epi8(h, h0, 0);
	u = mm_set1_epi8(0); // useless, but for compiler
	for(i=0;i<W;i++){
		u = mm_load(((xint*)us[0]) + i);
		// max(e, h)
		e = mm_adds_epi8(u, GapE);
		c = mm_cmpgt_epi8(e, z);
		d = mm_and(c, D); // bt
		h = mm_max_epi8(e, z);
		// max(f, h)
		c = mm_cmpgt_epi8(f, h);
		d = mm_blendv(d, I, c);
		h = mm_max_epi8(f, h);
		mm_store(((xint*)bs) + i, d);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate f(x, y)
		f = mm_adds_epi8(h, GapE);
		f = mm_subs_epi8(f, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W, ph, rh, bt0);
}

static inline int banded_striped_epi8_seqalign_piece1_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, f, u, v, c, d, z;
	xint I, IE, D, DE, GapOE, GapE;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(qs);
	UNUSED(gapo2);
	UNUSED(gape2);
	UNUSED(mov);
	I     = mm_set1_epi8(SEQALIGN_BT_I);
	IE    = mm_set1_epi8(SEQALIGN_BT1_IE);
	D     = mm_set1_epi8(SEQALIGN_BT_D);
	DE    = mm_set1_epi8(SEQALIGN_BT1_DE);
	GapOE = mm_set1_epi8(gapo1 + gape1);
	GapE  = mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + es[0][0]; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
		bt0 = SEQALIGN_BT_M;
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
		bt0 = SEQALIGN_BT_D;
	}
	// preparing initial f for each running block
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	h = mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;i++){
		u = mm_load(((xint*)us[0]) + i);
		e = mm_load(((xint*)es[0]) + i);
		// max h, e, f
		e = mm_adds_epi8(e, u);
		h = mm_max_epi8(e, h);
		h = mm_max_epi8(f, h);
		// preparing next f
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		f = mm_max_epi8(f, h);
		f = mm_subs_epi8(f, u);
		// preparing next h
		h = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	f = banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(W, f, fs, ubegs, gape1);
	// main loop
	// don't use the h from last W - 1 to calculate v = h - u, because this h maybe updated when F-penetration
	// will revise v after this loop
	u = mm_set1_epi8(0); // useless, but for compiler
	v = mm_set1_epi8(0);
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	z = mm_insert_epi8(h, h0, 0);
	for(i=0;__builtin_expect(i<W, 0);i++){
		u = mm_load(((xint*)us[0]) + i);
		e = mm_load(((xint*)es[0]) + i);
		// max(e, h)
		e = mm_adds_epi8(e, u);
		c = mm_cmpgt_epi8(e, z);
		d = mm_and(c, D); // bt
		h = mm_max_epi8(e, z);
		// max(f, h)
		c = mm_cmpgt_epi8(f, h);
		d = mm_blendv(d, I, c);
		h = mm_max_epi8(f, h);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate e(x, y)
		e = mm_adds_epi8(e, GapE);
		e = mm_subs_epi8(e, h);
		c = mm_cmpgt_epi8(e, GapOE);
		d = mm_or(d, mm_and(c, DE));
		e = mm_max_epi8(e, GapOE);
		mm_store(((xint*)es[1]) + i, e);
		// calculate f(x, y)
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		c = mm_cmpgt_epi8(f, h);
		d = mm_or(d, mm_and(c, IE));
		f = mm_max_epi8(f, h);
		mm_store(((xint*)bs) + i, d);
		f = mm_subs_epi8(f, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	h = mm_subs_epi8(h, GapOE);
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W, ph, rh, bt0);
}

static inline int banded_striped_epi8_seqalign_piece2_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, q, f, g, u, v, c, d, z;
	xint I1, I2, IE1, IE2, D1, D2, DE1, DE2, GapOE, GapE, GapQP, GapP, GapOQ;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(mov);
	I1     = mm_set1_epi8(SEQALIGN_BT2_I1);
	I2     = mm_set1_epi8(SEQALIGN_BT2_I2);
	IE1    = mm_set1_epi8(SEQALIGN_BT2_IE1);
	IE2    = mm_set1_epi8(SEQALIGN_BT2_IE2);
	D1     = mm_set1_epi8(SEQALIGN_BT2_D1);
	D2     = mm_set1_epi8(SEQALIGN_BT2_D2);
	DE1    = mm_set1_epi8(SEQALIGN_BT2_DE1);
	DE2    = mm_set1_epi8(SEQALIGN_BT2_DE2);
	GapOE = mm_set1_epi8(gapo1 + gape1);
	GapE  = mm_set1_epi8(gape1);
	GapQP = mm_set1_epi8(gapo2 + gape2);
	GapP  = mm_set1_epi8(gape2);
	GapOQ = mm_subs_epi8(GapOE, GapQP);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	g = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + num_max(es[0][0], qs[0][0]); // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
		bt0 = SEQALIGN_BT_M;
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
		bt0 = SEQALIGN_BT_D;
	}
	// preparing initial f for each running block
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	h = mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;i++){
		u = mm_load(((xint*)us[0]) + i);
		// max h, e, q, f, g
		e = mm_load(((xint*)es[0]) + i);
		q = mm_load(((xint*)qs[0]) + i);
		e = mm_adds_epi8(e, u);
		q = mm_adds_epi8(q, u);
		h = mm_max_epi8(e, h);
		h = mm_max_epi8(q, h);
		h = mm_max_epi8(f, h);
		h = mm_max_epi8(g, h);
		// preparing next f and g
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		f = mm_max_epi8(f, h);
		f = mm_subs_epi8(f, u);
		g = mm_adds_epi8(g, GapP);
		h = mm_subs_epi8(h, GapOQ); // gapo1 + gape1 - gapo2 - gape2
		g = mm_max_epi8(g, h);
		g = mm_subs_epi8(g, u);
		// preparing next h
		h = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	f = banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(W, f, fs, ubegs, gape1);
	g = banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(W, g, fs, ubegs, gape2);
	// main loop
	v = mm_set1_epi8(0);
	h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
	z = mm_insert_epi8(h, h0, 0);
	u = mm_set1_epi8(0); // useless, but for compiler
	for(i=0;i<W;i++){
		u = mm_load(((xint*)us[0]) + i);
		e = mm_load(((xint*)es[0]) + i);
		q = mm_load(((xint*)qs[0]) + i);
		// max(e, q, h)
		e = mm_adds_epi8(e, u);
		c = mm_cmpgt_epi8(e, z);
		d = mm_and(c, D1); // bt
		h = mm_max_epi8(e, z);
		q = mm_adds_epi8(q, u);
		c = mm_cmpgt_epi8(q, h);
		d = mm_blendv(d, D2, c);
		h = mm_max_epi8(q, h);
		// max(f, g, h)
		c = mm_cmpgt_epi8(f, h);
		d = mm_blendv(d, I1, c);
		h = mm_max_epi8(f, h);
		c = mm_cmpgt_epi8(g, h);
		d = mm_blendv(d, I2, c);
		h = mm_max_epi8(g, h);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate e(x, y) and q(x, y)
		e = mm_adds_epi8(e, GapE);
		e = mm_subs_epi8(e, h);
		c = mm_cmpgt_epi8(e, GapOE);
		d = mm_or(d, mm_and(c, DE1));
		e = mm_max_epi8(e, GapOE);
		mm_store(((xint*)es[1]) + i, e);
		q = mm_adds_epi8(q, GapP);
		q = mm_subs_epi8(q, h);
		c = mm_cmpgt_epi8(q, GapQP);
		d = mm_or(d, mm_and(c, DE2));
		q = mm_max_epi8(q, GapQP);
		mm_store(((xint*)qs[1]) + i, q);
		// calculate f(x, y) and g(x, y)
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		c = mm_cmpgt_epi8(f, h);
		d = mm_or(d, mm_and(c, IE1));
		f = mm_max_epi8(f, h);
		f = mm_subs_epi8(f, u);
		g = mm_adds_epi8(g, GapP);
		h = mm_subs_epi8(h, GapOQ); // (gapo1 + gape1 - gapo2 - gape2)
		c = mm_cmpgt_epi8(g, h);
		d = mm_or(d, mm_and(c, IE2));
		g = mm_max_epi8(g, h);
		mm_store(((xint*)bs) + i, d);
		g = mm_subs_epi8(g, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	h = mm_subs_epi8(h, GapQP);
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W, ph, rh, bt0);
}

static inline u4i banded_striped_epi8_seqalign_row_max(b1i *us, int *ubegs, u4i W, int *max_score){
	xint h, c, xi, Max[4], Scr[4], Idx[4], Num[2];
	u4i i, x;
	int ary[3][WORDSIZE/4];
	for(i=0;i<4;i++){
		Scr[i] = mm_load(((xint*)ubegs) + i);
		Max[i] = mm_set1_epi32(SEQALIGN_SCORE_MIN);
		Idx[i] = mm_set1_epi32(0);
	}
	for(i=0;i<W;i++){
		h = mm_load(((xint*)us) + i);
		Scr[0] = mm_add_epi32(Scr[0], mm_cvtepi8x0_epi32(h));
		Scr[1] = mm_add_epi32(Scr[1], mm_cvtepi8x1_epi32(h));
		Scr[2] = mm_add_epi32(Scr[2], mm_cvtepi8x2_epi32(h));
		Scr[3] = mm_add_epi32(Scr[3], mm_cvtepi8x3_epi32(h));
		xi = mm_set1_epi32(i);
		Idx[0] = mm_blendv(Idx[0], xi, mm_cmpgt_epi32(Scr[0], Max[0])); Max[0] = mm_max_epi32(Max[0], Scr[0]);
		Idx[1] = mm_blendv(Idx[1], xi, mm_cmpgt_epi32(Scr[1], Max[1])); Max[1] = mm_max_epi32(Max[1], Scr[1]);
		Idx[2] = mm_blendv(Idx[2], xi, mm_cmpgt_epi32(Scr[2], Max[2])); Max[2] = mm_max_epi32(Max[2], Scr[2]);
		Idx[3] = mm_blendv(Idx[3], xi, mm_cmpgt_epi32(Scr[3], Max[3])); Max[3] = mm_max_epi32(Max[3], Scr[3]);
	}
	if(0){
		int rs[WORDSIZE];
		mm_store(((xint*)rs) + 0, Scr[0]);
		mm_store(((xint*)rs) + 1, Scr[1]);
		mm_store(((xint*)rs) + 2, Scr[2]);
		mm_store(((xint*)rs) + 3, Scr[3]);
		for(i=1;i<=WORDSIZE;i++){
			if(ubegs[i] != rs[i-1]){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	c = mm_cmpgt_epi32(Max[1], Max[0]); Idx[0] = mm_blendv(Idx[0], Idx[1], c); Num[0] = mm_blendv(mm_set1_epi32(0), mm_set1_epi32(1), c); Max[0] = mm_max_epi32(Max[0], Max[1]);
	c = mm_cmpgt_epi32(Max[3], Max[2]); Idx[1] = mm_blendv(Idx[2], Idx[3], c); Num[1] = mm_blendv(mm_set1_epi32(2), mm_set1_epi32(3), c); Max[1] = mm_max_epi32(Max[2], Max[3]);
	c = mm_cmpgt_epi32(Max[1], Max[0]); Idx[0] = mm_blendv(Idx[0], Idx[1], c); Num[0] = mm_blendv(Num[0], Num[1], c); Max[0] = mm_max_epi32(Max[0], Max[1]);
	mm_store(((xint*)ary[0]), Max[0]);
	mm_store(((xint*)ary[1]), Idx[0]);
	mm_store(((xint*)ary[2]), Num[0]);
	*max_score = ary[0][0];
	x = 0;
	for(i=1;i<WORDSIZE/4;i++){
		if(ary[0][i] > *max_score){
			*max_score = ary[0][i];
			x = i;
		}
	}
	x = (ary[2][x] * WORDSIZE / 4 + x) * W + ary[1][x];
	return x;
}

static inline int banded_striped_epi8_seqalign_band_mov(b1i *us, int *ubegs, u4i W, u4i tidx, u4i qoff, u4i qlen){
	u4i i;
	int noisy;
	UNUSED(us);
	if(tidx <= W * WORDSIZE / 4) return 0;
	if(qoff + W * WORDSIZE >= qlen) return 0;
	noisy = 0;
	for(i=1;i<=WORDSIZE;i++){
		noisy += num_diff(ubegs[i], ubegs[i-1]);
	}
	noisy = num_max(2 * WORDSIZE / 2, noisy / WORDSIZE / W * WORDSIZE / 2);
	if(ubegs[0] + noisy < ubegs[WORDSIZE]){
		return 2;
	} else if(ubegs[0] > ubegs[WORDSIZE] + noisy){
		return 0;
	} else {
		return 1;
	}
}

static inline u4i banded_striped_epi8_seqalign_band_comb(u4i cnt, u4i *qoffs, int **ubegs, b1v *mempool, u4i W, u4i qlen){
	typedef struct { u4i qpos[2]; int sval; } run_t;
	run_t *runs, *r2;
	int *row, *ary, smax, s;
	u4i i, j, xlen, ri, rn, rx, ry, x, y, mpsz[3];
	if(W * WORDSIZE >= qlen) return 0; // full band
	xlen = roundup_times(qlen, WORDSIZE);
	mpsz[0] = xlen * sizeof(int); // row
	mpsz[1] = mpsz[0] + (rn = cnt * WORDSIZE) * sizeof(run_t); // run_t
	mpsz[2] = mpsz[1] + 2 * rn * sizeof(int);
	clear_and_encap_b1v(mempool, mpsz[2]);
	row   = (int*)(mempool->buffer + 0);
	runs  = (run_t*)(mempool->buffer + mpsz[0]);
	ary   = (int*)(mempool->buffer + mpsz[1]);
	memset(row, 0x80, xlen * sizeof(int)); // set to min int
	for(ri=i=0;i<cnt;i++){
		for(j=0;j<W;j++){
			runs[ri].qpos[0] = qoffs[i] + j * W;
			runs[ri].qpos[1] = qoffs[i] + j * W + W;
			runs[ri].sval = (ubegs[i][j] + ubegs[i][j]) / 2;
			ri ++;
		}
	}
	// set max score for cells of row
	sort_array(runs, rn, run_t, num_cmpgt(a.qpos[1],  b.qpos[1]));
	for(i=0;i<2*rn;i++) ary[i] = i;
	sort_array(ary, rn * 2, int, num_cmpgtx(runs[a >> 1].qpos[a & 0x1], runs[b >> 1].qpos[b & 0x1], b & 0x1, a & 0x1));
	rx = ry = 0;
	x  = runs[ary[0] >> 1].qpos[0];
	for(ri=0;ri<rn*2;ri++){
		if((ary[ri] & 0x01) == 0){
			ry = num_max(ary[ri] >> 1, Int(ry));
		}
		r2 = runs + (ary[ri] >> 1);
		y = r2->qpos[0] + (ary[ri] & 0x1) * W;
		if(y > x){
			smax = - MAX_B4;
			while(runs[rx].qpos[1] <= x) rx ++;
			for(i=rx;i<=ry;i++){
				if(runs[i].qpos[0] > x) continue;
				if(runs[i].sval > smax) smax = runs[i].sval;
			}
			for(i=x;i<y;i++) row[i] = smax;
		}
		x = y;
	}
	// find a band with max sums of scores, but there MUST NOT gap (== (int)0x80808080) inside the band
	smax = (int)0x80808080;
	rx = MAX_U4;
	x = 0;
	while(x + W * WORDSIZE <= qlen){
		while(x + W * WORDSIZE < qlen && row[x] == (int)0x80808080) x ++;
		s = 0;
		for(i=0;i<W*WORDSIZE;i++){
			if(row[x + i] == (int)0x80808080) break;
			s += row[x + i];
		}
		if(i < W * WORDSIZE){
			x = i + 1;
			continue;
		}
		if(s > smax){
			smax = s;
			rx = x;
		}
		while(x + W * WORDSIZE < qlen){
			if(row[x + W * WORDSIZE] == (int)0x80808080){
				x = x + W * WORDSIZE + 1;
				break;
			}
			s -= row[x];
			s += row[x + W * WORDSIZE];
			x ++;
			if(s > smax){
				smax = s;
				rx = x;
			}
		}
	}
	// prepare next move for max band
	s = 0;
	for(i=1;i<WORDSIZE;i++){
		s += num_diff(row[rx + i * W], row[rx + (i + 1) * W]);
	}
	s = num_max(2 * WORDSIZE / 2, s / (WORDSIZE - 1) / W * WORDSIZE / 2);
}

static inline int banded_striped_epi8_seqalign_getscore(b1i *us, int *ubegs, u4i W, u4i pos){
	u4i i, x, y;
	int s;
	x = pos % W;
	y = pos / W;
	s = ubegs[y];
	for(i=0;i<=x;i++){
		s += us[i * WORDSIZE + y];
	}
	return s;
}

static inline void banded_striped_epi8_seqalign_row_print(FILE *out, u1i *qseq, u4i tidx, u4i tpos, u1i tbase, u4i bandwidth, u4i mov, u4i rbeg, u4i rmax, int max_score, int *ubegs, b1i *us, b1i *bs, int detail){
	u4i i, x, W;
	int score;
	W = bandwidth / WORDSIZE;
	fprintf(out, "ROW[%d][%d][%c]\tMOV=%d\tBAND=%d,%d\tMAX=%d(%d),%d", tidx, tpos, "ACGTN-"[tbase], mov, rbeg, rbeg + bandwidth, rbeg + rmax, rmax, max_score);
	if(detail > 2){
		score = ubegs[0];
		for(i=0;i<bandwidth;i++){
			x = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
			fprintf(out, "\t%d:%c%d:%d:%d", i + rbeg, "ACGTN-"[qseq[rbeg + i]], score + us[x], us[x], bs[x] * 0xf);
			score += us[x];
		}
	}
	fprintf(out, "\n");
	if(detail > 1){
		for(i=0;i<=WORDSIZE;i++){
			fprintf(out, "\t%d", ubegs[i]);
		}
		fprintf(out, "\n");
	}
	fflush(out);
}

static inline int array_epi32_sum(int *rs, int n){
	xint v, s;
	int i, m, ary[WORDSIZE / 4], sum;
	s = mm_set1_epi32(0);
	m = n / (WORDSIZE / 4);
	for(i=0;i<m;i++){
		v = mm_load(((xint*)rs) + i);
		s = mm_add_epi32(s, v);
	}
	mm_store((xint*)ary, s);
	sum = 0;
	for(i=0;i<WORDSIZE/4;i++) sum += ary[i];
	for(i=m*(WORDSIZE/4);i<n;i++){
		sum += rs[i];
	}
	return sum;
}

static inline int banded_striped_epi8_seqalign_piecex_row_verify(int rbeg, int W, int ph, int rh, int hh, u1i tbase, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	int i, x, s1, s2, t, h, e, q, f, g, piecewise, h0;
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		piecewise = 2;
		t = us[0][0] + num_max(es[0][0], qs[0][0]);
	} else if(gapo1){
		piecewise = 1;
		t = us[0][0] + es[0][0];
	} else {
		piecewise = 0;
		t = us[0][0] + gape1;
	}
	f = g = SEQALIGN_SCORE_MIN;
	h0 = (rh - ph) + qprof[((rbeg + 0) * 4 + tbase) * WORDSIZE];
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX;
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
	}
	s1 = ph;
	s2 = hh;
	//fprintf(stdout, " -- %s -- %s -- piecewise=%d rbeg=%d\tph=%d\trh=%d\thh=%d\n", __FUNCTION__, __FILE__, piecewise, rbeg, ph, rh, hh); fflush(stdout);
	for(i=0;i<W*WORDSIZE;i++){
		x = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
		s2 += us[1][x];
		e = s1 + us[0][x] + (piecewise != 0? es[0][x] : gape1);
		q = s1 + us[0][x] + (piecewise == 2? qs[0][x] : SEQALIGN_SCORE_MIN);
		h = num_max(s1 + h0, e);
		h = num_max(h, q);
		h = num_max(h, f);
		h = num_max(h, g);
		if(h != s2){
			fprintf(stderr, " -- %c:%d i=%d h=%d(%d,%d,%d,%d) s2=%d s1=%d u0=%d u1=%d something wrong in %s -- %s:%d --\n", "ACGTN-"[tbase], h0, i, h, e, q, f, g, s2, s1, us[0][x], us[1][x], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			//abort();
			return 1;
		}
		f = num_max(f, h + gapo1) + gape1;
		g = (piecewise == 2)? num_max(g + gape2, h + gapo2 + gape2) : SEQALIGN_SCORE_MIN;
		s1 += us[0][x];
		h0 = banded_striped_epi8_seqalign_get_qprof_value(qprof, rbeg + i + 1, tbase);
	}
	return 0;
}

static inline u4i banded_striped_epi8_seqalign_piecex_backtrace(u1i *qseq, u1i *tseq, b1i *bs, int *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, u4v *cigars){
	u4i btx, bty, cg;
	u1i qbase, tbase, op, bt;
	rs->qb = rs->qe; rs->qe ++;
	rs->tb = rs->te; rs->te ++;
	rs->mat = rs->mis = rs->ins = rs->del = rs->aln = 0;
	bty = 0;
	cg = 0;
	if(cigars) clear_u4v(cigars);
	while(1){
		rs->aln ++;
		btx = bs[rs->tb * bandwidth + banded_striped_epi8_pos2idx(bandwidth, rs->qb - begs[rs->tb])];
		if(piecewise == 0){
			bt = bty = btx;
		} else if(piecewise == 1){
			if(bty == SEQALIGN_BT_M){
				bt = bty = btx & 0x3;
			} else {
				bt = bty = ((btx >> (bt + 1)) & 0x1)? bty : SEQALIGN_BT_M;
			}
		} else {
			if(bty == SEQALIGN_BT_M){
				bty = btx & 0xf;
			} else {
				bty = ((btx >> (bt + 2)) & 0x1)? bty : SEQALIGN_BT_M;
			}
			bt = (0b1001100100 >> (bty * 2)) & 0x3;
		}
		if(bt == SEQALIGN_BT_M){
			if(qseq && tseq){
				qbase = qseq[rs->qb];
				tbase = tseq[rs->tb];
				if(qbase == tbase){
					rs->mat ++;
					op = 7; // =
				} else {
					rs->mis ++;
					op = 8; // X
				}
			} else {
				rs->mat ++;
				op = 0; // M
			}
			rs->qb --;
			rs->tb --;
			if(rs->qb < 0 || rs->tb < 0) break;
		} else if(bt == SEQALIGN_BT_D){
			rs->del ++;
			rs->tb --;
			op = 2; // D
			if(rs->tb < 0) break;
		} else {
			rs->ins ++;
			rs->qb --;
			op = 1; // I
			if(rs->qb < 0) break;
		}
		if(op == (cg & 0xf)){
			cg += 0x10;
		} else {
			if(cg && cigars) push_u4v(cigars, cg);
			cg = 0x10 | op;
		}
	}
	if(op == (cg & 0xf)){
		cg += 0x10;
		if(cigars) push_u4v(cigars, cg);
	} else {
		if(cg && cigars) push_u4v(cigars, cg);
		cg = 0x10 | op;
		if(cigars) push_u4v(cigars, cg);
	}
	rs->qb ++;
	rs->tb ++;
	if(cigars) reverse_u4v(cigars);
	return rs->aln;
}

/*
 * backtrace by recurrence score
 * u and v are used to backtrace
 * if s(x, y) == v(x, y) + u(x, y - 1); then match
 * else if(gapo + n * gape == sum(u(x - n + 1 ... x, y)); then insertion
 * else if(gapo + n * gape == sum(v(x, y - n + 1 ... y)); then deletion
 * also special opertations on the boundrary of band
 *
 */
static inline u4i banded_striped_epi8_seqalign_piecex_trace(u1i *qseq, u1i *tseq, b1i ***ups, b1i ***vps, int *begs, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, seqalign_result_t *rs, u4v *cigars){
	u4i bt, cg, op, sz;
	int piecewise, x, y, s, h, hx, hy, u, v, gap1, gap2, W, state;
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		piecewise = 2;
	} else if(gapo1){
		piecewise = 1;
	} else {
		piecewise = 0;
	}
	rs->qb = rs->qe; rs->qe ++;
	rs->tb = rs->te; rs->te ++;
	rs->mat = rs->mis = rs->ins = rs->del = rs->aln = 0;
	cg = 0;
	W = bandwidth / WORDSIZE;
	if(cigars) clear_u4v(cigars);
	while(1){
#define fetch_ups(x, y) ups[(y) + 1][((x) - begs[(y)]) % W][((x) - begs[(y)]) / W]
#define fetch_vps(x, y) vps[(y) + 1][((x) - begs[(y)]) % W][((x) - begs[(y)]) / W]
		s = matrix[qseq[rs->qb] * 4 + tseq[rs->tb]];
		v = fetch_vps(rs->qb, rs->tb);
		if(rs->qb == Int(bandwidth + begs[rs->tb - 1])){ // here v indicates SEQALIGN_BT_M
			h = v? s : s + 1;
		} else if(rs->qb > Int(bandwidth + begs[rs->tb - 1])){ // SEQALIGN_BT_I only
			h = s + 1;
		} else if(rs->qb == begs[rs->tb]){ // here u indicates whether SEQALIGN_BT_M
			u = fetch_ups(rs->qb, rs->tb);
			h = u? s : s + 1;
		} else {
			u = fetch_ups(rs->qb, rs->tb - 1);
			h = u + v;
		}
		sz = 0;
		bt = 0;
		if(h == s){
			rs->mat ++;
			//printf("%c - %c\t%d,%d\n", "ACGTN-"[qseq[rs->qb]], "ACGTN-"[tseq[rs->tb]], rs->qb, rs->tb);
			//fflush(stdout);
			rs->qb --;
			rs->tb --;
			sz = 1;
		} else {
			hx = hy = 0;
			x = rs->qb; y = rs->tb;
			gap1 = gapo1 + gape1;
			gap2 = gapo2 + gape2;
			state = 0;
			while(1){
				if(x > (int)begs[rs->tb]){
					hx += fetch_ups(x, rs->tb);
					if(gap1 == hx){
						bt = 1;
					} else if(piecewise == 2 && gap2 == hx){
						bt = 1;
					}
					x --;
				} else state |= 1;
				if(y > 0 && rs->qb < Int(begs[y - 1] + bandwidth)){ // yes, begs[y - 1], not begs[y]
					hy += fetch_vps(rs->qb, y);
					if(gap1 == hy){
						bt = 2;
					} else if(piecewise == 2 && gap2 == hy){
						bt = 2;
					}
					y --;
				} else state |= 2;
				gap1 += gape1;
				gap2 += gape2;
				if(bt) break;
				if(state == 3){
					fprintf(stderr, " -- Backtrace Error at (%d,%d) (%d,%d) in %s -- %s:%d --\n", rs->qb, rs->tb, rs->qe, rs->te, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					return 0;
				}
			}
			//int i;
			if(bt == 1){
				sz = rs->qb - x;
				//for(i=0;i<sz;i++){
					//printf("%c | -\t%d,%d\n", "ACGTN-"[qseq[rs->qb - i]], rs->qb - i, rs->tb);
					//fflush(stdout);
				//}
				rs->qb = x;
				rs->ins += sz;
			} else if(bt == 2){
				sz = rs->tb - y;
				//for(i=0;i<sz;i++){
					//printf("- | %c\t%d,%d\n", "ACGTN-"[tseq[rs->tb - i]], rs->qb, rs->tb - i);
					//fflush(stdout);
				//}
				rs->tb = y;
				rs->del += sz;
			}
		}
		op = bt;
		//printf("%d:%d\n", op, sz);
		if(op == (cg & 0xf)){
			cg += 0x10 * sz;
		} else {
			if(cg && cigars){
				//printf("PUSH %d:%d\n", cg & 0xf, cg >> 4);
				push_u4v(cigars, cg);
			}
			cg = (0x10 * sz) | op;
		}
		rs->aln += sz;
		if(rs->qb < 0 || rs->tb < 0){
			if(cg && cigars){
				//printf("PUSH %d:%d\n", cg & 0xf, cg >> 4);
				push_u4v(cigars, cg);
			}
			break;
		}
	}
	rs->qb ++;
	rs->tb ++;
	if(cigars) reverse_u4v(cigars);
	return rs->aln;
}

static inline u4i seqalign_cigar2alnstr(u1i *qseq, u1i *tseq, seqalign_result_t *rs, u4v *cigars, char *alnstr[3], u4i length){
	u4i i, j, x, y, z, op, sz;
	if(length == 0) return 0;
	z = 0;
	x = rs->qb;
	y = rs->tb;
	for(i=0;i<cigars->size;i++){
		op = cigars->buffer[i] & 0xf;
		sz = cigars->buffer[i] >> 4;
		sz = num_min(sz, length - z);
		switch(op){
			case 0:
			case 7:
			case 8:
			for(j=0;j<sz;j++){
				alnstr[2][z] = (qseq[x] == tseq[y])? '|' : '*';
				alnstr[0][z] = "ACGTN-"[qseq[x ++]];
				alnstr[1][z] = "ACGTN-"[tseq[y ++]];
				z ++;
			}
			break;
			case 1:
			case 4:
			for(j=0;j<sz;j++){
				alnstr[2][z] = '-';
				alnstr[0][z] = "ACGTN-"[qseq[x ++]];
				alnstr[1][z] = '-';
				z ++;
			}
			break;
			case 2:
			case 3:
			for(j=0;j<sz;j++){
				alnstr[2][z] = '-';
				alnstr[0][z] = '-';
				alnstr[1][z] = "ACGTN-"[tseq[y ++]];
				z ++;
			}
			if(z == length) break;
		}
	}
	alnstr[0][z] = 0;
	alnstr[1][z] = 0;
	alnstr[2][z] = 0;
	return z;
}

static inline seqalign_result_t banded_striped_epi8_seqalign_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int mode, u4i bandwidth, b1i matrix[16], b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, int verbose){
	banded_striped_epi8_seqalign_set_query_prof_func    set_qprof;
	banded_striped_epi8_seqalign_piecex_row_init_func   row_init;
	banded_striped_epi8_seqalign_piecex_row_mov_func    row_mov;
	banded_striped_epi8_seqalign_piecex_row_cal_func    row_cal;
	banded_striped_epi8_seqalign_piecex_row_verify_func row_verify;
	banded_striped_epi8_seqalign_row_max_func           row_max;
	banded_striped_epi8_seqalign_band_mov_func          band_mov;
	banded_striped_epi8_seqalign_getscore_func          getscore;
	banded_striped_epi8_seqalign_piecex_backtrace_func  backtrace;
	seqalign_result_t rs;
	b1i *memp, *mempb, *qprof, *us[2], *es[2], *qs[2], *bs, *bs0;
	u8i mpsize;
	u4i i, W, rbeg, rmax, mov;
	int piecewise, smax, smin, lst_score, score, ph, rh, hh, rbx, rby, rbz, *begs, ubegs[WORDSIZE + 1];
	u1i tbase;
	set_qprof  = banded_striped_epi8_seqalign_set_query_prof;
	row_init   = banded_striped_epi8_seqalign_piecex_row_init;
	row_mov    = banded_striped_epi8_seqalign_piecex_row_mov;
	row_verify = banded_striped_epi8_seqalign_piecex_row_verify;
	row_max    = banded_striped_epi8_seqalign_row_max;
	band_mov   = banded_striped_epi8_seqalign_band_mov;
	getscore   = banded_striped_epi8_seqalign_getscore;
	backtrace  = banded_striped_epi8_seqalign_piecex_backtrace;
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		piecewise = 2;
		row_cal = banded_striped_epi8_seqalign_piece2_row_cal;
	} else if(gapo1){
		piecewise = 1;
		row_cal = banded_striped_epi8_seqalign_piece1_row_cal;
	} else {
		piecewise = 0;
		row_cal = banded_striped_epi8_seqalign_piece0_row_cal;
	}
	bandwidth = roundup_times(bandwidth, WORDSIZE);
	W = bandwidth / WORDSIZE;
	if(verbose){
		fprintf(stdout, "[%d,%d][%d,%d] PIECEWISE=%d\tW=%d\n", gapo1, gape1, gapo2, gape2, piecewise, W);
	}
	smax = - MAX_B1;;
	smin = MAX_B1;
	for(i=0;i<16;i++){
		smax = num_max(smax, matrix[i]);
		smin = num_min(smin, matrix[i]);
	}
	// allocate memory
	mpsize = WORDSIZE;
	mpsize += banded_striped_epi8_seqalign_qprof_size(qlen, bandwidth); // qprof[]
	mpsize += bandwidth * ((piecewise + 1) * 2 + tlen); // ups[][], vps[][], us[0/1][], es[0/1][], qs[0/1][], bs
	mpsize += roundup_times((tlen + 1) * sizeof(int), WORDSIZE); // row offset
	if(mempool){
		clear_and_encap_b1v(mempool, mpsize);
		memp = mempool->buffer + WORDSIZE;
		mempb = NULL;
	} else {
		mempb = malloc(mpsize);
		memp = mempb + WORDSIZE;
	}
	qprof = memp; memp += banded_striped_epi8_seqalign_qprof_size(qlen, bandwidth);
	us[0] = (b1i*)memp; memp += bandwidth;
	us[1] = (b1i*)memp; memp += bandwidth;
	if(piecewise){
		es[0] = (b1i*)memp; memp += bandwidth;
		es[1] = (b1i*)memp; memp += bandwidth;
	}
	if(piecewise == 2){
		qs[0] = (b1i*)memp; memp += bandwidth;
		qs[1] = (b1i*)memp; memp += bandwidth;
	}
	bs = bs0 = (b1i*)memp; memp += bandwidth * tlen;
	begs = ((int*)memp) + 1; memp += roundup_times((tlen + 1) * sizeof(u4i), WORDSIZE);
	// prepare
	set_qprof(qseq, qlen, qprof, bandwidth, matrix);
	memset(&rs, 0, sizeof(seqalign_result_t));
	rs.score = SEQALIGN_SCORE_MIN;
	row_init(us, es, qs, ubegs, begs, mode, bandwidth, smax, smin, gapo1, gape1, gapo2, gape2);
	// loop rows
	rbeg = rmax = 0;
	mov = 0;
	hh = 0;
	lst_score = 0;
	for(i=0;i<tlen;i++){
		ph = hh;
		tbase = tseq[i];
		if(mov && rbeg + bandwidth < qlen){
			mov = num_min(mov, num_max(0, Int(qlen) - Int(rbeg + bandwidth)));
			rbeg += mov;
			ph = getscore(us[1], ubegs, W, mov - 1);
			rh = ph + matrix[tbase * 4 + qseq[rbeg]]; // just guess
		} else {
			mov  = 0;
			rh = (rbeg == 0 && (mode == SEQALIGN_MODE_OVERLAP || i == 0))? 0 : SEQALIGN_SCORE_MIN;
		}
		row_mov(us, es, qs, ubegs, W, mov, piecewise); // mov and swap
		hh = row_cal(rbeg, tbase, us, es, qs, bs, ubegs, qprof, gapo1, gape1, gapo2, gape2, W, mov, ph, rh);
		if(verbose > 3){
			row_verify(rbeg, W, ph, rh, hh, tbase, us, es, qs, qprof, gapo1, gape1, gapo2, gape2);
		}
		if(mode == SEQALIGN_MODE_GLOBAL && verbose == 0){
			rmax = 0;
			lst_score = SEQALIGN_SCORE_MIN;
		} else {
			rmax = row_max(us[1], ubegs, W, &lst_score);
		}
		if(verbose){
			banded_striped_epi8_seqalign_row_print(stdout, qseq, 1, i, tbase, bandwidth, mov, rbeg, rmax, lst_score, ubegs, us[1], bs, verbose);
		}
		bs += bandwidth;
		// adaptive banded
		rbx = band_mov(us[1], ubegs, W, i, rbeg, qlen);
		if(mode == SEQALIGN_MODE_GLOBAL){
			rbz = 2 * num_max(Int(tlen / qlen), 1); // suggested max step
			rby = Int((1.0 * i / tlen) * qlen); // diagonal line
			if(rbeg + rbz * (tlen - i - 1) + bandwidth <= (qlen + rbz - 1)){ // be quick to move to end
				mov = 1 + ((qlen - (rbeg + bandwidth)) / num_max(1, tlen - i - 1));
			} else if(Int(rbeg) < rby - Int(bandwidth)){
				mov = rbx + 1;
			} else if(Int(rbeg) > rby){
				mov = num_max(0, rbx - 1);
			} else {
				mov = rbx;
			}
		} else {
			mov = rbx;
		}
		begs[i] = rbeg;
		if(mode != SEQALIGN_MODE_GLOBAL){
			if(rbeg + bandwidth >= qlen){
				score = getscore(us[1], ubegs, W, qlen - 1 - rbeg);
				if(score > rs.score){
					rs.score = score;
					rs.qe = qlen - 1;
					rs.te = i;
				}
			}
		}
	}
	if(mode == SEQALIGN_MODE_GLOBAL){
		rs.score = getscore(us[1], ubegs, W, qlen - 1 - rbeg);
		rs.qe = qlen - 1;
		rs.te = tlen - 1;
	} else if(lst_score > rs.score){
		rs.score = lst_score;
		rs.qe = rbeg + rmax;
		rs.te = tlen - 1;
	}
	// backtrace
	backtrace(qseq, tseq, bs0, begs, bandwidth, piecewise, &rs, cigars);
	if(mempb) free(mempb);
	return rs;
}

#endif
