/*
 * bsalign.h includes 'pairwise edit alignment' and 'pairwise gap-weighting alignment'
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
 * To use bsalign.h in your program, please copy bsalign.h, list.h, sort.h and mem_share.h together
 *
 */

#ifndef BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H
#define BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H

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
#define SEQALIGN_MODE_EDIT		8

#define SEQALIGN_BT_M	0
#define SEQALIGN_BT_I	1
#define SEQALIGN_BT_D	2

#define SEQALIGN_BT1_IE	4
#define SEQALIGN_BT1_DE	8

#define SEQALIGN_BT2_I1	1
#define SEQALIGN_BT2_D1	2
#define SEQALIGN_BT2_I2	3
#define SEQALIGN_BT2_D2	4
#define SEQALIGN_BT2_IE1	8
#define SEQALIGN_BT2_DE1	16
#define SEQALIGN_BT2_IE2	32
#define SEQALIGN_BT2_DE2	64

#define SEQALIGN_SCORE_EPI8_MIN	(-(MAX_B1 >> 1))
#define SEQALIGN_SCORE_EPI8_MAX	(MAX_B1 >> 1)
#define SEQALIGN_SCORE_MIN	(-(MAX_B4 >> 2))
#define SEQALIGN_SCORE_MAX	(MAX_B4 >> 2)

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
#define mm_setzero	_mm256_setzero_si256
#define mm_set1_epi8	_mm256_set1_epi8
#define mm_set1_epi16	_mm256_set1_epi16
#define mm_set1_epi32	_mm256_set1_epi32
#define mm_set1_epi64x	_mm256_set1_epi64x
#define mm_srli	_mm256_srli_si256
#define mm_slli	_mm256_slli_si256
#define mm_srli_epi32	_mm256_srli_epi32
#define mm_srai_epi32	_mm256_srai_epi32
#define mm_slli_epi32	_mm256_slli_epi32
#define mm_slai_epi32	_mm256_slai_epi32
#define mm_srli_epi64	_mm256_srli_epi64
#define mm_srai_epi64	_mm256_srai_epi64
#define mm_slli_epi64	_mm256_slli_epi64
#define mm_slai_epi64	_mm256_slai_epi64
#define mm_insert_epi8	_mm256_insert_epi8
#define mm_extract_epi16	_mm256_extract_epi16
#define mm_extract_epi32	_mm256_extract_epi32
#define mm_adds_epi8	_mm256_adds_epi8
#define mm_adds_epu8	_mm256_adds_epu8
#define mm_adds_epi16	_mm256_adds_epi16
#define mm_add_epi32	_mm256_add_epi32
#define mm_subs_epi8	_mm256_subs_epi8
#define mm_subs_epi16	_mm256_subs_epi16
#define mm_sub_epi32	_mm256_sub_epi32
#define mm_cmpeq_epi8	_mm256_cmpeq_epi8
#define mm_cmpgt_epi8	_mm256_cmpgt_epi8
#define mm_cmpeq_epi16	_mm256_cmpeq_epi16
#define mm_cmpgt_epi16	_mm256_cmpgt_epi16
#define mm_cmpeq_epi32	_mm256_cmpeq_epi32
#define mm_cmpgt_epi32	_mm256_cmpgt_epi32
#define mm_movemask_epi8	_mm256_movemask_epi8
#define mm_max_epi8	_mm256_max_epi8
#define mm_min_epi8	_mm256_min_epi8
#define mm_max_epi16	_mm256_max_epi16
#define mm_min_epi16	_mm256_min_epi16
#define mm_max_epi32	_mm256_max_epi32
#define mm_min_epi32	_mm256_min_epi32
#define mm_shuffle	_mm256_shuffle_epi8
#define mm_blendv	_mm256_blendv_epi8
#define mm_cvtepi8x0_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(a))
#define mm_cvtepi8x1_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
#define mm_cvtepi16x0_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(a))
#define mm_cvtepi16x1_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
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
#define mm_setzero	_mm_setzero_si128
#define mm_set1_epi8	_mm_set1_epi8
#define mm_set1_epi16	_mm_set1_epi16
#define mm_set1_epi32	_mm_set1_epi32
#define mm_set1_epi64x	_mm_set1_epi64x
#define mm_srli	_mm_srli_si128
#define mm_slli	_mm_slli_si128
#define mm_srli_epi32	_mm_srli_epi32
#define mm_srai_epi32	_mm_srai_epi32
#define mm_slli_epi32	_mm_slli_epi32
#define mm_slai_epi32	_mm_slai_epi32
#define mm_srli_epi64	_mm_srli_epi64
#define mm_srai_epi64	_mm_srai_epi64
#define mm_slli_epi64	_mm_slli_epi64
#define mm_slai_epi64	_mm_slai_epi64
#define mm_insert_epi8	_mm_insert_epi8
#define mm_extract_epi16	_mm_extract_epi16
#define mm_extract_epi32	_mm_extract_epi32
#define mm_adds_epi8	_mm_adds_epi8
#define mm_adds_epu8	_mm_adds_epu8
#define mm_adds_epi16	_mm_adds_epi16
#define mm_add_epi32	_mm_add_epi32
#define mm_subs_epi8	_mm_subs_epi8
#define mm_subs_epi16	_mm_subs_epi16
#define mm_sub_epi32	_mm_sub_epi32
#define mm_cmpeq_epi8	_mm_cmpeq_epi8
#define mm_cmpgt_epi8	_mm_cmpgt_epi8
#define mm_cmpeq_epi16	_mm_cmpeq_epi16
#define mm_cmpgt_epi16	_mm_cmpgt_epi16
#define mm_cmpeq_epi32	_mm_cmpeq_epi32
#define mm_cmpgt_epi32	_mm_cmpgt_epi32
#define mm_movemask_epi8	_mm_movemask_epi8
#define mm_max_epi8	_mm_max_epi8
#define mm_min_epi8	_mm_min_epi8
#define mm_max_epi16	_mm_max_epi16
#define mm_min_epi16	_mm_min_epi16
#define mm_max_epi32	_mm_max_epi32
#define mm_min_epi32	_mm_min_epi32
#define mm_shuffle	_mm_shuffle_epi8
#define mm_blendv	_mm_blendv_epi8
#define mm_cvtepi8x0_epi16(a)	_mm_cvtepi8_epi16(a)
#define mm_cvtepi8x1_epi16(a)	_mm_cvtepi8_epi16(_mm_srli_si128(a, 8))
#define mm_cvtepi16x0_epi32(a)	_mm_cvtepi16_epi32(a)
#define mm_cvtepi16x1_epi32(a)	_mm_cvtepi16_epi32(_mm_srli_si128(a, 8))
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

#define striped_epi2_seqedit_getval(xs, W, pos) (((xs)[((((pos) % (W)) * WORDSIZE) + (((pos) / (W)) >> 3))] >> (((((pos) / ((W)))) & 0x7))) & 0x1)
#define striped_epi2_seqedit_qprof_size(qlen) (roundup_times(qlen, WORDSIZE * 8) / 2)
typedef void (*striped_epi2_seqedit_set_query_prof_func)(u1i *qseq, u4i qlen, b1i *qprof);
typedef void (*striped_epi2_seqedit_row_init_func)(b1i *us[2], u4i W, int mode);
typedef void (*striped_epi2_seqedit_row_cal_func)(u4i rbeg, b1i *us[2][2], b1i *hs, b1i *qprof, u4i W, u1i base, int mode);
typedef seqalign_result_t (*striped_epi2_seqedit_backtrace_func)(b1i *uts[2], u4i W, int mode, u1i *qseq, int qend, u1i *tseq, int tend, u4v *cigars);
static inline seqalign_result_t striped_epi2_seqedit_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int mode, int verbose);

// only support short query (<=128*127=16256)
// edit-overlap mode
typedef void (*striped_epi2_seqedit_row_merge_func)(b1i *us[3][2], u2i *hs[2], u1i W);

/**
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
typedef void (*banded_striped_epi8_seqalign_piecex_row_init_func)(b1i *us, b1i *es, b1i *qs, int *ubegs, int *rbeg, int mode, u4i bandwidth, b1i max_nt_score, b1i min_nt_score, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// W = bandwidth / WORDSIZE, WORDSIZE = 16 (SSEx 128) or 32 (AVX2 512)
// converting from [1] to [0]
// make the two row aligned
// ubegs is the absolute scores for the first striped block
typedef void (*banded_striped_epi8_seqalign_piecex_row_mov_func)(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], u4i W, u4i mov, int piecewise, b1i nt_max, b1i nt_min, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// core func to update row scores, from us[0] to us[1]
// rh is the score of H(-1, y - 1), 0, or SCORE_MIN, please see banded_striped_epi8_seqalign_pairwise
// ubegs[i] - ubegs[i-1] is used in F-penetration
// ubegs will be updated
// @return: score of H(-1, y), note that H(-1, y) is useful to restore all scores of row in row_max
typedef int (*banded_striped_epi8_seqalign_piecex_row_btcal_func)(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh, int piecewise);
typedef int (*banded_striped_epi8_seqalign_piecex_row_cal_func)(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh, int piecewise);

// us[0] and us[1] is two row-offset-aligned progenitors
// cmp [0] and [1], merges into [0], records original in os[u, e, q], {us, es, qs, ubegs}
// os stores array of the origins, 0 or 1, respectively for u, e, q
// this function is used to merge multiple progenitors into current node in graph based alignment, program should call row_mov(i)+row_cal(i) before row_merge(1 .. n in pairwise)
typedef void (*banded_striped_epi8_seqalign_piecex_row_merge_func)(b1i *us[3], b1i *es[3], b1i *qs[3], int *ubegs[3], b1i *os[3], b1i osmarks[2], u4i W, int piecewise);

typedef int (*banded_striped_epi8_seqalign_piecex_row_verify_func)(int rowidx, int rowoff, int mode, int W, int mov, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

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
typedef u4i (*banded_striped_epi8_seqalign_piecex_backtrace_func)(u1i *qseq, u1i *tseq, b1i *bs, int *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, u4v *cigars);

// backcal
// restores the best path by revise calculating
typedef u4i(* banded_striped_epi8_seqalign_piecex_backcal_func)(u1i *qseq, u1i *tseq, b1i *ups, b1i *eps, b1i *qps, b1i *ubs, int *roffs, int mode, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, seqalign_result_t *rs, u4v *cigars);


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
						if(0){
							if(striped_epi2_seqedit_getval(qprof + i * W * WORDSIZE, W, pos) != 0x1){
								fprintf(stderr, " -- something wrong i=%d j=%d k=%d x=%d pos=%d in %s -- %s:%d --\n", i, j, k, x, pos, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
								abort();
							}
						}
					}
				}
				qp ++;
			}
		}
	}
}

static inline void striped_epi2_seqedit_row_init(b1i *us[2], u4i W, int mode){
	xint BIT0, BIT1;
	u4i i;
	UNUSED(mode);
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

static u1i *_DEBUG_QSEQ = NULL;
static u4i  _DEBUG_QLEN = 0;
static int _DEBUG_SCORE = 0;

static inline void print_epi2_row(u1i tbase, b1i *us[2], u4i W){
	u4i j, b1, b2;
	int score;
	score = _DEBUG_SCORE;
	fprintf(stdout, "[%c]\t", "ACGTN"[tbase]);
	for(j=0;j<W*WORDSIZE*8;j++){
		b1 = striped_epi2_seqedit_getval(us[0], W, j);
		b2 = striped_epi2_seqedit_getval(us[1], W, j);
		if(b1 == 0 && b2 == 1) score ++;
		else if(b1 == 1 && b2 == 0) score --;
		fprintf(stdout, "%3d%c%03d:%d%d ", j, j < _DEBUG_QLEN? "ACGTN"[_DEBUG_QSEQ[j]] : '*', score, b1, b2);
	}
	fprintf(stdout, "\n");
}

		//    s       u       v   =    h      v'
		// -1(10)  -1(10)  -1(10) = -1(10)  0(00)
		// -1(10)  -1(10)   0(00) = -1(10)  0(00)
		// -1(10)  -1(10)   1(01) =  0(00)  1(01)
		// -1(10)  -1(10)   2(11) =  1(01)  2(11)
		// -1(10)   0(00)  -1(10) = -1(10) -1(10)
		// -1(10)   0(00)   0(00) = -1(10) -1(10)
		// -1(10)   0(00)   1(01) =  0(00)
		// -1(10)   0(00)   2(11) =  1(01)
		// -1(10)   1(01)  -1(10) =  0(00)
		// -1(10)   1(01)   0(00) =  0(00)
		// -1(10)   1(01)   1(01) =  0(00)
		// -1(10)   1(01)   2(11) =  1(01)
		// -1(10)   2(11)  -1(10) =  1(01)
		// -1(10)   2(11)   0(00) =  1(01)
		// -1(10)   2(11)   1(01) =  1(01)
		// -1(10)   2(11)   2(11) =  1(01)
		//  1(00)  -1(10)  -1(10) =  0(00)
		//  1(00)  -1(10)   0(00) =  0(00)
		//  1(00)  -1(10)   1(01) =  0(00)
		//  1(00)  -1(10)   2(11) =  1(01)
		//  1(00)   0(00)  -1(10) =  0(00)
		//  1(00)   0(00)   0(00) =  1(01)
		//  1(00)   0(00)   1(01) =  1(01)
		//  1(00)   0(00)   2(11) =  1(01)
		//  1(00)   1(01)  -1(10) =  0(00)
		//  1(00)   1(01)   0(00) =  1(01)
		//  1(00)   1(01)   1(01) =  1(01)
		//  1(00)   1(01)   2(11) =  1(01)
		//  1(00)   2(11)  -1(10) =  1(01)
		//  1(00)   2(11)   0(00) =  1(01)
		//  1(00)   2(11)   1(01) =  1(01)
		//  1(00)   2(11)   2(11) =  1(01)
		//    ^^      ^^      ^^       ^^
		//    ab      cd      ef       xy




// global DNA edit distance
// rbeg == 0
// mode is useless
static void striped_epi2_seqedit_row_cal(u4i rbeg, b1i *us[2][2], b1i *hs, b1i *qprof, u4i W, u1i base, int mode){
	xint s, u1, u2, u3, u4, v1, v2, h, h2, BIT0, BIT1, BMSK;
	u4i i, running;
	UNUSED(rbeg);
	UNUSED(mode);
	BIT0 = mm_set1_epi8(0x00);
	BIT1 = mm_set1_epi8(0xFF);
	BMSK = mm_srli(mm_set1_epi8(0x1), WORDSIZE - 1);
	v1 = BIT0;
	v2 = BIT1;
	for(i=0;i<W;i++){
		// set s = 0/1 = 10/00
		//s  = mm_load(((xint*)qprof) + (rbeg + i) * 4 + base);
		s  = mm_load(((xint*)qprof) + base * W + i);
		// u = H(x, y - 1) - H(x - 1, y - 1) = -1/0/1 = 10/00/01
		// v = H(x - 1, y) - H(x - 1, y - 1) = 10/00/01
		u1  = mm_load(((xint*)us[0][0]) + i);
		u2  = mm_load(((xint*)us[0][1]) + i);
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
		if(0){
			fprintf(stdout, "i  = %d\n", i);
			fprintf(stdout, "s  = "); fprint_striped_epi1_word(stdout, s); fprintf(stdout, "\n");
			fprintf(stdout, "u1 = "); fprint_striped_epi1_word(stdout, u1); fprintf(stdout, "\n");
			fprintf(stdout, "u2 = "); fprint_striped_epi1_word(stdout, u2); fprintf(stdout, "\n");
			fprintf(stdout, "v1 = "); fprint_striped_epi1_word(stdout, v1); fprintf(stdout, "\n");
			fprintf(stdout, "v2 = "); fprint_striped_epi1_word(stdout, v2); fprintf(stdout, "\n");
			fprintf(stdout, "h  = "); fprint_striped_epi1_word(stdout,  h); fprintf(stdout, "\n");
			fprintf(stdout, "u3 = "); fprint_striped_epi1_word(stdout, u3); fprintf(stdout, "\n");
			fprintf(stdout, "u4 = "); fprint_striped_epi1_word(stdout, u4); fprintf(stdout, "\n");
		}
		v1 = mm_andnot(h, u2);
		v2 = mm_xor(u2, mm_or(h, mm_or(u1, u2)));
		mm_store(((xint*)hs) + i, h);
		mm_store(((xint*)us[1][0]) + i, u3);
		mm_store(((xint*)us[1][1]) + i, u4);
		//print_epi2_row(base, us[1], W);
	}
	running = 1;
	while(running){
		v1 = mm_sl1bit(v1);
		v2 = mm_sl1bit(v2);
		v2 = mm_or(v2, BMSK);
		for(i=0;i<W;i++){
			//s   = mm_load(((xint*)qprof) + (rbeg + i) * 4 + base);
			s   = mm_load(((xint*)qprof) + base * W + i);
			h2  = mm_load(((xint*)hs) + i);
			u1  = mm_load(((xint*)us[0][0]) + i);
			u2  = mm_load(((xint*)us[0][1]) + i);
			h   = mm_andnot(mm_or(s, mm_or(u1, v1)), BIT1);
			u3 = mm_andnot(h, v2);
			u4 = mm_xor(v2, mm_or(h, mm_or(v1, v2)));
			v1 = mm_andnot(h, u2);
			v2 = mm_xor(u2, mm_or(h, mm_or(u1, u2)));
			mm_store(((xint*)hs) + i, h);
			mm_store(((xint*)us[1][0]) + i, u3);
			mm_store(((xint*)us[1][1]) + i, u4);
			//print_epi2_row(base, us[1], W);
			//when there is nothing to update, break
			if(!mm_movemask_epi8(mm_andnot(mm_cmpeq_epi8(h2, h), BIT1))){
				running = 0;
				break;
			}
		}
	}
}

static inline seqalign_result_t striped_epi2_seqedit_backtrace(b1i *uts[2], u4i W, int mode, u1i *qseq, int x, u1i *tseq, int y, u4v *cigars){
	seqalign_result_t rs;
	u4i cg, op;
	int u1, u2, u3, u4;
	UNUSED(mode);
	ZEROS(&rs);
	rs.qe = x + 1;
	rs.te = y + 1;
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
	if(rs.qb){
		op = 1;
		if(op == (cg & 0xf)){
			cg += 0x10 * rs.qb;
		} else {
			if(cg && cigars) push_u4v(cigars, cg);
			cg = (0x10 * rs.qb) | op;
		}
		rs.ins += rs.qb;
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
		rs.del += rs.tb;
		rs.tb = 0;
	}
	rs.aln = rs.mat + rs.mis + rs.ins + rs.del;
	if(cg && cigars) push_u4v(cigars, cg);
	if(cigars) reverse_u4v(cigars);
	return rs;
}

// W <= 63
// the size of hs[0/1/2] is WORDSIZE * 2
// NOTE THAT THIS FUNCTION IS UN-FINISHED
static inline void striped_epi2_seqedit_row_merge(u2i sbegs[2], b1i *us[3][2], b1i *hs[3], u1i W){
	xint p1, p2, q1, q2, u1, u2, h1, h2, h3, h4;
	u4i i;
	b2i ts[3];
	u1i rs[4][WORDSIZE];
//#define CHECK_NO_SSE
#ifdef CHECK_NO_SSE
#define DEBUG_NO_SSE
#endif
#ifdef CHECK_NO_SSE
	int ds[2][128];
	memset(ds[0], 0, 128 * sizeof(int));
	memset(ds[1], 0, 128 * sizeof(int));
#endif
	if(W > 63){
		W = 63;
	}
	void striped_seqedit_sum_4bits(u1i *dst, b1i *src){
		u8i *s, t[WORDSIZE/8], *r, v;
		u4i i, k, l, j;
		memset(dst, 0, WORDSIZE);
		for(k=0;k<W;k+=3){ // 0xF/count_ones(0xF)=15/4=3
			l = num_min(k + 3, W);
			memset(t, 0, WORDSIZE);
			for(i=k;i<l;i++){
				s = (u8i*)(src + i * WORDSIZE);
				for(j=0;j<WORDSIZE/8;j++){ // xint to u8i
					// count ones per 4 bits
					v = s[j] - ((s[j] >> 1) & 0x5555555555555555LLU);
					v = (v & 0x3333333333333333LLU) + ((v >> 2) & 0x3333333333333333LLU);
					t[j] += v;
				}
			}
			for(j=0;j<WORDSIZE/8;j++){
				r = ((u8i*)dst) + j * 2;
				v = ((t[j] << 16) & 0x0000FFFF00000000LLU) | (t[j] & 0x000000000000FFFFLLU);
				v = ((v << 8) | v) & 0x00FF00FF00FF00FFLLU;
				v = ((v << 4) | v) & 0x0F0F0F0F0F0F0F0FLLU;
				r[0] += v;
				v = ((t[j] >> 16) & 0x0000FFFF00000000LLU) | ((t[j] >> 32) & 0x000000000000FFFFLLU);
				v = ((v << 8) | v) & 0x00FF00FF00FF00FFLLU;
				v = ((v << 4) | v) & 0x0F0F0F0F0F0F0F0FLLU;
				r[1] += v;
			}
		}
	}
	striped_seqedit_sum_4bits(rs[0], us[0][0]);
	striped_seqedit_sum_4bits(rs[1], us[0][1]);
	striped_seqedit_sum_4bits(rs[2], us[1][0]);
	striped_seqedit_sum_4bits(rs[3], us[1][1]);
	ts[0] = sbegs[0];
	ts[1] = sbegs[1];
	for(i=0;i<WORDSIZE*2;i++){
		ts[2] = (ts[0] + ts[1]) / 2;
		ts[0] -= ts[2];
		ts[1] -= ts[2];
		if(ts[0] <= ts[1]){
			hs[0][i] = ts[0] < -7? -7 : ts[0];
			hs[1][i] = ts[1] >  7?  7 : ts[1];
		} else {
			hs[0][i] = ts[0] >  7?  7 : ts[0];
			hs[1][i] = ts[1] < -7? -7 : ts[1];
		}
		ts[0] += rs[0][i];
		ts[1] += rs[1][i];
	}
	{
		h1 = mm_load(((xint*)hs[0]));
		h2 = mm_load(((xint*)hs[1]));
		h4 = mm_min_epi16(h1, h2);
		mm_store(((xint*)hs[2]), h4);
	}
#ifdef CHECK_NO_SSE
	void sub_print_3hs(int tag){
		UNUSED(tag);
#ifdef DEBUG_NO_SSE
		fprintf(stdout, "--- HS %d ---\n", tag);
		fprintf(stdout, "hs1");
		for(j=0;j<WORDSIZE*8;j++){
			k = (((j) & 0xf) << (WORDSHIFT - 1)) + ((j) >> 4);
			fprintf(stdout, "\t%d", hs[0][k]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, "hs2");
		for(j=0;j<WORDSIZE*8;j++){
			k = (((j) & 0xf) << (WORDSHIFT - 1)) + ((j) >> 4);
			fprintf(stdout, "\t%d", hs[1][k]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, "hs3");
		for(j=0;j<WORDSIZE*8;j++){
			k = (((j) & 0xf) << (WORDSHIFT - 1)) + ((j) >> 4);
			fprintf(stdout, "\t%d", hs[2][k]);
		}
		fprintf(stdout, "\n");
#endif
		for(j=0;j<WORDSIZE*8;j++){
			k = (((j) & 0xf) << (WORDSHIFT - 1)) + ((j) >> 4);
			if(hs[2][k] != num_min(hs[0][k], hs[1][k])){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	sub_print_3hs(-1);
#endif
	__attribute__((unused)) void striped_epi2_seqedit_row_merge_core_2(const int shift){
		xint MASK, ONES, A, B, S, C, D;
		MASK = mm_set1_epi8(0x11);
		ONES = mm_set1_epi8(0x0F);
		C    = mm_set1_epi8(0x10);
		D    = mm_set1_epi8(0x01);
		for(i=0;i<W;i++){ // TODO: blocks of 7
			q1 = mm_load(((xint*)us[0][0]) + i);
			p1 = mm_load(((xint*)us[0][1]) + i);
			h1 = mm_adds_epi8(h1, mm_and(mm_srli_epi64(p1, shift), MASK));
			h1 = mm_subs_epi8(h1, mm_and(mm_srli_epi64(q1, shift), MASK));
			q2 = mm_load(((xint*)us[1][0]) + i);
			p2 = mm_load(((xint*)us[1][1]) + i);
			h2 = mm_adds_epi8(h2, mm_and(mm_srli_epi64(p2, shift), MASK));
			h2 = mm_subs_epi8(h2, mm_and(mm_srli_epi64(q2, shift), MASK));
			u1 = mm_load(((xint*)us[2][0]) + i);
			u2 = mm_load(((xint*)us[2][1]) + i);
			// min epi4
			A = mm_and(h1, ONES);
			B = mm_and(h2, ONES);
			S = mm_and(mm_adds_epi8(mm_andnot(A, ONES), B), C); // ((~A) + B) & 010
			S = mm_srli_epi64(mm_subs_epi8(S, D), 4); // (S - 1) >> 4
			S = mm_adds_epi8(mm_andnot(S, A), mm_and(S, B)); // ((~S) & A) + (S & B) -> min_epi4(A, B)
			// (S - h3) -> {-1, 0, 1} -> {0xFF, 0x00, 0x01}
			A = mm_subs_epi8(S, h3);
			u1 = mm_or(u1, mm_slli_epi64(mm_and(mm_srli_epi64(A, 1), D), shift)); // bit2
			u2 = mm_or(u2, mm_slli_epi64(mm_and(mm_xor(mm_srli_epi64(A, 2), A), D), shift)); // bit1 ^ bit3
			h3 = S;

			A = mm_and(mm_srli_epi64(h1, 4), ONES);
			B = mm_and(mm_srli_epi64(h2, 4), ONES);
			S = mm_and(mm_adds_epi8(mm_andnot(A, ONES), B), C); // ((~A) + B) & 010
			S = mm_srli_epi64(mm_subs_epi8(S, D), 4); // (S - 1) >> 4 -> {0x00, 0x0F}
			S = mm_adds_epi8(mm_andnot(S, A), mm_and(S, B)); // ((~S) & A) + (S & B) -> min_epi4(A, B)
			// (S - h3) -> {-1, 0, 1} -> {0xFF, 0x00, 0x01}
			A = mm_subs_epi8(S, h4);
			u1 = mm_or(u1, mm_slli_epi64(mm_and(mm_srli_epi64(A, 1), D), shift + 4)); // bit2
			u2 = mm_or(u2, mm_slli_epi64(mm_and(mm_xor(mm_srli_epi64(A, 2), A), D), shift + 4)); // bit1 ^ bit3
			h4 = S;

			mm_store(((xint*)us[2][0]) + i, u1);
			mm_store(((xint*)us[2][1]) + i, u2);
		}
	}
	void striped_epi2_seqedit_row_merge_core(const int shift){
		xint MASK, SHFT;
		MASK = mm_set1_epi8(1);
		SHFT = mm_set1_epi8((short)(1 << shift));
		for(i=0;i<W;i++){
			q1 = mm_load(((xint*)us[0][0]) + i);
			p1 = mm_load(((xint*)us[0][1]) + i);
			h1 = mm_adds_epi8(h1, mm_and(mm_srli_epi64(p1, shift), MASK));
			h1 = mm_subs_epi8(h1, mm_and(mm_srli_epi64(q1, shift), MASK));
			q2 = mm_load(((xint*)us[1][0]) + i);
			p2 = mm_load(((xint*)us[1][1]) + i);
			h2 = mm_adds_epi8(h2, mm_and(mm_srli_epi64(p2, shift), MASK));
			h2 = mm_subs_epi8(h2, mm_and(mm_srli_epi64(q2, shift), MASK));
			h4 = mm_min_epi8(h1, h2);
			u1 = mm_load(((xint*)us[2][0]) + i);
			u2 = mm_load(((xint*)us[2][1]) + i);
			u1 = mm_or(u1, mm_and(mm_cmpgt_epi8(h3, h4), SHFT));
			u2 = mm_or(u2, mm_and(mm_cmpgt_epi8(h4, h3), SHFT));
			h3 = h4;
			mm_store(((xint*)us[2][0]) + i, u1);
			mm_store(((xint*)us[2][1]) + i, u2);
		}
	}
	striped_epi2_seqedit_row_merge_core(0);
	striped_epi2_seqedit_row_merge_core(1);
	striped_epi2_seqedit_row_merge_core(2);
	striped_epi2_seqedit_row_merge_core(3);
	striped_epi2_seqedit_row_merge_core(4);
	striped_epi2_seqedit_row_merge_core(5);
	striped_epi2_seqedit_row_merge_core(6);
	striped_epi2_seqedit_row_merge_core(7);
#ifdef CHECK_NO_SSE
	sub_print_3hs(i);
	int score[4], b1, b2;
	score[0] = sbegs[0];
	score[1] = sbegs[1];
	score[2] = num_min(score[0], score[1]);
	for(i=0;i<W*WORDSIZE;i++){
		b1 = striped_epi2_seqedit_getval(us[0][0], W, i);
		b2 = striped_epi2_seqedit_getval(us[0][1], W, i);
		score[0] = score[0] + b2 - b1;
		b1 = striped_epi2_seqedit_getval(us[1][0], W, i);
		b2 = striped_epi2_seqedit_getval(us[1][1], W, i);
		score[1] = score[1] + b2 - b1;
		b1 = striped_epi2_seqedit_getval(us[2][0], W, i);
		b2 = striped_epi2_seqedit_getval(us[2][1], W, i);
		score[2] = score[2] + b2 - b1;
		score[3] = num_min(score[0], score[1]);
		if(score[2] != score[3]){
			fflush(stdout); fprintf(stderr, " -- [%d] [%d,%d] %d != %d in %s -- %s:%d --\n", i, score[0], score[1], score[2], score[3], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
	}
#endif
}

static inline seqalign_result_t striped_epi2_seqedit_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int mode, int verbose){
	striped_epi2_seqedit_set_query_prof_func set_qprof;
	striped_epi2_seqedit_row_init_func       row_init;
	striped_epi2_seqedit_row_cal_func        row_cal;
	striped_epi2_seqedit_backtrace_func      backtrace;
	seqalign_result_t rs;
	b1i *memp, *mempb, *qprof, *uts[2], *us[2][2], *hs;
	u8i mpsize;
	u4i i, W;
	int rx, ry, sbeg, smin;
	set_qprof = striped_epi2_seqedit_set_query_prof;
	row_cal   = striped_epi2_seqedit_row_cal;
	row_init  = striped_epi2_seqedit_row_init;
	backtrace = striped_epi2_seqedit_backtrace;
	UNUSED(mode);
	W = roundup_times(qlen, WORDSIZE * 8) / (WORDSIZE * 8);
	//fprintf(stdout, "qlen=%d\tW=%d\n", qlen, W);
	mpsize = 0;
	mpsize += striped_epi2_seqedit_qprof_size(qlen); // qprof[]
	mpsize += 2 * W * WORDSIZE * ((tlen + 1)); // uts[]
	mpsize += W * WORDSIZE; // hs[]
	//mpsize += 2 * W * WORDSIZE + 3 * 2 * 8 * WORDSIZE; // TODO: remove it after debug
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
	_DEBUG_QSEQ = qseq;
	_DEBUG_QLEN = qlen;
	rx   = qlen - 1;
	ry   = tlen - 1;
	smin = MAX_B4;
	row_init(us[0], W, mode);
	for(i=0;i<tlen;i++){
		sbeg = i + 1;
		_DEBUG_SCORE = sbeg;
		us[0][0] = uts[0] + (i + 0) * W * WORDSIZE;
		us[0][1] = uts[1] + (i + 0) * W * WORDSIZE;
		us[1][0] = uts[0] + (i + 1) * W * WORDSIZE;
		us[1][1] = uts[1] + (i + 1) * W * WORDSIZE;
		row_cal(0, us, hs, qprof, W, tseq[i], mode);
		if(verbose){
			int vals[2][2] = {{0, 1}, {-1, 2}};
			int j, b1, b2, u, u2, v, v2, score, error;
			score = i + 1;
			v2    = 1;
			error = 0;
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
				if(error == 0 && u2 != vals[b1][b2]){
					fflush(stdout);
					fprintf(stderr, " -- j=%d s=%c%c u=%d v=%d u2=%d should be %d in %s -- %s:%d --\n", j, "ACGTN"[qseq[j]], "ACGTN"[tseq[i]], u, v, vals[b1][b2], u2, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					error = 1;
				}
				if(b1 == 0 && b2 == 1) score ++;
				else if(b1 == 1 && b2 == 0) score --;
				fprintf(stdout, "%c%03d:%c:%c ", "ACGTN"[qseq[j]], score, "-*+"[vals[b1][b2] + 1], "-*+"[v2 + 1]);
			}
			fprintf(stdout, "\n");
		}
	}
	rs = backtrace(uts, W, mode, qseq, rx, tseq, ry, cigars);
	if(mempb) free(mempb);
	return rs;
}

static inline int banded_striped_epi8_seqalign_get_piecewise(b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, int bandwidth){
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1 && (gapo1 - gapo2) / (gape1 - gape2) < bandwidth){
		return 2;
	} else if(gapo1){
		return 1;
	} else {
		return 0;
	}
}

static void banded_striped_epi8_seqalign_piecex_row_init(b1i *us, b1i *es, b1i *qs, int *ubegs, int *rbeg, int mode, u4i bandwidth, b1i max_nt, b1i min_nt, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
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
			for(k=0;k<W;k++) mm_store(((xint*)us) + k, GAP);
			for(k=0;k<WORDSIZE;k++) ubegs[k] = gape2 * W;
			us[0] = gapo1 + gape1 + min_nt - max_nt;
			ubegs[0] += us[0] - gape2;
			for(k=1;k<xp;k++){
				us[banded_striped_epi8_pos2idx(bandwidth, k)] = gape1;
				ubegs[k / W] += gape1 - gape2;
			}
		} else {
			GAP = mm_set1_epi8(gape1);
			for(k=0;k<W;k++) mm_store(((xint*)us) + k, GAP);
			us[0] = gapo1 + gape1 + min_nt - max_nt;
			for(k=0;k<WORDSIZE;k++) ubegs[k] = gape1 * W;
			ubegs[0] += us[0] - gape1;
		}
		s = max_nt - min_nt;
		for(k=0;k<WORDSIZE;k++){
			t = ubegs[k];
			ubegs[k] = s;
			s += t;
		}
		ubegs[k] = s;
		if(rbeg) *rbeg = 0;
	} else {
		for(k=0;k<W;k++) mm_store(((xint*)us) + k, ZERO);
		if(rbeg) *rbeg = 0;
		memset(ubegs, 0, (WORDSIZE + 1) * sizeof(int));
	}
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1 && (gapo1 - gapo2) / (gape1 - gape2) < Int(bandwidth)){
		for(k=0;k<W;k++) mm_store(((xint*)es) + k, MIN);
		for(k=0;k<W;k++) mm_store(((xint*)qs) + k, MIN);
	} else if(gapo1){
		for(k=0;k<W;k++) mm_store(((xint*)es) + k, MIN);
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

#if 0

static inline void banded_striped_epi8_seqalign_piecex_row_check_ubegs(b1i *us, int *ubegs, int W){
	int i, j, sc;
	sc = ubegs[0];
	for(j=0;j<WORDSIZE;j++){
		for(i=0;i<W;i++) sc += us[i*WORDSIZE + j];
		if(sc != ubegs[j+1]){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

#else

#define banded_striped_epi8_seqalign_piecex_row_check_ubegs(us, ubegs, W) 

#endif

// movx can be any value
static inline void banded_striped_epi8_seqalign_piecex_row_movx(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], u4i W, u4i movx, int piecewise, b1i nt_max, b1i nt_min, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	xint u, x, UBS[4], SHUFF1, SHUFF2;
	u4i i, div, mov, cyc;
	u1i shuff[WORDSIZE];
#ifdef __AVX2__
	fflush(stdout); fprintf(stderr, " -- Cannot work with AVX2 in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	abort();
#endif
	banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[1], ubegs[1], W);
	if(movx >= W * WORDSIZE){
		memset(us[0], 0, W * WORDSIZE);
		if(piecewise) memset(es[0], 0, W * WORDSIZE);
		if(piecewise == 2) memset(qs[0], 0, W * WORDSIZE);
		memset(us[0], 0, W * WORDSIZE);
		for(i=0;i<=WORDSIZE;i++) ubegs[0][i] = SEQALIGN_SCORE_MIN;
		return;
	} else if(movx == 0){
		memcpy(us[0], us[1], W * WORDSIZE);
		if(piecewise){
			memcpy(es[0], es[1], W * WORDSIZE);
		}
		if(piecewise == 2){
			memcpy(qs[0], qs[1], W * WORDSIZE);
		}
		memcpy(ubegs[0], ubegs[1], (WORDSIZE + 1) * sizeof(int));
		return;
	}
	cyc = movx / W;
	for(i=0;i+cyc<WORDSIZE;i++) shuff[i] = i + cyc;
	for(;i<WORDSIZE;i++) shuff[i] = 0xFF;
	SHUFF1 = mm_load((xint*)shuff);
	SHUFF2 = mm_srli(SHUFF1, 1);
	SHUFF2 = mm_insert_epi8(SHUFF2, 0xFF, WORDSIZE - 1);
	cyc = movx / W;
	mov = movx % W;
	div = W - mov;
	//memcpy(ubegs[0], ubegs[1] + cyc, (WORDSIZE + 1 - cyc) * sizeof(int));
	if(cyc){
		for(i=0;i<div;i++){
			u = mm_load(((xint*)us[1]) + i + mov);
			u = mm_shuffle(u, SHUFF1); // shift cyc epi8
			mm_store(((xint*)us[0]) + i, u);
		}
		if(piecewise){
			for(i=0;i<div;i++){
				x = mm_load(((xint*)es[1]) + i + mov);
				x = mm_shuffle(x, SHUFF1);
				mm_store(((xint*)es[0]) + i, x);
			}
		}
		if(piecewise == 2){
			for(i=0;i<div;i++){
				x = mm_load(((xint*)qs[1]) + i + mov);
				x = mm_shuffle(x, SHUFF1);
				mm_store(((xint*)qs[0]) + i, x);
			}
		}
	} else {
		memcpy(us[0], us[1] + mov * WORDSIZE, div * WORDSIZE);
		if(piecewise){
			memcpy(es[0], es[1] + mov * WORDSIZE, div * WORDSIZE);
		}
		if(piecewise == 2){
			memcpy(qs[0], qs[1] + mov * WORDSIZE, div * WORDSIZE);
		}
	}
	if(mov){
		{
			UBS[0] = mm_load(((xint*)ubegs[1]) + 0);
			UBS[1] = mm_load(((xint*)ubegs[1]) + 1);
			UBS[2] = mm_load(((xint*)ubegs[1]) + 2);
			UBS[3] = mm_load(((xint*)ubegs[1]) + 3);
		}
		for(i=div;i<W;i++){
			u = mm_load(((xint*)us[1]) + i - div);
			UBS[0] = mm_add_epi32(UBS[0], mm_cvtepi8x0_epi32(u));
			UBS[1] = mm_add_epi32(UBS[1], mm_cvtepi8x1_epi32(u));
			UBS[2] = mm_add_epi32(UBS[2], mm_cvtepi8x2_epi32(u));
			UBS[3] = mm_add_epi32(UBS[3], mm_cvtepi8x3_epi32(u));
			u = mm_shuffle(u, SHUFF2); // shift cyc+1 epi8
			mm_store(((xint*)us[0]) + i, u);
		}
		{
			mm_store(((xint*)ubegs[0]) + 0, UBS[0]);
			mm_store(((xint*)ubegs[0]) + 1, UBS[1]);
			mm_store(((xint*)ubegs[0]) + 2, UBS[2]);
			mm_store(((xint*)ubegs[0]) + 3, UBS[3]);
		}
		if(piecewise){
			for(i=div;i<W;i++){
				x = mm_load(((xint*)es[1]) + i - div);
				x = mm_shuffle(x, SHUFF2);
				mm_store(((xint*)es[0]) + i, x);
			}
		}
		if(piecewise == 2){
			for(i=div;i<W;i++){
				x = mm_load(((xint*)qs[1]) + i - div);
				x = mm_shuffle(x, SHUFF2);
				mm_store(((xint*)qs[0]) + i, x);
			}
		}
	}
	{
		u4i a, a2, b, b2, d;
		int c;
		if(mov){
			memmove(ubegs[0], ubegs[0] + cyc, (WORDSIZE - cyc) * sizeof(int));
		} else {
			memmove(ubegs[0], ubegs[1] + cyc, (WORDSIZE - cyc) * sizeof(int));
		}
		for(i=WORDSIZE-cyc;i<=WORDSIZE;i++) ubegs[0][i] = ubegs[1][WORDSIZE];
		banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[0], ubegs[0], W);
		// mimic insertions, make sure 1) the max score in next row never come from this overhang, 2) the difference value with next row never overflow EPI8
		if(piecewise == 2) d = (gapo1 - gapo2) / (gape2 - gape1);
		else d = W * WORDSIZE + 1;
		i = W * WORDSIZE - movx;
		a  = i % W;
		a2 = (i + d) % W;
		b  = i / W;
		b2 = (i + d) / W;
		if(piecewise == 2){
			c = num_min(nt_min, gapo2 + gape2) - 1 - nt_max + (gapo2 + gape2);
		} else {
			c = num_min(nt_min, gapo1 + gape1) - 1 - nt_max + (gapo1 + gape1);
		}
		us[0][(i % W) * WORDSIZE + (i / W)] = c;
		a ++;
		for(;b<WORDSIZE&&b<=b2;b++){
			if(b == b2){
				c += (a2 - a) * gape1;
				for(;a<a2;a++) us[0][a * WORDSIZE + b] = gape1;
				a = a2;
				if(a2 < W) break;
			}
			c += (W - a) * gape1;
			for(;a<W;a++) us[0][a * WORDSIZE + b] = gape1;
			ubegs[0][b + 1] += c;
			a = 0;
		}
		for(;b<WORDSIZE;b++){
			c += (W - a) * gape2;
			for(;a<W;a++) us[0][a * WORDSIZE + b] = gape2;
			ubegs[0][b + 1] += c;
			a = 0;
		}
		banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[0], ubegs[0], W);
	}
}

static inline void banded_striped_epi8_seqalign_piecex_row_mov(b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], u4i W, u4i mov, int piecewise, b1i nt_min, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	xint u, x, UBS[4];
	u4i i, div;
	UNUSED(nt_min);
	UNUSED(gapo1);
	UNUSED(gape1);
	UNUSED(gapo2);
	UNUSED(gape2);
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
	if(!mov){
		memcpy(ubegs[0], ubegs[1], (WORDSIZE + 1) * sizeof(int));
		return;
	}
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
		UBS[0] = mm_load(((xint*)ubegs[1]) + 0);
		UBS[1] = mm_load(((xint*)ubegs[1]) + 1);
		UBS[2] = mm_load(((xint*)ubegs[1]) + 2);
		UBS[3] = mm_load(((xint*)ubegs[1]) + 3);
	}
	for(i=0;i<mov;i++){
		u = mm_load(((xint*)us[1]) + i);
		UBS[0] = mm_add_epi32(UBS[0], mm_cvtepi8x0_epi32(u));
		UBS[1] = mm_add_epi32(UBS[1], mm_cvtepi8x1_epi32(u));
		UBS[2] = mm_add_epi32(UBS[2], mm_cvtepi8x2_epi32(u));
		UBS[3] = mm_add_epi32(UBS[3], mm_cvtepi8x3_epi32(u));
	}
	{
		mm_store(((xint*)ubegs[0]) + 0, UBS[0]);
		mm_store(((xint*)ubegs[0]) + 1, UBS[1]);
		mm_store(((xint*)ubegs[0]) + 2, UBS[2]);
		mm_store(((xint*)ubegs[0]) + 3, UBS[3]);
		ubegs[0][WORDSIZE] = ubegs[1][WORDSIZE] + SEQALIGN_SCORE_EPI8_MIN;
	}
	if(piecewise){
		for(i=div;i<W;i++){
			x = mm_load(((xint*)es[1]) + i + mov - W);
			x = mm_srli(x, 1);
			mm_store(((xint*)es[0]) + i, x);
		}
	}
	if(piecewise == 2){
		for(i=div;i<W;i++){
			x = mm_load(((xint*)qs[1]) + i + mov - W);
			x = mm_srli(x, 1);
			mm_store(((xint*)qs[0]) + i, x);
		}
	}
}

// us[2] can be the same with us[0/1], also es, qs, ubegs
static inline void banded_striped_epi8_seqalign_piecex_row_merge(b1i *us[3], b1i *es[3], b1i *qs[3], int *ubegs[3], u4i W, int piecewise){
	xint s[2][4], t[2][2], m[2][2], u[2], x[2][2];
	u4i i, ie, j, k, p;
	banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[0], ubegs[0], W);
	banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[1], ubegs[1], W);
	for(k=0;k<2;k++){
		for(j=0;j<4;j++){
			s[k][j] = mm_load(((xint*)ubegs[k]) + j);
		}
	}
	p = 0;
	for(j=0;j<4;j++){
		mm_store(((xint*)ubegs[2]) + j, mm_max_epi32(s[0][j], s[1][j]));
	}
	ubegs[2][WORDSIZE] = num_max(ubegs[0][WORDSIZE], ubegs[1][WORDSIZE]);
#if 1
	int rs[3][WORDSIZE];
	memcpy(rs[0], ubegs[0], WORDSIZE * sizeof(int));
	memcpy(rs[1], ubegs[1], WORDSIZE * sizeof(int));
	for(i=0;i<WORDSIZE;i++) rs[2][i] = num_max(rs[0][i], rs[1][i]);
#endif
	for(i=0;i<W;){
		ie = num_min(i + 256, W); // 256 * 127 = 32512 < 0x7FFFF(32767)
		{
			// epi32 -> saturated delta epi16
			u[0] = mm_sub_epi32(s[0][0], s[1][0]);
			u[0] = mm_max_epi32(u[0], mm_set1_epi32(-0x7FFF));
			u[0] = mm_min_epi32(u[0], mm_set1_epi32( 0x7FFF));
			x[0][0] = mm_srai_epi32(u[0], 1); // Please note: -1 srai 1 -> -1
			x[1][0] = mm_sub_epi32(x[0][0], u[0]);
			u[0] = mm_sub_epi32(s[0][1], s[1][1]);
			u[0] = mm_max_epi32(u[0], mm_set1_epi32(-0x7FFF));
			u[0] = mm_min_epi32(u[0], mm_set1_epi32( 0x7FFF));
			x[0][1] = mm_srai_epi32(u[0], 1);
			x[1][1] = mm_sub_epi32(x[0][1], u[0]);
			s[0][0] = mm_sub_epi32(s[0][0], x[0][0]);
			s[0][1] = mm_sub_epi32(s[0][1], x[0][1]);
			s[1][0] = mm_sub_epi32(s[1][0], x[1][0]);
			s[1][1] = mm_sub_epi32(s[1][1], x[1][1]);
			t[0][0] = mm_packs_epi32(x[0][0], x[0][1]);
			t[1][0] = mm_packs_epi32(x[1][0], x[1][1]);
			u[0] = mm_sub_epi32(s[0][2], s[1][2]);
			u[0] = mm_max_epi32(u[0], mm_set1_epi32(-0x7FFF));
			u[0] = mm_min_epi32(u[0], mm_set1_epi32( 0x7FFF));
			x[0][0] = mm_srai_epi32(u[0], 1);
			x[1][0] = mm_sub_epi32(x[0][0], u[0]);
			u[0] = mm_sub_epi32(s[0][3], s[1][3]);
			u[0] = mm_max_epi32(u[0], mm_set1_epi32(-0x7FFF));
			u[0] = mm_min_epi32(u[0], mm_set1_epi32( 0x7FFF));
			x[0][1] = mm_srai_epi32(u[0], 1);
			x[1][1] = mm_sub_epi32(x[0][1], u[0]);
			s[0][2] = mm_sub_epi32(s[0][2], x[0][0]);
			s[0][3] = mm_sub_epi32(s[0][3], x[0][1]);
			s[1][2] = mm_sub_epi32(s[1][2], x[1][0]);
			s[1][3] = mm_sub_epi32(s[1][3], x[1][1]);
			t[0][1] = mm_packs_epi32(x[0][0], x[0][1]);
			t[1][1] = mm_packs_epi32(x[1][0], x[1][1]);
			p = 0;
			m[p][0] = mm_max_epi16(t[0][0], t[1][0]);
			m[p][1] = mm_max_epi16(t[0][1], t[1][1]);
		}
		for(;i<ie;i++){
			p = !p;
			u[0]    = mm_load(((xint*)us[0]) + i);
			t[0][0] = mm_adds_epi16(t[0][0], mm_cvtepi8x0_epi16(u[0]));
			t[0][1] = mm_adds_epi16(t[0][1], mm_cvtepi8x1_epi16(u[0]));
			u[1]    = mm_load(((xint*)us[1]) + i);
			t[1][0] = mm_adds_epi16(t[1][0], mm_cvtepi8x0_epi16(u[1]));
			t[1][1] = mm_adds_epi16(t[1][1], mm_cvtepi8x1_epi16(u[1]));
			//c[0]     = mm_cmpgt_epi16(t[1][0], t[0][0]);
			m[p][0]  = mm_max_epi16(t[0][0], t[1][0]);
			m[!p][0] = mm_subs_epi16(m[p][0], m[!p][0]);
			//c[1]     = mm_cmpgt_epi16(t[1][1], t[0][1]);
			m[p][1]  = mm_max_epi16(t[0][1], t[1][1]);
			m[!p][1] = mm_subs_epi16(m[p][1], m[!p][1]);
			u[0] = mm_packs_epi16(m[!p][0], m[!p][1]);
			mm_store(((xint*)us[2]) + i, u[0]);
#if 1
			for(j=0;j<WORDSIZE;j++){
				rs[0][j] += us[0][i * WORDSIZE + j];
				rs[1][j] += us[1][i * WORDSIZE + j];
				int rm = num_max(rs[0][j], rs[1][j]);
				if(rm - rs[2][j] != us[2][i * WORDSIZE + j]){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				rs[2][j] = rm;
			}
#endif
			//c[0] = mm_packs_epi16(c[0], c[1]);
			//c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
			//mm_store(((xint*)os[0]) + i, c[0]);
			if(piecewise == 0) continue;
			u[0]    = mm_load(((xint*)es[0]) + i);
			x[0][0] = mm_adds_epi16(t[0][0], mm_cvtepi8x0_epi16(u[0]));
			x[0][1] = mm_adds_epi16(t[0][1], mm_cvtepi8x1_epi16(u[0]));
			u[1]    = mm_load(((xint*)es[1]) + i);
			x[1][0] = mm_adds_epi16(t[1][0], mm_cvtepi8x0_epi16(u[1]));
			x[1][1] = mm_adds_epi16(t[1][1], mm_cvtepi8x1_epi16(u[1]));
			//c[0]     = mm_cmpgt_epi16(x[1][0], x[0][0]);
			m[!p][0] = mm_max_epi16(x[0][0], x[1][0]);
			m[!p][0] = mm_subs_epi16(m[!p][0], m[p][0]);
			//c[1]     = mm_cmpgt_epi16(x[1][1], x[0][1]);
			m[!p][1] = mm_max_epi16(x[0][1], x[1][1]);
			m[!p][1] = mm_subs_epi16(m[!p][1], m[p][1]);
			u[0] = mm_packs_epi16(m[!p][0], m[!p][1]);
			mm_store(((xint*)es[2]) + i, u[0]);
			//c[0] = mm_packs_epi16(c[0], c[1]);
			//c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
			//mm_store(((xint*)os[1]) + i, c[0]);
			if(piecewise == 1) continue;
			u[0]    = mm_load(((xint*)qs[0]) + i);
			x[0][0] = mm_adds_epi16(t[0][0], mm_cvtepi8x0_epi16(u[0]));
			x[0][1] = mm_adds_epi16(t[0][1], mm_cvtepi8x1_epi16(u[0]));
			u[1]    = mm_load(((xint*)qs[1]) + i);
			x[1][0] = mm_adds_epi16(t[1][0], mm_cvtepi8x0_epi16(u[1]));
			x[1][1] = mm_adds_epi16(t[1][1], mm_cvtepi8x1_epi16(u[1]));
			//c[0]     = mm_cmpgt_epi16(x[1][0], x[0][0]);
			m[!p][0] = mm_max_epi16(x[0][0], x[1][0]);
			m[!p][0] = mm_subs_epi16(m[!p][0], m[p][0]);
			//c[1]     = mm_cmpgt_epi16(x[1][1], x[0][1]);
			m[!p][1] = mm_max_epi16(x[0][1], x[1][1]);
			m[!p][1] = mm_subs_epi16(m[!p][1], m[p][1]);
			u[0] = mm_packs_epi16(m[!p][0], m[!p][1]);
			mm_store(((xint*)qs[2]) + i, u[0]);
			//c[0] = mm_packs_epi16(c[0], c[1]);
			//c[0] = mm_blendv(mrk[0], mrk[1], c[0]);
			//mm_store(((xint*)os[2]) + i, c[0]);
		}
		{
			// delta epi16 -> epi32
			s[0][0] = mm_add_epi32(s[0][0], mm_cvtepi16x0_epi32(t[0][0]));
			s[0][1] = mm_add_epi32(s[0][1], mm_cvtepi16x1_epi32(t[0][0]));
			s[0][2] = mm_add_epi32(s[0][2], mm_cvtepi16x0_epi32(t[0][1]));
			s[0][3] = mm_add_epi32(s[0][3], mm_cvtepi16x1_epi32(t[0][1]));
			s[1][0] = mm_add_epi32(s[1][0], mm_cvtepi16x0_epi32(t[1][0]));
			s[1][1] = mm_add_epi32(s[1][1], mm_cvtepi16x1_epi32(t[1][0]));
			s[1][2] = mm_add_epi32(s[1][2], mm_cvtepi16x0_epi32(t[1][1]));
			s[1][3] = mm_add_epi32(s[1][3], mm_cvtepi16x1_epi32(t[1][1]));
		}
	}
	banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[2], ubegs[2], W);
}

static inline int banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(xint h, xint u, xint v, b1i *us[2], b1i *fs, int *ubegs[2], u4i W){
	u4i i;
	UNUSED(W);
	// revise the first striped block of u
	v = mm_subs_epi8(h, u); // v(x - 1, y) = h(x - 1, y) - u(x - 1, y - 1)
	mm_store((xint*)fs, v);
	for(i=1;i<=WORDSIZE;i++){
		ubegs[1][i] = ubegs[0][i] + fs[i - 1];
	}
	v = mm_slli(v, 1);
	u = mm_load(((xint*)us[1]) + 0);
	u = mm_subs_epi8(u, v); // because I previously set v to zero, now update it, u(x, y) = h(x, y) - v(x - 1, y)
	mm_store(((xint*)us[1]) + 0, u);
	// shift score to fit EPI8
	ubegs[1][0] = ubegs[0][0] + us[1][0];
	us[1][0] = 0;
	banded_striped_epi8_seqalign_piecex_row_check_ubegs(us[1], ubegs[1], W);
	return ubegs[1][0];
}

//__attribute__((always_inline))
static inline xint banded_striped_epi8_seqalign_piecex_row_cal_FPenetration_codes(u4i W, xint f, b1i *fs, int *ubegs[2], b1i gape){
	int i, s, t;
	f = mm_slli(f, 1);
	mm_store((xint*)fs, f);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape;
	s = t + fs[0] - (ubegs[0][1] - ubegs[0][0]);
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s) fs[i] = s;
		s = t + fs[i] - (ubegs[0][i + 1] - ubegs[0][i]);
	}
	f = mm_load((xint*)fs);
	return f;
}

static inline int banded_striped_epi8_seqalign_piece0_row_btcal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
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
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + gape1; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ubegs[0][0] >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
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
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piece0_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
	xint h, e, f, u, v, z;
	xint GapE;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(es);
	UNUSED(qs);
	UNUSED(gapo1);
	UNUSED(gapo2);
	UNUSED(gape2);
	UNUSED(mov);
	//I     = mm_set1_epi8(SEQALIGN_BT_I);
	//D     = mm_set1_epi8(SEQALIGN_BT_D);
	GapE  = mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + gape1; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ubegs[0][0] >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
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
		h = mm_max_epi8(e, z);
		// max(f, h)
		h = mm_max_epi8(f, h);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate f(x, y)
		f = mm_adds_epi8(h, GapE);
		f = mm_subs_epi8(f, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piece1_row_btcal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
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
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
	t = us[0][0] + es[0][0]; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ubegs[0][0] >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
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
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piece1_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
	xint h, e, f, u, v, z;
	xint GapOE, GapE;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(qs);
	UNUSED(gapo2);
	UNUSED(gape2);
	UNUSED(mov);
	GapOE = mm_set1_epi8(gapo1 + gape1);
	GapE  = mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
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
		h = mm_max_epi8(e, z);
		// max(f, h)
		h = mm_max_epi8(f, h);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate e(x, y)
		e = mm_adds_epi8(e, GapE);
		e = mm_subs_epi8(e, h);
		e = mm_max_epi8(e, GapOE);
		mm_store(((xint*)es[1]) + i, e);
		// calculate f(x, y)
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		f = mm_max_epi8(f, h);
		f = mm_subs_epi8(f, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	h = mm_subs_epi8(h, GapOE);
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piece2_row_btcal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
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
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
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
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piece2_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh){
	xint h, e, q, f, g, u, v, z;
	xint GapOE, GapE, GapQP, GapP, GapOQ;
	b1i fs[WORDSIZE];
	u4i i;
	int t, h0, bt0;
	UNUSED(mov);
	GapOE = mm_set1_epi8(gapo1 + gape1);
	GapE  = mm_set1_epi8(gape1);
	GapQP = mm_set1_epi8(gapo2 + gape2);
	GapP  = mm_set1_epi8(gape2);
	GapOQ = mm_subs_epi8(GapOE, GapQP);
	// ::: max(h, e, f)
	f = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	g = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ubegs[0][0]) + qprof[((rbeg + 0) * 4 + base) * WORDSIZE]; // h
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
		h = mm_max_epi8(e, z);
		q = mm_adds_epi8(q, u);
		h = mm_max_epi8(q, h);
		// max(f, g, h)
		h = mm_max_epi8(f, h);
		h = mm_max_epi8(g, h);
		// calculate u(x, y)
		v = mm_subs_epi8(h, v);
		mm_store(((xint*)us[1]) + i, v);
		v = mm_subs_epi8(h, u);
		// calculate e(x, y) and q(x, y)
		e = mm_adds_epi8(e, GapE);
		e = mm_subs_epi8(e, h);
		e = mm_max_epi8(e, GapOE);
		mm_store(((xint*)es[1]) + i, e);
		q = mm_adds_epi8(q, GapP);
		q = mm_subs_epi8(q, h);
		q = mm_max_epi8(q, GapQP);
		mm_store(((xint*)qs[1]) + i, q);
		// calculate f(x, y) and g(x, y)
		f = mm_adds_epi8(f, GapE);
		h = mm_adds_epi8(h, GapOE);
		f = mm_max_epi8(f, h);
		f = mm_subs_epi8(f, u);
		g = mm_adds_epi8(g, GapP);
		h = mm_subs_epi8(h, GapOQ); // (gapo1 + gape1 - gapo2 - gape2)
		g = mm_max_epi8(g, h);
		g = mm_subs_epi8(g, u);
		z = mm_load(((xint*)qprof) + (rbeg + i + 1) * 4 + base);
	}
	h = mm_subs_epi8(h, GapQP);
	return banded_striped_epi8_seqalign_piecex_row_cal_tail_codes(h, u, v, us, fs, ubegs, W);
}

static inline int banded_striped_epi8_seqalign_piecex_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int rh, int piecewise){
	if(piecewise == 2)      return banded_striped_epi8_seqalign_piece2_row_cal(rbeg, base, us, es, qs, ubegs, qprof, gapo1, gape1, gapo2, gape2, W, mov, rh);
	else if(piecewise == 1) return banded_striped_epi8_seqalign_piece1_row_cal(rbeg, base, us, es, qs, ubegs, qprof, gapo1, gape1, gapo2, gape2, W, mov, rh);
	else                    return banded_striped_epi8_seqalign_piece0_row_cal(rbeg, base, us, es, qs, ubegs, qprof, gapo1, gape1, gapo2, gape2, W, mov, rh);
}

static inline u4i banded_striped_epi8_seqalign_row_max(b1i *us, int *ubegs, u4i W, int *max_score){
	xint h, c, xi, Max[4], max[2], Scr[4], scr[2], Idx[4], Pos[4];
	u4i i, j, x, y, STEP;
	int ary[2][WORDSIZE/4], uscr, umax;
	for(i=0;i<4;i++){
		Scr[i] = mm_load(((xint*)ubegs) + i);
		Max[i] = mm_set1_epi32(SEQALIGN_SCORE_MIN);
		Idx[i] = mm_set1_epi32(0);
	}
	for(i=0;i<WORDSIZE/4;i++){ ary[0][i] = i; } Idx[0] =  Pos[0] = mm_load((xint*)ary[0]);
	for(i=0;i<WORDSIZE/4;i++){ ary[0][i] = i + WORDSIZE / 4; } Idx[1] = Pos[1] = mm_load((xint*)ary[0]);
	for(i=0;i<WORDSIZE/4;i++){ ary[0][i] = i + 2 * WORDSIZE / 4; } Idx[2] = Pos[2] = mm_load((xint*)ary[0]);
	for(i=0;i<WORDSIZE/4;i++){ ary[0][i] = i + 3 * WORDSIZE / 4; } Idx[3] = Pos[3] = mm_load((xint*)ary[0]);
	STEP = 32;
	for(i=0;i<W;i+=x){
		x = num_min(i + STEP, W) - i;
		{
			max[0] = max[1] = mm_set1_epi16(-(MAX_B2));
			scr[0] = scr[1] = mm_set1_epi16(0);
		}
		for(j=0;j<x;j++){
			h = mm_load(((xint*)us) + i + j);
			scr[0] = mm_adds_epi16(scr[0], mm_cvtepi8x0_epi16(h));
			scr[1] = mm_adds_epi16(scr[1], mm_cvtepi8x1_epi16(h));
			max[0] = mm_max_epi16(max[0], scr[0]);
			max[1] = mm_max_epi16(max[1], scr[1]);
		}
		{
			xi = mm_set1_epi32(x << 8);
			Pos[0] = mm_add_epi32(Pos[0], xi);
			Pos[1] = mm_add_epi32(Pos[1], xi);
			Pos[2] = mm_add_epi32(Pos[2], xi);
			Pos[3] = mm_add_epi32(Pos[3], xi);
			h = mm_add_epi32(Scr[0], mm_cvtepi16x0_epi32(max[0])); Idx[0] = mm_blendv(Idx[0], Pos[0], mm_cmpgt_epi32(h, Max[0])); Max[0] = mm_max_epi32(Max[0], h);
			h = mm_add_epi32(Scr[1], mm_cvtepi16x1_epi32(max[0])); Idx[1] = mm_blendv(Idx[1], Pos[1], mm_cmpgt_epi32(h, Max[1])); Max[1] = mm_max_epi32(Max[1], h);
			h = mm_add_epi32(Scr[2], mm_cvtepi16x0_epi32(max[1])); Idx[2] = mm_blendv(Idx[2], Pos[2], mm_cmpgt_epi32(h, Max[2])); Max[2] = mm_max_epi32(Max[2], h);
			h = mm_add_epi32(Scr[3], mm_cvtepi16x1_epi32(max[1])); Idx[3] = mm_blendv(Idx[3], Pos[3], mm_cmpgt_epi32(h, Max[3])); Max[3] = mm_max_epi32(Max[3], h);
			Scr[0] = mm_add_epi32(Scr[0], mm_cvtepi16x0_epi32(scr[0]));
			Scr[1] = mm_add_epi32(Scr[1], mm_cvtepi16x1_epi32(scr[0]));
			Scr[2] = mm_add_epi32(Scr[2], mm_cvtepi16x0_epi32(scr[1]));
			Scr[3] = mm_add_epi32(Scr[3], mm_cvtepi16x1_epi32(scr[1]));
		}
	}
	c = mm_cmpgt_epi32(Max[1], Max[0]); Idx[0] = mm_blendv(Idx[0], Idx[1], c); Max[0] = mm_max_epi32(Max[0], Max[1]);
	c = mm_cmpgt_epi32(Max[3], Max[2]); Idx[1] = mm_blendv(Idx[2], Idx[3], c); Max[1] = mm_max_epi32(Max[2], Max[3]);
	c = mm_cmpgt_epi32(Max[1], Max[0]); Idx[0] = mm_blendv(Idx[0], Idx[1], c); Max[0] = mm_max_epi32(Max[0], Max[1]);
	mm_store(((xint*)ary[0]), Max[0]);
	mm_store(((xint*)ary[1]), Idx[0]);
	*max_score = ary[0][0];
	x = 0;
	for(i=1;i<WORDSIZE/4;i++){
		if(ary[0][i] > *max_score){
			*max_score = ary[0][i];
			x = i;
		}
	}
	x = ary[1][x];
	i = x & 0xFF;
	x = x >> 8;
	y = num_min(x + STEP, W);
	j = x;
	umax = uscr = 0;
	for(;x<y;x++){
		uscr += us[x * WORDSIZE + i];
		if(uscr > umax){
			j = x;
			umax = uscr;
		}
	}
	x = i * W + j;
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

static inline void banded_striped_epi8_seqalign_row_print(FILE *out, u1i *qseq, u4i qlen, u4i tidx, u4i tpos, u1i tbase, u4i bandwidth, u4i mov, u4i rbeg, u4i rmax, int max_score, int *ubegs, b1i *us, b1i *es, int detail){
	u4i i, x, W;
	int score;
	W = bandwidth / WORDSIZE;
	fprintf(out, "ROW[%d][%d][%c]\tMOV=%d\tBAND=%d,%d\tMAX=%d(%d),%d", tidx, tpos, "ACGTN-"[tbase], mov, rbeg, rbeg + bandwidth, rbeg + rmax, rmax, max_score);
	if(detail > 2){
		score = ubegs[0];
		for(i=0;i<bandwidth;i++){
			x = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
			fprintf(out, "\t%d:%c%d:%d:%d", i + rbeg, "ACGTN-"[rbeg + i < qlen? qseq[rbeg + i] : 4], score + us[x], us[x], es[x]);
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

static inline int banded_striped_epi8_seqalign_piecex_row_verify(int rowidx, int rowoff, int mode, int W, int mov, u1i tbase, b1i *us[2], b1i *es[2], b1i *qs[2], int *ubegs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	int i, x1, x2, s1, s2, s, h, u, e, q, f, g, piecewise;
	piecewise = banded_striped_epi8_seqalign_get_piecewise(gapo1, gape1, gapo2, gape2, W * WORDSIZE);
	s1 = ubegs[0][0];
	s2 = ubegs[1][0];
	f = g = SEQALIGN_SCORE_MIN;
	if(mov == 0){
		if(rowoff) s = SEQALIGN_SCORE_MIN;
		else if(mode == SEQALIGN_MODE_OVERLAP || rowidx == 0) s = 0;
		else if(piecewise < 2) s = gapo1 + gape1 * rowidx;
		else s = num_max(gapo1 + gape1 * rowidx, gapo2 + gape2 * rowidx);
		s += banded_striped_epi8_seqalign_get_qprof_value(qprof, rowoff, tbase);
		s -= s1;
		//if(rowidx == 0 || (rowoff == 0 && mode == SEQALIGN_MODE_OVERLAP)){
			//s = banded_striped_epi8_seqalign_get_qprof_value(qprof, rowoff, tbase);
		//} else {
			//s = SEQALIGN_SCORE_MIN;
		//}
	} else if(mov >= W * WORDSIZE){
		return 1;
	} else {
		for(i=0;i<mov&&i<W*WORDSIZE;i++){
			x1 = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
			u = us[0][x1];
			s1 += u;
			if(((i + 1) % W) == 0){
				if(ubegs[0][(i + 1) / W] != s1){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
		}
		s = banded_striped_epi8_seqalign_get_qprof_value(qprof, rowoff + mov, tbase);
	}
	for(i=mov;i<W*WORDSIZE;i++){
		x1 = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
		u = us[0][x1];
		e = (piecewise == 0)? gape1 : es[0][x1];
		q = (piecewise == 2)? qs[0][x1] : SEQALIGN_SCORE_MIN;
		h = num_max(s, u + e);
		h = num_max(h, u + q);
		h += s1;
		h = num_max(h, f);
		h = num_max(h, g);
		x2 = banded_striped_epi8_pos2idx(W << WORDSHIFT, i - mov);
		if(h > SEQALIGN_SCORE_MIN / 2 && h != s2 + us[1][x2]){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		f = (piecewise == 0)? h + gape1 : num_max(f + gape1, h + gapo1 + gape1);
		g = (piecewise == 2)? num_max(g + gape2, h + gapo2 + gape2) : SEQALIGN_SCORE_MIN;
		if(piecewise > 0 && h > SEQALIGN_SCORE_MIN / 2 && h + es[1][x2] != num_max(s1 + es[0][x1] + us[0][x1] + gape1, h + gapo1 + gape1)){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		if(piecewise == 2 && h > SEQALIGN_SCORE_MIN / 2 && h + qs[1][x2] != num_max(s1 + qs[0][x1] + us[0][x1] + gape2, h + gapo2 + gape2)){
			fflush(stdout); fprintf(stderr, " -- i=%d q:%d!=%d in %s -- %s:%d --\n", i, num_max(s1 + qs[0][x1] + us[0][x1] + gape2, h + gapo2 + gape2) - h, qs[1][x2], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		s1 += us[0][x1];
		s2 += us[1][x2];
		if(((i + 1) % W) == 0 && ubegs[0][(i + 1) / W] != s1){
			fflush(stdout); fprintf(stderr, " -- something wrong i = %d in %s -- %s:%d --\n", i, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		if(((i - mov + 1) % W) == 0 && ubegs[1][(i - mov + 1) / W] != s2){
			fflush(stdout); fprintf(stderr, " -- something wrong i = %d in %s -- %s:%d --\n", i, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		s = banded_striped_epi8_seqalign_get_qprof_value(qprof, rowoff + i + 1, tbase);
	}
	u = SEQALIGN_SCORE_EPI8_MIN;
	for(;i<W*WORDSIZE+mov;i++){
		h = num_max(s1 + s, f);
		h = num_max(h, g);
		//h = num_max(h, SEQALIGN_SCORE_MIN);
		x2 = banded_striped_epi8_pos2idx(W << WORDSHIFT, i - mov);
		if(h > SEQALIGN_SCORE_MIN / 2 && h != s2 + us[1][x2]){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		f = (piecewise == 0)? h + gape1 : num_max(f + gape1, h + gapo1 + gape1);
		g = (piecewise == 2)? num_max(g + gape2, h + gapo2 + gape2) : SEQALIGN_SCORE_MIN;
		if(piecewise > 0 && h > SEQALIGN_SCORE_MIN / 2 && h + es[1][x2] != h + gapo1 + gape1){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		if(piecewise == 2 && h > SEQALIGN_SCORE_MIN / 2 && h + qs[1][x2] != h + gapo2 + gape2){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		s1 += u;
		s2 += us[1][x2];
		if(((i - mov + 1) % W) == 0 && ubegs[1][(i - mov + 1) / W] != s2){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		u = 0;
		s = SEQALIGN_SCORE_MIN;
	}
	return 0;
}

static inline u4i banded_striped_epi8_seqalign_piecex_backtrace(u1i *qseq, u1i *tseq, b1i *bs, int *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, u4v *cigars){
	u4i btx, bty, cg;
	u1i qbase, tbase, op, bt;
	rs->qb = rs->qe; rs->qe ++;
	rs->tb = rs->te; rs->te ++;
	rs->mat = rs->mis = rs->ins = rs->del = rs->aln = 0;
	bty = bt = cg = 0;
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

static inline u4i _push_cigar_bsalign(u4v *cigars, u4i cg, u1i op, u4i sz){
	if(op == (cg & 0xf)){
		cg += sz << 4;
	} else {
		if(cigars && cg) push_u4v(cigars, cg);
		cg = sz << 4 | op;
	}
	return cg;
}

// x is the real_pos in row[tb] - row_offset[tb - 1]
// u, e, q come from (x, tb - 1)
static inline int banded_striped_epi8_seqalign_piecex_backcal_cell(int x, u1i qbase, u1i tbase, int Hs[2], b1i u, b1i e, b1i q, u4i W, b1i *matrix, int piecewise, int prior_match){
	int s, h;
	s = matrix[qbase * 4 + tbase];
	h = Hs[1] - Hs[0];
	if(x > Int(W * WORDSIZE)){ // SEQALIGN_BT_I only
		return SEQALIGN_BT_I;
	} else if(x == Int(W * WORDSIZE)){
		if(h == s){
			return SEQALIGN_BT_M;
		}
	} else if(prior_match){
		if(h == s){
			return SEQALIGN_BT_M;
		} else {
			if(h == u + e){
				return SEQALIGN_BT_D;
			} else if(piecewise == 2 && h == u + q){
				return SEQALIGN_BT2_D2;
			}
			return SEQALIGN_BT_I;
		}
	} else {
		if(h == u + e){
			return SEQALIGN_BT_D;
		} else if(piecewise == 2 && h == u + q){
			return SEQALIGN_BT2_D2;
		} else if(h == s){
			return SEQALIGN_BT_M;
		}
	}
	return SEQALIGN_BT_I;
}

static inline u4i banded_striped_epi8_seqalign_piecex_backcal(u1i *qseq, u1i *tseq, b1i *ups, b1i *eps, b1i *qps, b1i *ubs, int *roffs, int mode, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, seqalign_result_t *rs, u4v *cigars){
	u4i bt, cg, op, sz;
	int piecewise, prior_match, W, t, Hs[3];
	piecewise = banded_striped_epi8_seqalign_get_piecewise(gapo1, gape1, gapo2, gape2, bandwidth);
	rs->qb = rs->qe; rs->qe ++;
	rs->tb = rs->te; rs->te ++;
	rs->mat = rs->mis = rs->ins = rs->del = rs->aln = 0;
	cg = 0;
	W = bandwidth / WORDSIZE;
	if(cigars) clear_u4v(cigars);
	Hs[0] = 0; // H score to be checked
	Hs[1] = banded_striped_epi8_seqalign_getscore(ups + rs->tb * ((b8i)bandwidth), (int*)(ubs + rs->tb * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)), W, (rs->qb - roffs[rs->tb]));
	Hs[2] = 0; // (length << 4) | op; op <- 0:three_ways, 1:D, 3:D2
	prior_match = 0;
	while(1){
		if((Hs[2] & 0xf) == SEQALIGN_BT_D){ // E
			Hs[0] = banded_striped_epi8_seqalign_getscore(ups + (rs->tb) * ((b8i)bandwidth), (int*)(ubs + (rs->tb) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)), W, rs->qb - roffs[rs->tb]);
			t = gapo1 + (Hs[2] >> 4) * gape1;
			if(Hs[0] + t == Hs[1]){
				cg = _push_cigar_bsalign(cigars, cg, SEQALIGN_BT_D, Hs[2] >> 4);
				rs->ins += Hs[2] >> 4;
				rs->aln += Hs[2] >> 4;
				Hs[1] = Hs[0];
				Hs[2] = 0;
			} else {
				Hs[2] += 1 << 4;
				rs->tb --;
				continue;
			}
		} else if((Hs[2] & 0xf) == SEQALIGN_BT2_D2){ // Q
			Hs[0] = banded_striped_epi8_seqalign_getscore(ups + (rs->tb) * ((b8i)bandwidth), (int*)(ubs + (rs->tb) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)), W, rs->qb - roffs[rs->tb]);
			t = gapo2 + (Hs[2] >> 4) * gape2;
			if(Hs[0] + t == Hs[1]){
				cg = _push_cigar_bsalign(cigars, cg, SEQALIGN_BT_D, Hs[2] >> 4);
				rs->del += Hs[2] >> 4;
				rs->aln += Hs[2] >> 4;
				Hs[1] = Hs[0];
				Hs[2] = 0;
			} else {
				Hs[2] += 1 << 4;
				rs->tb --;
				continue;
			}
		}
		if(rs->qb < 0 || rs->tb < 0){
			break;
		}
		if(rs->qb == roffs[rs->tb - 1]){
			if(rs->qb){
				Hs[0] = SEQALIGN_SCORE_MIN;
				//Hs[0] = ((int*)(ubs + (rs->tb - 1) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)))[0];
			} else {
				if(mode == SEQALIGN_MODE_OVERLAP || rs->tb == 0) Hs[0] = 0;
				else if(piecewise < 2) Hs[0] = gapo1 + gape1 * rs->tb;
				else Hs[0] = num_max(gapo1 + gape1 * rs->tb, gapo2 + gape2 * rs->tb);
			}
		} else {
			Hs[0] = banded_striped_epi8_seqalign_getscore(ups + (rs->tb - 1) * ((b8i)bandwidth), (int*)(ubs + (rs->tb - 1) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)), W, rs->qb - 1 - roffs[rs->tb - 1]);
		}
		t = rs->qb - roffs[rs->tb - 1];
		t = (rs->tb - 1) * bandwidth + (t % W) * WORDSIZE + (t / W);
		bt = banded_striped_epi8_seqalign_piecex_backcal_cell(rs->qb - roffs[rs->tb - 1], qseq[rs->qb], tseq[rs->tb], Hs, ups[t], eps? eps[t] : gapo1 + gape1, qps? qps[t] : 0, W, matrix, piecewise, prior_match);
		prior_match = 1;
		if(bt == SEQALIGN_BT_M){
			if(qseq[rs->qb] == tseq[rs->tb]){
				rs->mat ++;
			} else {
				rs->mis ++;
			}
			rs->qb --;
			rs->tb --;
			rs->aln ++;
			cg = _push_cigar_bsalign(cigars, cg, 0, 1);
			Hs[1] = Hs[0];
		} else if(bt == SEQALIGN_BT_I){
			if(rs->qb <= 0){
				cg = _push_cigar_bsalign(cigars, cg, 1, sz);
				Hs[1] = Hs[0];
				rs->qb -= sz;
				rs->ins += sz;
				rs->aln += sz;
			} else {
				for(sz=1;Int(sz)+roffs[rs->tb]<=rs->qb;sz++){
					if(piecewise == 2){
						t = num_max(gapo1 + sz * gape1, gapo2 + sz * gape2);
					} else {
						t = gapo1 + sz * gape1;
					}
					Hs[0] = banded_striped_epi8_seqalign_getscore(ups + (rs->tb) * ((b8i)bandwidth), (int*)(ubs + (rs->tb) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE)), W, rs->qb - sz - roffs[rs->tb]);
					if(Hs[0] + t == Hs[1]){
						cg = _push_cigar_bsalign(cigars, cg, 1, sz);
						Hs[1] = Hs[0];
						rs->qb -= sz;
						rs->ins += sz;
						rs->aln += sz;
						break;
					}
				}
			}
		} else {
			Hs[2] = (1 << 4) | bt;
			rs->tb --;
			continue;
		}
	}
	if(mode == SEQALIGN_MODE_OVERLAP){
		if(cg && cigars){
			//printf("PUSH %d:%d\n", cg & 0xf, cg >> 4);
			push_u4v(cigars, cg);
		}
	} else {
		if(rs->qb >= 0){
			op = 1;
			sz = rs->qb + 1;
			rs->ins += sz;
			rs->qb = -1;
		} else if(rs->tb >= 0){
			op = 2;
			sz = rs->tb + 1;
			rs->del += sz;
			rs->tb = -1;
		} else {
			op = sz = 0;
		}
		rs->aln += sz;
		cg = _push_cigar_bsalign(cigars, cg, op, sz);
		if(cg && cigars){
			//printf("PUSH %d:%d\n", cg & 0xf, cg >> 4);
			push_u4v(cigars, cg);
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
	banded_striped_epi8_seqalign_piecex_row_mov_func    row_movx;
	//banded_striped_epi8_seqalign_piecex_row_mov_func    row_mov;
	banded_striped_epi8_seqalign_piecex_row_cal_func    row_cal;
	banded_striped_epi8_seqalign_piecex_row_verify_func row_verify;
	banded_striped_epi8_seqalign_row_max_func           row_max;
	banded_striped_epi8_seqalign_band_mov_func          band_mov;
	banded_striped_epi8_seqalign_getscore_func          getscore;
	banded_striped_epi8_seqalign_piecex_backtrace_func  backtrace;
	banded_striped_epi8_seqalign_piecex_backcal_func    backcal;
	seqalign_result_t rs;
	b1i *memp, *mempb, *qprof, *us[2], *es[2], *qs[2], *ups, *eps, *qps, *ubs;
	u8i mpsize;
	u4i i, W, rbeg, rmax, mov;
	int piecewise, smax, smin, max_score, score, rh, rbx, rby, rbz, *begs, *ubegs[2];
	u1i tbase;
	set_qprof  = banded_striped_epi8_seqalign_set_query_prof;
	row_init   = banded_striped_epi8_seqalign_piecex_row_init;
	//row_mov    = banded_striped_epi8_seqalign_piecex_row_mov;
	row_movx   = banded_striped_epi8_seqalign_piecex_row_movx;
	row_cal    = banded_striped_epi8_seqalign_piecex_row_cal;
	row_verify = banded_striped_epi8_seqalign_piecex_row_verify;
	row_max    = banded_striped_epi8_seqalign_row_max;
	band_mov   = banded_striped_epi8_seqalign_band_mov;
	getscore   = banded_striped_epi8_seqalign_getscore;
	backtrace  = banded_striped_epi8_seqalign_piecex_backtrace;
	backcal    = banded_striped_epi8_seqalign_piecex_backcal;
	bandwidth = roundup_times(bandwidth, WORDSIZE);
	W = bandwidth / WORDSIZE;
	piecewise = banded_striped_epi8_seqalign_get_piecewise(gapo1, gape1, gapo2, gape2, bandwidth);
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
	mpsize += ((u8i)bandwidth) * ((tlen + 1) * (piecewise + 1)); // ups[], eps[], qps[]
	mpsize += (tlen + 2) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE); // ubegs
	mpsize += bandwidth * ((piecewise + 1)); // us[0][], es[0][], qs[0][]
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
	ups = memp + bandwidth; memp += bandwidth * (tlen + 1);
	if(piecewise){
		eps = memp + bandwidth; memp += bandwidth * (tlen + 1);
	} else eps = NULL;
	if(piecewise == 2){
		qps = memp + bandwidth; memp += bandwidth * (tlen + 1);
	} else qps = NULL;
	ubs = memp + roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE); memp += (tlen + 1) * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE);
	ubegs[0] = (int*)memp; memp += 1 * roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE);
	us[0] = memp; memp += bandwidth;
	if(piecewise){
		es[0] = memp; memp += bandwidth;
	} else es[0] = NULL;
	if(piecewise == 2){
		qs[0] = memp; memp += bandwidth;
	} else qs[0] = NULL;
	begs = ((int*)memp) + 1; memp += roundup_times((tlen + 1) * sizeof(u4i), WORDSIZE);
	// prepare
	set_qprof(qseq, qlen, qprof, bandwidth, matrix);
	memset(&rs, 0, sizeof(seqalign_result_t));
	rs.score = SEQALIGN_SCORE_MIN;
	us[1] = ups - bandwidth;
	ubegs[1] = (int*)(ubs - roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE));
	es[1] = eps? eps - bandwidth : NULL;
	qs[1] = qps? qps - bandwidth : NULL;
	row_init(us[1], es[1], qs[1], ubegs[1], begs - 1, mode, bandwidth, smax, smin, gapo1, gape1, gapo2, gape2);
	// loop rows
	rbeg = rmax = 0;
	mov = 0;
	for(i=0;i<tlen;i++){
		tbase = tseq[i];
		if(mov && rbeg + bandwidth < qlen){
			mov = num_min(mov, num_max(0, Int(qlen) - Int(rbeg + bandwidth)));
			rbeg += mov;
			rh = getscore(us[1], ubegs[1], W, mov - 1);
		} else {
			mov = 0;
			//rh = (rbeg == 0 && (mode == SEQALIGN_MODE_OVERLAP || i == 0))? 0 : SEQALIGN_SCORE_MIN;
			if(rbeg){
				rh = SEQALIGN_SCORE_MIN;
			} else {
				if(mode == SEQALIGN_MODE_OVERLAP || i == 0) rh = 0;
				else if(piecewise < 2) rh = gapo1 + gape1 * i;
				else rh = num_max(gapo1 + gape1 * i, gapo2 + gape2 * i);
			}
		}
		row_movx(us, es, qs, ubegs, W, mov, piecewise, smax, smin, gapo1, gape1, gapo2, gape2); // mov and swap
#if 0
		{
			b1i *tus[2], *tes[2], *tqs[2];
			int *tubegs[2];
			u4i j;
			tus[0] = (b1i*)malloc(bandwidth);
			tes[0] = eps? (b1i*)malloc(bandwidth) : NULL;
			tqs[0] = qps? (b1i*)malloc(bandwidth) : NULL;
			tubegs[0] = (int*)malloc((WORDSIZE + 1) * sizeof(int));
			tus[1] = us[1];
			tes[1] = es[1];
			tqs[1] = qs[1];
			tubegs[1] = ubegs[1];
			row_mov(tus, tes, tqs, tubegs, W, mov, piecewise); // mov and swap
			if(memcmp(tus[0], us[0], bandwidth)){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				for(j=0;j<bandwidth;j++){
					if(tus[0][j] != us[0][j]){
						fflush(stdout); fprintf(stderr, " -- j=%d %d!=%d in %s -- %s:%d --\n", j, tus[0][j], us[0][j], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					}
				}
				abort();
			}
			if(memcmp(tubegs[0], ubegs[0], (WORDSIZE + 1) * sizeof(int))){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				for(j=0;j<=WORDSIZE;j++){
					if(tubegs[0][j] != ubegs[0][j]){
						fflush(stdout); fprintf(stderr, " -- j=%d %d!=%d in %s -- %s:%d --\n", j, tubegs[0][j], ubegs[0][j], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					}
				}
				abort();
			}
		}
#endif
		us[1] = ups + ((u8i)bandwidth) * i;
		es[1] = eps? eps + ((u8i)bandwidth) * i : NULL;
		qs[1] = qps? qps + ((u8i)bandwidth) * i : NULL;
		ubegs[1] = (int*)(ubs + roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE) * i);
		__builtin_prefetch(us[1], 1);
		row_cal(rbeg, tbase, us, es, qs, ubegs, qprof, gapo1, gape1, gapo2, gape2, W, mov, rh, piecewise);
		if(verbose){
			rmax = row_max(us[1], ubegs[1], W, &max_score);
			banded_striped_epi8_seqalign_row_print(stdout, qseq, qlen, 1, i, tbase, bandwidth, mov, rbeg, rmax, max_score, ubegs[1], us[1], es[1], verbose);
		}
		if(verbose > 3){
			b1i *tus[2], *tes[2], *tqs[2];
			int *tubegs[2];
			tus[0] = ups + ((u8i)bandwidth) * (Int(i) - 1);
			tes[0] = eps? eps + ((u8i)bandwidth) * (Int(i) - 1) : NULL;
			tqs[0] = qps? qps + ((u8i)bandwidth) * (Int(i) - 1) : NULL;
			tubegs[0] = (int*)(ubs + roundup_times((WORDSIZE + 1) * sizeof(int), WORDSIZE) * (Int(i) - 1));
			tus[1] = us[1];
			tes[1] = es[1];
			tqs[1] = qs[1];
			tubegs[1] = ubegs[1];
			row_verify(i, rbeg - mov, mode, W, mov, tbase, tus, tes, tqs, tubegs, qprof, gapo1, gape1, gapo2, gape2);
		}
		// adaptive banded
		rbx = band_mov(us[1], ubegs[1], W, i, rbeg, qlen);
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
				score = getscore(us[1], ubegs[1], W, qlen - 1 - rbeg);
				if(score > rs.score){
					rs.score = score;
					rs.qe = qlen - 1;
					rs.te = i;
				}
			}
		}
	}
	if(mode == SEQALIGN_MODE_GLOBAL){
		rs.score = getscore(us[1], ubegs[1], W, qlen - 1 - rbeg);
		rs.qe = qlen - 1;
		rs.te = tlen - 1;
	} else {
		rmax = row_max(us[1], ubegs[1], W, &max_score);
		if(max_score > rs.score){
			rs.score = max_score;
			rs.qe = rbeg + rmax;
			rs.te = tlen - 1;
		}
	}
	// backtrace
	backcal(qseq, tseq, ups, eps, qps, ubs, begs, mode, bandwidth, matrix, gapo1, gape1, gapo2, gape2, &rs, cigars);
	if(mempb) free(mempb);
	return rs;
}

#endif
