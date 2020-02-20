#ifndef BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H
#define BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H

/**
 *
 * bsalign.h
 *
 * Jue Ruan <ruanjue@gmail.com>
 *
 * References:
 * Farrar, Michael. 2007. "Striped Smith-Waterman Speeds Database Searches Six Times over Other SIMD Implementations." Bioinformatics 23 (2): 156¨C61. https://doi.org/10.1093/bioinformatics/btl5
 * Suzuki, Hajime, and Masahiro Kasahara. 2017. "Acceleration of Nucleotide Semi-Global Alignment with Adaptive Banded Dynamic Programming" BioRxiv, September, 130633. https://doi.org/10.1101/13063
 * Suzuki, Hajime, and Masahiro Kasahara. 2018. "Introducing Difference Recurrence Relations for Faster Semi-Global Alignment of Long Sequences." BMC Bioinformatics 19 (Suppl 1). https://doi.org/10.1186/s12859-018-2018.
 * Li, Heng. 2018. "Minimap2: Pairwise Alignment for Nucleotide Sequences." Bioinformatics 34 (18): 3094¨C3100. https://doi.org/10.1093/bioinformatics/b
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
 *  By spliting it into pieces of 256 epi8, I can use epi16 instead of epi32 to speed up the sum and max
 *  Also, scores of the first and last running blocks can be stored during this procedure
 *
 * 6, adaptive band
 * if sum(H[0 .. W - 1]) > sum(H[15 * W .. 16 * W - 1])  row_offset = row_offset + 0 // in normal coordinate
 * if sum(H[0 .. W - 1]) == sum(H[15 * W .. 16 * W - 1]) row_offset = row_offset + 1
 * if sum(H[0 .. W - 1]) < sum(H[15 * W .. 16 * W - 1])  row_offset = row_offset + 2
 * keep the bandwidth, but shift the offset of band in this row for next call, see 1) shifting
 *
 */

/*
 * To use bsalign.h in your program, please copy bsalign.h, list.h and mem_share.h together
 */

#include "list.h"
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
#define mm_store	_mm256_store_si256
#define mm_or	_mm256_or_si256
#define mm_and	_mm256_and_si256
#define mm_set1_epi8	_mm256_set1_epi8
#define mm_set1_epi16	_mm256_set1_epi16
#define mm_set1_epi32	_mm256_set1_epi32
#define mm_srli	_mm256_srli_si256
#define mm_slli	_mm256_slli_si256
#define mm_insert_epi8	_mm256_insert_epi8
#define mm_extract_epi16	_mm256_extract_epi16
#define mm_extract_epi32	_mm256_extract_epi32
#define mm_adds_epi8	_mm256_adds_epi8
#define mm_adds_epi16	_mm256_adds_epi16
#define mm_add_epi32	_mm256_add_epi32
#define mm_subs_epi8	_mm256_subs_epi8
#define mm_subs_epi16	_mm256_subs_epi16
#define mm_cmpgt_epi8	_mm256_cmpgt_epi8
#define mm_cmpgt_epi32	_mm256_cmpgt_epi32
#define mm_max_epi8	_mm256_max_epi8
#define mm_max_epi16	_mm256_max_epi16
#define mm_blendv_epi8	_mm256_blendv_epi8
#define mm_cvtepi8lo_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(a))
#define mm_cvtepi8hi_epi16(a)	_mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))
#define mm_cvtepi16lo_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(a))
#define mm_cvtepi16hi_epi32(a)	_mm256_cvtepi16_epi32(_mm256_castsi256_si128(_mm256_srli_si256(a, 16)))

#else

//#pragma message("Choose SSE4.2 in " __FILE__ ". Just a message, ignore it")
#define WORDSIZE	16
#define WORDSHIFT	4
typedef __m128i	xint;
#define mm_load	_mm_load_si128
#define mm_store	_mm_store_si128
#define mm_or	_mm_or_si128
#define mm_and	_mm_and_si128
#define mm_set1_epi8	_mm_set1_epi8
#define mm_set1_epi16	_mm_set1_epi16
#define mm_set1_epi32	_mm_set1_epi32
#define mm_srli	_mm_srli_si128
#define mm_slli	_mm_slli_si128
#define mm_insert_epi8	_mm_insert_epi8
#define mm_extract_epi16	_mm_extract_epi16
#define mm_extract_epi32	_mm_extract_epi32
#define mm_adds_epi8	_mm_adds_epi8
#define mm_adds_epi16	_mm_adds_epi16
#define mm_add_epi32	_mm_add_epi32
#define mm_subs_epi8	_mm_subs_epi8
#define mm_subs_epi16	_mm_subs_epi16
#define mm_cmpgt_epi8	_mm_cmpgt_epi8
#define mm_cmpgt_epi32	_mm_cmpgt_epi32
#define mm_max_epi8	_mm_max_epi8
#define mm_max_epi16	_mm_max_epi16
#define mm_blendv_epi8	_mm_blendv_epi8
#define mm_cvtepi8lo_epi16(a)	_mm_cvtepi8_epi16(a)
#define mm_cvtepi8hi_epi16(a)	_mm_cvtepi8_epi16(_mm_srli_si128(a, 8))
#define mm_cvtepi16lo_epi32(a)	_mm_cvtepi16_epi32(a)
#define mm_cvtepi16hi_epi32(a)	_mm_cvtepi16_epi32(_mm_srli_si128(a, 8))

#endif

typedef struct {
	int score;
	int qb, qe;
	int tb, te;
	int mat, mis, ins, del, aln;
} seqalign_result_t;

/**
 * Basic function referings for global DNA sequence alignment, there are implementations in this file
 */

#define banded_striped_epi8_pos2idx(bandwidth, pos) ((((pos) % (bandwidth >> WORDSHIFT)) << WORDSHIFT) + ((pos) / ((bandwidth) >> WORDSHIFT)))

static inline void banded_striped_epi8_seqalign_set_score_matrix(b1i matrix[16], b1i mat, b1i mis){ u4i i; for(i=0;i<16;i++) matrix[i] = ((i ^ (i >> 2)) & 0x3)? mis : mat; }

#define banded_striped_epi8_seqalign_qprof_size(qlen, bandwidth) (((num_max(qlen, bandwidth) + 1) * 4) * WORDSIZE)
#define banded_striped_epi8_seqalign_get_qprof_value(qprof, pos, base) (qprof)[((pos) * 4 + (base)) * WORDSIZE]

// prepare query profile
typedef void (*banded_striped_epi8_seqalign_set_query_prof_func)(u1i *qseq, u4i qlen, b1i *qprof, u4i bandwidth, b1i score_matrix[16]);

// prepare the rows[-1]
// mode = SEQALIGN_MODE_GLOBAL/SEQALIGN_MODE_EXTEND/SEQALIGN_MODE_OVERLAP
typedef void (*banded_striped_epi8_seqalign_piecex_row_init_func)(b1i *us[2], b1i *es[2], b1i *qs[2], int *begs, int mode, u4i bandwidth, b1i max_nt_score, b1i min_nt_score, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// W = bandwidth / WORDSIZE, WORDSIZE = 16 (SSEx 128) or 32 (AVX2 512)
// mov <= W
// converting from us[1] to us[0]
// make the two row aligned
typedef void (*banded_striped_epi8_seqalign_piecex_row_mov_func)(b1i *us[2], b1i *es[2], b1i *qs[2], u4i W, u4i mov);

// core func to update row scores, from us[0] to us[1]
// ph is the score of H(-1, y - 1)
// rh is the revised score of H(-1, y - 1)
// @return: score of H(-1, y), note that H(-1, y) is useful to restore all scores of row in row_max
typedef int (*banded_striped_epi8_seqalign_piecex_row_cal_func)(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh);

typedef int (*banded_striped_epi8_seqalign_piecex_row_verify_func)(int rbeg, int W, int ph, int rh, int hh, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// find max score and return the real offset in band
// wscores[0] store the first running block
// wscores[1] store the last  running block
typedef u4i (*banded_striped_epi8_seqalign_piecex_row_max_func)(b1i *us, u4i W, int hh, int *wscores[2], int *max_score);

// wscores[0] store the first running block
// wscores[1] store the last  running block
typedef void (*banded_striped_epi8_seqalign_piecex_row_look_func)(b1i *us, u4i W, int hh, int *wscores[2]);

// backtrace
// rs->qe, rs->te and rs->score MUST be set before call this function
// begs provides the band's offset of rows, it is continously suming of mov in banded_striped_epi8_seqalign_piecex_row_mov_func
//typedef u4i (*banded_striped_epi8_seqalign_piecex_trace_func)(u1i *qseq, u1i *tseq, b1i **ups, b1i **vps, int *begs, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, seqalign_result_t *rs, u4v *cigars);
typedef u4i (*banded_striped_epi8_seqalign_piecex_backtrace_func)(u1i *qseq, u1i *tseq, b1i *bs, int *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, u4v *cigars);


static inline u4i banded_striped_epi8_seqalign_cigar2alnstr(u1i *qseq, u1i *tseq, seqalign_result_t *rs, u4v *cigars, char *alnstr[3], u4i length);

// implementation of overlap alignment for two sequences
// bandwidth should be times of WORDSIZE
static inline seqalign_result_t banded_striped_epi8_seqalign_pairwise(u1i *qseq, u4i qlen, u1i *tseq, u4i tlen, b1v *mempool, u4v *cigars, int mode, u4i bandwidth, b1i matrix[16], b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, int verbose);

static void banded_striped_epi8_seqalign_piecex_row_init(b1i *us[2], b1i *es[2], b1i *qs[2], int *begs, int mode, u4i bandwidth, b1i max_nt, b1i min_nt, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	xint ZERO, MIN, GAP;
	u4i k, W, xp;
	W = bandwidth / WORDSIZE;
	ZERO = mm_set1_epi8(0);
	MIN = mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	if(mode == SEQALIGN_MODE_GLOBAL || mode == SEQALIGN_MODE_EXTEND){
		if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1 && (gapo1 - gapo2) / (gape1 - gape2) < Int(bandwidth)){
			xp = (gapo2 - gapo1) / (gape1 - gape2);
			GAP = mm_set1_epi8(gape2);
			for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, GAP);
			us[1][0] = gapo1 + gape1 + min_nt - max_nt;
			for(k=2;k<=xp;k++){
				us[1][banded_striped_epi8_pos2idx(bandwidth, k)] = gape1;
			}
		} else {
			GAP = mm_set1_epi8(gape1);
			for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, GAP);
			us[1][0] = gapo1 + gape1 + min_nt - max_nt;
		}
		begs[-1] = - Int(bandwidth);
	} else {
		for(k=0;k<W;k++) mm_store(((xint*)us[1]) + k, ZERO);
		begs[-1] = 0;
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

static inline void banded_striped_epi8_seqalign_piece0_row_mov(b1i *us[2], b1i *es[2], b1i *qs[2], u4i W, u4i mov){
	xint u;
	u4i i, div;
	UNUSED(es);
	UNUSED(qs);
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = mm_load(((xint*)us[1]) + i + mov);
		mm_store(((xint*)us[0]) + i, u);
	}
	if(!mov) return;
	{
		u = mm_load(((xint*)us[1]) + i + mov - W);
		u = mm_srli(u, 1);
		u = mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, WORDSIZE - 1);
	}
	for(;i<W;i++){
		mm_store(((xint*)us[0]) + i, u);
		u = mm_load(((xint*)us[1]) + (i + mov - W + 1) % W);
		u = mm_srli(u, 1);
	}
}

static inline void banded_striped_epi8_seqalign_piece1_row_mov(b1i *us[2], b1i *es[2], b1i *qs[2], u4i W, u4i mov){
	xint u, e;
	u4i i, div;
	UNUSED(qs);
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = mm_load(((xint*)us[1]) + i + mov);
		e = mm_load(((xint*)es[1]) + i + mov);
		mm_store(((xint*)us[0]) + i, u);
		mm_store(((xint*)es[0]) + i, e);
	}
	if(!mov) return;
	{
		u = mm_load(((xint*)us[1]) + i + mov - W);
		u = mm_srli(u, 1);
		u = mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, WORDSIZE - 1);
	}
	for(;i<W;i++){
		e = mm_load(((xint*)es[1]) + i + mov - W);
		e = mm_srli(e, 1);
		mm_store(((xint*)us[0]) + i, u);
		mm_store(((xint*)es[0]) + i, e);
		u = mm_load(((xint*)us[1]) + (i + mov - W + 1) % W);
		u = mm_srli(u, 1);
	}
}

static inline void banded_striped_epi8_seqalign_piece2_row_mov(b1i *us[2], b1i *es[2], b1i *qs[2], u4i W, u4i mov){
	xint u, e, x;
	u4i i, div;
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = mm_load(((xint*)us[1]) + i + mov);
		e = mm_load(((xint*)es[1]) + i + mov);
		x = mm_load(((xint*)qs[1]) + i + mov);
		mm_store(((xint*)us[0]) + i, u);
		mm_store(((xint*)es[0]) + i, e);
		mm_store(((xint*)qs[0]) + i, x);
	}
	if(!mov) return;
	{
		u = mm_load(((xint*)us[1]) + i + mov - W);
		u = mm_srli(u, 1);
		u = mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, WORDSIZE - 1);
	}
	for(;i<W;i++){
		e = mm_load(((xint*)es[1]) + i + mov - W);
		x = mm_load(((xint*)qs[1]) + i + mov - W);
		e = mm_srli(e, 1);
		x = mm_srli(x, 1);
		mm_store(((xint*)us[0]) + i, u);
		mm_store(((xint*)es[0]) + i, e);
		mm_store(((xint*)qs[0]) + i, x);
		u = mm_load(((xint*)us[1]) + (i + mov - W + 1) % W);
		u = mm_srli(u, 1);
	}
}

static inline int banded_striped_epi8_seqalign_piece0_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, f, u, v, c, d, z, m2[2], m4[4];
	xint I, D, GapE;
	int ms[WORDSIZE];
	b1i fs[WORDSIZE];
	u4i i, k;
	int s, t, h0, bt0;
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
	m2[0] = mm_set1_epi16(0);
	m2[1] = mm_set1_epi16(0);
	m4[0] = mm_set1_epi32(0);
	m4[1] = mm_set1_epi32(0);
	m4[2] = mm_set1_epi32(0);
	m4[3] = mm_set1_epi32(0);
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
	for(i=0;i<W;){
		k = num_min(i + 256, W);
		for(;i<k;i++){
			u = mm_load(((xint*)us[0]) + i);
			// sums us[0] within each running block
			m2[0] = mm_subs_epi16(m2[0], mm_cvtepi8lo_epi16(u));
			m2[1] = mm_subs_epi16(m2[1], mm_cvtepi8hi_epi16(u));
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
		{
			m4[0] = mm_add_epi32(m4[0], mm_cvtepi16lo_epi32(m2[0]));
			m4[1] = mm_add_epi32(m4[1], mm_cvtepi16hi_epi32(m2[0]));
			m4[2] = mm_add_epi32(m4[2], mm_cvtepi16lo_epi32(m2[1]));
			m4[3] = mm_add_epi32(m4[3], mm_cvtepi16hi_epi32(m2[1]));
			m2[0] = mm_set1_epi16(0);
			m2[1] = mm_set1_epi16(0);
		}
	}
	mm_store(((xint*)ms) + 0, m4[0]);
	mm_store(((xint*)ms) + 1, m4[1]);
	mm_store(((xint*)ms) + 2, m4[2]);
	mm_store(((xint*)ms) + 3, m4[3]);
	f = mm_slli(f, 1);
	mm_store((xint*)fs, f);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + fs[0] + ms[0];
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s) fs[i] = s;
		s = t + fs[i] + ms[i];
	}
	f = mm_load((xint*)fs);
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
		d = mm_blendv_epi8(d, I, c);
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
	// revise the first striped block of u
	v = mm_subs_epi8(h, u);
	v = mm_slli(v, 1);
	u = mm_load(((xint*)us[1]) + 0);
	u = mm_subs_epi8(u, v);
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
	return rh;
}

static inline int banded_striped_epi8_seqalign_piece1_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, f, u, v, c, d, z, m2[2], m4[4];
	xint I, IE, D, DE, GapOE, GapE;
	int ms[WORDSIZE];
	b1i fs[WORDSIZE];
	u4i i, k;
	int s, t, h0, bt0;
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
	m2[0] = mm_set1_epi16(0);
	m2[1] = mm_set1_epi16(0);
	m4[0] = mm_set1_epi32(0);
	m4[1] = mm_set1_epi32(0);
	m4[2] = mm_set1_epi32(0);
	m4[3] = mm_set1_epi32(0);
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
	inline void banded_striped_epi8_seqalign_piece1_row_cal_sub1(){
		h = mm_load(((xint*)qprof) + (rbeg + 0) * 4 + base);
		h = mm_insert_epi8(h, h0, 0);
		for(i=0;i<W;){
			k = num_min(i + 256, W);
			for(;__builtin_expect(i<k, 0);i++){
				u = mm_load(((xint*)us[0]) + i);
				e = mm_load(((xint*)es[0]) + i);
				// sums us[0] within each running block
				m2[0] = mm_subs_epi16(m2[0], mm_cvtepi8lo_epi16(u));
				m2[1] = mm_subs_epi16(m2[1], mm_cvtepi8hi_epi16(u));
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
			{
				m4[0] = mm_add_epi32(m4[0], mm_cvtepi16lo_epi32(m2[0]));
				m4[1] = mm_add_epi32(m4[1], mm_cvtepi16hi_epi32(m2[0]));
				m4[2] = mm_add_epi32(m4[2], mm_cvtepi16lo_epi32(m2[1]));
				m4[3] = mm_add_epi32(m4[3], mm_cvtepi16hi_epi32(m2[1]));
				m2[0] = mm_set1_epi16(0);
				m2[1] = mm_set1_epi16(0);
			}
		}
	}
	banded_striped_epi8_seqalign_piece1_row_cal_sub1();
	mm_store(((xint*)ms) + 0, m4[0]);
	mm_store(((xint*)ms) + 1, m4[1]);
	mm_store(((xint*)ms) + 2, m4[2]);
	mm_store(((xint*)ms) + 3, m4[3]);
	f = mm_slli(f, 1);
	mm_store((xint*)fs, f);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + fs[0] + ms[0];
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s) fs[i] = s;
		s = t + fs[i] + ms[i];
	}
	f = mm_load((xint*)fs);
	// main loop
	// don't use the h from last W - 1 to calculate v = h - u, because this h maybe updated when F-penetration
	// will revise v after this loop
	u = mm_set1_epi8(0); // useless, but for compiler
	inline void banded_striped_epi8_seqalign_piece1_row_cal_sub2(){ // for gprof
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
			d = mm_blendv_epi8(d, I, c);
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
	}
	banded_striped_epi8_seqalign_piece1_row_cal_sub2();
	// revise the first striped block of u and v
	v = mm_subs_epi8(h, u);
	v = mm_slli(v, 1);
	u = mm_load(((xint*)us[1]) + 0);
	u = mm_subs_epi8(u, v);
	mm_store(((xint*)us[1]) + 0, u);
	// shift score to fit EPI8
	if(bt0 == SEQALIGN_BT_M){ // set us[1][0][0] to non-zero
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
	return rh;
}

static inline int banded_striped_epi8_seqalign_piece2_row_cal(u4i rbeg, u1i base, b1i *us[2], b1i *es[2], b1i *qs[2], b1i *bs, b1i *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, u4i mov, int ph, int rh){
	xint h, e, q, f, g, u, v, c, d, z, m2[2], m4[4];
	xint I1, I2, IE1, IE2, D1, D2, DE1, DE2, GapOE, GapE, GapQP, GapP, GapOQ;
	int ms[WORDSIZE];
	b1i fs[WORDSIZE];
	u4i i, k;
	int s, t, h0, bt0;
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
	m2[0] = mm_set1_epi16(0);
	m2[1] = mm_set1_epi16(0);
	m4[0] = mm_set1_epi32(0);
	m4[1] = mm_set1_epi32(0);
	m4[2] = mm_set1_epi32(0);
	m4[3] = mm_set1_epi32(0);
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
	for(i=0;i<W;){
		k = num_min(i + 256, W);
		for(;i<k;i++){
			u = mm_load(((xint*)us[0]) + i);
			// sums us[0] within each running block
			m2[0] = mm_subs_epi16(m2[0], mm_cvtepi8lo_epi16(u));
			m2[1] = mm_subs_epi16(m2[1], mm_cvtepi8hi_epi16(u));
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
		{
			m4[0] = mm_add_epi32(m4[0], mm_cvtepi16lo_epi32(m2[0]));
			m4[1] = mm_add_epi32(m4[1], mm_cvtepi16hi_epi32(m2[0]));
			m4[2] = mm_add_epi32(m4[2], mm_cvtepi16lo_epi32(m2[1]));
			m4[3] = mm_add_epi32(m4[3], mm_cvtepi16hi_epi32(m2[1]));
			m2[0] = mm_set1_epi16(0);
			m2[1] = mm_set1_epi16(0);
		}
	}
	mm_store(((xint*)ms) + 0, m4[0]);
	mm_store(((xint*)ms) + 1, m4[1]);
	mm_store(((xint*)ms) + 2, m4[2]);
	mm_store(((xint*)ms) + 3, m4[3]);
	f = mm_slli(f, 1);
	mm_store((xint*)fs, f);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + fs[0] + ms[0];
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s) fs[i] = s;
		s = t + fs[i] + ms[i];
	}
	f = mm_load((xint*)fs);
	g = mm_slli(g, 1);
	mm_store((xint*)fs, g);
	fs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape2;
	s = t + fs[0] + ms[0];
	for(i=1;i<WORDSIZE;i++){
		if(fs[i] < s){
			fs[i] = s;
		}
		s = t + fs[i] + ms[i];
	}
	g = mm_load((xint*)fs);
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
		d = mm_blendv_epi8(d, D2, c);
		h = mm_max_epi8(q, h);
		// max(f, g, h)
		c = mm_cmpgt_epi8(f, h);
		d = mm_blendv_epi8(d, I1, c);
		h = mm_max_epi8(f, h);
		c = mm_cmpgt_epi8(g, h);
		d = mm_blendv_epi8(d, I2, c);
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
	// revise the first striped block of u
	v = mm_subs_epi8(h, u);
	v = mm_slli(v, 1);
	u = mm_load(((xint*)us[1]) + 0);
	u = mm_subs_epi8(u, v);
	mm_store(((xint*)us[1]) + 0, u);
	// set last mov elements of v to SEQALIGN_BT
	// shift score to fit EPI8
	if(bt0 == SEQALIGN_BT_M){ // set us[1][0][0] to non-zero
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
	return rh;
}

static inline u4i banded_striped_epi8_seqalign_row_max(b1i *us, u4i W, int shift, int *wscores[2], int *max_score){
	xint h, c, m, Max[4], max[2], Scr[4], scr[2], Idx[4], num[2], ONE;
	u4i k, i, j, x, y, WX;
	int score, tmp, ary[WORDSIZE];
	// find max
	for(i=0;i<4;i++){
		Idx[i] = Scr[i] = Max[i] = mm_set1_epi32(0);
	}
	ONE = mm_set1_epi16(1);
	WX = roundup_times(W, WORDSIZE / 4) / (WORDSIZE / 4);
	if(wscores[0] == NULL){
		wscores[0] = calloc(WX * WORDSIZE / 4, sizeof(int));
	}
	if(wscores[1] == NULL){
		wscores[1] = calloc(WX * WORDSIZE / 4, sizeof(int));
	}
	for(k=0;k<W;k+=256){
		x = k;
		y = num_min(x + 256, W);
		for(i=0;i<2;i++){
			scr[i] = max[i] = mm_set1_epi16(0);
		}
		for(i=x;i<y;i++){
			h = mm_load(((xint*)us) + i);
			{
				m = mm_cvtepi8lo_epi16(h);
				scr[0] = mm_adds_epi16(scr[0], m);
				max[0] = mm_max_epi16(max[0], scr[0]);
				m = mm_cvtepi8hi_epi16(h);
				scr[1] = mm_adds_epi16(scr[1], m);
				max[1] = mm_max_epi16(max[1], scr[1]);
			}
			wscores[0][i] = (b2i)mm_extract_epi16(scr[0], 0);
			wscores[1][i] = (b2i)mm_extract_epi16(scr[1], WORDSIZE / 2 - 1);
		}
		num[0] = mm_set1_epi32(mm_extract_epi32(Scr[0], 0));
		num[1] = mm_set1_epi32(mm_extract_epi32(Scr[3], WORDSIZE / 4 - 1));
		for(i=x;i<y;i+=WORDSIZE/4){
			h = mm_load((xint*)(wscores[0] + i));
			h = mm_add_epi32(h, num[0]);
			mm_store((xint*)(wscores[0] + i), h);
			h = mm_load((xint*)(wscores[1] + i));
			h = mm_add_epi32(h, num[1]);
			mm_store((xint*)(wscores[1] + i), h);
		}
		h = mm_set1_epi32(k);
		for(i=0;i<4;i++){
			j = i / 2;
			// value of max score
			m = mm_cvtepi16lo_epi32(max[j]);
			max[j] = mm_srli(max[j], WORDSIZE / 2);
			m = mm_add_epi32(m, Scr[i]);
			c = mm_cmpgt_epi32(m, Max[i]);
			Max[i] = mm_blendv_epi8(Max[i], m, c);
			Idx[i] = mm_blendv_epi8(Idx[i], h, c);
			// add scores
			m = mm_cvtepi16lo_epi32(scr[j]);
			Scr[i] = mm_add_epi32(Scr[i], m);
			scr[j] = mm_srli(scr[j], WORDSIZE / 2);
		}
	}
	mm_store(((xint*)ary) + 0, Scr[0]);
	mm_store(((xint*)ary) + 1, Scr[1]);
	mm_store(((xint*)ary) + 2, Scr[2]);
	mm_store(((xint*)ary) + 3, Scr[3]);
	score = shift;
	for(i=0;i<WORDSIZE;i++){
		tmp = ary[i];
		ary[i] = score;
		score += tmp;
	}
	num[0] = mm_set1_epi32(ary[0]);
	num[1] = mm_set1_epi32(ary[WORDSIZE - 1]);
	for(i=0;i<WX;i++){
		h = mm_load(((xint*)wscores[0]) + i);
		h = mm_add_epi32(h, num[0]);
		mm_store(((xint*)wscores[0]) + i, h);
		h = mm_load(((xint*)wscores[1]) + i);
		h = mm_add_epi32(h, num[1]);
		mm_store(((xint*)wscores[1]) + i, h);
	}
	Max[0] = mm_add_epi32(mm_load(((xint*)ary) + 0), Max[0]);
	Max[1] = mm_add_epi32(mm_load(((xint*)ary) + 1), Max[1]);
	Max[2] = mm_add_epi32(mm_load(((xint*)ary) + 2), Max[2]);
	Max[3] = mm_add_epi32(mm_load(((xint*)ary) + 3), Max[3]);
	mm_store(((xint*)ary) + 0, Max[0]);
	mm_store(((xint*)ary) + 1, Max[1]);
	mm_store(((xint*)ary) + 2, Max[2]);
	mm_store(((xint*)ary) + 3, Max[3]);
	x = 0;
	for(i=1;i<WORDSIZE;i++){
		if(ary[i] > ary[x]){
			x = i;
		}
	}
	*max_score = ary[x];
	mm_store(((xint*)ary) + 0, Idx[x / (WORDSIZE / 4)]);
	y = ary[x % (WORDSIZE / 4)];
	// find the index of max score in <= 256 cells
	j = num_min(y + 256, W);
	tmp = score = us[y * WORDSIZE + x];
	k = y;
	for(i=y+1;i<j;i++){
		tmp += us[i * WORDSIZE + x];
		if(tmp > score){
			score = tmp;
			k = i;
		}
	}
	return x * W + k; // convert striped into normal coordinate
}

static inline void banded_striped_epi8_seqalign_row_look(b1i *us, u4i W, int shift, int *wscores[2]){
	xint h, m, Scr[4], scr[2], num[2];
	u4i k, i, j, x, y, WX;
	int score, tmp, ary[WORDSIZE];
	// find max
	for(i=0;i<4;i++){
		Scr[i] = mm_set1_epi32(0);
	}
	WX = roundup_times(W, WORDSIZE / 4) / (WORDSIZE / 4);
	if(wscores[0] == NULL){
		wscores[0] = calloc(WX * WORDSIZE / 4, sizeof(int));
	}
	if(wscores[1] == NULL){
		wscores[1] = calloc(WX * WORDSIZE / 4, sizeof(int));
	}
	for(k=0;k<W;k+=256){
		x = k;
		y = num_min(x + 256, W);
		for(i=0;i<2;i++){
			scr[i] = mm_set1_epi16(0);
		}
		for(i=x;i<y;i++){
			h = mm_load(((xint*)us) + i);
			{
				m = mm_cvtepi8lo_epi16(h);
				scr[0] = mm_adds_epi16(scr[0], m);
				m = mm_cvtepi8hi_epi16(h);
				scr[1] = mm_adds_epi16(scr[1], m);
			}
			wscores[0][i] = (b2i)mm_extract_epi16(scr[0], 0);
			wscores[1][i] = (b2i)mm_extract_epi16(scr[1], WORDSIZE / 2 - 1);
		}
		num[0] = mm_set1_epi32(mm_extract_epi32(Scr[0], 0));
		num[1] = mm_set1_epi32(mm_extract_epi32(Scr[3], WORDSIZE / 4 - 1));
		for(i=x;i<y;i+=WORDSIZE/4){
			h = mm_load((xint*)(wscores[0] + i));
			h = mm_add_epi32(h, num[0]);
			mm_store((xint*)(wscores[0] + i), h);
			h = mm_load((xint*)(wscores[1] + i));
			h = mm_add_epi32(h, num[1]);
			mm_store((xint*)(wscores[1] + i), h);
		}
		h = mm_set1_epi32(k);
		for(i=0;i<4;i++){
			j = i / 2;
			// add scores
			m = mm_cvtepi16lo_epi32(scr[j]);
			Scr[i] = mm_add_epi32(Scr[i], m);
			scr[j] = mm_srli(scr[j], WORDSIZE / 2);
		}
	}
	mm_store(((xint*)ary) + 0, Scr[0]);
	mm_store(((xint*)ary) + 1, Scr[1]);
	mm_store(((xint*)ary) + 2, Scr[2]);
	mm_store(((xint*)ary) + 3, Scr[3]);
	score = shift;
	for(i=0;i<WORDSIZE;i++){
		tmp = ary[i];
		ary[i] = score;
		score += tmp;
	}
	num[0] = mm_set1_epi32(ary[0]);
	num[1] = mm_set1_epi32(ary[WORDSIZE - 1]);
	for(i=0;i<WX;i++){
		h = mm_load(((xint*)wscores[0]) + i);
		h = mm_add_epi32(h, num[0]);
		mm_store(((xint*)wscores[0]) + i, h);
		h = mm_load(((xint*)wscores[1]) + i);
		h = mm_add_epi32(h, num[1]);
		mm_store(((xint*)wscores[1]) + i, h);
	}
}

static inline void banded_striped_epi8_seqalign_row_print(FILE *out, u1i *qseq, u4i tidx, u4i tpos, u1i tbase, u4i bandwidth, u4i mov, u4i rbeg, u4i rmax, int max_score, int shift, b1i *us, b1i *bs, int wl, int wr, int detail){
	u4i i, x, W;
	int score;
	W = bandwidth / WORDSIZE;
	fprintf(out, "ROW[%d][%d][%c]\tMOV=%d\tBAND=%d,%d\tMAX=%d(%d),%d\tSHIFT=%d\t\tWLR=%d,%d", tidx, tpos, "ACGTN-"[tbase], mov, rbeg, rbeg + bandwidth, rbeg + rmax, rmax, max_score, shift, wl, wr);
	if(detail){
		score = shift;
		for(i=0;i<bandwidth;i++){
			x = banded_striped_epi8_pos2idx(W << WORDSHIFT, i);
			fprintf(out, "\t%d:%c%d:%d:%d", i + rbeg, "ACGTN-"[qseq[rbeg + i]], score + us[x], us[x], bs[x] * 0xf);
			score += us[x];
		}
	}
	fprintf(out, "\n");
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

static inline u4i banded_striped_epi8_seqalign_cigar2alnstr(u1i *qseq, u1i *tseq, seqalign_result_t *rs, u4v *cigars, char *alnstr[3], u4i length){
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
	banded_striped_epi8_seqalign_piecex_row_max_func    row_max;
	banded_striped_epi8_seqalign_piecex_row_look_func   row_look;
	banded_striped_epi8_seqalign_piecex_backtrace_func  backtrace;
	seqalign_result_t rs;
	b1i *memp, *mempb, *qprof, *us[2], *es[2], *qs[2], *bs, *bs0;
	u8i mpsize;
	u4i i, W, rbeg, rmax, mov;
	int piecewise, smax, smin, lst_score, score, ph, rh, hh, wl, wr, wt, rbx, rby, rbz, *begs, *wscores[2];
	u1i tbase;
	set_qprof  = banded_striped_epi8_seqalign_set_query_prof;
	row_init   = banded_striped_epi8_seqalign_piecex_row_init;
	row_verify = banded_striped_epi8_seqalign_piecex_row_verify;
	row_max    = banded_striped_epi8_seqalign_row_max;
	row_look   = banded_striped_epi8_seqalign_row_look;
	backtrace  = banded_striped_epi8_seqalign_piecex_backtrace;
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		piecewise = 2;
		row_mov = banded_striped_epi8_seqalign_piece2_row_mov;
		row_cal = banded_striped_epi8_seqalign_piece2_row_cal;
	} else if(gapo1){
		piecewise = 1;
		row_mov = banded_striped_epi8_seqalign_piece1_row_mov;
		row_cal = banded_striped_epi8_seqalign_piece1_row_cal;
	} else {
		piecewise = 0;
		row_mov = banded_striped_epi8_seqalign_piece0_row_mov;
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
	mpsize += roundup_times(W * sizeof(int), WORDSIZE) * 2; // wscores[0], wscores[1]
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
	wscores[0] = (int*)memp; memp += roundup_times(W * sizeof(int), WORDSIZE);
	wscores[1] = (int*)memp; memp += roundup_times(W * sizeof(int), WORDSIZE);
	// prepare
	set_qprof(qseq, qlen, qprof, bandwidth, matrix);
	memset(&rs, 0, sizeof(seqalign_result_t));
	rs.score = SEQALIGN_SCORE_MIN;
	row_init(us, es, qs, begs, mode, bandwidth, smax, smin, gapo1, gape1, gapo2, gape2);
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
			ph = wscores[0][mov - 1];
			rh = ph + matrix[tbase * 4 + qseq[rbeg]]; // just guess
		} else {
			mov  = 0;
			rh = (rbeg == 0 && (mode == SEQALIGN_MODE_OVERLAP || i == 0))? 0 : SEQALIGN_SCORE_MIN;
		}
		row_mov(us, es, qs, W, mov); // mov and swap
		hh = row_cal(rbeg, tbase, us, es, qs, bs, qprof, gapo1, gape1, gapo2, gape2, W, mov, ph, rh);
		if(verbose > 2){
			row_verify(rbeg, W, ph, rh, hh, tbase, us, es, qs, qprof, gapo1, gape1, gapo2, gape2);
		}
		if(mode == SEQALIGN_MODE_GLOBAL && verbose == 0){
			rmax = 0;
			lst_score = SEQALIGN_SCORE_MIN;
			row_look(us[1], W, hh, wscores);
		} else {
			rmax = row_max(us[1], W, hh, wscores, &lst_score);
		}
		// adaptive banded
		wl = array_epi32_sum(wscores[0], W);
		wr = array_epi32_sum(wscores[1], W);
		wt = W * smax;
		if(verbose){
			banded_striped_epi8_seqalign_row_print(stdout, qseq, 1, i, tbase, bandwidth, mov, rbeg, rmax, lst_score, hh, us[1], bs, wl, wr, verbose > 1);
		}
		bs += bandwidth;
		if(i < (WORDSIZE / 2) * W){
			rbx = 0;
		} else if(wl + wt < wr){
			rbx = 2;
		} else if(wl > wr + wt){
			rbx = 0;
		} else {
			rbx = 1;
		}
		if(mode == SEQALIGN_MODE_GLOBAL){
			rbz = 2 * num_max(Int(tlen / qlen), 1); // suggested max step
			rby = Int((1.0 * i / tlen) * qlen); // diagonal line
			if(rbeg + rbz * (tlen - i - 1) + bandwidth <= (qlen + rbz - 1)){ // be quick to move to end
				mov = 1 + ((qlen - (rbeg + bandwidth)) / num_max(1, tlen - i - 1));
			} else if(rbeg < rby - bandwidth){
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
				score = wscores[1][qlen - 1 - rbeg - (WORDSIZE - 1) * W];
				if(score > rs.score){
					rs.score = score;
					rs.qe = qlen - 1;
					rs.te = i;
				}
			}
		}
	}
	if(mode == SEQALIGN_MODE_GLOBAL){
		rs.score = wscores[1][qlen - 1 - rbeg - (WORDSIZE - 1) * W];
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
