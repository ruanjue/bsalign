#ifndef BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H
#define BAND_STRIPED_DNA_SEQ_ALIGNMENT_RJ_H

/**
 * References:
 * Farrar, Michael. 2007. ¡°Striped Smith-Waterman Speeds Database Searches Six Times over Other SIMD Implementations.¡± Bioinformatics 23 (2): 156¨C61. https://doi.org/10.1093/bioinformatics/btl5
 * Suzuki, Hajime, and Masahiro Kasahara. 2017. ¡°Acceleration of Nucleotide Semi-Global Alignment with Adaptive Banded Dynamic Programming.¡± BioRxiv, September, 130633. https://doi.org/10.1101/13063
 * Suzuki, Hajime, and Masahiro Kasahara. 2018. ¡°Introducing Difference Recurrence Relations for Faster Semi-Global Alignment of Long Sequences.¡± BMC Bioinformatics 19 (Suppl 1). https://doi.org/10.1186/s12859-018-2018.
 * Li, Heng. 2018. ¡°Minimap2: Pairwise Alignment for Nucleotide Sequences.¡± Bioinformatics 34 (18): 3094¨C3100. https://doi.org/10.1093/bioinformatics/b
 *
 * My Algorithm = global overlap alignment + striped vectorization + difference recurrence relation + adaptive banded + active F-loop
 * 
 * I compute the score matrix in the way of one row by one row, that is y always increases 1 in next call.
 * x is resorted into striped vector based on W and B.
 * B: number of values in a SIMD word, W = bandwidth / B.
 * When B = 4 and W = 2, there have
 * normal array     [0, 1, 2, 3, 4, 5, 6, 7]
 * striped blocks   [0, 2, 4, 6; 1, 3, 5, 7]
 * running blocks [0, 1; 2, 3; 4, 5; 6, 7]
 * In implementation, I use __m128i to store 16 int8_t, B = 16.
 * H, E, F, Q, G are the absolute scores
 * e, f, u, q, g are the relative scores
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


#include "chararray.h"
#include "dna.h"
#include "list.h"
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>

#define SEQALIGN_BT_M	0
#define SEQALIGN_BT_D1	1
#define SEQALIGN_BT_D2	2
#define SEQALIGN_BT_I1	3
#define SEQALIGN_BT_I2	4
#define SEQALIGN_BT_DE1	8
#define SEQALIGN_BT_DE2	16
#define SEQALIGN_BT_IE1	32
#define SEQALIGN_BT_IE2	64

// Please DO NOT set M/X/O/E to be bigger than 63
#define SEQALIGN_SCORE_EPI8_MIN	(-(MAX_B1 >> 1))
#define SEQALIGN_SCORE_EPI8_MAX	(MAX_B1 >> 1)
#define SEQALIGN_SCORE_MIN	(-(MAX_B4 >> 1))
#define SEQALIGN_SCORE_MAX	(MAX_B4 >> 1)

typedef struct {
	int score;
	int qb, qe;
	int tb, te;
	int mat, mis, ins, del;
} seqalign_result_t;

/**
 * Basic function referings for global DNA sequence alignment, please their implementations within this file
 */

#define banded_striped_epi8_pos2idx(bandwidth, pos) ((((pos) % (bandwidth >> 4)) << 4) + ((pos) / ((bandwidth) >> 4)))

static inline void banded_striped_epi8_seqalign_set_score_matrix(b1i matrix[16], b1i mat, b1i mis){ u4i i; for(i=0;i<16;i++) matrix[i] = ((i ^ (i >> 2)) & 0x3)? mis : mat; }

// prepare query profile
typedef void (*banded_striped_epi8_seqalign_set_query_prof_func)(BaseBank *seqs, u8i qoff, u4i qlen, b1v *qprof, u4i bandwidth, b1i score_matrix[16]);

// mov <= W
// converting from [1] to [0]
// make the two row aligned
// here, us, es, and qs are array of pointers, aims to facilitate the memory management for applications
typedef void (*banded_striped_epi8_seqalign_piecex_row_mov_func)(b1i **us[2], b1i **es[2], b1i **qs[2], u4i W, u4i mov);

// core func to update row scores, from [0] to [1]
// ph is the score of H(-1, y - 1)
// rh is the revised score of H(-1, y - 1)
// @return: score of H(-1, y), note that H(-1, y) is useful to restore all scores of row
typedef int (*banded_striped_epi8_seqalign_piecex_row_cal_func)(u4i rbeg, u1i base, b1i **us[2], b1i **es[2], b1i **qs[2], b1i *bs, b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, int ph, int rh);

typedef void (*banded_striped_epi8_seqalign_piecex_row_verify_func)(int rbeg, int W, int ph, int rh, int hh, u1i base, b1i **us[2], b1i **es[2], b1i **qs[2], b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2);

// find max score and return the real offset in row
// wscores[0] store the first W real scores in normal cell order
// wscores[1] store the last  W real scores in normal cell order
typedef u4i (*banded_striped_epi8_seqalign_piecex_row_max_func)(b1i **us, u4i W, int hh, int *wscores[2], int *max_score);

// backtrace
// rs->qe, rs->te and rs->score MUST be set before call this function
// begs provides the band's offset of rows, it is continously suming of mov in banded_striped_epi8_seqalign_piecex_row_mov_func
// piecewise = 0/1/2
typedef void (*banded_striped_epi8_seqalign_piecex_backtrace_func)(BaseBank *seqs, u8i qoff, u4i qlen, u8i toff, u4i tlen, b1v *btds, u4v *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, String **alnstr);

// implementation of overlap alignment for two sequences
// bandwidth should be times of 16
static inline seqalign_result_t banded_striped_epi8_seqalign_pairwise_overlap(BaseBank *seqs, u8i qoff, u4i qlen, u8i toff, u4i tlen, b1v *qprof, b1v *rows, b1v *btds, u4v *begs, u4i bandwidth, b1i matrix[16], b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, String **alnstr, int verbose);


#define my_print_epi32X4(title, v)	\
{	\
	int vals[16], iter;	\
	_mm_store_si128(((__m128i*)vals) + 0, (v)[0]);	\
	_mm_store_si128(((__m128i*)vals) + 1, (v)[1]);	\
	_mm_store_si128(((__m128i*)vals) + 2, (v)[2]);	\
	_mm_store_si128(((__m128i*)vals) + 3, (v)[3]);	\
	fprintf(stdout, TOSTR(title));	\
	for(iter=0;iter<16;iter++){	\
		fprintf(stdout, "\t%d", vals[iter]);	\
	}	\
	fprintf(stdout, "\n");	\
}

#define my_print_epi32(title, v)	\
{	\
	int vals[4], iter;	\
	_mm_store_si128(((__m128i*)vals), v);	\
	fprintf(stdout, TOSTR(title));	\
	for(iter=0;iter<4;iter++){	\
		fprintf(stdout, "\t%d", vals[iter]);	\
	}	\
	fprintf(stdout, "\n");	\
}

#define my_print_epi16X2(title, v)	\
{	\
	b2i vals[16], iter;	\
	_mm_store_si128(((__m128i*)vals) + 0, (v)[0]);	\
	_mm_store_si128(((__m128i*)vals) + 1, (v)[1]);	\
	fprintf(stdout, TOSTR(title));	\
	for(iter=0;iter<16;iter++){	\
		fprintf(stdout, "\t%d", vals[iter]);	\
	}	\
	fprintf(stdout, "\n");	\
}

#define my_print_epi16(title, v)	\
{	\
	b2i vals[8], iter;	\
	_mm_store_si128(((__m128i*)vals), v);	\
	fprintf(stdout, TOSTR(title));	\
	for(iter=0;iter<16;iter++){	\
		fprintf(stdout, "\t%d", vals[iter]);	\
	}	\
	fprintf(stdout, "\n");	\
}

#define my_print_epi8(title, v)	\
{	\
	b1i vals[16], iter;	\
	_mm_store_si128(((__m128i*)vals), v);	\
	fprintf(stdout, TOSTR(title));	\
	for(iter=0;iter<16;iter++){	\
		fprintf(stdout, "\t%d", vals[iter]);	\
	}	\
	fprintf(stdout, "\n");	\
}

static inline void banded_striped_epi8_seqalign_piece0_row_mov(b1i **us[2], b1i **es[2], b1i **qs[2], u4i W, u4i mov){
	__m128i u;
	u4i i, div;
	UNUSED(es);
	UNUSED(qs);
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = _mm_load_si128((__m128i*)us[1][i + mov]);
		_mm_store_si128((__m128i*)us[0][i], u);
	}
	if(!mov) return;
	{
		u = _mm_load_si128((__m128i*)us[1][i + mov - W]);
		u = _mm_srli_si128(u, 1);
		u = _mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, 15);
	}
	for(;i<W;i++){
		_mm_store_si128((__m128i*)us[0][i], u);
		u = _mm_load_si128((__m128i*)us[1][(i + mov - W + 1) % W]);
		u = _mm_srli_si128(u, 1);
	}
}

static inline void banded_striped_epi8_seqalign_piece1_row_mov(b1i **us[2], b1i **es[2], b1i **qs[2], u4i W, u4i mov){
	__m128i u, e;
	u4i i, div;
	UNUSED(qs);
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = _mm_load_si128((__m128i*)us[1][i + mov]);
		e = _mm_load_si128((__m128i*)es[1][i + mov]);
		_mm_store_si128((__m128i*)us[0][i], u);
		_mm_store_si128((__m128i*)es[0][i], e);
	}
	if(!mov) return;
	{
		u = _mm_load_si128((__m128i*)us[1][i + mov - W]);
		u = _mm_srli_si128(u, 1);
		u = _mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, 15);
	}
	for(;i<W;i++){
		e = _mm_load_si128((__m128i*)es[1][i + mov - W]);
		e = _mm_srli_si128(e, 1);
		_mm_store_si128((__m128i*)us[0][i], u);
		_mm_store_si128((__m128i*)es[0][i], e);
		u = _mm_load_si128((__m128i*)us[1][(i + mov - W + 1) % W]);
		u = _mm_srli_si128(u, 1);
	}
}

static inline void banded_striped_epi8_seqalign_piece2_row_mov(b1i **us[2], b1i **es[2], b1i **qs[2], u4i W, u4i mov){
	__m128i u, e, x;
	u4i i, div;
	if(mov > W){
		fprintf(stderr, " -- mov(%d) > W(%d) in %s -- %s:%d --\n", mov, W, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	div = W - mov;
	for(i=0;i<div;i++){
		u = _mm_load_si128((__m128i*)us[1][i + mov]);
		e = _mm_load_si128((__m128i*)es[1][i + mov]);
		x = _mm_load_si128((__m128i*)qs[1][i + mov]);
		_mm_store_si128((__m128i*)us[0][i], u);
		_mm_store_si128((__m128i*)es[0][i], e);
		_mm_store_si128((__m128i*)qs[0][i], x);
	}
	if(!mov) return;
	{
		u = _mm_load_si128((__m128i*)us[1][i + mov - W]);
		u = _mm_srli_si128(u, 1);
		u = _mm_insert_epi8(u, SEQALIGN_SCORE_EPI8_MIN, 15);
	}
	for(;i<W;i++){
		e = _mm_load_si128((__m128i*)es[1][i + mov - W]);
		x = _mm_load_si128((__m128i*)qs[1][i + mov - W]);
		e = _mm_srli_si128(e, 1);
		x = _mm_srli_si128(x, 1);
		_mm_store_si128((__m128i*)us[0][i], u);
		_mm_store_si128((__m128i*)es[0][i], e);
		_mm_store_si128((__m128i*)qs[0][i], x);
		u = _mm_load_si128((__m128i*)us[1][(i + mov - W + 1) % W]);
		u = _mm_srli_si128(u, 1);
	}
}

static inline int banded_striped_epi8_seqalign_piece0_row_cal(u4i rbeg, u1i base, b1i **us[2], b1i **es[2], b1i **qs[2], b1i *bs, b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, int ph, int rh){
	__m128i h, e, f, u, v, c, d, m2[2], m4[4];
	__m128i I, D, GapE;
	int ms[16];
	b1i vs[16];
	u4i i, k;
	int s, t, h0;
	UNUSED(es);
	UNUSED(qs);
	UNUSED(gapo1);
	UNUSED(gapo2);
	UNUSED(gape2);
	I     = _mm_set1_epi8(SEQALIGN_BT_I1);
	D     = _mm_set1_epi8(SEQALIGN_BT_D1);
	GapE  = _mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	m2[0] = _mm_set1_epi16(0);
	m2[1] = _mm_set1_epi16(0);
	m4[0] = _mm_set1_epi32(0);
	m4[1] = _mm_set1_epi32(0);
	m4[2] = _mm_set1_epi32(0);
	m4[3] = _mm_set1_epi32(0);
	f = _mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof->buffer[((rbeg + 0) * 4 + base) * 16]; // h
	t = us[0][0][0] + gape1; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
	}
	// preparing initial f for each running block
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;){
		k = num_min(i + 256, W);
		for(;i<k;i++){
			u = _mm_load_si128((__m128i*)us[0][i]);
			// sums us[0] within each running block
			m2[0] = _mm_subs_epi16(m2[0], _mm_cvtepi8_epi16(u));
			m2[1] = _mm_subs_epi16(m2[1], _mm_cvtepi8_epi16(_mm_srli_si128(u, 8)));
			// max h, e, f
			e = _mm_adds_epi8(u, GapE);
			h = _mm_max_epi8(e, h);
			h = _mm_max_epi8(f, h);
			// preparing next f
			f = _mm_adds_epi8(h, GapE);
			f = _mm_subs_epi8(f, u);
			// preparing next h
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
		{
			m4[0] = _mm_add_epi32(m4[0], _mm_cvtepi16_epi32(m2[0]));
			m4[1] = _mm_add_epi32(m4[1], _mm_cvtepi16_epi32(_mm_srli_si128(m2[0], 8)));
			m4[2] = _mm_add_epi32(m4[2], _mm_cvtepi16_epi32(m2[1]));
			m4[3] = _mm_add_epi32(m4[3], _mm_cvtepi16_epi32(_mm_srli_si128(m2[1], 8)));
			m2[0] = _mm_set1_epi16(0);
			m2[1] = _mm_set1_epi16(0);
		}
	}
	_mm_store_si128(((__m128i*)ms) + 0, m4[0]);
	_mm_store_si128(((__m128i*)ms) + 1, m4[1]);
	_mm_store_si128(((__m128i*)ms) + 2, m4[2]);
	_mm_store_si128(((__m128i*)ms) + 3, m4[3]);
	f = _mm_slli_si128(f, 1);
	_mm_store_si128((__m128i*)vs, f);
	vs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + vs[0] + ms[0];
	for(i=1;i<16;i++){
		if(vs[i] < s) vs[i] = s;
		s = t + vs[i] + ms[i];
	}
	f = _mm_load_si128((__m128i*)vs);
	// main loop
	// don't use the h from last W - 1 to calculate v = h - u, because this h maybe updated when F-penetration
	// will revise v after this loop
	v = _mm_set1_epi8(0);
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	u = _mm_set1_epi8(0); // useless, but for compiler
	for(i=0;i<W;i++){
		u = _mm_load_si128((__m128i*)us[0][i]);
		// max(e, h)
		e = _mm_adds_epi8(u, GapE);
		c = _mm_cmpgt_epi8(e, h);
		d = _mm_and_si128(c, D); // bt
		h = _mm_max_epi8(e, h);
		// max(f, h)
		c = _mm_cmpgt_epi8(f, h);
		d = _mm_blendv_epi8(d, I, c);
		h = _mm_max_epi8(f, h);
		_mm_store_si128(((__m128i*)bs) + i, d);
		// calculate u(x, y)
		v = _mm_subs_epi8(h, v);
		_mm_store_si128((__m128i*)us[1][i], v);
		v = _mm_subs_epi8(h, u);
		// calculate f(x, y)
		f = _mm_adds_epi8(h, GapE);
		f = _mm_subs_epi8(f, u);
		if(__builtin_expect(i + 1 == W, 0)){
		} else {
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
	}
	// revise the first striped block of u
	v = _mm_subs_epi8(h, u);
	v = _mm_slli_si128(v, 1);
	u = _mm_load_si128((__m128i*)us[1][0]);
	u = _mm_subs_epi8(u, v);
	_mm_store_si128((__m128i*)us[1][0], u);
	// shift score to fit EPI8
	rh = ph + us[1][0][0];
	us[1][0][0] = 0;
	return rh;
}

static inline int banded_striped_epi8_seqalign_piece1_row_cal(u4i rbeg, u1i base, b1i **us[2], b1i **es[2], b1i **qs[2], b1i *bs, b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, int ph, int rh){
	__m128i h, e, f, u, v, c, d, m2[2], m4[4];
	__m128i I, IE, D, DE, GapOE, GapE;
	int ms[16];
	b1i vs[16];
	u4i i, k;
	int s, t, h0;
	UNUSED(qs);
	UNUSED(gapo2);
	UNUSED(gape2);
	I     = _mm_set1_epi8(SEQALIGN_BT_I1);
	IE    = _mm_set1_epi8(SEQALIGN_BT_IE1);
	D     = _mm_set1_epi8(SEQALIGN_BT_D1);
	DE    = _mm_set1_epi8(SEQALIGN_BT_DE1);
	GapOE = _mm_set1_epi8(gapo1 + gape1);
	GapE  = _mm_set1_epi8(gape1);
	// ::: max(h, e, f)
	m2[0] = _mm_set1_epi16(0);
	m2[1] = _mm_set1_epi16(0);
	m4[0] = _mm_set1_epi32(0);
	m4[1] = _mm_set1_epi32(0);
	m4[2] = _mm_set1_epi32(0);
	m4[3] = _mm_set1_epi32(0);
	f = _mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof->buffer[((rbeg + 0) * 4 + base) * 16]; // h
	t = us[0][0][0] + es[0][0][0]; // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
	}
	// preparing initial f for each running block
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;){
		k = num_min(i + 256, W);
		for(;i<k;i++){
			u = _mm_load_si128((__m128i*)us[0][i]);
			// sums us[0] within each running block
			m2[0] = _mm_subs_epi16(m2[0], _mm_cvtepi8_epi16(u));
			m2[1] = _mm_subs_epi16(m2[1], _mm_cvtepi8_epi16(_mm_srli_si128(u, 8)));
			// max h, e, f
			e = _mm_load_si128((__m128i*)es[0][i]);
			e = _mm_adds_epi8(e, u);
			h = _mm_max_epi8(e, h);
			h = _mm_max_epi8(f, h);
			// preparing next f
			f = _mm_adds_epi8(f, GapE);
			h = _mm_adds_epi8(h, GapOE);
			f = _mm_max_epi8(f, h);
			f = _mm_subs_epi8(f, u);
			// preparing next h
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
		{
			m4[0] = _mm_add_epi32(m4[0], _mm_cvtepi16_epi32(m2[0]));
			m4[1] = _mm_add_epi32(m4[1], _mm_cvtepi16_epi32(_mm_srli_si128(m2[0], 8)));
			m4[2] = _mm_add_epi32(m4[2], _mm_cvtepi16_epi32(m2[1]));
			m4[3] = _mm_add_epi32(m4[3], _mm_cvtepi16_epi32(_mm_srli_si128(m2[1], 8)));
			m2[0] = _mm_set1_epi16(0);
			m2[1] = _mm_set1_epi16(0);
		}
	}
	_mm_store_si128(((__m128i*)ms) + 0, m4[0]);
	_mm_store_si128(((__m128i*)ms) + 1, m4[1]);
	_mm_store_si128(((__m128i*)ms) + 2, m4[2]);
	_mm_store_si128(((__m128i*)ms) + 3, m4[3]);
	f = _mm_slli_si128(f, 1);
	_mm_store_si128((__m128i*)vs, f);
	vs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + vs[0] + ms[0];
	for(i=1;i<16;i++){
		if(vs[i] < s) vs[i] = s;
		s = t + vs[i] + ms[i];
	}
	f = _mm_load_si128((__m128i*)vs);
	// main loop
	// don't use the h from last W - 1 to calculate v = h - u, because this h maybe updated when F-penetration
	// will revise v after this loop
	v = _mm_set1_epi8(0);
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	u = _mm_set1_epi8(0); // useless, but for compiler
	for(i=0;i<W;i++){
		u = _mm_load_si128((__m128i*)us[0][i]);
		e = _mm_load_si128((__m128i*)es[0][i]);
		// max(e, h)
		e = _mm_adds_epi8(e, u);
		c = _mm_cmpgt_epi8(e, h);
		d = _mm_and_si128(c, D); // bt
		h = _mm_max_epi8(e, h);
		// max(f, h)
		c = _mm_cmpgt_epi8(f, h);
		d = _mm_blendv_epi8(d, I, c);
		h = _mm_max_epi8(f, h);
		// calculate u(x, y)
		v = _mm_subs_epi8(h, v);
		_mm_store_si128((__m128i*)us[1][i], v);
		v = _mm_subs_epi8(h, u);
		// calculate e(x, y)
		e = _mm_adds_epi8(e, GapE);
		e = _mm_subs_epi8(e, h);
		c = _mm_cmpgt_epi8(e, GapOE);
		d = _mm_or_si128(d, _mm_and_si128(c, DE));
		e = _mm_max_epi8(e, GapOE);
		_mm_store_si128((__m128i*)es[1][i], e);
		// calculate f(x, y)
		f = _mm_adds_epi8(f, GapE);
		h = _mm_adds_epi8(h, GapOE);
		c = _mm_cmpgt_epi8(f, h);
		d = _mm_or_si128(d, _mm_and_si128(c, IE));
		_mm_store_si128(((__m128i*)bs) + i, d);
		f = _mm_max_epi8(f, h);
		f = _mm_subs_epi8(f, u);
		if(__builtin_expect(i + 1 == W, 0)){
			h = _mm_subs_epi8(h, GapOE);
		} else {
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
	}
	// revise the first striped block of u
	v = _mm_subs_epi8(h, u);
	v = _mm_slli_si128(v, 1);
	u = _mm_load_si128((__m128i*)us[1][0]);
	u = _mm_subs_epi8(u, v);
	_mm_store_si128((__m128i*)us[1][0], u);
	// shift score to fit EPI8
	rh = ph + us[1][0][0];
	us[1][0][0] = 0;
	return rh;
}

static inline int banded_striped_epi8_seqalign_piece2_row_cal(u4i rbeg, u1i base, b1i **us[2], b1i **es[2], b1i **qs[2], b1i *bs, b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, u4i W, int ph, int rh){
	__m128i h, e, q, f, g, u, v, c, d, m2[2], m4[4];
	__m128i I1, I2, IE1, IE2, D1, D2, DE1, DE2, GapOE, GapE, GapQP, GapP, GapOQ;
	int ms[16];
	b1i vs[16];
	u4i i, k;
	int s, t, h0;
	I1     = _mm_set1_epi8(SEQALIGN_BT_I1);
	I2     = _mm_set1_epi8(SEQALIGN_BT_I2);
	IE1    = _mm_set1_epi8(SEQALIGN_BT_IE1);
	IE2    = _mm_set1_epi8(SEQALIGN_BT_IE2);
	D1     = _mm_set1_epi8(SEQALIGN_BT_D1);
	D2     = _mm_set1_epi8(SEQALIGN_BT_D2);
	DE1    = _mm_set1_epi8(SEQALIGN_BT_DE1);
	DE2    = _mm_set1_epi8(SEQALIGN_BT_DE2);
	GapOE = _mm_set1_epi8(gapo1 + gape1);
	GapE  = _mm_set1_epi8(gape1);
	GapQP = _mm_set1_epi8(gapo2 + gape2);
	GapP  = _mm_set1_epi8(gape2);
	GapOQ = _mm_subs_epi8(GapOE, GapQP);
	// ::: max(h, e, f)
	m2[0] = _mm_set1_epi16(0);
	m2[1] = _mm_set1_epi16(0);
	m4[0] = _mm_set1_epi32(0);
	m4[1] = _mm_set1_epi32(0);
	m4[2] = _mm_set1_epi32(0);
	m4[3] = _mm_set1_epi32(0);
	f = _mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	g = _mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	h0 = (rh - ph) + qprof->buffer[((rbeg + 0) * 4 + base) * 16]; // h
	t = us[0][0][0] + num_max(es[0][0][0], qs[0][0][0]); // e
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX; // score will loss, please never set rh - ph >= SEQALIGN_SCORE_EPI8_MAX - base_match_score
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
	}
	// preparing initial f for each running block
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	for(i=0;i<W;){
		k = num_min(i + 256, W);
		for(;i<k;i++){
			u = _mm_load_si128((__m128i*)us[0][i]);
			// sums us[0] within each running block
			m2[0] = _mm_subs_epi16(m2[0], _mm_cvtepi8_epi16(u));
			m2[1] = _mm_subs_epi16(m2[1], _mm_cvtepi8_epi16(_mm_srli_si128(u, 8)));
			// max h, e, q, f, g
			e = _mm_load_si128((__m128i*)es[0][i]);
			q = _mm_load_si128((__m128i*)qs[0][i]);
			e = _mm_adds_epi8(e, u);
			q = _mm_adds_epi8(q, u);
			h = _mm_max_epi8(e, h);
			h = _mm_max_epi8(q, h);
			h = _mm_max_epi8(f, h);
			h = _mm_max_epi8(g, h);
			// preparing next f and g
			f = _mm_adds_epi8(f, GapE);
			h = _mm_adds_epi8(h, GapOE);
			f = _mm_max_epi8(f, h);
			f = _mm_subs_epi8(f, u);
			g = _mm_adds_epi8(g, GapP);
			h = _mm_subs_epi8(h, GapOQ); // gapo1 + gape1 - gapo2 - gape2
			g = _mm_max_epi8(g, h);
			g = _mm_subs_epi8(g, u);
			// preparing next h
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
		{
			m4[0] = _mm_add_epi32(m4[0], _mm_cvtepi16_epi32(m2[0]));
			m4[1] = _mm_add_epi32(m4[1], _mm_cvtepi16_epi32(_mm_srli_si128(m2[0], 8)));
			m4[2] = _mm_add_epi32(m4[2], _mm_cvtepi16_epi32(m2[1]));
			m4[3] = _mm_add_epi32(m4[3], _mm_cvtepi16_epi32(_mm_srli_si128(m2[1], 8)));
			m2[0] = _mm_set1_epi16(0);
			m2[1] = _mm_set1_epi16(0);
		}
	}
	_mm_store_si128(((__m128i*)ms) + 0, m4[0]);
	_mm_store_si128(((__m128i*)ms) + 1, m4[1]);
	_mm_store_si128(((__m128i*)ms) + 2, m4[2]);
	_mm_store_si128(((__m128i*)ms) + 3, m4[3]);
	f = _mm_slli_si128(f, 1);
	_mm_store_si128((__m128i*)vs, f);
	vs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape1;
	s = t + vs[0] + ms[0];
	for(i=1;i<16;i++){
		if(vs[i] < s) vs[i] = s;
		s = t + vs[i] + ms[i];
	}
	f = _mm_load_si128((__m128i*)vs);
	g = _mm_slli_si128(g, 1);
	_mm_store_si128((__m128i*)vs, g);
	vs[0] = SEQALIGN_SCORE_EPI8_MIN;
	t = W * gape2;
	s = t + vs[0] + ms[0];
	for(i=1;i<16;i++){
		if(vs[i] < s){
			vs[i] = s;
		}
		s = t + vs[i] + ms[i];
	}
	g = _mm_load_si128((__m128i*)vs);
	// main loop
	v = _mm_set1_epi8(0);
	h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + 0) * 4 + base);
	h = _mm_insert_epi8(h, h0, 0);
	u = _mm_set1_epi8(0); // useless, but for compiler
	for(i=0;i<W;i++){
		u = _mm_load_si128((__m128i*)us[0][i]);
		e = _mm_load_si128((__m128i*)es[0][i]);
		q = _mm_load_si128((__m128i*)qs[0][i]);
		// max(e, q, h)
		e = _mm_adds_epi8(e, u);
		c = _mm_cmpgt_epi8(e, h);
		d = _mm_and_si128(c, D1); // bt
		h = _mm_max_epi8(e, h);
		q = _mm_adds_epi8(q, u);
		c = _mm_cmpgt_epi8(q, h);
		d = _mm_blendv_epi8(d, D2, c);
		h = _mm_max_epi8(q, h);
		// max(f, g, h)
		c = _mm_cmpgt_epi8(f, h);
		d = _mm_blendv_epi8(d, I1, c);
		h = _mm_max_epi8(f, h);
		c = _mm_cmpgt_epi8(g, h);
		d = _mm_blendv_epi8(d, I2, c);
		h = _mm_max_epi8(g, h);
		// calculate u(x, y)
		v = _mm_subs_epi8(h, v);
		_mm_store_si128((__m128i*)us[1][i], v);
		v = _mm_subs_epi8(h, u);
		// calculate e(x, y) and q(x, y)
		e = _mm_adds_epi8(e, GapE);
		e = _mm_subs_epi8(e, h);
		c = _mm_cmpgt_epi8(e, GapOE);
		d = _mm_or_si128(d, _mm_and_si128(c, DE1));
		e = _mm_max_epi8(e, GapOE);
		_mm_store_si128((__m128i*)es[1][i], e);
		q = _mm_adds_epi8(q, GapP);
		q = _mm_subs_epi8(q, h);
		c = _mm_cmpgt_epi8(q, GapQP);
		d = _mm_or_si128(d, _mm_and_si128(c, DE2));
		q = _mm_max_epi8(q, GapQP);
		_mm_store_si128((__m128i*)qs[1][i], q);
		// calculate f(x, y) and g(x, y)
		f = _mm_adds_epi8(f, GapE);
		h = _mm_adds_epi8(h, GapOE);
		c = _mm_cmpgt_epi8(f, h);
		d = _mm_or_si128(d, _mm_and_si128(c, IE1));
		f = _mm_max_epi8(f, h);
		f = _mm_subs_epi8(f, u);
		g = _mm_adds_epi8(g, GapP);
		h = _mm_subs_epi8(h, GapOQ); // (gapo1 + gape1 - gapo2 - gape2)
		c = _mm_cmpgt_epi8(g, h);
		d = _mm_or_si128(d, _mm_and_si128(c, IE2));
		g = _mm_max_epi8(g, h);
		g = _mm_subs_epi8(g, u);
		_mm_store_si128(((__m128i*)bs) + i, d);
		if(__builtin_expect(i + 1 == W, 0)){
			h = _mm_subs_epi8(h, GapQP);
		} else {
			h = _mm_load_si128(((__m128i*)qprof->buffer) + (rbeg + i + 1) * 4 + base);
		}
	}
	// revise the first striped block of u
	v = _mm_subs_epi8(h, u);
	v = _mm_slli_si128(v, 1);
	u = _mm_load_si128((__m128i*)us[1][0]);
	u = _mm_subs_epi8(u, v);
	_mm_store_si128((__m128i*)us[1][0], u);
	// shift score to fit EPI8
	rh = ph + us[1][0][0];
	us[1][0][0] = 0;
	return rh;
}

// array Wscores records absolute values of the first and last W cells in natural order
static inline u4i banded_striped_epi8_seqalign_row_max(b1i **us, u4i W, int shift, int *wscores[2], int *max_score){
	__m128i h, c, m, Max[4], max[2], Scr[4], scr[2], Idx[4], num[2], ONE;
	u4i k, i, j, x, y, WX;
	int score, tmp, ary[16];
	// find max
	for(i=0;i<4;i++){
		Idx[i] = Scr[i] = Max[i] = _mm_set1_epi32(0);
	}
	ONE = _mm_set1_epi16(1);
	WX = roundup_times(W, 4) / 4;
	if(wscores[0] == NULL){
		wscores[0] = calloc(WX * 4, sizeof(int));
	} else {
		//memset(wscores[0], 0, WX * 4 * sizeof(int));
	}
	if(wscores[1] == NULL){
		wscores[1] = calloc(WX * 4, sizeof(int));
	} else {
		//memset(wscores[1], 0, WX * 4 * sizeof(int));
	}
	for(k=0;k<W;k+=256){
		x = k;
		y = num_min(x + 256, W);
		for(i=0;i<2;i++){
			scr[i] = max[i] = _mm_set1_epi16(0);
		}
		for(i=x;i<y;i++){
			h = _mm_load_si128((__m128i*)us[i]);
			{
				m = _mm_cvtepi8_epi16(h);
				h = _mm_srli_si128(h, 8);
				scr[0] = _mm_adds_epi16(scr[0], m);
				max[0] = _mm_max_epi16(max[0], scr[0]);
				m = _mm_cvtepi8_epi16(h);
				scr[1] = _mm_adds_epi16(scr[1], m);
				max[1] = _mm_max_epi16(max[1], scr[1]);
			}
			wscores[0][i] = (b2i)_mm_extract_epi16(scr[0], 0);
			wscores[1][i] = (b2i)_mm_extract_epi16(scr[1], 7);
		}
		num[0] = _mm_set1_epi32(_mm_extract_epi32(Scr[0], 0));
		num[1] = _mm_set1_epi32(_mm_extract_epi32(Scr[3], 3));
		for(i=x;i<y;i+=4){
			h = _mm_load_si128((__m128i*)(wscores[0] + i));
			h = _mm_add_epi32(h, num[0]);
			_mm_store_si128((__m128i*)(wscores[0] + i), h);
			h = _mm_load_si128((__m128i*)(wscores[1] + i));
			h = _mm_add_epi32(h, num[1]);
			_mm_store_si128((__m128i*)(wscores[1] + i), h);
		}
		h = _mm_set1_epi32(k);
		for(i=0;i<4;i++){
			j = i / 2;
			// value of max score
			m = _mm_cvtepi16_epi32(max[j]);
			max[j] = _mm_srli_si128(max[j], 8);
			m = _mm_add_epi32(m, Scr[i]);
			c = _mm_cmpgt_epi32(m, Max[i]);
			Max[i] = _mm_blendv_epi8(Max[i], m, c);
			Idx[i] = _mm_blendv_epi8(Idx[i], h, c);
			// add scores
			m = _mm_cvtepi16_epi32(scr[j]);
			scr[j] = _mm_srli_si128(scr[j], 8);
			Scr[i] = _mm_add_epi32(Scr[i], m);
		}
	}
	_mm_store_si128(((__m128i*)ary) + 0, Scr[0]);
	_mm_store_si128(((__m128i*)ary) + 1, Scr[1]);
	_mm_store_si128(((__m128i*)ary) + 2, Scr[2]);
	_mm_store_si128(((__m128i*)ary) + 3, Scr[3]);
	score = shift;
	for(i=0;i<16;i++){
		tmp = ary[i];
		ary[i] = score;
		score += tmp;
	}
	num[0] = _mm_set1_epi32(ary[0]);
	num[1] = _mm_set1_epi32(ary[15]);
	for(i=0;i<WX;i++){
		h = _mm_load_si128(((__m128i*)wscores[0]) + i);
		h = _mm_add_epi32(h, num[0]);
		_mm_store_si128(((__m128i*)wscores[0]) + i, h);
		h = _mm_load_si128(((__m128i*)wscores[1]) + i);
		h = _mm_add_epi32(h, num[1]);
		_mm_store_si128(((__m128i*)wscores[1]) + i, h);
	}
	Scr[0] = _mm_load_si128(((__m128i*)ary) + 0);
	Scr[1] = _mm_load_si128(((__m128i*)ary) + 1);
	Scr[2] = _mm_load_si128(((__m128i*)ary) + 2);
	Scr[3] = _mm_load_si128(((__m128i*)ary) + 3);
	Max[0] = _mm_add_epi32(Scr[0], Max[0]);
	Max[1] = _mm_add_epi32(Scr[1], Max[1]);
	Max[2] = _mm_add_epi32(Scr[2], Max[2]);
	Max[3] = _mm_add_epi32(Scr[3], Max[3]);
	_mm_store_si128(((__m128i*)ary) + 0, Max[0]);
	_mm_store_si128(((__m128i*)ary) + 1, Max[1]);
	_mm_store_si128(((__m128i*)ary) + 2, Max[2]);
	_mm_store_si128(((__m128i*)ary) + 3, Max[3]);
	x = 0;
	for(i=1;i<16;i++){
		if(ary[i] > ary[x]){
			x = i;
		}
	}
	*max_score = ary[x];
	_mm_store_si128(((__m128i*)ary) + 0, Idx[x / 4]);
	y = ary[x % 4];
	// find the index of max score in <= 256 cells
	j = num_min(y + 256, W);
	tmp = score = us[y][x];
	k = y;
	for(i=y+1;i<j;i++){
		tmp += us[i][x];
		if(tmp > score){
			score = tmp;
			k = i;
		}
	}
	return x * W + k; // convert stripped into normal coordinate
}

static inline void banded_striped_epi8_seqalign_set_query_prof(BaseBank *seqs, u8i qoff, u4i qlen, b1v *qprof, u4i bandwidth, b1i mtx[16]){
	b1i *qp;
	u4i xlen, x, pos, i, j, W, b;
	W = bandwidth / 16;
	xlen = num_max(qlen, bandwidth); // In case of bandwidth > qlen
	clear_and_encap_b1v(qprof, ((xlen + 1) * 4) * 16);
	for(x=0;x<=xlen;x++){ // leave the last block all SEQALIGN_SCORE_EPI8_MIN for accessing [W]
		qp = qprof->buffer + x * 4 * 16;
		for(j=0;j<16;j++){
			pos = x + j * W;
			if(pos < qlen){
				b = get_basebank(seqs, qoff + pos);
				b *= 4;
				for(i=0;i<4;i++){
					qp[i * 16 + j] = mtx[b + i];
				}
			} else {
				for(i=0;i<4;i++){
					qp[i * 16 + j] = SEQALIGN_SCORE_EPI8_MIN;
				}
			}
		}
	}
}

#define banded_striped_epi8_seqalign_get_rd_query_prof(qprof, pos, base) (((base) < 4)? (qprof)->buffer + ((pos) * 4 + base) * 16 : (qprof)->buffer - 16)

static inline void banded_striped_epi8_seqalign_row_print(FILE *out, BaseBank *seqs, u4i tidx, u4i tpos, u1i tbase, u4i bandwidth, u4i mov, u8i qoff, u4i rbeg, u4i rmax, int max_score, int shift, b1i **us, b1i *bs, int wl, int wr, int detail){
	u4i i, x, y, W;
	int score;
	W = bandwidth / 16;
	fprintf(out, "ROW[%d][%d][%c]\tMOV=%d\tBAND=%d,%d\tMAX=%d(%d),%d\tSHIFT=%d\t\tWLR=%d,%d", tidx, tpos, bit_base_table[tbase], mov, rbeg, rbeg + bandwidth, rbeg + rmax, rmax, max_score, shift, wl, wr);
	if(detail){
		score = shift;
		for(i=0;i<bandwidth;i++){
			x = i % W;
			y = i / W;
			fprintf(out, "\t%d:%c%d:%d", i + rbeg, bit_base_table[get_basebank(seqs, qoff + rbeg + i)], score + us[x][y], bs[x * 16 + y] & 0x07);
			score += us[x][y];
		}
	}
	fprintf(out, "\n");
}

static inline int array_epi32_sum(int *rs, int n){
	__m128i v, s;
	int i, m, ary[4], sum;
	s = _mm_set1_epi32(0);
	m = n / 4;
	for(i=0;i<m;i++){
		v = _mm_load_si128(((__m128i*)rs) + i);
		s = _mm_add_epi32(s, v);
	}
	_mm_store_si128((__m128i*)ary, s);
	sum = ary[0] + ary[1] + ary[2] + ary[3];
	for(i=m*4;i<n;i++){
		sum += rs[i];
	}
	return sum;
}

static inline void banded_striped_epi8_seqalign_piecex_row_verify(int rbeg, int W, int ph, int rh, int hh, u1i tbase, b1i **us[2], b1i **es[2], b1i **qs[2], b1v *qprof, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2){
	int i, x, y, s1, s2, t, h, e, q, f, g, piecewise, h0;
	if(gapo2 < gapo1 && gape2 > gape1 && gapo2 + gape2 < gapo1 + gape1){
		piecewise = 2;
		t = us[0][0][0] + num_max(es[0][0][0], qs[0][0][0]);
	} else if(gapo1){
		piecewise = 1;
		t = us[0][0][0] + es[0][0][0];
	} else {
		piecewise = 0;
		t = us[0][0][0] + gape1;
	}
	f = g = SEQALIGN_SCORE_MIN;
	h0 = (rh - ph) + qprof->buffer[((rbeg + 0) * 4 + tbase) * 16];
	if(h0 >= t){
		if(h0 > SEQALIGN_SCORE_EPI8_MAX) h0 = SEQALIGN_SCORE_EPI8_MAX;
	} else {
		h0 = SEQALIGN_SCORE_EPI8_MIN;
	}
	s1 = ph;
	s2 = hh;
	//fprintf(stdout, " -- %s -- %s -- piecewise=%d rbeg=%d\tph=%d\trh=%d\thh=%d\n", __FUNCTION__, __FILE__, piecewise, rbeg, ph, rh, hh); fflush(stdout);
	for(i=0;i<W*16;i++){
		x = i % W;
		y = i / W;
		s2 += us[1][x][y];
		e = s1 + us[0][x][y] + (piecewise != 0? es[0][x][y] : gape1);
		q = s1 + us[0][x][y] + (piecewise == 2? qs[0][x][y] : SEQALIGN_SCORE_MIN);
		h = num_max(s1 + h0, e);
		h = num_max(h, q);
		h = num_max(h, f);
		h = num_max(h, g);
		if(h != s2){
			fprintf(stdout, " -- %c:%d i=%d h=%d(%d,%d,%d,%d) s2=%d s1=%d u0=%d u1=%d something wrong in %s -- %s:%d --\n", bit_base_table[tbase], h0, i, h, e, q, f, g, s2, s1, us[0][x][y], us[1][x][y], __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
			//abort();
			return;
		}
		f = num_max(f, h + gapo1) + gape1;
		g = (piecewise == 2)? num_max(g + gape2, h + gapo2 + gape2) : SEQALIGN_SCORE_MIN;
		s1 += us[0][x][y];
		h0 = banded_striped_epi8_seqalign_get_rd_query_prof(qprof, rbeg + i + 1, tbase)[0];
	}
}

static inline void banded_striped_epi8_seqalign_piecex_backtrace(BaseBank *seqs, u8i qoff, u4i qlen, u8i toff, u4i tlen, b1v *btds, u4v *begs, u4i bandwidth, int piecewise, seqalign_result_t *rs, String **alnstr){
	u4i i, bt, btx;
	u1i qbase, tbase;
	rs->qb = rs->qe; rs->qe ++;
	rs->tb = rs->te; rs->te ++;
	rs->mat = rs->mis = rs->ins = rs->del = 0;
	if(alnstr){
		clear_string(alnstr[0]);
		clear_string(alnstr[1]);
		clear_string(alnstr[2]);
		if(0){
			for(i=qlen;(int)i>rs->qe;i--){
				push_string(alnstr[0], bit_base_table[get_basebank(seqs, qoff + i - 1)]);
				push_string(alnstr[1], '-');
				push_string(alnstr[2], '-');
			}
			for(i=tlen;(int)i>rs->te;i--){
				push_string(alnstr[0], '-');
				push_string(alnstr[1], bit_base_table[get_basebank(seqs, toff + i - 1)]);
				push_string(alnstr[2], '-');
			}
		}
	}
	bt = 0;
	while(1){
		btx = get_b1v(btds, rs->tb * bandwidth + banded_striped_epi8_pos2idx(bandwidth, rs->qb - get_u4v(begs, rs->tb)));
		if(piecewise == 0 || bt == SEQALIGN_BT_M){
			bt = btx & 0x7;
		} else {
			bt = ((btx >> (bt + 2)) & 0x1)? bt : SEQALIGN_BT_M;
		}
		if(bt == SEQALIGN_BT_M){
			qbase = get_basebank(seqs, qoff + rs->qb);
			tbase = get_basebank(seqs, toff + rs->tb);
			if(qbase == tbase){
				rs->mat ++;
				if(alnstr){
					push_string(alnstr[0], bit_base_table[qbase]);
					push_string(alnstr[1], bit_base_table[tbase]);
					push_string(alnstr[2], '|');
				}
			} else {
				rs->mis ++;
				if(alnstr){
					push_string(alnstr[0], bit_base_table[qbase]);
					push_string(alnstr[1], bit_base_table[tbase]);
					push_string(alnstr[2], '*');
				}
			}
			rs->qb --;
			rs->tb --;
			if(rs->qb < 0 || rs->tb < 0) break;
		} else if(bt == SEQALIGN_BT_D1 || bt == SEQALIGN_BT_D2){
			tbase = get_basebank(seqs, toff + rs->tb);
			rs->del ++;
			if(alnstr){
				push_string(alnstr[0], '-');
				push_string(alnstr[1], bit_base_table[tbase]);
				push_string(alnstr[2], '-');
			}
			rs->tb --;
			if(rs->tb < 0) break;
		} else {
			qbase = get_basebank(seqs, qoff + rs->qb);
			rs->ins ++;
			if(alnstr){
				push_string(alnstr[0], bit_base_table[qbase]);
				push_string(alnstr[1], '-');
				push_string(alnstr[2], '-');
			}
			rs->qb --;
			if(rs->qb < 0) break;
		}
	}
	rs->qb ++;
	rs->tb ++;
	if(alnstr){
		if(0){
			for(i=rs->qb;i>0;i--){
				push_string(alnstr[0], bit_base_table[get_basebank(seqs, qoff + i - 1)]);
				push_string(alnstr[1], '-');
				push_string(alnstr[2], '-');
			}
			for(i=rs->tb;i>0;i--){
				push_string(alnstr[0], '-');
				push_string(alnstr[1], bit_base_table[get_basebank(seqs, toff + i - 1)]);
				push_string(alnstr[2], '-');
			}
		}
		reverse_string(alnstr[0]);
		reverse_string(alnstr[1]);
		reverse_string(alnstr[2]);
	}
}

static inline seqalign_result_t banded_striped_epi8_seqalign_pairwise_overlap(BaseBank *seqs, u8i qoff, u4i qlen, u8i toff, u4i tlen, b1v *qprof, b1v *rows, b1v *btds, u4v *begs, u4i bandwidth, b1i *matrix, b1i gapo1, b1i gape1, b1i gapo2, b1i gape2, String **alnstr, int verbose){
	banded_striped_epi8_seqalign_set_query_prof_func    set_qprof;
	banded_striped_epi8_seqalign_piecex_row_mov_func    row_mov;
	banded_striped_epi8_seqalign_piecex_row_cal_func    row_cal;
	banded_striped_epi8_seqalign_piecex_row_verify_func row_verify;
	banded_striped_epi8_seqalign_piecex_row_max_func    row_max;
	banded_striped_epi8_seqalign_piecex_backtrace_func  backtrace;
	seqalign_result_t rs;
	b1i **us[2], **es[2], **qs[2], *bs;
	__m128i ZERO, MIN;
	u4i i, k, W, rbeg, rmax, mov;
	int piecewise, smax, lst_score, score, ph, rh, hh, wl, wr, wt, *wscores[2];
	u1i tbase;
	set_qprof  = banded_striped_epi8_seqalign_set_query_prof;
	row_verify = banded_striped_epi8_seqalign_piecex_row_verify;
	row_max    = banded_striped_epi8_seqalign_row_max;
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
	bandwidth = roundup_times(bandwidth, 16);
	W = bandwidth / 16;
	if(verbose){
		fprintf(stdout, "[%d,%d][%d,%d] PIECEWISE=%d\tW=%d\n", gapo1, gape1, gapo2, gape2, piecewise, W);
	}
	smax = 0;
	for(i=0;i<16;i++) smax = num_max(smax, num_abs(matrix[i]));
	set_qprof(seqs, qoff, qlen, qprof, bandwidth, matrix);
	clear_and_encap_b1v(rows, bandwidth * (piecewise + 1) * 2);
	clear_and_encap_b1v(btds, bandwidth * (tlen + 1));
	clear_and_encap_u4v(begs, tlen);
	us[0] = malloc(W * sizeof(b1i*));
	us[1] = malloc(W * sizeof(b1i*));
	if(piecewise){
		es[0] = malloc(W * sizeof(b1i*));
		es[1] = malloc(W * sizeof(b1i*));
	} else {
		es[0] = es[1] = NULL;
	}
	if(piecewise == 2){
		qs[0] = malloc(W * sizeof(b1i*));
		qs[1] = malloc(W * sizeof(b1i*));
	} else {
		qs[0] = qs[1] = NULL;
	}
	i = 0;
	for(k=0;k<W;k++) us[1][k] = rows->buffer + 16 * k + i * bandwidth;
	i ++;
	for(k=0;k<W;k++) us[0][k] = rows->buffer + 16 * k + i * bandwidth;
	i ++;
	if(piecewise){
		for(k=0;k<W;k++) es[1][k] = rows->buffer + 16 * k + i * bandwidth;
		i ++;
		for(k=0;k<W;k++) es[0][k] = rows->buffer + 16 * k + i * bandwidth;
		i ++;
	}
	if(piecewise == 2){
		for(k=0;k<W;k++) qs[1][k] = rows->buffer + 16 * k + i * bandwidth;
		i ++;
		for(k=0;k<W;k++) qs[0][k] = rows->buffer + 16 * k + i * bandwidth;
		i ++;
	}
	bs = btds->buffer;
	wscores[0] = calloc(roundup_times(W, 4), sizeof(int));
	wscores[1] = calloc(roundup_times(W, 4), sizeof(int));
	// overlap alignment
	ZERO = _mm_set1_epi8(0);
	MIN = _mm_set1_epi8(SEQALIGN_SCORE_EPI8_MIN);
	for(k=0;k<W;k++) _mm_store_si128((__m128i*)us[1][k], ZERO);
	if(piecewise){
		for(k=0;k<W;k++) _mm_store_si128((__m128i*)es[1][k], MIN);
	}
	if(piecewise == 2){
		for(k=0;k<W;k++) _mm_store_si128((__m128i*)qs[1][k], MIN);
	}
	for(k=0;k<W;k++) _mm_store_si128(((__m128i*)bs) + k, ZERO);
	memset(&rs, 0, sizeof(seqalign_result_t));
	rs.score = SEQALIGN_SCORE_MIN;
	rbeg = rmax = 0;
	mov = 0;
	hh = 0;
	lst_score = 0;
	for(i=0;i<tlen;i++){
		ph = hh;
		tbase = get_basebank(seqs, toff + i);
		if(mov && rbeg + bandwidth < qlen){
			mov = num_min(mov, num_max(0, Int(qlen) - Int(rbeg + bandwidth)));
			rbeg += mov;
			ph = wscores[0][mov - 1];
			rh = ph + matrix[tbase * 4 + get_basebank(seqs, qoff + rbeg)]; // just guess
		} else {
			mov  = 0;
			rh = rbeg? SEQALIGN_SCORE_MIN : 0; // overlap alignment mode
		}
		row_mov(us, es, qs, W, mov); // mov and swap
		hh = row_cal(rbeg, tbase, us, es, qs, bs, qprof, gapo1, gape1, gapo2, gape2, W, ph, rh);
		if(verbose > 2){
			row_verify(rbeg, W, ph, rh, hh, tbase, us, es, qs, qprof, gapo1, gape1, gapo2, gape2);
		}
		rmax = row_max(us[1], W, hh, wscores, &lst_score);
		// adaptive banded
		wl = array_epi32_sum(wscores[0], W);
		wr = array_epi32_sum(wscores[1], W);
		wt = W * smax;
		if(verbose){
			banded_striped_epi8_seqalign_row_print(stdout, seqs, 1, i, tbase, bandwidth, mov, qoff, rbeg, rmax, lst_score, hh, us[1], bs, wl, wr, verbose > 1);
		}
		if(i < 4 * W){
			mov = 0;
		} else if(wl + wt < wr){
			mov = 2;
		} else if(wl > wr + wt){
			mov = 0;
		} else {
			mov = 1;
		}
		bs += bandwidth;
		push_u4v(begs, rbeg);
		if(rbeg + bandwidth >= qlen){
			score = wscores[1][qlen - 1 - rbeg - 15 * W];
			if(score > rs.score){
				rs.score = score;
				rs.qe = qlen - 1;
				rs.te = i;
			}
		}
	}
	if(lst_score > rs.score){
		rs.score = lst_score;
		rs.qe = rbeg + rmax;
		rs.te = tlen - 1;
	}
	// backtrace
	backtrace(seqs, qoff, qlen, toff, tlen, btds, begs, bandwidth, piecewise, &rs, alnstr);
	free(us[0]);
	free(us[1]);
	if(es[0]) free(es[0]);
	if(es[1]) free(es[1]);
	if(qs[0]) free(qs[0]);
	if(qs[1]) free(qs[1]);
	free(wscores[0]);
	free(wscores[1]);
	return rs;
}

#endif
