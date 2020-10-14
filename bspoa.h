#ifndef BANDED_STRIPED_SIMD_RECURRENT_ALIGNMENT_GRAPH_MSA_CNS_RJ_H
#define BANDED_STRIPED_SIMD_RECURRENT_ALIGNMENT_GRAPH_MSA_CNS_RJ_H

#include <math.h>
#include "mem_share.h"
#include "sort.h"
#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "chararray.h"
#include "bsalign.h"

#if __BYTE_ORDER == 1234
//#pragma message(" ** " __FILE__ " has been tested in LITTLE_ENDIAN **\n")
#else
#pragma message(" ** " __FILE__ " hasn't been tested in BIG_ENDIAN **\n")
#endif

//#define BSPOACNS_NO_RECUR	1

#define SEQALIGN_MODE_POA	16

#define BSPOA_RDLEN_MAX	0x0FFFFFFFU
#define BSPOA_RDCNT_MAX	0x3FFF
#define BSPOA_HEAD_NODE	0
#define BSPOA_TAIL_NODE	1

#define BSPOA_VST_MAX	MAX_U2

typedef struct {
	u2i rid;
	u1i base:7, ref:1, state, aux;
	u2i vst;
	u2i nin, nou;
	u2i nct, cov;
	u4i pos, mpos, cpos, rpos;
	u4i edge, erev;
	u4i aligned, header;
	u4i mmidx;
} bspoanode_t;
define_list(bspoanodev, bspoanode_t);

// put pair of edges together, idx % 2 = 0, 1
typedef struct {
	u4i node;
	u4i cov:29, is_aux:1, mask:1, vst:1;
	u4i next;
} bspoaedge_t;
define_list(bspoaedgev, bspoaedge_t);

typedef struct {
	int refmode; // 0: no
	int alnmode; // SEQALIGN_MODE_OVERLAP
	int remsa; // 1
	int nrec; // 20, new read will align against previous nrec reads on POA
	int ksz; // 15
	int bwtrigger; // 1, min aligned reads to tigger bandwidth for next read
	int bandwidth; // 128, when tiggering bandwidth, first generate cns, then align read against cns to get the offset of each node' band
	int bwmax; // 1024 * 64
	int M, X, O, E, Q, P; // 2, -6, -3, -2, -8, -1
	int refbonus; // 1
	int rma_win; // 5, min length of flinking high quality cns bases
	int qltlo, qlthi; // 30, 35. trigger local remsa when cns quality <= qltlo
	float psub, pins, pdel, piex, pdex, hins, hdel; // 0.05, 0.05, 0.10, 0.25, 0.30, 0.10, 0.20
} BSPOAPar;

static const BSPOAPar DEFAULT_BSPOA_PAR = (BSPOAPar){0, SEQALIGN_MODE_OVERLAP, 1, 20, 15, 1, 128, 64 * 1024, 2, -6, -3, -2, -8, -1, 1, 5, 30, 35, 0.05, 0.05, 0.10, 0.25, 0.30, 0.10, 0.20};

typedef struct {
	u4i coff:29, bt:3;
	float max;
} bspoacns_t;

typedef struct {
	u4i cpos, mpos;
	u2i refn, altn;
	u1i refb, altb;
} bspoavar_t;
define_list(bspoavarv, bspoavar_t);

typedef struct {
	String   *mtag;
	SeqBank  *seqs;
	u4v      *ndoffs;
	u4v      *cigars, *cgbs, *cges; // refmode, concatenate all cigars together
	bspoanodev *nodes;
	bspoaedgev *edges;
	u4v      *ecycs;
	BSPOAPar *par;
	int piecewise, bwtrigger;
	u4i  bandwidth, qlen, slen, qb, qe, nrds; // real bandwidth used in alignment
	u1v *qseq;
	b1i matrix[2][16];
	b1v *qprof[2];
	b1v *memp; // memory pool
	u8i mmblk, mmcnt;
	int maxscr, maxidx, maxoff;
	u8v *heap, *todels;
	u4v *stack;
	u4i backbone;
	u1v *msacols;
	u4v *msaidxs, *msacycs;
	float dpvals[8], dporis[8];
	u1i *dptable;
	u1v *cns, *qlt, *alt;
	bspoavarv *var;
	String *strs;
	u8i ncall, obj_uid;
} BSPOA;

static inline void gen_cns_aln_event_table_bspoa(BSPOAPar *par, float ps[8], float os[8], u1i *table){
	u4i i, a, b, c, d;
	os[0] = 1 - par->psub; // M
	os[1] = par->psub; // X
	os[2] = par->pins; // I
	os[3] = par->pdel; // D
	os[4] = par->piex; // IE
	os[5] = par->pdex; // DE
	os[6] = par->hins; // HI
	os[7] = par->hdel; // HD
	for(i=0;i<8;i++) ps[i] = log(os[i]);
	for(i=0;i<5*5*5*5;i++){
		a = (i % 5); // cur cns base
		b = (i % (5 * 5)) / 5; // cur read base
		c = (i % (5 * 5 * 5)) / (5 * 5); // last cns non-N base; if b < 4, c = b; else = last last ...
		d = (i % (5 * 5 * 5 * 5)) / (5 * 5 * 5); // last state M/I/D
		if(a < 4){
			if(b < 4){
				if(a == b){
					table[i] = (0 << 3) | 0; // ps[M] << 3 | M
				} else {
					if(b == c && ps[6] > ps[1]){
						table[i] = (6 << 3) | 0; // ps[HI] << 3 | M
					} else {
						table[i] = (1 << 3) | 0; // ps[X] << 3 | M
					}
				}
			} else {
				if(d == 2){
					if(a == c && ps[7] > ps[5]){
						table[i] = (7 << 3) | 2; // ps[HD] << 3 | D
					} else {
						table[i] = (5 << 3) | 2; // ps[DE] << 3 | D
					}
				} else {
					if(a == c && ps[7] > ps[3]){
						table[i] = (7 << 3) | 2; // ps[HD] << 3 | D
					} else {
						table[i] = (3 << 3) | 2; // ps[D] << 3 | D
					}
				}
			}
		} else {
			if(b < 4){
				if(d == 1){
					if(b == c && ps[6] > ps[4]){
						table[i] = (6 << 3) | 1; // ps[HI] << 3 | I
					} else {
						table[i] = (4 << 3) | 1; // ps[IE] << 3 | I
					}
				} else {
					if(b == c && ps[6] > ps[2]){
						table[i] = (6 << 3) | 1; // ps[HI] << 3 | I
					} else {
						table[i] = (2 << 3) | 1; // ps[I] << 3 | I
					}
				}
			} else {
				table[i] = (0 << 3) | d; // ps[M] | last_state
			}
		}
	}
}

static inline BSPOA* init_bspoa(BSPOAPar par){
	BSPOA *g;
	g = malloc(sizeof(BSPOA));
	g->mtag = init_string(8);
	g->seqs = init_seqbank();
	g->ndoffs = init_u4v(1024);
	g->cigars = init_u4v(1024);
	g->cgbs   = init_u4v(32);
	g->cges   = init_u4v(32);
	g->nodes = init_bspoanodev(16 * 1024);
	g->edges = init_bspoaedgev(16 * 1024);
	g->ecycs = init_u4v(32);
	g->par   = malloc(sizeof(BSPOAPar));
	memcpy(g->par, &par, sizeof(BSPOAPar));
	g->par->bandwidth = roundup_times(g->par->bandwidth, WORDSIZE);
	g->piecewise = 1;
	g->nrds      = 0;
	g->bwtrigger  = 0;
	g->bandwidth = 0;
	g->qseq  = init_u1v(1024);
	g->qlen  = 0;
	g->qprof[0] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->qprof[1] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->memp  = adv_init_b1v(1024, 0, WORDSIZE, 0);
	g->mmblk = 0;
	g->mmcnt = 0;
	g->heap = init_u8v(32);
	g->todels = init_u8v(4);
	g->stack = init_u4v(32);
	g->backbone = 0;
	g->msacols = init_u1v(16 * 1024);
	g->msaidxs = init_u4v(1024);
	g->msacycs = init_u4v(8);
	g->dptable = malloc(5 * 5 * 5 * 5);
	gen_cns_aln_event_table_bspoa(g->par, g->dpvals, g->dporis, g->dptable);
	g->strs = init_string(1024);
	g->cns = init_u1v(1024);
	g->qlt = init_u1v(1024);
	g->alt = init_u1v(1024);
	g->var = init_bspoavarv(8);
	g->ncall = 0;
	g->obj_uid = 0;
	return g;
}

static inline void renew_bspoa(BSPOA *g){
	free_seqbank(g->seqs); g->seqs = init_seqbank();
	renew_u4v(g->ndoffs, 1024);
	renew_u4v(g->cigars, 1024);
	renew_u4v(g->cgbs, 32);
	renew_u4v(g->cges, 32);
	renew_bspoanodev(g->nodes, 16 * 1024);
	renew_bspoaedgev(g->edges, 16 * 1024);
	renew_u4v(g->ecycs, 32);
	renew_u1v(g->qseq, 1024);
	renew_b1v(g->qprof[0], 4 * 1024);
	renew_b1v(g->qprof[1], 4 * 1024);
	renew_b1v(g->memp, 4 * 1024);
	g->nrds  = 0;
	g->mmblk = 0;
	g->mmcnt = 0;
	renew_u8v(g->heap, 32);
	renew_u4v(g->stack, 32);
	renew_u1v(g->msacols, 16 * 1024);
	renew_u4v(g->msaidxs, 1024);
	renew_u4v(g->msacycs, 1024);
	renew_u1v(g->cns, 1024);
	renew_u1v(g->qlt, 1024);
	renew_u1v(g->alt, 1024);
	recap_string(g->strs, 1024);
}

static inline void free_bspoa(BSPOA *g){
	free_string(g->mtag);
	free_seqbank(g->seqs);
	free_u4v(g->ndoffs);
	free_u4v(g->cigars);
	free_u4v(g->cgbs);
	free_u4v(g->cges);
	free_bspoanodev(g->nodes);
	free_bspoaedgev(g->edges);
	free_u4v(g->ecycs);
	free_u1v(g->qseq);
	free_b1v(g->qprof[0]);
	free_b1v(g->qprof[1]);
	free_b1v(g->memp);
	free_u8v(g->heap);
	free_u8v(g->todels);
	free_u4v(g->stack);
	free_u1v(g->msacols);
	free_u4v(g->msaidxs);
	free_u4v(g->msacycs);
	free(g->dptable);
	free_u1v(g->cns);
	free_u1v(g->qlt);
	free_u1v(g->alt);
	free_bspoavarv(g->var);
	free_string(g->strs);
	free(g->par);
	free(g);
}

static inline void push_bspoacore(BSPOA *g, char *seq, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX && len){
		len = num_min(len, BSPOA_RDLEN_MAX);
		push_seqbank(g->seqs, NULL, 0, seq, len);
		push_u4v(g->cgbs, g->cigars->size);
		append_array_u4v(g->cigars, cgs, ncg);
		push_u4v(g->cges, g->cigars->size);
	}
}

static inline void push_bspoa(BSPOA *g, char *seq, u4i len){ push_bspoacore(g, seq, len, NULL, 0); }

static inline void fwdbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		fwdbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u4v(g->cgbs, g->cigars->size);
		append_array_u4v(g->cigars, cgs, ncg);
		push_u4v(g->cges, g->cigars->size);
	}
}

static inline void fwdbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ fwdbitpush_bspoacore(g, bits, off, len, NULL, 0); }

static inline void revbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		revbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u4v(g->cgbs, g->cigars->size);
		append_array_u4v(g->cigars, cgs, ncg);
		push_u4v(g->cges, g->cigars->size);
	}
}

static inline void revbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ revbitpush_bspoacore(g, bits, off, len, NULL, 0); }

static inline void print_dot_bspoa(BSPOA *g, u4i posbeg, u4i posend, u4i mincnt, FILE *out){
	bspoanode_t *n;
	bspoaedge_t *e;
	u4i nidx, eidx;
	fprintf(out, "digraph {\n");
	fprintf(out, "rankdir=LR\n");
	fprintf(out, "N0 [label=\"BEG\"]\n");
	fprintf(out, "N1 [label=\"END\"]\n");
	for(nidx=BSPOA_TAIL_NODE+1;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->mpos < posbeg || n->mpos >= posend) continue;
#if 1
		if(n->nin == 0 && n->nou == 0) continue;
		if(n->cov >= mincnt){
			fprintf(out, "N%u [label=%c%u_%d_%d_N%u color=blue]\n", nidx, bit_base_table[(n->base) & 0x03], n->mpos, n->cpos, n->cov, nidx);
		} else {
			fprintf(out, "N%u [label=%c%u_%d_%d_N%u]\n", nidx, bit_base_table[(n->base) & 0x03], n->mpos, n->cpos, n->cov, nidx);
		}
#else
		fprintf(out, "N%u [label=R%u_%u_%c]\n", nidx, n->rid, n->pos, bit_base_table[(n->base) & 0x03]);
#endif
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->mpos < posbeg || n->mpos >= posend) continue;
#if 1
		if(n->nin == 0 && n->nou == 0) continue;
#else
		if(n->aligned != nidx){
			fprintf(out, "N%u -> N%u [color=magenta style=dashed]\n", nidx, n->aligned);
		}
#endif
		eidx = n->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
#if DEBUG
			if(e->next == eidx){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			eidx = e->next;
			if(e->is_aux) continue;
			fprintf(out, "N%u -> N%u [label=%u%s]\n", nidx, e->node, e->cov, e->cov >= mincnt? " color=blue" : "");
		}
	}
	fprintf(out, "}\n");
}

static inline void fprint_dot_bspoa(BSPOA *g, u4i msabeg, u4i msaend, u4i msacnt, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_dot_bspoa(g, msabeg, msaend, msacnt, out);
	fclose(out);
}

static inline void print_vstdot_bspoa(BSPOA *g, char *fname){
	FILE *out;
	bspoanode_t *n;
	bspoaedge_t *e;
	u4i nidx, eidx;
	out = open_file_for_write(fname, NULL, 1);
	fprintf(out, "digraph {\n");
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->state == 0) continue;
		//if(nidx && n->vst == 0) continue;
		fprintf(out, "N%u [label=\"N%u:%u:%u:%d\"]\n", nidx, nidx, n->nin, n->nct, n->vst);
		eidx = n->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(e->is_aux) continue;
			if(ref_bspoanodev(g->nodes, e->node)->state == 0) continue;
			fprintf(out, "N%u -> N%u [label=%u]\n", nidx, e->node, e->cov);
		}
	}
	fprintf(out, "}\n");
	fclose(out);
}

static inline void print_local_dot_bspoa(BSPOA *g, u4i ori_nidx, u4i step, FILE *out){
	bspoanode_t *n;
	bspoaedge_t *e;
	u4i nidx, eidx, k, *auxs;
	fprintf(out, "digraph {\n");
	fprintf(out, "rankdir=LR\n");
	auxs = (u4i*)calloc(g->nodes->size, sizeof(u4i));
	for(k=0;k<2;k++){
		clear_u4v(g->stack);
		push_u4v(g->stack, ori_nidx);
		auxs[ori_nidx] = 1;
		while(pop_u4v(g->stack, &nidx)){
			n = ref_bspoanodev(g->nodes, nidx);
			eidx = k? n->erev : n->edge;
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				if(auxs[e->node]) continue;
				if(auxs[nidx] > step) continue;
				push_u4v(g->stack, e->node);
				auxs[e->node] = auxs[nidx] + 1;
			}
		}
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		if(auxs[nidx] == 0) continue;
		n = ref_bspoanodev(g->nodes, nidx);
		fprintf(out, "N%u [label=%c%u_%d_%d_N%u%s]\n", nidx, bit_base_table[(n->base) & 0x03], n->rpos, n->cpos, n->cov, nidx, (nidx == ori_nidx)? " style=filled fillcolor=yellow" : "");
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		if(auxs[nidx] == 0) continue;
		n = ref_bspoanodev(g->nodes, nidx);
		eidx = n->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(e->is_aux) continue;
			fprintf(out, "N%u -> N%u [label=%u]\n", nidx, e->node, e->cov);
		}
	}
	fprintf(out, "}\n");
	free(auxs);
}

static inline void fprint_local_dot_bspoa(BSPOA *g, u4i ori_nidx, u4i step, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_local_dot_bspoa(g, ori_nidx, step, out);
	fclose(out);
}

static inline int print_node_edges_bspoa(BSPOA *g, u4i nidx, int rev, FILE *out){
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	u4i eidx, ret;
	v = ref_bspoanodev(g->nodes, nidx);
	eidx = rev? v->erev : v->edge;
	ret = 0;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		w = ref_bspoanodev(g->nodes, e->node);
#if DEBUG
		fprintf(out, "E%u\t%c\t%d\tN%u -> N%u\t%d\tR%u_%u_%u\n", eidx, "+-"[rev], e->vst, nidx, e->node, e->cov, w->rid, w->pos, w->base);
#else
		fprintf(out, "E%u\t%c\tN%u -> N%u\t%d\tR%u_%u_%u\n", eidx, "+-"[rev], nidx, e->node, e->cov, w->rid, w->pos, w->base);
#endif
		eidx = e->next;
		ret ++;
	}
	return ret;
}

static inline void print_seqs_bspoa(BSPOA *g, char *prefix, char *suffix){
	FILE *out;
	u4i i;
	out = open_file_for_write(prefix, suffix, 1);
	for(i=0;i<g->seqs->nseq;i++){
		fprintf(out, ">S%u len=%u\n", i, g->seqs->rdlens->buffer[i]);
		println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[i], g->seqs->rdlens->buffer[i], out);
	}
	fclose(out);
}

static inline char* str_msa_bspoa(BSPOA *g, String *strp){
	u1i *col;
	u4i i, j, b, nseq, mrow, mlen;
	char *str, *cp;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	if(strp){
		clear_string(strp);
		encap_string(strp, (mrow + 2) * (mlen + 1));
		str = strp->string;
	} else {
		str = malloc((mrow + 2) * (mlen + 1));
	}
	for(i=0;i<mlen;i++){
		col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
		for(j=0;j<=nseq;j++){
			//if(col[nseq] < 4 && col[j] < 4 && col[j] != col[nseq]){
			if(col[j] < 4 && col[j] != col[nseq]){
				str[j * (mlen + 1) + i] = "acgt-.*"[col[j]];
			} else {
				str[j * (mlen + 1) + i] = "ACGT-.*"[col[j]];
			}
		}
		for(;j<mrow;j++){
			str[j * (mlen + 1) + i] = '!' + col[j]; // phred qualty score
		}
	}
	for(j=0;j<mrow;j++) str[j * (mlen + 1) + mlen] = '\0';
	cp = str + mrow * (mlen + 1);
	for(i=j=0;i<mlen;i++){
		if((i % 10) == 0 && j + 6 <= mlen){
			sprintf(cp, "|%05u", i + 1);
			j  += 6;
			cp += 6;
		} else if(i >= j){
			cp[0] = ' ';
			cp ++;
			j ++;
		}
	}
	cp[0] = '\0';
	cp ++;
	for(i=j=b=0;i<mlen;i++){
		if(g->msacols->buffer[i * mrow + nseq] < 4){
			j ++;
			if((j % 10) == 1){
				while(b < i){
					cp[0] = ' ';
					cp ++;
					b  ++;
				}
				if(b + 6 < mlen){
					sprintf(cp, "|%05u", j);
					cp += 6;
					b  += 6;
				}
			}
		}
	}
	while(b < mlen){ cp[0] = ' '; cp ++; b ++; }
	cp[0] = '\0';
	return str;
}

static inline void print_msa_sline_bspoa(BSPOA *g, FILE *out){
	char *str;
	u4i i, nseq, mrow, mlen;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	str = str_msa_bspoa(g, g->strs);
	fprintf(out, "MSA [POS] %s\n", str + mrow * (mlen + 1));
	for(i=0;i<mrow;i++){
		if(i == nseq){
			fprintf(out, "MSA [CNS] ");
		} else if(i == nseq + 1){
			fprintf(out, "MSA [QLT] ");
		} else {
			fprintf(out, "MSA [%03u] ", i);
		}
		fprintf(out, "%s\n", str + i * (mlen + 1));
	}
	fprintf(out, "MSA [POS] %s\n", str + (mrow + 1) * (mlen + 1));
	for(i=0;i<g->cns->size;i++){
		str[i] = bit_base_table[g->cns->buffer[i]];
	}
	str[i] = 0;
	fprintf(out, "CNS\t%d\t%s\n", Int(g->cns->size), str);
	str += g->cns->size + 1;
	for(i=0;i<g->qlt->size;i++){
		str[i] = '!' +  g->qlt->buffer[i]; // Phred + 33
	}
	str[i] = 0;
	fprintf(out, "QLT\t%d\t%s\n", Int(g->qlt->size), str);
	str += g->cns->size + 1;
	for(i=0;i<g->alt->size;i++){
		str[i] = '!' +  g->alt->buffer[i]; // Phred + 33
	}
	str[i] = 0;
	fprintf(out, "ALT\t%d\t%s\n", Int(g->alt->size), str);
}

static inline void print_msa_mline_bspoa(BSPOA *g, FILE *out){
	char *str, *cp, c;
	u4i i, j, b, nseq, mrow, mlen;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	str = str_msa_bspoa(g, g->strs);
	for(i=0;i<mlen;i+=100){
		b = num_min(i + 100, mlen);
		cp = str + mrow * (mlen + 1);
		c = cp[b]; cp[b] = '\0';
		fprintf(out, "MSA [POS] %s\n", cp + i);
		cp[b] = c;
		for(j=0;j<mrow;j++){
			cp = str + j * (mlen + 1);
			c = cp[b]; cp[b] = '\0';
			if(j == nseq){
				fprintf(out, "MSA [CNS] ");
			} else if(j == nseq + 1){
				fprintf(out, "MSA [QLT] ");
			} else if(j == nseq + 2){
				fprintf(out, "MSA [ALT] ");
			} else {
				fprintf(out, "MSA [%03u] ", j);
			}
			fprintf(out, "%s\n", cp + i);
			cp[b] = c;
		}
		cp = str + (mrow + 1) * (mlen + 1);
		c = cp[b]; cp[b] = '\0';
		fprintf(out, "MSA [POS] %s\n", cp + i);
		cp[b] = c;
		fprintf(out, "MSA\n");
	}
	for(i=0;i<g->cns->size;i++){
		str[i] = bit_base_table[g->cns->buffer[i]];
	}
	str[i] = 0;
	fprintf(out, "CNS\t%d\t%s\n", Int(g->cns->size), str);
	str += g->cns->size + 1;
	for(i=0;i<g->qlt->size;i++){
		str[i] = '!' +  g->qlt->buffer[i]; // Phred + 33
	}
	str[i] = 0;
	fprintf(out, "QLT\t%d\t%s\n", Int(g->qlt->size), str);
	str += g->cns->size + 1;
	for(i=0;i<g->alt->size;i++){
		str[i] = '!' +  g->alt->buffer[i]; // Phred + 33
	}
	str[i] = 0;
	fprintf(out, "ALT\t%d\t%s\n", Int(g->alt->size), str);
}

static inline void print_snp_bspoa(BSPOA *g, FILE *out){
	bspoavar_t *var;
	u4i i, j, fsz, fct;
	char flanks[2][6];
	fsz = 5;
	for(i=0;i<g->var->size;i++){
		var = ref_bspoavarv(g->var, i);
		fct = num_min(var->cpos, fsz);
		memcpy(flanks[0], g->cns->buffer + var->cpos - fct, fct);
		for(j=0;j<fct;j++) flanks[0][j] = bit_base_table[(int)flanks[0][j]];
		flanks[0][fct] = 0;
		fct = num_min(g->cns->size - var->cpos - 1, fsz);
		memcpy(flanks[1], g->cns->buffer + var->cpos + 1, fct);
		for(j=0;j<fct;j++) flanks[1][j] = bit_base_table[(int)flanks[1][j]];
		flanks[1][fct] = 0;
		if(g->mtag->size){
			fprintf(out, "SNP %s\t", g->mtag->string);
		} else {
			fprintf(out, "SNP %llu\t", g->ncall);
		}
		fprintf(out, "%d\t%d\t%s\t%c\t%d\t%c\t%d\t%s\t%d\n", var->mpos, var->cpos, flanks[0], bit_base_table[var->refb], var->refn, bit_base_table[var->altb], var->altn, flanks[1], g->alt->buffer[var->cpos]);
	}
}

static inline void check_node_edges_bspoa(BSPOA *g, u4i nidx, int rev){
	bspoanode_t *v, *w;
	bspoaedge_t *e, *r;
	u4i eidx, ridx;
	v = ref_bspoanodev(g->nodes, nidx);
	eidx = rev? v->erev : v->edge;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		r = rev? e - 1 : e + 1;
		if(r->node != nidx){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		w = ref_bspoanodev(g->nodes, e->node);
		ridx = rev? w->edge : w->erev;
		while(ridx){
			r = ref_bspoaedgev(g->edges, ridx);
			if(r->node == nidx){
				break;
			}
			ridx = r->next;
		}
		if(ridx == 0){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

static inline void check_all_node_edges_bspoa(BSPOA *g){
	u4i nidx;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		check_node_edges_bspoa(g, nidx, 0);
		check_node_edges_bspoa(g, nidx, 1);
	}
}

static inline bspoanode_t* add_node_bspoa(BSPOA *g, u2i rid, u4i pos, u1i base, int node_exists){
	bspoanode_t *u;
	if(node_exists){
		u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + pos);
	} else {
		u = next_ref_bspoanodev(g->nodes);
	}
	ZEROS(u);
	u->rid  = rid;
	u->pos  = pos;
	u->base = base;
	u->aux  = 0;
	u->aligned = offset_bspoanodev(g->nodes, u);
	u->header  = offset_bspoanodev(g->nodes, u);
	//check_all_node_edges_bspoa(g);
	return u;
}

static inline bspoaedge_t* new_edge_bspoa(BSPOA *g, bspoanode_t *u, bspoanode_t *v, int cov, int is_aux){
	bspoaedge_t *e, *r;
	if(g->ecycs->size){
		e = ref_bspoaedgev(g->edges, g->ecycs->buffer[--g->ecycs->size]);
		r = e + 1;
	} else {
		encap_bspoaedgev(g->edges, 2);
		e = next_ref_bspoaedgev(g->edges);
		r = next_ref_bspoaedgev(g->edges);
	}
	e->node = offset_bspoanodev(g->nodes, v);
	r->node = offset_bspoanodev(g->nodes, u);
	e->cov = cov;
	r->cov = cov;
	e->is_aux = is_aux;
	r->is_aux = is_aux;
	e->next = 0;
	r->next = 0;
	return e;
}

static inline void add_edge_bspoacore(BSPOA *g, bspoanode_t *v, u4i eidx){
	bspoaedge_t *e, *f, *p;
	u4i *eptr;
	if(eidx&1){
		v->nin ++;
		eptr = &v->erev;
	} else {
		v->nou ++;
		eptr = &v->edge;
	}
	if((*eptr) == 0){
		*eptr = eidx;
		return;
	}
	e = ref_bspoaedgev(g->edges, eidx);
	p = ref_bspoaedgev(g->edges, *eptr);
	if(e->cov > p->cov){
		e->next = *eptr;
		*eptr = eidx;
		return;
	}
	while(p->next){
		f = ref_bspoaedgev(g->edges, p->next);
		if(e->cov > f->cov){
			break;
		}
		p = f;
	}
	e->next = p->next;
	p->next = eidx;
}

static inline bspoaedge_t* add_edge_bspoa(BSPOA *g, bspoanode_t *u, bspoanode_t *v, int cov, int is_aux, int *exists){
	bspoaedge_t *e, *f, *r;
	u4i *eptr;
	int exs;
#if DEBUG
	if(u == v){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	//if(offset_bspoanodev(g->nodes, u) == 0 && offset_bspoanodev(g->nodes, v) == 154901){
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
	encap_bspoaedgev(g->edges, 2);
	e = f = NULL;
	eptr = &u->edge;
	exs = 0;
	while(*eptr){
		f = ref_bspoaedgev(g->edges, *eptr);
		if(f->node == offset_bspoanodev(g->nodes, v)){
			e = f;
			break;
		}
		eptr = &f->next;
	}
	if(e == NULL){
		exs = 0;
		e = new_edge_bspoa(g, u, v, cov, is_aux);
		r = e + 1;
	} else {
		exs = 1;
		*eptr = e->next;
		e->next = 0;
		v->nin --;
		u->nou --;
		eptr = &v->erev;
		r = f = NULL;
		while(*eptr){
			f = ref_bspoaedgev(g->edges, *eptr);
			if(f->node == offset_bspoanodev(g->nodes, u)){
				r = f;
				break;
			}
			eptr = &f->next;
		}
		if(r == NULL || r != e + 1){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		*eptr = r->next;
		r->next = 0;
		e->is_aux |= is_aux;
		r->is_aux |= is_aux;
		e->cov += cov;
		r->cov += cov;
	}
	add_edge_bspoacore(g, u, offset_bspoaedgev(g->edges, e));
	add_edge_bspoacore(g, v, offset_bspoaedgev(g->edges, r));
#if DEBUG
	if(e == NULL){ // for debug
		check_all_node_edges_bspoa(g);
	}
#endif
	if(exists) *exists = exs;
	return e;
}

static inline bspoaedge_t* hadd_edge_bspoa(BSPOA *g, bspoanode_t *_u, bspoanode_t *_v, int cov, int is_aux, int *exists){
	bspoanode_t *u, *v;
	u = ref_bspoanodev(g->nodes, _u->header);
	v = ref_bspoanodev(g->nodes, _v->header);
	return add_edge_bspoa(g, u, v, cov, is_aux, exists);
}

static inline void del_edge_bspoacore(BSPOA *g, bspoanode_t *v, u4i eidx){
	bspoaedge_t *e, *p;
	u4i *eptr;
	eptr = (eidx&1)? &v->erev : &v->edge;
	e = NULL;
	while(*eptr){
		if(*eptr == eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			*eptr = e->next;
			e->next = 0;
			break;
		}
		p = ref_bspoaedgev(g->edges, *eptr);
		eptr = &p->next;
	}
	if(e == NULL){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(eidx&1){
		v->nin --;
	} else {
		v->nou --;
		push_u4v(g->ecycs, eidx);
	}
}

static inline void del_edge_bspoa(BSPOA *g, bspoanode_t *v, u4i eidx){
	bspoanode_t *w;
	bspoaedge_t *e;
	e = ref_bspoaedgev(g->edges, eidx);
	w = ref_bspoanodev(g->nodes, e->node);
	del_edge_bspoacore(g, v, eidx);
	del_edge_bspoacore(g, w, eidx ^ 0x1);
#if DEBUG
	if(eidx == 0){ // for debug
		check_all_node_edges_bspoa(g);
	}
#endif
}

static inline void mov_node_edges_bspoa(BSPOA *g, u4i src, u4i dst){
	bspoanode_t *u, *v, *w;
	bspoaedge_t *e;
	u4i eidx;
	u = ref_bspoanodev(g->nodes, src);
	v = ref_bspoanodev(g->nodes, dst);
	eidx = u->edge;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		w = ref_bspoanodev(g->nodes, e->node);
		del_edge_bspoa(g, u, offset_bspoaedgev(g->edges, e));
		add_edge_bspoa(g, v, w, e->cov, e->is_aux, NULL);
	}
	eidx = u->erev;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		w = ref_bspoanodev(g->nodes, e->node);
		del_edge_bspoa(g, w, offset_bspoaedgev(g->edges, e) - 1);
		add_edge_bspoa(g, w, v, e->cov, e->is_aux, NULL);
	}
}

static inline u4i print_aligned_nodes_bspoa(BSPOA *g, u4i nidx, FILE *out){
	bspoanode_t *v;
	u4i xidx, ret;
	ret = 0;
	xidx = nidx;
	do {
		v = ref_bspoanodev(g->nodes, xidx);
		fprintf(out, "N%u rid=%u, pos=%u, base=%u, state=%d, vst=%d, nct=%d, nin=%d, edge=%u(N%u)\n", xidx, v->rid, v->pos, v->base, v->state, v->vst, v->nct, v->nin, v->edge, v->edge? ref_bspoaedgev(g->edges, v->edge)->node : MAX_U4);
		ret ++;
		xidx = v->aligned;
	} while(xidx != nidx);
	return ret;
}

static inline bspoanode_t* merge_node_bspoa(BSPOA *g, u2i rid, bspoanode_t *x, bspoanode_t *u){
	bspoanode_t *v, *w, *y, *t;
	UNUSED(rid);
	v = x;
	y = u;
	// find y = the header of x
	do {
		if(v->base == u->base){
			y = ref_bspoanodev(g->nodes, v->header);
			break;
		}
		v = ref_bspoanodev(g->nodes, v->aligned);
	} while(v != x);
	if(y != u){ // found
		// mov edges from u to y
		mov_node_edges_bspoa(g, offset_bspoanodev(g->nodes, u), offset_bspoanodev(g->nodes, y));
		v = y;
	} else {
		v = x;
		w = ref_bspoanodev(g->nodes, x->aligned);
		while(w != x && w->base == x->base){
			v = w;
			w = ref_bspoanodev(g->nodes, w->aligned);
		}
	}
	t = u;
	while(t->aligned != offset_bspoanodev(g->nodes, u)){
		t->header = offset_bspoanodev(g->nodes, y);
		t = ref_bspoanodev(g->nodes, t->aligned);
	}
	t->aligned = v->aligned;
	v->aligned = offset_bspoanodev(g->nodes, u);
	u->header  = offset_bspoanodev(g->nodes, y);
	return y;
}

static inline void mov_nodes_header_bspoa(BSPOA *g, u4i nidx, u4i uidx){
	bspoanode_t *u, *v, *x;
	u4i xidx, found;
	if(nidx == uidx) return;
	found = 0;
	xidx = nidx;
	do {
		if(xidx == uidx) found ++;
		x = ref_bspoanodev(g->nodes, xidx);
		xidx = x->aligned;
		x->header = uidx;
	} while(xidx != nidx);
	if(found != 1){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	u = ref_bspoanodev(g->nodes, uidx);
	v = ref_bspoanodev(g->nodes, nidx);
#if DEBUG
	if(u->edge || u->erev){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	mov_node_edges_bspoa(g, nidx, uidx);
}

static inline bspoanode_t* seperate_node_bspoa(BSPOA *g, u2i rid, u4i pos){
	bspoanode_t *u, *v, *x;
	bspoaedge_t *e, *f;
	u4i nidx, xidx, eidx, found;
	nidx = g->ndoffs->buffer[rid] + pos;
	v = ref_bspoanodev(g->nodes, nidx);
	if(v->aligned == nidx){
		return v;
	}
	if(v->header == nidx){
		mov_nodes_header_bspoa(g, nidx, v->aligned);
	}
	u = ref_bspoanodev(g->nodes, v->header);
	found = 0;
	for(x=ref_bspoanodev(g->nodes, v->aligned);x->aligned!=v->aligned;x=ref_bspoanodev(g->nodes, x->aligned)){
		if(x->aligned == nidx){
			x->aligned = v->aligned;
			v->aligned = nidx;
			v->header  = nidx;
			found = 1;
			break;
		}
	}
	if(!found){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	{
		if(pos + 1U >= g->seqs->rdlens->buffer[rid]){
			xidx = BSPOA_TAIL_NODE;
		} else {
			xidx = g->ndoffs->buffer[rid] + pos + 1;
			x = ref_bspoanodev(g->nodes, xidx);
			xidx = x->header;
		}
		found = 0;
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			if(e->node == xidx){
				f = e + 1;
				e->cov --;
				f->cov --;
				if(e->cov == 0 && e->is_aux == 0){
					del_edge_bspoa(g, u, eidx);
				}
				found = 1;
				break;
			}
			eidx = e->next;
		}
		if(!found){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		x = ref_bspoanodev(g->nodes, xidx);
		add_edge_bspoa(g, v, x, 1, 0, NULL);
	}
	{
		if(pos == 0){
			xidx = BSPOA_HEAD_NODE;
		} else {
			xidx = g->ndoffs->buffer[rid] + pos - 1;
			x = ref_bspoanodev(g->nodes, xidx);
			xidx = x->header;
		}
		found = 0;
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			if(e->node == xidx){
				f = e - 1;
				e->cov --;
				f->cov --;
				if(e->cov == 0 && e->is_aux == 0){
					del_edge_bspoa(g, u, eidx);
				}
				found = 1;
				break;
			}
			eidx = e->next;
		}
		if(!found){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		x = ref_bspoanodev(g->nodes, xidx);
		add_edge_bspoa(g, x, v, 1, 0, NULL);
	}
	return v;
}

static inline void del_read_nodes_bspoa(BSPOA *g, u2i rid, u4i rb, u4i re){
	bspoanode_t *u, *v, *x;
	bspoaedge_t *e;
	u4i nidx, nide, xidx, rdlen, eidx, brk, dir;
	rdlen = g->seqs->rdlens->buffer[rid];
	nidx = g->ndoffs->buffer[rid];
	if(re >= rdlen){
		u = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
		re = rdlen;
		nidx = nidx + re - 1;
	} else {
		nidx = nidx + re - 1;
		u = ref_bspoanodev(g->nodes, nidx + 1);
		u = ref_bspoanodev(g->nodes, u->header);
	}
	if(rb >= re) return;
	nide = g->ndoffs->buffer[rid] + rb - 1;
	for(;;nidx--){
		if(nidx == nide){
			if(rb == 0){
				nidx = BSPOA_HEAD_NODE;
			}
			brk = 1;
		} else {
			brk = 0;
		}
		v = ref_bspoanodev(g->nodes, nidx);
		xidx = v->header;
		x = ref_bspoanodev(g->nodes, xidx);
		// find edge from u -> x
		if(u->erev == 0){
			if(brk){
				break;
			} else {
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
		eidx = x->edge;
		e = NULL;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			if(e->node == offset_bspoanodev(g->nodes, u)){
				break;
			}
			eidx = e->next;
		}
#if DEBUG
		if(u->nin == 0 || eidx == 0 || e->cov == 0){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
#endif
		e->cov --;
		(e + 1)->cov --;
		//if(e->cov == 0 && e->is_aux == 0){ // keep aux edge
		if(e->cov == 0){
			del_edge_bspoa(g, x, eidx);
		}
		u = x;
		if(brk) break;
		// remove v from aligned loop
		while(x->aligned != nidx){
			x = ref_bspoanodev(g->nodes, x->aligned);
		}
		x->aligned = v->aligned;
		if(v->header == nidx){
			x = ref_bspoanodev(g->nodes, v->aligned);
			if(x != v && x->base == v->base){
				u = x;
				mov_node_edges_bspoa(g, nidx, v->aligned);
				while(x->base == v->base){
					x->header = v->aligned;
					if(x->aligned == v->aligned) break;
					x = ref_bspoanodev(g->nodes, x->aligned);
				}
			} else {
				for(dir=0;dir<2;dir++){
					eidx = dir? v->erev : v->edge;
					while(eidx){
						e = ref_bspoaedgev(g->edges, eidx);
						eidx = e->next;
						if(e->is_aux){
							del_edge_bspoa(g, v, offset_bspoaedgev(g->edges, e));
						}
					}
				}
			}
		}
		v->aligned = nidx;
#if DEBUG
		if(v->edge){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
#endif
	}
}

static inline void beg_bspoacore(BSPOA *g, u1i *cns, u4i len, int clear_all){
	bspoanode_t *head, *tail, *u, *v;
	bspoaedge_t *e;
	u4i i, bb;
	g->bwtrigger = 0;
	g->ncall ++;
	clear_u4v(g->ndoffs);
	clear_and_encap_bspoanodev(g->nodes, 2 + len);
	clear_bspoaedgev(g->edges);
	ZEROS(next_ref_bspoaedgev(g->edges));
	ZEROS(next_ref_bspoaedgev(g->edges)); // even idx for edge
	clear_u4v(g->ecycs);
	head = next_ref_bspoanodev(g->nodes);
	ZEROS(head);
	head->base = 4;
	head->header  = BSPOA_HEAD_NODE;
	head->aligned = BSPOA_HEAD_NODE;
	head->cpos = 0;
	tail = next_ref_bspoanodev(g->nodes);
	ZEROS(tail);
	tail->base = 4;
	tail->header  = BSPOA_TAIL_NODE;
	tail->aligned = BSPOA_TAIL_NODE;
	tail->cpos = len;
	u = head;
	// add backbone nodes
	g->backbone = len;
	for(i=0;i<len;i++){
		bb = cns[i];
		v = add_node_bspoa(g, BSPOA_RDCNT_MAX, i, bb, 0);
		v->ref  = 1;
		v->aux  = 1;
		v->cpos = i;
		e = add_edge_bspoa(g, u, v, 0, 1, NULL);
		u = v;
	}
	e = add_edge_bspoa(g, u, tail, 0, 1, NULL);
	if(g->par->refmode == 0){
		clear_u4v(g->cigars);
		clear_u4v(g->cgbs);
		clear_u4v(g->cges);
	}
	if(clear_all){
		clear_seqbank(g->seqs);
		clear_u1v(g->cns);
		clear_u1v(g->qlt);
		clear_u1v(g->alt);
	}
}

static inline void beg_bspoa(BSPOA *g){
	if((g->ncall % 16) == 0){
		renew_bspoa(g);
	}
	beg_bspoacore(g, NULL, 0, 1);
}

static inline b1i* dpalign_row_prepare_data(BSPOA *g, u4i mmidx, b1i **us, b1i **es, b1i **qs, int **ubegs){
	us[0] = g->memp->buffer + mmidx * g->mmblk;
	es[0] = (g->piecewise > 0)? us[0] + g->bandwidth : NULL;
	qs[0] = (g->piecewise > 1)? es[0] + g->bandwidth : NULL;
	ubegs[0] = (int*)(us[0] + g->bandwidth * (g->piecewise + 1));
	return ((b1i*)ubegs[0]) + (WORDSIZE + 1) * sizeof(int);
}

static inline void prepare_rd_align_bspoa(BSPOA *g, u2i rid){
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	seqalign_result_t rs;
	u4i seqoff, seqlen, *cgs, ncg;
	u4i nidx, eidx, i, j, x, y, op, sz, tb, te;
	u4i *rmap;
	int rpos, exists;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	if(seqlen == 0) return;
	g->qlen = g->slen = seqlen;
	g->qb = 0;
	g->qe = g->qlen;
	clear_and_encap_u1v(g->qseq, g->qlen);
	bitseq_basebank(g->seqs->rdseqs, seqoff, g->qlen, g->qseq->buffer);
	g->qseq->size = g->qlen;
	rmap = NULL;
	cgs  = NULL;
	ncg  = 0;
	tb = 0; te = g->cns->size;
	x = y = 0;
	if(g->par->refmode && g->backbone && g->cges->buffer[rid] > g->cgbs->buffer[rid]){
		if(g->par->bandwidth == 0){
			g->bandwidth = roundup_times(seqlen, WORDSIZE);
		} else {
			g->bandwidth = num_min(g->par->bandwidth, Int(seqlen));
			g->bandwidth = roundup_times(g->bandwidth, WORDSIZE);
		}
		cgs = g->cigars->buffer + g->cgbs->buffer[rid];
		ncg = g->cges->buffer[rid] - g->cgbs->buffer[rid];
		// find tailing margin, see 'D/N/H' cigar in SAM
		x = y = 0;
		for(i=0;i<ncg;i++){
			if((cgs[i]&0xf) == 2 || (cgs[i]&0xf) == 3 || (cgs[i]&0xf) == 5){ // D/N/H
				y += cgs[i] >> 4;
			} else if((cgs[i]&0xf) == 1 || (cgs[i]&0xf) == 4){ // I/S
				x += cgs[i] >> 4;
			} else {
				break;
			}
		}
		cgs += i; ncg -= i;
		g->qb = x;
		tb = y;
		x = y = 0;
		for(i=ncg;i;i--){
			if((cgs[i - 1]&0xf) == 2 || (cgs[i - 1]&0xf) == 3 || (cgs[i - 1]&0xf) == 5){ // D/N/H
				y += cgs[i - 1] >> 4;
			} else if((cgs[i - 1]&0xf) == 1 || (cgs[i - 1]&0xf) == 4){ // I/S
				x += cgs[i] >> 4;
			} else {
				break;
			}
		}
		ncg = i;
		g->qe = g->qlen - x;
		g->slen = g->qe - g->qb;
		te = g->backbone - y;
		x = 0;
		y = tb;
		tb = (tb >= (g->bandwidth / 2))? tb - g->bandwidth / 4 : 0;
		te = (g->cns->size - te >= g->bandwidth / 2)? te + g->bandwidth / 4 : (g->cns->size);
		clear_and_encap_b1v(g->memp, (g->cns->size + 1) * sizeof(u4i));
		rmap = (u4i*)(g->memp->buffer);
	} else if(g->bwtrigger && Int(roundup_times(seqlen, WORDSIZE)) > g->par->bandwidth){
		if(g->par->bandwidth == 0){
			g->bandwidth = roundup_times(seqlen, WORDSIZE);
		} else {
			g->bandwidth = num_min(g->par->bandwidth, Int(seqlen));
			g->bandwidth = roundup_times(g->bandwidth, WORDSIZE);
		}
		if(g->par->ksz){
			rs = kmer_striped_seqedit_pairwise(g->par->ksz, g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->memp, g->stack, 0);
		} else {
			rs = striped_seqedit_pairwise(g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->par->alnmode, 0, g->memp, g->stack, 0);
		}
		//rs = banded_striped_epi8_seqalign_pairwise(g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->memp, g->stack, SEQALIGN_MODE_GLOBAL, g->qseq->size, g->matrix[0], 0, g->par->E, 0, 0, 0);
		if(_DEBUG_LOG_){
			char *alnstr[3];
			alnstr[0] = malloc(rs.aln + 1);
			alnstr[1] = malloc(rs.aln + 1);
			alnstr[2] = malloc(rs.aln + 1);
			seqalign_cigar2alnstr(g->qseq->buffer, g->cns->buffer, &rs, g->stack, alnstr, rs.aln);
			fprintf(stderr, "#RID%d\t%d\t%d\t%d\tCNS\t%d\t%d\t%d\tmat=%d\taln=%d\n", rid, Int(g->qseq->size), rs.qb, rs.qe, Int(g->cns->size), rs.tb, rs.te, rs.mat, rs.aln);
			fprintf(stderr, "#%s\n#%s\n#%s\n", alnstr[0], alnstr[2], alnstr[1]);
			free(alnstr[0]);
			free(alnstr[1]);
			free(alnstr[2]);
		}
		g->qb = rs.qb;
		g->qe = rs.qe;
		g->slen = g->qe - g->qb;
		tb = (rs.tb >= Int(g->bandwidth / 2))? rs.tb - g->bandwidth / 4 : 0;
		te = (g->cns->size - rs.te >= g->bandwidth / 2)? Int(rs.te + g->bandwidth / 4) : Int(g->cns->size);
		// cigar2rpos
		clear_and_encap_b1v(g->memp, (g->cns->size + 1) * sizeof(u4i));
		rmap = (u4i*)(g->memp->buffer);
		cgs = g->stack->buffer;
		ncg = g->stack->size;
		//x = rs.qb;
		x = 0;
		y = rs.tb;
	} else {
		g->bandwidth = roundup_times(seqlen, WORDSIZE);
	}
	g->bandwidth = num_min(Int(g->bandwidth), g->par->bwmax);
	if(rmap && cgs && ncg){
		rmap[0] = 0;
		for(i=1;i<y;i++) rmap[i] = i * g->qb / (y + 1);
		for(i=0;i<ncg;i++){
			op = cgs[i] & 0xf;
			sz = cgs[i] >> 4;
			switch(op){
				case 0:
				case 7:
				case 8: // match/mistmatch
				for(j=0;j<sz;j++){
					rmap[y++] = x ++;
				}
				break;
				case 1:
				case 4: // insertion
				x += sz;
				break;
				case 2:
				case 3: // deletion
				for(j=0;j<sz;j++){
					rmap[y++] = x;
				}
				break;
				case 5: // H, used to locate the begin and end
				for(j=0;j<sz;j++){
					rmap[y++] = x;
				}
				break;
			}
		}
		for(i=y;i<g->cns->size;i++){
			rmap[i] = x + (i - y + 1) * (g->slen - x) / (g->cns->size - y + 1);
		}
#if DEBUG && 0
		for(i=1;i<g->cns->size;i++){
			if(rmap[i] < rmap[i-1]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
#endif
		rmap[g->cns->size] = g->slen;
		if(_DEBUG_LOG_ > 1){
			fprintf(stderr, "RMAP:");
			for(i=0;i<=g->cns->size;i++){
				fprintf(stderr, "\t[%d]%d", i, rmap[i]);
			}
			fprintf(stderr, "\n");
		}
		// set rpos
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			rpos = rmap[u->cpos] - g->bandwidth / 2;
			if(rpos < 0) rpos = 0;
			else if(rpos + g->bandwidth > g->qlen){
				rpos = g->qlen - g->bandwidth;
			}
			u->rpos = rpos;
			if(u->cpos == tb && tb){
				e = hadd_edge_bspoa(g, ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE), u, 0, 1, &exists);
				if(!exists) push_u8v(g->todels, (((u8i)BSPOA_HEAD_NODE) << 32) | offset_bspoaedgev(g->edges, e));
				tb = 0;
			}
			if(u->cpos == te && te != g->cns->size){
				e = hadd_edge_bspoa(g, u, ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE), 0, 1, &exists);
				if(!exists) push_u8v(g->todels, (((u8i)u->header) << 32) | offset_bspoaedgev(g->edges, e));
				te = g->cns->size;
			}
		}
	} else {
		g->bandwidth = roundup_times(seqlen, WORDSIZE);
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			u->rpos = 0;
		}
	}
#if DEBUG
	if(rid == 10000){
		fprint_local_dot_bspoa(g, BSPOA_TAIL_NODE, 50, "1.dot", NULL);
	}
#endif
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[0], g->par->M, g->par->X);
	clear_and_encap_b1v(g->qprof[0], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[0]->buffer, g->bandwidth, g->matrix[0]);
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[1], g->par->M + g->par->refbonus, g->par->X);
	clear_and_encap_b1v(g->qprof[1], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[1]->buffer, g->bandwidth, g->matrix[1]);
	// calculate nct
	for(i=0;i<g->nodes->size;i++){
		v = ref_bspoanodev(g->nodes, i);
		v->mmidx = 0;
		v->vst  = 0;
		v->cov  = 0;
		if(g->par->nrec == 0){
			v->state = 1;
			v->nct = v->nin;
		} else if(i < UInt(2 + g->backbone)){
			v->state = 1;
			v->nct = 0;
		} else {
			v->state = 0;
			v->nct = 0;
		}
	}
	for(i=0;i<g->nodes->size;i++){
		v = ref_bspoanodev(g->nodes, i);
		if(v->cov) continue;
		j = 1;
		u = v;
		while(u->aligned != i){
			u = ref_bspoanodev(g->nodes, u->aligned);
			j ++;
		}
		u = v;
		u->cov = j;
		while(u->aligned != i){
			u = ref_bspoanodev(g->nodes, u->aligned);
			u->cov = j;
		}
	}
	if(g->par->nrec){
		ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE)->state = 1;
		ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE)->state = 1;
		// latest nrec reads
		for(i=0;i<g->nodes->size;i++){
			v = ref_bspoanodev(g->nodes, i);
			if(v->rid + g->par->nrec >= rid){
				u = ref_bspoanodev(g->nodes, v->header);
				u->state = 1;
			}
		}
		// connect read's beg and end with HEAD and TAIL
		for(i=0;i+g->par->nrec<rid;i++){
			for(j=0;j<g->seqs->rdlens->buffer[i];j++){
				v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[i] + j);
				u = ref_bspoanodev(g->nodes, v->header);
				if(u->state) break;
				u->state = 1;
			}
			for(j=1;j<=g->seqs->rdlens->buffer[i];j++){
				v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[i] + g->seqs->rdlens->buffer[i] - j);
				u = ref_bspoanodev(g->nodes, v->header);
				if(u->state) break;
				u->state = 1;
			}
		}
		// find nodes can be visited from HEAD
		clear_u4v(g->stack);
		ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE)->vst = 1;
		push_u4v(g->stack, BSPOA_HEAD_NODE);
		while(pop_u4v(g->stack, &nidx)){
			u = ref_bspoanodev(g->nodes, nidx);
			eidx = u->edge;
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				v = ref_bspoanodev(g->nodes, e->node);
				if(v->state == 0) continue;
				if(v->vst) continue;
				v->vst = 1;
				push_u4v(g->stack, e->node);
			}
		}
		// reverse calculate nct(number of living input edges, visited by BSPOAPar->nrec latest reads)
		for(i=0;i<g->nodes->size;i++){
			u = ref_bspoanodev(g->nodes, i);
			if(u->vst == 0) continue;
			u->vst = 0; // set all vst to 0
			if(u->state == 0) continue;
			eidx = u->edge;
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				v = ref_bspoanodev(g->nodes, e->node);
				if(v->state){
					v->nct ++;
				}
			}
		}
	}
	// prepare space for nodes and edges
	encap_bspoanodev(g->nodes, seqlen + 2);
	encap_bspoaedgev(g->edges, seqlen + 2);
	for(i=0;i<seqlen;i++){
		v = add_node_bspoa(g, rid, i, g->qseq->buffer[i], 0);
		v->cpos = i;
	}
	g->piecewise = banded_striped_epi8_seqalign_get_piecewise(g->par->O, g->par->E, g->par->Q, g->par->P, g->bandwidth);
	g->mmblk = roundup_times(g->bandwidth * (g->piecewise + 1) + (WORDSIZE + 1) * sizeof(int), WORDSIZE); // us, es, qs, ubegs, qoff, max_score, max_index
	g->mmcnt = 2; // reserve two blocks
	for(i=0;i<g->nodes->size;i++){
		v = ref_bspoanodev(g->nodes, i);
		if(v->state){
			v->mmidx = g->mmcnt ++;
		}
	}
	clear_and_encap_b1v(g->memp, g->mmcnt * g->mmblk);
	u = ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE);
	{
		b1i *us, *es, *qs;
		int *ubegs;
		dpalign_row_prepare_data(g, u->mmidx, &us, &es, &qs, &ubegs);
		banded_striped_epi8_seqalign_piecex_row_init(us, es, qs, ubegs, NULL, g->par->alnmode, g->bandwidth, g->par->M, g->par->X, g->par->O, g->par->E, g->par->Q, g->par->P);
	}
	g->maxscr = SEQALIGN_SCORE_MIN;
	g->maxidx = -1;
	g->maxoff = -1;
}

static inline void dpalign_row_update_bspoa(BSPOA *g, b1i *qp, u4i mmidx1, u4i mmidx2, u4i toff, u4i qoff1, u4i qoff2, u1i next_base){
	b1i *us[2], *es[2], *qs[2];
	int W, rh, *ubegs[2];
	W = g->bandwidth / WORDSIZE;
	dpalign_row_prepare_data(g,      0, us + 0, es + 0, qs + 0, ubegs + 0);
	__builtin_prefetch(us[0], 0);
	dpalign_row_prepare_data(g, mmidx1, us + 1, es + 1, qs + 1, ubegs + 1);
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_movx(us + 0, es + 0, qs + 0, ubegs + 0, W, qoff2 - qoff1, g->piecewise, g->par->M, g->par->X, g->par->O, g->par->E, g->par->Q, g->par->P);
	dpalign_row_prepare_data(g, mmidx2, us + 1, es + 1, qs + 1, ubegs + 1);
	if(qoff1 == qoff2){
		if(qoff1){
			rh = SEQALIGN_SCORE_MIN;
		} else {
			if(seqalign_mode_type(g->par->alnmode) == SEQALIGN_MODE_OVERLAP || toff == 0) rh = 0;
			else if(g->piecewise < 2) rh = g->par->O + g->par->E * toff;
			else rh = num_max(g->par->O + g->par->E * toff, g->par->Q + g->par->P * toff);
		}
	} else if(qoff1 + W * WORDSIZE >= qoff2){
		rh = ubegs[0][0]; // movx -> aligned
	} else {
		rh = SEQALIGN_SCORE_MIN;
	}
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_cal(qoff2, next_base, us, es, qs, ubegs, qp, g->par->O, g->par->E, g->par->Q, g->par->P, W, qoff2 - qoff1, rh, g->piecewise);
#if 0
	dpalign_row_prepare_data(g, mmidx1, us + 0, es + 0, qs + 0, ubegs + 0);
	banded_striped_epi8_seqalign_piecex_row_verify(toff, qoff1, g->par->alnmode, W, qoff2 - qoff1, next_base, us, es, qs, ubegs, qp, g->par->O, g->par->E, g->par->Q, g->par->P);
#endif
}

static inline void dpalign_row_merge_bspoa(BSPOA *g, u4i mmidx1, u4i mmidx2){
	b1i *us[3], *es[3], *qs[3];
	int *ubegs[3];
	dpalign_row_prepare_data(g, mmidx1, us + 0, es + 0, qs + 0, ubegs + 0);
	__builtin_prefetch(us[0], 0);
	dpalign_row_prepare_data(g, mmidx2, us + 1, es + 1, qs + 1, ubegs + 1);
	__builtin_prefetch(us[1], 1);
	dpalign_row_prepare_data(g, mmidx2, us + 2, es + 2, qs + 2, ubegs + 2);
	banded_striped_epi8_seqalign_piecex_row_merge(us, es, qs, ubegs, g->bandwidth / WORDSIZE, g->piecewise);
}

static inline u1i get_rdbase_bspoa(BSPOA *g, u4i rid, u4i pos){
	u4i seqoff;
	seqoff = g->seqs->rdoffs->buffer[rid];
	return get_basebank(g->seqs->rdseqs, seqoff + pos);
}

// Hs[1] stores the beg H score, which comes from E or Q
// Hs[0] return the end H score
static inline u4i backcal_del_bspoa(BSPOA *g, u4i node, int qb, int Hs[2], int bt, u8v *heap){
	typedef struct { u8i nidx:24, step:16; int qscr:24; } bcal_t;
	bspoanode_t *v;
	bspoaedge_t *e;
	bcal_t BC, CC;
	b1i *us, *es, *qs;
	int *ubegs, q;
	u4i nidx, eidx, bsoff;
	bsoff = banded_striped_epi8_pos2idx(g->bandwidth, qb);
	clear_u8v(heap);
	BC.nidx = node;
	BC.qscr = Hs[1];
	BC.step = 1;
	array_heap_push(heap->buffer, heap->size, heap->cap, bcal_t, BC, num_cmp(a.step, b.step));
	while(heap->size){
		BC = ((bcal_t*)heap->buffer)[0];
		array_heap_remove((bcal_t*)heap->buffer, heap->size, heap->cap, bcal_t, 0, num_cmp(a.step, b.step));
		v = ref_bspoanodev(g->nodes, nidx);
		dpalign_row_prepare_data(g, v->mmidx, &us, &es, &qs, &ubegs);
		Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, qb);
		if(bt == SEQALIGN_BT_D){
			q = es[bsoff];
			if(Hs[0] + q != BC.qscr) continue;
			if(q == g->par->O + g->par->E){
				return BC.nidx;
			}
		} else {
			q = qs[bsoff];
			if(Hs[0] + q != BC.qscr) continue;
			if(q == g->par->P + g->par->Q){
				return BC.nidx;
			}
		}
		eidx = v->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			CC.nidx = e->node;
			CC.qscr = BC.qscr - ((bt == SEQALIGN_BT_D)? g->par->E : g->par->Q);
			CC.step = BC.step + 1;
			array_heap_push(heap->buffer, heap->size, heap->cap, bcal_t, CC, num_cmp(a.step, b.step));
		}
	}
	fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	abort();
	return BSPOA_HEAD_NODE;
}

static inline int _alignment2graph_bspoa(BSPOA *g, u4i rid, u4i midx, int xe, int *_mat, int *_tail, String **alnstrs){
	bspoanode_t *u, *v, *w, *n;
	bspoaedge_t *e;
	b1i *us, *es, *qs, q;
	u4i nidx, eidx, i, cpos, W, bt, ft;
	int x, xb, *ubegs, Hs[3], t, s, realn, mat, tail;
	realn = 0;
	W = g->bandwidth / WORDSIZE;
	x = xe;
	if(alnstrs){
		clear_string(alnstrs[0]);
		clear_string(alnstrs[1]);
		clear_string(alnstrs[2]);
	}
	v = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
	cpos = v->cpos;
	mat = 0;
	tail = g->qlen - x - g->qb - 1;
	if(x + g->qb + 1 < g->qlen){
		for(i=g->qlen-1;i>x+g->qb;i--){
			u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + i);
			u->cpos = cpos;
			add_edge_bspoa(g, u, v, 1, 0, NULL);
			if(0 && alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], '-');
				add_char_string(alnstrs[2], '-');
			}
			v = u;
		}
	}
	xb = x;
	nidx = midx;
	bt = MAX_U4; // to be calculated
	n = ref_bspoanodev(g->nodes, nidx);
	dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
	Hs[0] = 0;
	Hs[1] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - n->rpos);
	Hs[2] = 0;
	while(1){
		if(nidx == BSPOA_HEAD_NODE || x < 0){
			xb = x;
			x += g->qb;
			while(x >= 0){
				w = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + x);
				w->cpos = cpos;
				e = add_edge_bspoa(g, w, v, 1, 0, NULL);
				if(0 && alnstrs){
					add_char_string(alnstrs[0], bit_base_table[w->base&0x03]);
					add_char_string(alnstrs[1], '-');
					add_char_string(alnstrs[2], '-');
				}
				x --;
				tail ++;
				v = w;
			}
			{
				w = ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE);
				e = add_edge_bspoa(g, w, v, 1, 0, NULL);
				v = w;
			}
			break;
		}
		if(bt == SEQALIGN_BT_D || bt == SEQALIGN_BT2_D2){
			if(alnstrs){
				add_char_string(alnstrs[0], '-');
				add_char_string(alnstrs[1], bit_base_table[n->base&0x03]);
				add_char_string(alnstrs[2], '-');
			}
			eidx = n->erev;
#if DEBUG
			int found = 0;
#endif
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				w = ref_bspoanodev(g->nodes, e->node);
				if(w->state == 0) continue;
				if(x < Int(w->rpos) || x >= Int(w->rpos + g->bandwidth)) continue;
				dpalign_row_prepare_data(g, w->mmidx, &us, &es, &qs, &ubegs);
				Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos);
				if(bt == SEQALIGN_BT_D){
					q = g->piecewise? es[banded_striped_epi8_pos2idx(g->bandwidth, x - w->rpos)] : g->par->O + g->par->E;
				} else {
					q = qs[banded_striped_epi8_pos2idx(g->bandwidth, x - w->rpos)];
				}
				if(Hs[0] + q != Hs[1]) continue;
				nidx = e->node;
				n = ref_bspoanodev(g->nodes, nidx);
				if(q == ((bt == SEQALIGN_BT_D)? g->par->O + g->par->E : g->par->Q + g->par->P)){
					bt = MAX_U4;
					Hs[1] = Hs[0];
					Hs[2] = 0;
				} else {
					Hs[1] -= ((bt == SEQALIGN_BT_D)? g->par->E : g->par->P);
					Hs[2] ++;
				}
#if DEBUG
				found = 1;
#endif
				break;
			}
#if DEBUG
			if(!found){ // Not found
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			continue;
		} else if(bt == SEQALIGN_BT_I || bt == SEQALIGN_BT2_I2){
			if(g->piecewise == 2){
				t = num_max(g->par->O + g->par->E * Hs[2], g->par->Q + g->par->P * Hs[2]);
			} else {
				t = g->par->O + g->par->E * Hs[2];
			}
#if DEBUG
			if(x < Int(n->rpos)){ // should never happen
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + x + g->qb);
			u->cpos = cpos;
			e = add_edge_bspoa(g, u, v, 1, 0, NULL);
			v = u;
			if(alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], '-');
				add_char_string(alnstrs[2], '-');
			}
			x --;
			if(Hs[0] + t == Hs[1]){
				bt = MAX_U4;
				Hs[1] = Hs[0];
				Hs[2] = 0;
#if DEBUG
			} else if(Hs[0] + t > Hs[1]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
#endif
			} else if(x >= 0){
				Hs[0] -= us[banded_striped_epi8_pos2idx(g->bandwidth, x - n->rpos)];
				Hs[2] ++;
			} else {
				// the alignment starts with insertion
			}
			continue;
		} else if(bt == SEQALIGN_BT_M){
			u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + x + g->qb);
			u->cpos = cpos = n->cpos;
			e = add_edge_bspoa(g, u, v, 1, 0, NULL);
			if(alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], offset_bspoanodev(g->nodes, n)? bit_base_table[n->base&0x03] : '^');
				add_char_string(alnstrs[2], "*|"[(u->base&0x03) == (n->base&0x03)]);
			}
			x --;
			if(offset_bspoanodev(g->nodes, n) > BSPOA_TAIL_NODE && u->base == n->base){
				u = merge_node_bspoa(g, rid, n, u);
				mat ++;
			}
			v = u;
			n = ref_bspoanodev(g->nodes, nidx);
			bt = MAX_U4; // to be calculated
		} else {
#if DEBUG
			if(0){
				if(offset_bspoanodev(g->nodes, n) != nidx){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
				Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - n->rpos);
				if(Hs[0] != Hs[1]){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
#endif
			s = g->matrix[n->ref][g->qseq->buffer[x + g->qb] * 4 + n->base];
			eidx = n->erev;
			clear_u4v(g->stack);
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				w = ref_bspoanodev(g->nodes, e->node);
				if(w->state == 0) continue;
				dpalign_row_prepare_data(g, w->mmidx, &us, &es, &qs, &ubegs);
				ft = 0;
				if(x < Int(w->rpos) || x > Int(g->bandwidth + w->rpos)){
					continue;
				} else if(x == Int(g->bandwidth + w->rpos)){
					Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos - 1);
					ft |= 1 << SEQALIGN_BT_D;
					ft |= 1 << SEQALIGN_BT2_D2;
				} else if(x == Int(w->rpos)){
					if(w->rpos == 0 && (seqalign_mode_type(g->par->alnmode) == SEQALIGN_MODE_OVERLAP || offset_bspoanodev(g->nodes, w) == BSPOA_HEAD_NODE)){
						Hs[0] = 0;
					} else {
						Hs[0] = ubegs[0]; // useless
						ft |= 1 << SEQALIGN_BT_M; // forbid M
					}
				} else {
					Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos - 1);
				}
				t = ((x - w->rpos) % W) * WORDSIZE + ((x - w->rpos) / W);
				push_u4v(g->stack, e->node);
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT_M))? UInt(SEQALIGN_SCORE_MIN) : UInt(Hs[0] + s));
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT_D))? UInt(SEQALIGN_SCORE_MIN) : UInt(Hs[0] + us[t] + (es? es[t] : g->par->E)));
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT2_D2))? UInt(SEQALIGN_SCORE_MIN) : UInt(Hs[0] + (qs? us[t] + qs[t] : SEQALIGN_SCORE_MAX)));
			}
			bt = SEQALIGN_BT_I;
			for(i=0;i<g->stack->size;i+=4){
				if(Hs[1] == Int(g->stack->buffer[i + 1])){
					bt = SEQALIGN_BT_M;
					nidx = g->stack->buffer[i + 0];
					Hs[1] -= s;
					Hs[2] = 0;
					break;
				}
			}
			if(bt == SEQALIGN_BT_I){
				for(i=0;i<g->stack->size;i+=4){
					if(Hs[1] == Int(g->stack->buffer[i + 2])){
						bt = SEQALIGN_BT_D;
						//nidx = g->stack->buffer[i + 0];
						//n = ref_bspoanodev(g->nodes, nidx);
						//dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
						Hs[2] = 1;
						break;
					} else if(Hs[1] == Int(g->stack->buffer[i + 3])){
						bt = SEQALIGN_BT2_D2;
						//nidx = g->stack->buffer[i + 0];
						//n = ref_bspoanodev(g->nodes, nidx);
						//dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
						Hs[2] = 1;
						break;
					}
				}
			}
			if(bt == SEQALIGN_BT_I){
				Hs[2] = 1;
				dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
				Hs[0] = Hs[1] - us[banded_striped_epi8_pos2idx(g->bandwidth, x - n->rpos)];
			}
		}
	}
	*_mat  = mat;
	*_tail = tail;
	return xb;
}

static inline int align_rd_bspoacore(BSPOA *g, u2i rid){
	bspoanode_t *u, *v;
	bspoaedge_t *e, *f;
	u4i seqlen, nidx, eidx, mmidx;
	b1i *us, *es, *qs;
	int *ubegs, smax, rmax, maxoff;
	seqlen = g->seqs->rdlens->buffer[rid];
	if(seqlen == 0) return 0;
	prepare_rd_align_bspoa(g, rid);
	clear_u4v(g->stack);
	push_u4v(g->stack, BSPOA_HEAD_NODE);
	g->nodes->buffer[BSPOA_HEAD_NODE].mpos = 0;
	for(nidx=1;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->mpos = MAX_U4;
	}
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			f = e + 1; // v -> u
			eidx = e->next;
			//if(e->node == BSPOA_TAIL_NODE) tail = 1;
			v = ref_bspoanodev(g->nodes, e->node);
			if(u->mpos + 1 < v->mpos){
				v->mpos = u->mpos + 1;
			}
			if(v->state == 0) continue;
			if(e->node == BSPOA_TAIL_NODE){
				dpalign_row_prepare_data(g, u->mmidx, &us, &es, &qs, &ubegs);
				if(seqalign_mode_type(g->par->alnmode) == SEQALIGN_MODE_GLOBAL){
					maxoff = num_min(g->slen, u->rpos + g->bandwidth) - 1;
					smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, maxoff - u->rpos);
					if(smax > g->maxscr){
						g->maxscr = smax;
						g->maxidx = offset_bspoanodev(g->nodes, u);
						g->maxoff = maxoff;
					}
				} else {
					rmax = banded_striped_epi8_seqalign_row_max(us, ubegs, g->bandwidth / WORDSIZE, &smax);
					if(smax > g->maxscr){
						g->maxscr = smax;
						g->maxidx = offset_bspoanodev(g->nodes, u);
						g->maxoff = rmax + u->rpos;
					}
				}
				v->vst ++;
			} else {
				mmidx = v->vst? 1 : v->mmidx;
#if DEBUG
				if(v->rpos < u->rpos){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
#endif
				dpalign_row_update_bspoa(g, g->qprof[v->ref]->buffer, u->mmidx, mmidx, v->mpos, u->rpos, v->rpos, v->base);
				if(v->vst){
					dpalign_row_merge_bspoa(g, mmidx, v->mmidx);
				}
				v->vst ++;
				if(v->vst == v->nct){
					if(seqalign_mode_type(g->par->alnmode) != SEQALIGN_MODE_GLOBAL && v->rpos + g->bandwidth >= g->slen){
						dpalign_row_prepare_data(g, v->mmidx, &us, &es, &qs, &ubegs);
						smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, g->slen - 1 - v->rpos);
						if(smax > g->maxscr){
							g->maxscr = smax;
							g->maxidx = offset_bspoanodev(g->nodes, v);
							g->maxoff = g->slen - 1;
						}
					}
					push_u4v(g->stack, e->node);
				}
			}
		}
	}
	v = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
#if DEBUG
	if(v->vst != v->nct){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	return g->maxscr;
}

static inline int align_rd_bspoa(BSPOA *g, u2i rid){
	bspoaedge_t *e;
	String **alnstrs;
	u4i i;
	int xb, xe, rlen, score, mat, bad;
	clear_u8v(g->todels);
	score = align_rd_bspoacore(g, rid);
	rlen = g->seqs->rdlens->buffer[rid];
	if(rlen == 0) return 0;
	xe = g->maxoff;
	xb = 0;
	alnstrs = NULL;
	if(_DEBUG_LOG_){
		alnstrs = malloc(3 * sizeof(String*));
		alnstrs[0] = init_string(rlen * 2);
		alnstrs[1] = init_string(rlen * 2);
		alnstrs[2] = init_string(rlen * 2);
	}
	xb = _alignment2graph_bspoa(g, rid, g->maxidx, xe, &mat, &bad, alnstrs);
	xb += g->qb;
	xe += g->qb;
	for(i=0;i<g->todels->size;i++){
		e = ref_bspoaedgev(g->edges, g->todels->buffer[i] & MAX_U4);
		if(e->is_aux == 0) continue;
		del_edge_bspoa(g, ref_bspoanodev(g->nodes, g->todels->buffer[i] >> 32), g->todels->buffer[i] & MAX_U4);
	}
	clear_u8v(g->todels);
	if(alnstrs){
		reverse_string(alnstrs[0]);
		reverse_string(alnstrs[1]);
		reverse_string(alnstrs[2]);
		fprintf(stderr, "ALIGN[%03d] len=%u band=%d aligned=%d,%d mat=%d,%0.3f tail=%d score=%d\n", rid, rlen, g->bandwidth, xb + 1, xe + 1, mat, 1.0 * mat / rlen, bad, score);
		fprintf(stderr, "%s\n%s\n%s\n", alnstrs[0]->string, alnstrs[2]->string, alnstrs[1]->string);
		free_string(alnstrs[0]);
		free_string(alnstrs[1]);
		free_string(alnstrs[2]);
		free(alnstrs);
	}
	return score;
}

static inline void check_dup_edges_bspoa(BSPOA *g){
	bspoanode_t *u;
	bspoaedge_t *e;
	u32hash *hash;
	u4i nidx, eidx, *t;
	int exists;
	hash = init_u32hash(13);
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		clear_u32hash(hash);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			t = prepare_u32hash(hash, e->node, &exists);
			if(exists){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			} else {
				*t = e->node;
			}
		}
	}
	free_u32hash(hash);
}

static inline u4i sort_nodes_bspoa(BSPOA *g){
	bspoanode_t *v, *u, *x;
	bspoaedge_t *e;
	u4i nseq, mrow, nidx, eidx, xidx, moff, ready, nou;
	nseq = g->nrds;
	mrow = nseq + 3;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst   = 0;
		u->state = 0;
		u->mpos  = 0;
	}
#if DEBUG
	for(eidx=0;eidx<g->edges->size;eidx++){
		e = ref_bspoaedgev(g->edges, eidx);
		e->vst   = 0;
	}
#endif
	clear_u4v(g->stack);
	nidx = BSPOA_TAIL_NODE;
	push_u4v(g->stack, nidx);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
#if DEBUG
			e->vst = 1;
			(e-1)->vst = 1;
#endif
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			if(u->mpos + 1 > v->mpos){
				v->mpos = u->mpos + 1;
			}
			v->vst ++;
			if(v->vst > v->nou){
				print_vstdot_bspoa(g, "1.dot");
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				check_dup_edges_bspoa(g);
				abort();
			}
		}
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			if(v->state){
				continue; // already pushed
			}
			if(v->vst > v->nou){
				print_vstdot_bspoa(g, "1.dot");
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				check_dup_edges_bspoa(g);
				abort();
			}
			if(v->vst == v->nou){
				ready = 1;
				{
					xidx = v->aligned;
					moff = v->mpos;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						if(x->nou > x->vst){
							ready = 0;
							break;
						}
						if(x->mpos > moff){
							moff = x->mpos;
						}
						xidx = x->aligned;
					}
				}
				if(ready){
					v->mpos  = moff;
					v->state = 1;
					push_u4v(g->stack, e->node);
					xidx = v->aligned;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						x->mpos = moff;
						if(x->edge){
							push_u4v(g->stack, xidx);
							x->state = 1;
						}
						xidx = x->aligned;
					}
				}
			}
		}
	}
	if(nidx != BSPOA_HEAD_NODE){
		fprint_dot_bspoa(g, 0, MAX_U4, 1, "1.dot", NULL);
		print_seqs_bspoa(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	u = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
	clear_u4v(g->stack);
	eidx = u->erev;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		if(e->node == BSPOA_HEAD_NODE) continue;
		x = u;
		v = ref_bspoanodev(g->nodes, e->node);
		while(1){
			nou = 0;
			xidx = v->edge;
			while(xidx){
				if(g->edges->buffer[xidx].node != offset_bspoanodev(g->nodes, x) && g->edges->buffer[xidx].node != BSPOA_TAIL_NODE) nou ++;
				xidx = g->edges->buffer[xidx].next;
			}
			if(nou) break;
			if(v->nin != 1) break;
			x = v;
			v = ref_bspoanodev(g->nodes, g->edges->buffer[v->erev].node);
		}
		if(x == u) continue;
		moff = v->mpos - 1;
		v = x;
		if(v->mpos == moff) continue;
		while(v != u){
			xidx = v->aligned;
#if DEBUG
			if(v->mpos > moff){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			do {
				x = ref_bspoanodev(g->nodes, xidx);
				x->mpos = moff;
				xidx = x->aligned;
			} while(x != v);
			moff --;
			xidx = v->edge;
			v = NULL;
			while(xidx){
				if(g->edges->buffer[xidx].node != BSPOA_TAIL_NODE){
					if(v){
						fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						fprint_local_dot_bspoa(g, g->edges->buffer[xidx].node, 50, "1.dot", NULL);
						abort();
					} else {
						v = ref_bspoanodev(g->nodes, g->edges->buffer[xidx].node);
					}
				}
				xidx = g->edges->buffer[xidx].next;
			}
			if(v == NULL) break;
		}
	}
#if DEBUG && 0
	for(xidx=0;xidx<nseq;xidx++){
		x = ref_bspoanodev(g->nodes, g->ndoffs->buffer[xidx]);
		v = x + g->seqs->rdlens->buffer[xidx] - 1;
		while(x + 1 < v){
			if(x->mpos <= (x+1)->mpos){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			x ++;
		}
	}
#endif
	v = ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE);
	clear_u4v(g->msaidxs);
	clear_u4v(g->msacycs);
	for(nidx=0;nidx<v->mpos;nidx++){
		push_u4v(g->msaidxs, nidx);
	}
	clear_and_encap_u1v(g->msacols, g->msaidxs->size * mrow);
	g->msacols->size = g->msaidxs->size * mrow;
	memset(g->msacols->buffer, 4, g->msacols->size);
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst = 0;
		u->mpos = g->msaidxs->size - 1 - u->mpos;
	}
	return g->msaidxs->size;
}

// TODO: find sibling, then backtrace to make sure no cycle will be generated, merge them. I am not sure which direction, or both
static inline u4i merge_siblings_bspoa(BSPOA *g){
	bspoanode_t *u, *v, *w, *x;
	bspoaedge_t *e;
	u4i nidx, eidx, hsz, ssz, tidx, i, j, flag, maxpos, minpos, ret;
	ret = 0;
	// assign mpos for all nodes, to backtrace on DAG
	sort_nodes_bspoa(g);
	for(nidx=0;nidx<g->nodes->size;nidx++){
		v = ref_bspoanodev(g->nodes, nidx);
		v->vst = 0;
		v->state = 0;
	}
	clear_u4v(g->stack);
	push_u4v(g->stack, BSPOA_HEAD_NODE);
	while(g->stack->size){
		nidx = g->stack->buffer[-- g->stack->size];
		u = ref_bspoanodev(g->nodes, nidx);
		if(u->state) continue;
		u->state = 1;
		if(u->vst != u->nin){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		hsz = g->stack->size;
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			push_u4v(g->stack, e->node);
		}
		sort_array(g->stack->buffer + hsz, g->stack->size - hsz, u4i, num_cmpgtx(g->nodes->buffer[a].base, g->nodes->buffer[b].base, a, b)); // base and node_idx
		v = ref_bspoanodev(g->nodes, g->stack->buffer[hsz]);
		for(i=hsz+1;i<g->stack->size;i++){
			w = ref_bspoanodev(g->nodes, g->stack->buffer[i]);
			if(w->base == v->base){
				flag = 0;
				ssz = g->stack->size;
				for(j=0;flag==0&&j<2;j++){
					if(j & 0x1){
						push_u4v(g->stack, offset_bspoanodev(g->nodes, w));
						tidx = offset_bspoanodev(g->nodes, v);
					} else {
						push_u4v(g->stack, offset_bspoanodev(g->nodes, v));
						tidx = offset_bspoanodev(g->nodes, w);
					}
					while(g->stack->size > ssz){
						g->stack->size --;
						if(g->stack->buffer[g->stack->size] == tidx){
							flag = 1;
							break;
						}
						x = ref_bspoanodev(g->nodes, g->stack->buffer[g->stack->size]);
						if(x->vst == x->nin) continue;
						eidx = x->erev;
						while(eidx){
							e = ref_bspoaedgev(g->edges, eidx);
							eidx = e->next;
							push_u4v(g->stack, e->node);
						}
					}
				}
				g->stack->size = ssz;
				if(flag == 0){ // v and w not connected by nodes in one direction
					minpos = num_min(v->mpos, w->mpos);
					maxpos = num_max(v->mpos, w->mpos);
					fflush(stdout); fprintf(stderr, " -- merge %c%u %c%u from %c%u in %s -- %s:%d --\n", "ACGTN"[w->base], w->mpos, "ACGTN"[v->base], v->mpos, "ACGTN"[u->base], u->mpos, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					fprint_dot_bspoa(g, u->mpos > 2? u->mpos - 2 : 0, maxpos + 2, 1, "1.dot", NULL);
					w = merge_node_bspoa(g, 0, w, v);
					w->mpos = minpos;
					w->vst = 0;
					eidx = w->erev;
					while(eidx){
						e = ref_bspoaedgev(g->edges, eidx);
						eidx = e->next;
						x = ref_bspoanodev(g->nodes, e->node);
						if(x->state && x != u){
							w->vst ++;
						}
					}
					fprint_dot_bspoa(g, u->mpos > 2? u->mpos - 2 : 0, maxpos + 2, 1, "2.dot", NULL);
					//if(u->mpos == 251 && w->mpos == 256 && v->mpos == 259){
						//return ret;
					//}
					ret ++;
				}
			}
			v = w;
		}
		g->stack->size = hsz;
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			v->vst ++;
			if(v->vst == v->nin){
				push_u4v(g->stack, e->node);
			} else if(v->vst > v->nin){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				print_node_edges_bspoa(g, e->node, 1, stderr);
				abort();
			}
		}
	}
	return ret;
}

static inline u4i msa_bspoa(BSPOA *g){
	bspoanode_t *v, *u, *x;
	bspoaedge_t *e;
	u1i *qs;
	u4i nseq, mrow, mlen, nidx, eidx, xidx, ready, rid, pos;
	nseq = g->nrds;
	mrow = nseq + 3;
	sort_nodes_bspoa(g);
	mlen = g->msaidxs->size;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst = 0;
	}
	clear_u4v(g->stack);
	push_u4v(g->stack, BSPOA_HEAD_NODE);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			v->vst ++;
			if(v->vst == v->nin){
				ready = 1;
				xidx = v->aligned;
				while(xidx != e->node){
					x = ref_bspoanodev(g->nodes, xidx);
					if(x->nin > x->vst){
						ready = 0;
						break;
					}
					xidx = x->aligned;
				}
				if(ready){
					xidx = e->node;
					do {
						x = ref_bspoanodev(g->nodes, xidx);
						if(xidx != BSPOA_TAIL_NODE && x->aux == 0){
							set_u1v(g->msacols, x->mpos * mrow + x->rid, (x->base) & 0x03);
						}
						if(x->erev){
							push_u4v(g->stack, xidx);
						}
						xidx = x->aligned;
					} while(xidx != e->node);
				}
			} else if(v->vst > v->nin){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	for(rid=0;rid<nseq;rid++){
		for(pos=0;pos<mlen;pos++){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(qs[rid] < 4){
				break;
			} else if(qs[rid] == 4){
				qs[rid] = 5;
			}
		}
	}
	for(rid=0;rid<nseq;rid++){
		for(pos=mlen-1;pos>0;pos--){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(qs[rid] < 4){
				break;
			} else if(qs[rid] == 4){
				qs[rid] = 5;
			}
		}
	}
	if(nidx != BSPOA_TAIL_NODE){
		fprint_dot_bspoa(g, 0, MAX_U4, 1, "1.dot", NULL);
		print_seqs_bspoa(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		print_node_edges_bspoa(g, nidx, 1, stderr);
		print_node_edges_bspoa(g, nidx, 1, stderr);
		abort();
	}
	return g->msaidxs->size;
}

static inline u4i mask_free_ends_msa_bspoa(BSPOA *g){
	bspoanode_t *v, *u;
	bspoaedge_t *e;
	u4i dir, nseq, mrow, mlen, i, nidx, eidx, ridx, rlen, max, midx, ret;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	ret = 0;
	for(dir=0;dir<2;dir++){
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			u->vst = 0;
		}
		max = 0; midx = dir? BSPOA_TAIL_NODE : BSPOA_HEAD_NODE;
		clear_u4v(g->stack);
		push_u4v(g->stack, midx);
		while(pop_u4v(g->stack, &nidx)){
			u = ref_bspoanodev(g->nodes, nidx);
			if(dir){
				eidx = u->erev;
			} else {
				eidx = u->edge;
			}
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				v = ref_bspoanodev(g->nodes, e->node);
				if(v->vst) continue;
				if(dir){
					if(v->nou > 1) continue;
				} else {
					if(v->nin > 1) continue;
				}
				if(v->aligned != e->node) continue;
				v->vst = u->vst + 1;
				if(v->vst > max){
					max = v->vst;
					midx = e->node;
				}
				push_u4v(g->stack, e->node);
			}
		}
		if(max == 0) continue;
		midx = ref_bspoanodev(g->nodes, midx)->rid;
		for(ridx=0;ridx<nseq;ridx++){
			if(ridx == midx) continue;
			rlen = g->seqs->rdlens->buffer[ridx];
			for(i=0;i<rlen;i++){
				u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[ridx] + (dir? rlen - 1 - i : i));
				if(u->vst == 0) break;
				ret ++;
				set_u1v(g->msacols, g->msaidxs->buffer[u->mpos] * mrow + u->rid, 6); // 6: masked
			}
		}
		if(_DEBUG_LOG_){
			fflush(stdout); fprintf(stderr, " -- NSEQ[%d] dir=%d max=%d ridx=%d mask=%d in %s -- %s:%d --\n", nseq, dir, max, midx, ret, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
	}
	return ret;
}

static inline void fix_tenon_mortise_msa_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i rid, nseq, mrow, mlen, i, pos, *rps;
	u2i pairs[2][4], cut;
	u1i *cols[2], pmaxs[2];
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	if(mlen == 0) return;
	clear_and_encap_b1v(g->memp, nseq * sizeof(u4i));
	rps = (u4i*)g->memp->buffer;
	for(rid=0;rid<nseq;rid++){
		rps[rid] = g->seqs->rdlens->buffer[rid];
	}
	cut = num_max(2, nseq / 3);
	cols[0] = cols[1] = NULL;
	for(pos=mlen-1;pos!=MAX_U4;pos--){
		cols[1] = cols[0];
		cols[0] = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(rid=0;rid<nseq;rid++){
			if(cols[0][rid] < 4){
				rps[rid] --;
			}
		}
		if(cols[1] == NULL) continue;
		memset(pairs[0], 0, 4 * sizeof(u2i));
		memset(pairs[1], 0, 4 * sizeof(u2i));
		for(rid=0;rid<nseq;rid++){
			if(cols[0][rid] == cols[1][rid]) continue;
			if(cols[1][rid] == 4) pairs[0][cols[0][rid]] ++;
			else if(cols[0][rid] == 4) pairs[1][cols[1][rid]] ++;
		}
		pmaxs[0] = pmaxs[1] = 0;
		for(i=1;i<4;i++){
			if(pairs[0][i] > pairs[0][pmaxs[0]]) pmaxs[0] = i;
			if(pairs[1][i] > pairs[1][pmaxs[1]]) pmaxs[1] = i;
		}
		if(pairs[0][pmaxs[0]] < cut || pairs[1][pmaxs[1]] < cut) continue;
		if(pmaxs[0] == pmaxs[1]) continue; // TODO: hp_adjust
		if(_DEBUG_LOG_ > 2){
			fflush(stdout);
			fprintf(stderr, "[tenon_mortise:%d/%d] %c\t%c\t%d\t%d\n", pos, mlen, "ACGT-"[pmaxs[0]], "ACGT-"[pmaxs[1]], pairs[0][pmaxs[0]], pairs[1][pmaxs[1]]);
		}
		for(rid=0;rid<nseq;rid++){
			if(cols[0][rid] == pmaxs[0] && cols[1][rid] == 4){
				cols[1][rid] = pmaxs[0];
				cols[0][rid] = 4;
				v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + rps[rid]);
				v->mpos ++;
			}
		}
	}
}

// quick and dirty, major bases
static inline void simple_cns_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i nseq, mrow, mlen, pos, cpos, rid, i, b;
	u2i bcnts[7], brank[7];
	u1i *qs;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	if(mlen == 0) return;
	//mask_free_ends_msa_bspoa(g);
	clear_u1v(g->cns);
	clear_u1v(g->qlt);
	clear_u1v(g->alt);
	for(rid=0;rid<nseq;rid++){
		for(pos=0;pos<mlen;pos++){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(qs[rid] < 4){
				break;
			} else if(qs[rid] == 4){
				qs[rid] = 5;
			}
		}
	}
	for(rid=0;rid<nseq;rid++){
		for(pos=mlen-1;pos>0;pos--){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(qs[rid] < 4){
				break;
			} else if(qs[rid] == 4){
				qs[rid] = 5;
			}
		}
	}
	for(pos=0;pos<mlen;pos++){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		memset(bcnts, 0, 7 * sizeof(u2i));
		memset(brank, 0xff, 7 * sizeof(u2i));
		for(rid=0;rid<nseq;rid++){
			bcnts[qs[rid]] ++;
			if(brank[qs[rid]] == 0xFFFFU){
				brank[qs[rid]] = rid;
			}
		}
		b = 4;
		for(i=0;i<4;i++){
			if(bcnts[i] > bcnts[b]){
				b = i;
			} else if(bcnts[i] == bcnts[b]){
				if(brank[i] < brank[b]){
					b = i;
				}
			}
		}
		qs[nseq] = b;
		qs[nseq + 1] = 0;
		if(b < 4){
			push_u1v(g->cns, b);
			push_u1v(g->qlt, 0);
			push_u1v(g->alt, 0);
		}
	}
	for(rid=0;rid<nseq;rid++){
		cpos = 0;
		v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid]);
		for(pos=0;pos<mlen;pos++){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(qs[rid] != 4 && qs[rid] != 5){
				v->cpos = cpos;
				v ++;
			}
			if(qs[nseq] < 4) cpos ++;
		}
	}
	ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE)->cpos = 0;
	ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE)->cpos = g->cns->size;
}

#define MAX_LOG_CACHE	1000
static u4i _log_caches_n = 1;
static long double _log_caches[MAX_LOG_CACHE + 1];

static inline long double cal_permutation_bspoa(u4i n, u4i m){
	if(n > MAX_LOG_CACHE) return 1;
	_log_caches[0] = 0;
	while(_log_caches_n <= n){
		_log_caches[_log_caches_n] = _log_caches[_log_caches_n - 1] + logl(_log_caches_n);
		_log_caches_n ++;
	}
	return _log_caches[n] - _log_caches[m] - _log_caches[n - m];
}

static inline long double cal_binomial_bspoa(u4i n, u4i m, long double p){
	return logl(p) * m + logl(1 - p) * (n - m) + cal_permutation_bspoa(n, m);
}

static inline long double sum_log_nums(int cnt, long double *vals){
	long double delta, sum;
	int i;
	if(cnt <= 0) return 0;
	sort_array(vals, cnt, long double, num_cmpgt(b, a));
	sum = vals[0];
	for(i=1;i<cnt;i++){
		delta = vals[i] - sum;
		if(delta <= -40) break; // overflow
		delta = expl(delta);
		delta = logl(1 + delta);
		sum += delta;
	}
	return sum;
}

static inline long double cns_bspoa(BSPOA *g){
	typedef struct {long double sc[6]; u1i bt, lb;} dp_t;
	bspoanode_t *v;
	bspoavar_t *var;
	dp_t *dps[5], *dp[5], *lp;
	long double ret, errs[5], errd, erre, p, log10;
	u1i a, b, c, d, e, f, *r, *qs, *ts, *bs[10], *col;
	u4i nseq, mrow, mlen, mpos, cpos, rid, cnt[3], i, mpsize, cnts[7];
	int pos;
	nseq = g->nrds;
	mrow = nseq + 3;
	log10 = logl(10);
	mlen = g->msaidxs->size;
	mpsize = (mlen + 1) * 5 * sizeof(dp_t);
	mpsize += nseq * 10; // bs[0/1/2/3/4] is the state of the last read non-N base vs cns[A/C/G/T/-], M/I/D = 0/1/2; 5-9 is a backup
	clear_and_encap_b1v(g->memp, mpsize);
	r = (u1i*)(g->memp->buffer);
	for(i=0;i<5;i++){
		dps[i] = (dp_t*)r;
		r += (mlen + 1) * sizeof(dp_t);
		memset(dps[i][mlen].sc, 0, 6 * sizeof(long double));
		dps[i][mlen].bt = 4;
		dps[i][mlen].lb = 4;
	}
	for(i=0;i<10;i++){
		bs[i] = (u1i*)r; r += nseq;
		memset(bs[i], 0, nseq);
	}
	for(pos=mlen-1;pos>=0;pos--){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(a=0;a<=4;a++){
			dp[a] = dps[a] + pos;
			memset(dp[a]->sc, 0, 6 * sizeof(long double));
		}
		for(e=0;e<=4;e++){
			lp = dps[e] + pos + 1;
			for(rid=0;rid<nseq;rid++){
				b = qs[rid];
				if(b > 4) continue;
				c = lp->lb;
				d = bs[e][rid];
				ts = g->dptable + b * 5 + c * 25 + d * 125;
				for(a=0;a<=4;a++){
					dp[a]->sc[e] += g->dpvals[ts[a] >> 3];
				}
			}
		}
		for(a=0;a<=4;a++){
			errs[0] = dp[a]->sc[0] + dps[0][pos+1].sc[5];
			errs[1] = dp[a]->sc[1] + dps[1][pos+1].sc[5];
			errs[2] = dp[a]->sc[2] + dps[2][pos+1].sc[5];
			errs[3] = dp[a]->sc[3] + dps[3][pos+1].sc[5];
			errs[4] = dp[a]->sc[4] + dps[4][pos+1].sc[5];
			dp[a]->sc[5] = sum_log_nums(5, errs);
			dp[a]->bt = 4;
			for(e=0;e<4;e++){
				if(dp[a]->sc[e] + dps[e][pos+1].sc[5] > dp[a]->sc[dp[a]->bt] + dps[dp[a]->bt][pos+1].sc[5]) dp[a]->bt = e;
			}
			lp = dps[dp[a]->bt] + pos + 1;
			dps[a][pos].lb = (a < 4)? a : lp->lb;
			for(rid=0;rid<nseq;rid++){
				b = qs[rid];
				if(b > 4){
					bs[a + 5][rid] = 4;
					continue;
				}
				f = g->dptable[a + b * 5 + lp->lb * 25 + bs[dp[a]->bt][rid] * 125];
				bs[a + 5][rid] = f & 0x7;
			}
		}
		for(a=0;a<5;a++){
			memcpy(bs[a], bs[a + 5], nseq);
		}
	}
	c = 0;
	for(a=1;a<=4;a++){
		if(dps[a][0].sc[5] > dps[c][0].sc[5]){
			c = a;
		}
	}
	ret = dps[c][0].sc[5];
	clear_u1v(g->cns);
	clear_u1v(g->qlt);
	clear_u1v(g->alt);
	for(i=0;i<mlen;i++){
		r = g->msacols->buffer + g->msaidxs->buffer[i] * mrow + nseq;
		r[0] = c;
		c = dps[c][i].bt;
	}
	for(i=0;i<10;i++){
		memset(bs[i], 0, nseq);
	}
	for(pos=mlen-1;pos>=0;pos--){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		c = qs[nseq];
		errs[0] = dps[0][pos].sc[5];
		errs[1] = dps[1][pos].sc[5];
		errs[2] = dps[2][pos].sc[5];
		errs[3] = dps[3][pos].sc[5];
		errs[4] = dps[4][pos].sc[5];
		erre = sum_log_nums(5, errs);
		errd = dps[c][pos].sc[5];
		erre = - (10 * log(1 - exp(errd - erre)) / log10);
		qs[nseq + 1] = (int)num_min(erre, 40.0);
		if(1){
			a = (c + 1) % 5;
			for(e=0;e<=4;e++){
				if(e == c) continue;
				if(dps[e][pos].sc[5] > dps[a][pos].sc[5]) a = e;
			}
			// assumes number of alt base equals cns base
			cnt[0] = cnt[1] = 0;
			for(rid=0;rid<nseq;rid++){
				if(qs[rid] == a) cnt[0] ++;
				else if(qs[rid] == c) cnt[1] ++;
			}
			e = 4;
			for(i=pos+1;i<mlen;i++){
				e = g->msacols->buffer[g->msaidxs->buffer[i] * mrow + nseq];
				if(e < 4) break;
			}
			d = 4;
			for(i=pos;i>0;i--){
				d = g->msacols->buffer[g->msaidxs->buffer[i - 1] * mrow + nseq];
				if(d < 4) break;
			}
			p = num_max(g->dporis[g->dptable[c + a * 5 + e * 25 + e * 125] >> 3], g->dporis[g->dptable[c + a * 5 + d * 25 + d * 125] >> 3]);
			erre = 0;
			for(cnt[2]=0;cnt[2]<cnt[0];cnt[2]++){
				erre += expl(cal_binomial_bspoa(cnt[0] + cnt[1], cnt[2], p));
			}
			if(erre == 0){
				errd = 0;
			} else {
				errd = - (10 * logl(1 - erre) / log10);
			}
			if(_DEBUG_LOG_ > 2){
				char flanks[2][3];
				for(f=2,i=pos;i>0&&f>0;i--){
					d = g->msacols->buffer[g->msaidxs->buffer[i - 1] * mrow + nseq];
					if(d < 4){
						flanks[0][f - 1] = "ACGT-"[d];
						f --;
					}
				}
				for(f=0,i=pos+1;i<mlen&&f<2;i++){
					d = g->msacols->buffer[g->msaidxs->buffer[i] * mrow + nseq];
					if(d < 4) flanks[1][f++] = "ACGT-"[d];
				}
				flanks[0][2] = flanks[1][2] = 0;
				fprintf(stderr, "[ALTQ] %d/%d\t%s\t%c\t%c\t%s\t%d\t%d\t%Lf\t%Lf\t%Lf\n", pos, mlen, flanks[0], "ACGT- "[c], "ACGT-"[a], flanks[1], cnt[1], cnt[0], p, erre, errd);
			}
			qs[nseq + 2] = num_min(errd, 40.0);
		} else {
			qs[nseq + 2] = 0;
		}
		if(qs[nseq + 0] < 4){
			push_u1v(g->cns, qs[nseq + 0]);
			push_u1v(g->qlt, qs[nseq + 1]);
			push_u1v(g->alt, qs[nseq + 2]);
		}
	}
	reverse_u1v(g->cns);
	reverse_u1v(g->qlt);
	reverse_u1v(g->alt);
	if(_DEBUG_LOG_ > 2){
		for(pos=0;pos<Int(mlen);pos++){
			long double minp = dps[0][pos].sc[5];
			for(a=1;a<=4;a++){
				if(dps[a][pos].sc[5] < minp) minp = dps[a][pos].sc[5];
			}
			fprintf(stderr, "[CNSQ] %d\t%c\t%0.2Lf", pos, "ACGT- "[g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + nseq]], minp);
			for(a=0;a<=4;a++){
				fprintf(stderr, "\t%c:%0.2Lf[", "ACGT-"[a], dps[a][pos].sc[5] - minp);
				for(e=0;e<=4;e++){
					if(e == dps[a][pos].bt) fprintf(stderr, "*");
					fprintf(stderr, "%c:%0.2Lf", "acgt-"[e], dps[a][pos].sc[e] + dps[e][pos+1].sc[5] - minp);
					if(e < 4) fprintf(stderr, ",");
				}
				fprintf(stderr, "]");
			}
			fprintf(stderr, "\n");
		}
	}
	// set node_t->cpos to the position on cns, useful in bandwidth
	for(rid=0;rid<nseq;rid++){
		cpos = 0;
		//v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + g->seqs->rdlens->buffer[rid] - 1);
		v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid]);
		for(pos=0;pos<Int(mlen);pos++){
			if(g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + rid] < 4){
				v->cpos = cpos;
				v ++;
			}
			if(g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + nseq] < 4) cpos ++;
		}
	}
	ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE)->cpos = 0;
	ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE)->cpos = g->cns->size;
	clear_bspoavarv(g->var);
	for(mpos=cpos=0;mpos<mlen;mpos++){
		col = g->msacols->buffer + g->msaidxs->buffer[mpos] * mrow;
		if(col[nseq] >= 4) continue;
		while(g->alt->buffer[cpos] >= g->par->qltlo){
			memset(cnts, 0, 7 * sizeof(u4i));
			for(i=0;i<nseq;i++){
				cnts[col[i]] ++;
			}
			a = g->cns->buffer[cpos];
			b = (a + 1) % 5;
			for(i=0;i<=4;i++){
				if(i == a) continue;
				if(cnts[i] > cnts[b]) b = i;
			}
			if(b == 4) break; // Only supports SNP
			var = next_ref_bspoavarv(g->var);
			var->cpos = cpos;
			var->mpos = mpos;
			var->refn = cnts[a];
			var->refb = a;
			var->altn = cnts[b];
			var->altb = b;
			break;
		}
		cpos ++;
	}
	return ret;
}

static inline void hp_adjust_bspoa(BSPOA *g){
	bspoanode_t *u, *v;
	u4i rid, nseq, mrow, mlen, i, p, pos, uidx, *rps;
	u1i *col;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	if(mlen == 0) return;
	clear_and_encap_b1v(g->memp, nseq * sizeof(u4i));
	rps = (u4i*)g->memp->buffer;
	for(rid=0;rid<nseq;rid++){
		rps[rid] = g->seqs->rdlens->buffer[rid];
	}
	for(pos=mlen-1;pos!=MAX_U4;pos--){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		uidx = MAX_U4;
		for(rid=0;rid<nseq;rid++){
			if(col[rid] < 4){
				rps[rid] --;
				if(uidx == MAX_U4 && col[rid] == col[nseq]){
					uidx = g->ndoffs->buffer[rid] + rps[rid];
				}
			}
		}
		if(col[nseq] >= 4) continue;
		if(uidx == MAX_U4){ // all MSA base in this col were moved away
			continue;
		}
		for(rid=0;rid<nseq;rid++){
			if(col[rid] == 4 && rps[rid] && get_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[rid] + rps[rid] - 1) == col[nseq]){
				rps[rid] --;
				for(i=pos-1;;i--){
					p = g->msaidxs->buffer[i] * mrow + rid;
					if(g->msacols->buffer[p] >= 4) continue;
#if DEBUG
					if(g->msacols->buffer[p] != col[nseq]){
						fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					g->msacols->buffer[p] = rps[rid]? 4 : 5;
					col[rid] = col[nseq];
					u = ref_bspoanodev(g->nodes, uidx);
					v = seperate_node_bspoa(g, rid, rps[rid]);
					merge_node_bspoa(g, rid, u, v);
					break;
				}
			}
		}
	}
}

static inline void end_bspoa(BSPOA *g){
	u4i nidx, tigger;
	int score;
	if(g->seqs->nseq == 0){
		clear_u1v(g->cns);
		clear_u1v(g->qlt);
		clear_u1v(g->alt);
		return;
	}
	clear_u4v(g->ndoffs);
	tigger = g->par->bwtrigger;
	nidx = g->nodes->size;
	for(g->nrds=0;g->nrds<g->seqs->nseq;g->nrds++){
		//if(g->par->refmode == 0 && tigger > 0 && ((g->nrds + 1) % tigger) == 0){
		if(g->par->refmode == 0 && (g->bwtrigger || (tigger > 0 && g->nrds == tigger))){
			// generate MSA+CNS per bwtrigger reads, the cpos of node is needed to perform banded alignment
			g->bwtrigger = 1;
			msa_bspoa(g);
			simple_cns_bspoa(g);
			if(_DEBUG_LOG_ > 1){
				print_msa_sline_bspoa(g, stderr);
			}
		}
		push_u4v(g->ndoffs, nidx);
		nidx += g->seqs->rdlens->buffer[g->nrds];
		score = align_rd_bspoa(g, g->nrds);
	}
	push_u4v(g->ndoffs, nidx);
	// generate MSA and CNS
	msa_bspoa(g);
	cns_bspoa(g);
	hp_adjust_bspoa(g);
	fix_tenon_mortise_msa_bspoa(g);
	cns_bspoa(g);
}

static inline void check_msa_bspoa(BSPOA *g){
	u4i nseq, mrow, mlen, rid, pos, len;
	u1i c;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	for(rid=0;rid<nseq;rid++){
		len = 0;
		for(pos=0;pos<mlen;pos++){
			c = g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + rid];
			if(c < 4){
				if(c != get_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[rid] + len)){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				len ++;
			}
		}
		if(len != g->seqs->rdlens->buffer[rid]){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

static inline void check_graph_cov_bspoa(BSPOA *g){
	bspoanode_t *v, *x;
	bspoaedge_t *e;
	u4i node, eidx, xidx, ncov, ecov;
	for(node=2;node<g->nodes->size;node++){
		v = ref_bspoanodev(g->nodes, node);
		if(v->header != node) continue;
		ncov = 1;
		xidx = v->aligned;
		while(xidx != node){
			x = ref_bspoanodev(g->nodes, xidx);
			ncov ++;
			xidx = x->aligned;
		}
		ecov = 0;
		eidx = v->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			ecov += e->cov;
		}
		if(ecov != ncov){
			print_aligned_nodes_bspoa(g, node, stderr);
			print_node_edges_bspoa(g, node, 0, stderr);
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		ecov = 0;
		eidx = v->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			ecov += e->cov;
		}
		if(ecov != ncov){
			print_aligned_nodes_bspoa(g, node, stderr);
			print_node_edges_bspoa(g, node, 1, stderr);
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

// return updated msaend on whole MSA
static inline u4i local_remsa_bspoa(BSPOA *g, u4i msabeg, u4i msaend, u4i *_rbs, u4i *_res, u4i *_ridxs, float *_scores, u1i *_ld, BSPOA *lg){
	float *scores;
	u4i nseq, mrow, mlen, i, rid, pos;
	u4i *rbs, *res, *ridxs;
	u1i *qs, q, f, lc, *ld;
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	msaend = num_min(msaend, mlen);
	if(_rbs) rbs = _rbs;
	else {
		rbs = (u4i*)alloca(nseq * sizeof(u4i));
		memset(rbs, 0, nseq * sizeof(u4i));
		for(pos=0;pos<msabeg;pos++){
			qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			for(rid=0;rid<nseq;rid++){
				if(qs[rid] != 4 && qs[rid] != 5) rbs[rid] ++;
			}
		}
	}
	if(_res) res = _res;
	else res = (u4i*)alloca(nseq * sizeof(u4i));
	if(_ridxs) ridxs = _ridxs;
	else ridxs = (u4i*)alloca(nseq * sizeof(u4i));
	for(rid=0;rid<nseq;rid++) ridxs[rid] = rid;
	if(_scores) scores = _scores;
	else scores = (float*)alloca(nseq * sizeof(float));
	memset(scores, 0, nseq * sizeof(float));
	if(_ld) ld = _ld;
	else ld = (u1i*)alloca(nseq * sizeof(u1i));
	memset(ld, 0, nseq);
	lc = 4;
	memcpy(res, rbs, nseq * sizeof(u4i));
	for(pos=msabeg;pos<msaend;pos++){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(rid=0;rid<nseq;rid++){
			if(qs[rid] < 4){
				res[rid] ++;
				q = qs[rid];
			} else {
				if(qs[rid] == 6){
					res[rid] ++;
				}
				q = 4;
			}
			f = g->dptable[qs[nseq] + qs[rid] * 5 + lc * 25 + ld[rid] * 125];
			scores[rid] += g->dpvals[f >> 3] * qs[nseq + 1]; // probs * cns_phred_quality
			ld[rid] = f & 0x7;
		}
	}
	sort_array(ridxs, nseq, u4i, num_cmpgt(scores[b], scores[a]));
	beg_bspoa(lg);
	lg->par->refmode   = 0;
	lg->par->remsa     = 0;
	lg->par->nrec      = 0;
	lg->par->bwtrigger = 0;
	for(i=0;i<nseq;i++){
		rid = ridxs[i];
		if(_DEBUG_LOG_){
			fprintf(stderr, "LRMSA:RD[%d]\t%d\t%0.4f\t%d\t%d\n", i, rid, scores[rid], rbs[rid], res[rid]); fflush(stderr);
		}
		fwdbitpush_bspoa(lg, g->seqs->rdseqs->bits, g->seqs->rdoffs->buffer[rid] + rbs[rid], res[rid] - rbs[rid]);
		del_read_nodes_bspoa(g, rid, rbs[rid], res[rid]);
	}
	end_bspoa(lg);
	if(_DEBUG_LOG_){
		print_msa_sline_bspoa(lg, stderr);
	}
	// copy the local graph
	{
		bspoanode_t *v, *u, *x;
		bspoaedge_t *e;
		u4i nidx, eidx, xidx, vidxs[3];
		ref_bspoanodev(lg->nodes, BSPOA_HEAD_NODE)->vst = 0;
		ref_bspoanodev(lg->nodes, BSPOA_TAIL_NODE)->vst = 0;
		for(nidx=2;nidx<lg->nodes->size;nidx++){
			u = ref_bspoanodev(lg->nodes, nidx);
			u->vst = 0;
			if(u->header == nidx){
				rid = ridxs[u->rid];
				vidxs[0] = g->ndoffs->buffer[rid] + rbs[rid] + u->pos;
				vidxs[1] = vidxs[0];
				xidx = u->aligned;
				while(xidx != nidx){
					x = ref_bspoanodev(lg->nodes, xidx);
					rid = ridxs[x->rid];
					vidxs[2] = g->ndoffs->buffer[rid] + rbs[rid] + x->pos;
					ref_bspoanodev(g->nodes, vidxs[1])->aligned = vidxs[2];
					ref_bspoanodev(g->nodes, vidxs[1])->header  = vidxs[0];
					vidxs[1] = vidxs[2];
					xidx = x->aligned;
				}
				ref_bspoanodev(g->nodes, vidxs[1])->aligned = vidxs[0];
				ref_bspoanodev(g->nodes, vidxs[1])->header  = vidxs[0];
			}
			u->vst = 0;
		}
		clear_u4v(lg->stack);
		push_u4v(lg->stack, BSPOA_HEAD_NODE);
		while(pop_u4v(lg->stack, &nidx)){
			u = ref_bspoanodev(lg->nodes, nidx);
			rid = ridxs[u->rid];
			vidxs[0] = (nidx < 2)? 0 : g->ndoffs->buffer[rid] + rbs[rid] + u->pos;
			eidx = u->edge;
			while(eidx){
				e = ref_bspoaedgev(lg->edges, eidx);
				eidx = e->next;
				v = ref_bspoanodev(lg->nodes, e->node);
				rid = ridxs[v->rid];
				vidxs[1] = (e->node < 2)? 0 : g->ndoffs->buffer[rid] + rbs[rid] + v->pos;
				if(vidxs[0] && vidxs[1]){
					hadd_edge_bspoa(g, ref_bspoanodev(g->nodes, vidxs[0]), ref_bspoanodev(g->nodes, vidxs[1]), e->cov, e->is_aux, NULL);
				}
				v->vst ++;
				if(v->vst == v->nin){
					push_u4v(lg->stack, e->node);
#if DEBUG
				} else if(v->vst > v->nin){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
#endif
				}
			}
		}
		// join the head and tail
		u = ref_bspoanodev(lg->nodes, BSPOA_HEAD_NODE);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(lg->edges, eidx);
			eidx = e->next;
			if(e->is_aux) continue;
			//if(e->node < 2) continue; // HEAD or TAIL
			xidx = e->node;
			do {
				v = ref_bspoanodev(lg->nodes, xidx);
				xidx = v->aligned;
				if(offset_bspoanodev(lg->nodes, v) == BSPOA_TAIL_NODE) continue;
				rid = ridxs[v->rid];
				if(rbs[rid] == res[rid]) continue;
				vidxs[1] = g->ndoffs->buffer[rid] + rbs[rid] + v->pos;
				if(v->pos){ // within lg
					continue;
				} else if(rbs[rid]){
					vidxs[0] = g->ndoffs->buffer[rid] + rbs[rid] + v->pos - 1;
				} else {
					vidxs[0] = BSPOA_HEAD_NODE;
				}
				hadd_edge_bspoa(g, ref_bspoanodev(g->nodes, vidxs[0]), ref_bspoanodev(g->nodes, vidxs[1]), 1, 0, NULL);
			} while(xidx != e->node);
		}
		u = ref_bspoanodev(lg->nodes, BSPOA_TAIL_NODE);
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(lg->edges, eidx);
			eidx = e->next;
			if(e->is_aux) continue;
			//if(e->node < 2) continue; // HEAD or TAIL
			xidx = e->node;
			do {
				v = ref_bspoanodev(lg->nodes, xidx);
				xidx = v->aligned;
				if(offset_bspoanodev(lg->nodes, v) == BSPOA_HEAD_NODE) continue;
				rid = ridxs[v->rid];
				if(rbs[rid] == res[rid]) continue;
				vidxs[0] = g->ndoffs->buffer[rid] + rbs[rid] + v->pos;
				if(rbs[rid] + v->pos + 1 < res[rid]){ // within lg
					continue;
				} else if(res[rid] < g->seqs->rdlens->buffer[rid]){
					vidxs[1] = g->ndoffs->buffer[rid] + rbs[rid] + v->pos + 1;
				} else {
					vidxs[1] = BSPOA_TAIL_NODE;
				}
				hadd_edge_bspoa(g, ref_bspoanodev(g->nodes, vidxs[0]), ref_bspoanodev(g->nodes, vidxs[1]), 1, 0, NULL);
			} while(xidx != e->node);
		}
#if 0
		check_all_node_edges_bspoa(g); // TODO: debug
		check_dup_edges_bspoa(g);
		check_graph_cov_bspoa(g);
#endif
	}
	// copy the local msa
	{
		u4i msalen, colmax;
		u1i *src, *dst;
		msalen = num_min(msaend - msabeg, lg->msaidxs->size);
		for(pos=0;pos<msalen;pos++){
			src = lg->msacols->buffer + lg->msaidxs->buffer[pos] * mrow;
			dst = g->msacols->buffer + g->msaidxs->buffer[msabeg + pos] * mrow;
			for(i=0;i<nseq;i++) dst[ridxs[i]] = src[i];
			dst[nseq + 0] = src[nseq + 0];
			dst[nseq + 1] = src[nseq + 1];
			dst[nseq + 2] = src[nseq + 2];
		}
		if(msalen == lg->msaidxs->size){
			for(i=msabeg+lg->msaidxs->size;i<msaend;i++) push_u4v(g->msacycs, g->msaidxs->buffer[i]);
			memmove(g->msaidxs->buffer + msabeg + lg->msaidxs->size, g->msaidxs->buffer + msaend, (g->msaidxs->size - msaend) * sizeof(u4i));
			g->msaidxs->size -= msaend - msabeg - lg->msaidxs->size;
		} else {
			colmax = g->msacycs->size + g->msaidxs->size;
			if(g->msacycs->size < lg->msaidxs->size - (msaend - msabeg)){
				encap_u1v(g->msacols, (lg->msaidxs->size - (msaend - msabeg) - g->msacycs->size) * mrow);
				g->msacols->size += (lg->msaidxs->size - (msaend - msabeg) - g->msacycs->size) * mrow;
			}
			encap_u4v(g->msaidxs, lg->msaidxs->size - (msaend - msabeg));
			memmove(g->msaidxs->buffer + msabeg + lg->msaidxs->size, g->msaidxs->buffer + msaend, (g->msaidxs->size - msaend) * sizeof(u4i));
			g->msaidxs->size += lg->msaidxs->size - (msaend - msabeg);
			while(msalen < lg->msaidxs->size){
				if(g->msacycs->size){
					g->msaidxs->buffer[msabeg + msalen] = g->msacycs->buffer[--g->msacycs->size];
				} else {
					g->msaidxs->buffer[msabeg + msalen] = colmax ++;
				}
				src = lg->msacols->buffer + lg->msaidxs->buffer[msalen] * mrow;
				dst = g->msacols->buffer + g->msaidxs->buffer[msabeg + msalen] * mrow;
				for(i=0;i<nseq;i++) dst[ridxs[i]] = src[i];
				dst[nseq + 0] = src[nseq + 0];
				dst[nseq + 1] = src[nseq + 1];
				dst[nseq + 2] = src[nseq + 2];
				msalen ++;
			}
		}
#if 0
		check_msa_bspoa(g);
#endif
	}
	return msabeg + lg->msaidxs->size;
}

static inline u4i cal_complexity_col_bspoa(u1i *col, u2i nseq){
	u4i x;
	u2i i, j,  bs[7];
	memset(bs, 0, 7 * sizeof(u2i));
	for(i=0;i<nseq;i++){
		bs[col[i]] ++;
	}
	if(bs[5] + bs[6] == nseq){
		return MAX_U4;
	}
	x = 0;
	for(i=0;i<5;i++){
		for(j=i+1;j<5;j++){
			x += bs[i] * bs[j];
		}
	}
	return x;
}

static inline u4i remsa_bspoa(BSPOA *g, BSPOA *lg){
	u8i mpsize;
	float *scores;
	u4i ret, nseq, mrow, i, ac, bc, cc, wsz, msabeg, msaend, cnscnt, rb, *rbs, *res, *ridxs;
	u1i *col, *sts, *ld, *mem;
	int moffset;
	lg->par->bwtrigger = 0; // region of remsa is usually small
	wsz = g->par->rma_win;
	if(g->msaidxs->size < wsz) return 0;
	ret = 0;
	nseq = g->nrds;
	mrow = nseq + 3;
	mpsize  = roundup_times(g->msaidxs->size * sizeof(u1i), WORDSIZE);
	mpsize += roundup_times(3 * nseq * sizeof(u4i), WORDSIZE);
	mpsize += roundup_times(nseq * sizeof(float), WORDSIZE);
	mpsize += roundup_times(nseq, WORDSIZE);
	clear_and_encap_b1v(g->memp, mpsize);
	mem = (u1i*)g->memp->buffer;
	sts = (u1i*)mem; mem += roundup_times(g->msaidxs->size * sizeof(u1i), WORDSIZE);
	rbs = (u4i*)mem; mem += roundup_times(3 * nseq * sizeof(u4i), WORDSIZE);
	res = rbs + nseq;
	ridxs = res + nseq;
	scores = (float*)mem; mem += roundup_times(nseq * sizeof(float), WORDSIZE);
	ld = (u1i*)mem;
	memset(rbs, 0, nseq * sizeof(u4i));
	memset(res, 0, nseq * sizeof(u4i));
	for(i=0;i<g->msaidxs->size;i++){
		col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
		cc  = col[nseq + 1]; // quality score
		if(Int(cc) <= g->par->qltlo) sts[i] = 2;
		else if(Int(cc) <= g->par->qlthi || col[nseq] == 4) sts[i] = 1;
		else sts[i] = 0;
	}
#if 0
	check_graph_cov_bspoa(g);
#endif
	rb = 0;
	msabeg = msaend = 0;
	moffset = 0;
	while(msaend < g->msaidxs->size){
		cnscnt = 0;
		ac = msabeg;
		while(ac < g->msaidxs->size && sts[ac+moffset] != 2) ac ++;
		bc = ac;
		while(bc < g->msaidxs->size && sts[bc+moffset] == 2) bc ++;
		msaend = bc;
		cnscnt = 0;
		while(ac > msabeg && cnscnt < 2 * wsz){
			if(sts[ac+moffset] == 0) cnscnt ++;
			ac --;
		}
		if(cnscnt < wsz){
			msabeg = msaend;
			continue;
		}
		cc = 0;
		while(1){
			cnscnt = 0;
			while(bc < g->msaidxs->size && cnscnt < 2 * wsz){
				if(sts[bc+moffset] == 0) cnscnt ++;
				else if(sts[bc+moffset] == 2) break;
				bc ++;
			}
			cc += cnscnt;
			if(bc >= g->msaidxs->size || cnscnt >= wsz || cc >= wsz * 8) break;
			bc ++;
		}
		if(cnscnt < wsz){
			msabeg = msaend;
			continue;
		}
		msabeg = ac;
		cc = bc;
		while(rb < msabeg){
			col = g->msacols->buffer + g->msaidxs->buffer[rb] * mrow;
			for(i=0;i<nseq;i++){
				if(col[i] < 4) rbs[i] ++;
			}
			rb ++;
		}
		if(_DEBUG_LOG_){
			fflush(stdout); fprintf(stderr, "LRMSA[%d]: %d(%d) - %d(%d)\n", ret, msabeg + 1, msabeg + moffset + 1, cc + 1, cc + moffset + 1);
		}
		bc = local_remsa_bspoa(g, msabeg, cc, rbs, res, ridxs, scores, ld, lg);
		moffset = moffset + cc - bc;
		msabeg = bc;
		//memcpy(rbs, res, nseq * sizeof(u4i));
		if(_DEBUG_LOG_){
			fflush(stderr);
			print_msa_sline_bspoa(g, stdout);
			fflush(stdout);
		}
		ret ++;
	}
	msa_bspoa(g);
	hp_adjust_bspoa(g);
	fix_tenon_mortise_msa_bspoa(g);
	cns_bspoa(g);
	return ret;
}

#endif
