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

#define BSPOA_VST_MAX	MAX_U2

typedef struct {
	u2i rid;
	u1i base:7, ref:1, aux:2, colorful:1, inuse:1, bonus:1, bless:1, rdc:1, rdd:1; // rdc: connect previous node on read, rdd: next node
	u2i vst;
	u2i nin, nou;
	u2i nct, cov;
	int pos, mpos, cpos, rpos;
	u4i edge, erev;
	u4i next, prev, header;
	u4i mmidx;
} bspoanode_t;
define_list(bspoanodev, bspoanode_t);

// put pair of edges together, idx % 2 = 0, 1
typedef struct {
	u4i node;
	u4i cov:31, vst:1;
	u4i next;
} bspoaedge_t;
define_list(bspoaedgev, bspoaedge_t);

typedef struct {
	u4i node1, node2;
	int cov;
} bspoaedge_chg_t;
define_list(bspoaedgechgv, bspoaedge_chg_t);

typedef struct {
	int refmode; // 0
	int shuffle; // 1
	int alnmode; // SEQALIGN_MODE_OVERLAP
	int realn; // 3
	u4i seqcore; // 40, core seqs in constructing primary MSA, more reads will be fast map to primary MSA
	int nrec; // 20, new read will align against previous nrec reads on POA
	int ksz; // 15
	int bwtrigger; // 1
	int bandwidth; // 128, when tiggering bandwidth, first generate cns, then align read against cns to get the offset of each node' band
	int M, X, O, E, Q, P, T; // 2, -6, -3, -2, -8, -1, 20
	int refbonus; // 1
	int editbw; // 64
	int althi, qlthi; // 5: 1e-0.5, 70:1e-7
	float psub, pins, pdel, piex, pdex, hins, hdel; // 0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.40
	int min_varcnt; // 3
	float min_covfrq; // 0.5
	u4i min_snvqlt; // 5: 1e-5
} BSPOAPar;

static const BSPOAPar DEFAULT_BSPOA_PAR = (BSPOAPar){0, 1, SEQALIGN_MODE_OVERLAP, 3, 40, 20, 15, 1, 128,
	2, -6, -3, -2, -8, -1, 20, 1, 64, 5, 70,
	0.10, 0.10, 0.15, 0.15, 0.20, 0.20, 0.40, 3, 0.5, 5};

typedef struct {
	u4i coff:29, bt:3;
	float max;
} bspoacns_t;

#define BSPOA_MIN_LOGVAL	-1000000000
#define BSPOA_QLT_MAX	90

typedef struct {
	u2i rid;
	u4i rbeg, mbeg;
	u2i rlen, mlen;
	float scr;
} bspoalsp_t;
define_list(bspoalspv, bspoalsp_t);

typedef struct {
	u4i cpos, mpos;
	u2i refn, altn, covn;
	u1i refb, altb;
	u4i qual;
} bspoavar_t;
define_list(bspoavarv, bspoavar_t);

typedef struct {
	String     *mtag;
	SeqBank    *seqs;
	u4v        *ords;
	u4v        *ndoffs;
	int        keep_seqs; // if ture, don't clear seqs in beg_bspoa; then set to false when used
	u4v        *cigars, *cgbs, *cges; // refmode, concatenate all cigars together
	u4i         HEAD, TAIL;
	bspoanodev *nodes;
	bspoaedgev *edges;
	u4v        *ecycs;
	uuhash     *hash;
	BSPOAPar   *par;
	int piecewise;
	u4i  bandwidth, qlen, slen, qb, qe, nmsa, nrds; // real bandwidth used in alignment
	u1v *qseq;
	b1i matrix[4][16];
	b1v *qprof[4];
	b1v *memp; // memory pool
	u8i mmblk, mmcnt;
	int maxscr, maxidx, maxoff;
	u4v    *sels; // selected node idx
	b4v    *rdregs[2];
	BitVec *rdbits, *states;
	u8v *heap, *todels;
	u4v *stack;
	bspoaedgechgv *echgs;
	u4i backbone;
	u1v *msacols;
	u4v *msaidxs;
	double dpvals[8], dporis[8];
	u1i *dptable;
	u1v *cns, *qlt, *alt;
	bspoalspv *lsp;
	bspoavarv *var;
	String *strs;
	u8i ncall, obj_uid;
} BSPOA;

static inline void gen_cns_aln_event_table_bspoa(BSPOAPar *par, double ps[8], double os[8], u1i *table){
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
					//if(b == c && ps[6] > ps[1]){
						//table[i] = (6 << 3) | 0; // ps[HI] << 3 | M
					//} else {
						table[i] = (1 << 3) | 0; // ps[X] << 3 | M
					//}
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
	g->ords = init_u4v(32);
	g->ndoffs = init_u4v(1024);
	g->keep_seqs = 0;
	g->cigars = init_u4v(1024);
	g->cgbs   = init_u4v(32);
	g->cges   = init_u4v(32);
	g->HEAD   = 0;
	g->TAIL   = 1;
	g->nodes = init_bspoanodev(16 * 1024);
	g->edges = init_bspoaedgev(16 * 1024);
	g->ecycs = init_u4v(32);
	g->hash  = init_uuhash(13);
	g->par   = malloc(sizeof(BSPOAPar));
	memcpy(g->par, &par, sizeof(BSPOAPar));
	g->par->bandwidth = roundup_times(g->par->bandwidth, WORDSIZE);
	g->piecewise = 1;
	g->nmsa      = 0;
	g->nrds      = 0;
	g->bandwidth = 0;
	g->qseq  = init_u1v(1024);
	g->qlen  = 0;
	g->qprof[0] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->qprof[1] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->qprof[2] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->qprof[3] = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->memp  = adv_init_b1v(1024, 0, WORDSIZE, 0);
	g->mmblk = 0;
	g->mmcnt = 0;
	g->sels  = init_u4v(1024);
	g->rdregs[0] = init_b4v(32);
	g->rdregs[1] = init_b4v(32);
	g->rdbits = init_bitvec(1024);
	g->states = init_bitvec(1024);
	g->heap = init_u8v(32);
	g->todels = init_u8v(4);
	g->stack = init_u4v(32);
	g->echgs = init_bspoaedgechgv(8);
	g->backbone = 0;
	g->msacols = init_u1v(16 * 1024);
	g->msaidxs = init_u4v(1024);
	g->dptable = malloc(5 * 5 * 5 * 5);
	gen_cns_aln_event_table_bspoa(g->par, g->dpvals, g->dporis, g->dptable);
	g->strs = init_string(1024);
	g->cns = init_u1v(1024);
	g->qlt = init_u1v(1024);
	g->alt = init_u1v(1024);
	g->lsp = init_bspoalspv(1024);
	g->var = init_bspoavarv(8);
	g->ncall = 0;
	g->obj_uid = 0;
	return g;
}

static inline void renew_bspoa(BSPOA *g){
	if(g->keep_seqs){
		g->keep_seqs = 0;
	} else {
		free_seqbank(g->seqs); g->seqs = init_seqbank();
		renew_u4v(g->cigars, 1024);
		renew_u4v(g->cgbs, 32);
		renew_u4v(g->cges, 32);
		renew_u1v(g->cns, 1024);
		renew_u1v(g->qlt, 1024);
		renew_u1v(g->alt, 1024);
	}
	renew_u4v(g->ords, 1024);
	renew_u4v(g->ndoffs, 1024);
	renew_bspoanodev(g->nodes, 16 * 1024);
	renew_bspoaedgev(g->edges, 16 * 1024);
	ZEROS(next_ref_bspoaedgev(g->edges));
	ZEROS(next_ref_bspoaedgev(g->edges)); // even idx for edge
	renew_u4v(g->ecycs, 32);
	free_uuhash(g->hash); g->hash = init_uuhash(13);
	renew_u1v(g->qseq, 1024);
	renew_b1v(g->qprof[0], 4 * 1024);
	renew_b1v(g->qprof[1], 4 * 1024);
	renew_b1v(g->qprof[2], 4 * 1024);
	renew_b1v(g->qprof[3], 4 * 1024);
	renew_b1v(g->memp, 4 * 1024);
	g->nmsa  = 0;
	g->nrds  = 0;
	g->mmblk = 0;
	g->mmcnt = 0;
	renew_u4v(g->sels, 1024);
	renew_b4v(g->rdregs[0], 32);
	renew_b4v(g->rdregs[1], 32);
	renew_bitvec(g->rdbits, 1024);
	renew_bitvec(g->states, 1024);
	renew_u8v(g->heap, 32);
	renew_u4v(g->stack, 32);
	renew_bspoaedgechgv(g->echgs, 32);
	renew_u1v(g->msacols, 16 * 1024);
	renew_u4v(g->msaidxs, 1024);
	renew_bspoalspv(g->lsp, 1024);
	renew_bspoavarv(g->var, 32);
	recap_string(g->strs, 1024);
}

static inline void clear_bspoa(BSPOA *g){
	if(g->keep_seqs){
		g->keep_seqs = 0;
	} else {
		clear_seqbank(g->seqs);
		clear_u4v(g->cigars);
		clear_u4v(g->cgbs);
		clear_u4v(g->cges);
		clear_u1v(g->cns);
		clear_u1v(g->qlt);
		clear_u1v(g->alt);
	}
	clear_u4v(g->ords);
	clear_u4v(g->ndoffs);
	clear_bspoanodev(g->nodes);
	clear_bspoaedgev(g->edges);
	ZEROS(next_ref_bspoaedgev(g->edges));
	ZEROS(next_ref_bspoaedgev(g->edges)); // even idx for edge
	clear_u4v(g->ecycs);
	clear_uuhash(g->hash);
	clear_u1v(g->qseq);
	clear_b1v(g->qprof[0]);
	clear_b1v(g->qprof[1]);
	clear_b1v(g->qprof[2]);
	clear_b1v(g->qprof[3]);
	clear_b1v(g->memp);
	g->nrds  = 0;
	g->nmsa = 0;
	g->mmblk = 0;
	g->mmcnt = 0;
	clear_u4v(g->sels);
	clear_b4v(g->rdregs[0]);
	clear_b4v(g->rdregs[1]);
	clear_bitvec(g->rdbits);
	clear_bitvec(g->states);
	clear_u8v(g->heap);
	clear_u4v(g->stack);
	clear_bspoaedgechgv(g->echgs);
	clear_u1v(g->msacols);
	clear_u4v(g->msaidxs);
	clear_bspoalspv(g->lsp);
	clear_bspoavarv(g->var);
	clear_string(g->strs);
}

static inline void free_bspoa(BSPOA *g){
	free_string(g->mtag);
	free_seqbank(g->seqs);
	free_u4v(g->ords);
	free_u4v(g->ndoffs);
	free_u4v(g->cigars);
	free_u4v(g->cgbs);
	free_u4v(g->cges);
	free_bspoanodev(g->nodes);
	free_bspoaedgev(g->edges);
	free_u4v(g->ecycs);
	free_uuhash(g->hash);
	free_u1v(g->qseq);
	free_b1v(g->qprof[0]);
	free_b1v(g->qprof[1]);
	free_b1v(g->qprof[2]);
	free_b1v(g->qprof[3]);
	free_b1v(g->memp);
	free_u4v(g->sels);
	free_b4v(g->rdregs[0]);
	free_b4v(g->rdregs[1]);
	free_bitvec(g->rdbits);
	free_bitvec(g->states);
	free_u8v(g->heap);
	free_u8v(g->todels);
	free_u4v(g->stack);
	free_bspoaedgechgv(g->echgs);
	free_u1v(g->msacols);
	free_u4v(g->msaidxs);
	free(g->dptable);
	free_u1v(g->cns);
	free_u1v(g->qlt);
	free_u1v(g->alt);
	free_bspoalspv(g->lsp);
	free_bspoavarv(g->var);
	free_string(g->strs);
	free(g->par);
	free(g);
}

static inline bspoanode_t* new_node_bspoa(BSPOA *g, u2i rid, int pos, u1i base){
	bspoanode_t *u;
	u = next_ref_bspoanodev(g->nodes);
	ZEROS(u);
	u->rid  = rid;
	u->pos  = pos;
	u->base = base;
	u->cov  = 1;
	u->next    = offset_bspoanodev(g->nodes, u);
	u->prev    = offset_bspoanodev(g->nodes, u);
	u->header  = offset_bspoanodev(g->nodes, u);
	return u;
}

static inline bspoanode_t* get_rdnode_bspoa(BSPOA *g, u2i rid, int pos){
	return ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + pos);
}

static inline bspoanode_t* header_node_bspoa(BSPOA *g, bspoanode_t *v){
	return ref_bspoanodev(g->nodes, v->header);
}

static inline void set_nodecov_bspoa(BSPOA *g, bspoanode_t *u, u2i nodecov){
	bspoanode_t *v;
	v = ref_bspoanodev(g->nodes, u->header);
	v->cov = nodecov;
}

static inline u2i get_nodecov_bspoa(BSPOA *g, bspoanode_t *u){
	bspoanode_t *v;
	v = ref_bspoanodev(g->nodes, u->header);
	return v->cov;
}

static inline bspoaedge_t* get_edge_bspoa(BSPOA *g, bspoanode_t *_u, bspoanode_t *_v){
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	u4i eidx;
	u = ref_bspoanodev(g->nodes, _u->header);
	v = ref_bspoanodev(g->nodes, _v->header);
	eidx = u->edge;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		if(e->node == offset_bspoanodev(g->nodes, v)){
			return e;
		}
	}
	return NULL;
}

static inline bspoaedge_t* new_edge_bspoa(BSPOA *g, bspoanode_t *u, bspoanode_t *v, int cov){
	bspoaedge_t *e, *r;
	if(g->ecycs->size){
		e = ref_bspoaedgev(g->edges, g->ecycs->buffer[--g->ecycs->size]);
		r = e + 1;
	} else {
		encap_bspoaedgev(g->edges, 2);
		e = next_ref_bspoaedgev(g->edges);
		r = next_ref_bspoaedgev(g->edges);
	}
	ZEROS(e);
	ZEROS(r);
	e->node = offset_bspoanodev(g->nodes, v);
	r->node = offset_bspoanodev(g->nodes, u);
	e->cov = cov;
	r->cov = cov;
	return e;
}

static inline void _add_edge_bspoacore(BSPOA *g, bspoanode_t *v, u4i eidx){
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

static inline void _del_edge_bspoacore(BSPOA *g, bspoanode_t *v, u4i eidx){
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

static inline bspoaedge_t* chg_edge_bspoa(BSPOA *g, bspoanode_t *_u, bspoanode_t *_v, int cov, int *exists){
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	u4i eidx;
	int ncov;
	if(cov == 0 || _u == NULL || _v == NULL){
		return NULL;
	}
	u = ref_bspoanodev(g->nodes, _u->header);
	v = ref_bspoanodev(g->nodes, _v->header);
#if DEBUG
	if(v->header == g->HEAD){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	} else if(u->header == g->TAIL){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	if(u == v) return NULL;
	e = get_edge_bspoa(g, u, v);
	if(e == NULL){
		if(exists) exists[0] = 0;
		ncov = cov;
	} else {
		if(exists) exists[0] = 1;
		eidx = offset_bspoaedgev(g->edges, e);
		ncov = e->cov + cov;
		_del_edge_bspoacore(g, u, eidx);
		_del_edge_bspoacore(g, v, eidx + 1);
	}
	if(ncov > 0){
		e = new_edge_bspoa(g, u, v, ncov);
		eidx = offset_bspoaedgev(g->edges, e);
		_add_edge_bspoacore(g, u, eidx);
		_add_edge_bspoacore(g, v, eidx + 1);
#if DEBUG
		if(e == NULL){ // for debug
			check_all_node_edges_bspoa(g);
		}
#endif
		return e;
	} else {
		return NULL;
	}
}

static inline int is_all_sep_rdnodes_bspoa(BSPOA *g, u2i rid, int rbeg, int rlen){
	bspoanode_t *v;
	u4i nidx;
	int i;
	for(i=0;i<rlen;i++){
		v = get_rdnode_bspoa(g, rid, rbeg + i);
		nidx = offset_bspoanodev(g->nodes, v);
		if(v->next != nidx){
			return 0;
		}
	}
	return 1;
}

static inline void connect_rdnode_bspoa(BSPOA *g, u2i rid, int pos){
	bspoanode_t *u, *v;
	u = get_rdnode_bspoa(g, rid, pos - 1);
	v = get_rdnode_bspoa(g, rid, pos);
	if(v->rdc){
		return;
	}
	chg_edge_bspoa(g, u, v, 1, NULL);
	u->rdd = v->rdc = 1;
}

static inline void connect_rdnodes_bspoa(BSPOA *g, u2i rid){
	u4i i, rdlen;
	rdlen = g->seqs->rdlens->buffer[rid];
	for(i=0;i<=rdlen;i++){
		connect_rdnode_bspoa(g, rid, i);
	}
}

static inline void disconnect_rdnode_bspoa(BSPOA *g, u2i rid, int pos){
	bspoanode_t *u, *v;
	u = get_rdnode_bspoa(g, rid, pos - 1);
	if(u->rdd == 0){
		return;
	}
	v = get_rdnode_bspoa(g, rid, pos);
	chg_edge_bspoa(g, u, v, -1, NULL);
	u->rdd = v->rdc = 0;
}

static inline void disconnect_rdnodes_bspoa(BSPOA *g, u2i rid){
	u4i i, rdlen;
	rdlen = g->seqs->rdlens->buffer[rid];
	for(i=0;i<=rdlen;i++){
		disconnect_rdnode_bspoa(g, rid, i);
	}
}

static inline int isconnect_rdnode_bspoa(BSPOA *g, u2i rid, int pos){
	bspoanode_t *v;
	v = get_rdnode_bspoa(g, rid, pos);
	return v->rdc;
}

static inline void check_rdc_bspoa(BSPOA *g){
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	u4i nidx;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		v = ref_bspoanodev(g->nodes, nidx);
		if(v->rdc){
			w = get_rdnode_bspoa(g, v->rid, v->pos - 1);
			if((e = get_edge_bspoa(g, w, v)) == NULL){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
}

#define BSPOA_EMOVTYPE_MOVALL	0x0F0F // all edges
#define BSPOA_EMOVTYPE_KPTONE	0x1E0F // all edges but retaining spec one
#define BSPOA_EMOVTYPE_MOVONE	0xE1F0 // mov spec one edges but retianing rest
// Not safe
//#define BSPOA_EMOVTYPE_DELMOV	0x0E0F // del spec one and move others
//#define BSPOA_EMOVTYPE_DELKPT	0xE0F0 // del spec one and keep others
// u and v might be aligned
static inline void _mov_node_edges_bspoacore(BSPOA *g, bspoanode_t *u, bspoanode_t *v, u4i _spec_node, int dir, int movtype){
	bspoanode_t *w, *x;
	bspoaedge_t *e;
	bspoaedge_chg_t *chg;
	u4i i, j, eidx, spec_node;
	int ecov, covs[4];
	if(_spec_node < g->nodes->size){
		w = ref_bspoanodev(g->nodes, _spec_node);
		spec_node = w->header;
	} else spec_node = _spec_node;
	eidx = dir? u->erev : u->edge;
	clear_bspoaedgechgv(g->echgs);
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		ecov = e->cov;
		w = ref_bspoanodev(g->nodes, e->node);
		if(e->node == spec_node){
			covs[0] = 0; covs[1] = ecov;
		} else {
			covs[0] = ecov; covs[1] = 0;
		}
		covs[2] = covs[3] = 0;
		for(i=0;i<2;i++){
			for(j=0;j<2;j++){
				switch((movtype >> (4 * (i * 2 + j))) & 0xF){
					case 0xF: covs[3 - j] += covs[i]; break;
					case 0xE: covs[3 - j] += num_max(covs[i] - 1, 0); break;
					case 0x1: covs[3 - j] += num_min(covs[i], 1); break;
				}
			}
		}
		if(dir){
			push_bspoaedgechgv(g->echgs, (bspoaedge_chg_t){w? offset_bspoanodev(g->nodes, w) : MAX_U4, u? offset_bspoanodev(g->nodes, u) : MAX_U4, covs[2] - ecov});
			push_bspoaedgechgv(g->echgs, (bspoaedge_chg_t){w? offset_bspoanodev(g->nodes, w) : MAX_U4, v? offset_bspoanodev(g->nodes, v) : MAX_U4, covs[3]});
		} else {
			push_bspoaedgechgv(g->echgs, (bspoaedge_chg_t){u? offset_bspoanodev(g->nodes, u) : MAX_U4, w? offset_bspoanodev(g->nodes, w) : MAX_U4, covs[2] - ecov});
			push_bspoaedgechgv(g->echgs, (bspoaedge_chg_t){v? offset_bspoanodev(g->nodes, v) : MAX_U4, w? offset_bspoanodev(g->nodes, w) : MAX_U4, covs[3]});
		}
	}
	for(i=0;i<g->echgs->size;i++){
		chg = ref_bspoaedgechgv(g->echgs, i);
		w = chg->node1 == MAX_U4? NULL : ref_bspoanodev(g->nodes, chg->node1);
		x = chg->node2 == MAX_U4? NULL : ref_bspoanodev(g->nodes, chg->node2);
		chg_edge_bspoa(g, w, x, chg->cov, NULL);
	}
}

#define BSPOA_RDNODE_CUTEDGE	1
#define BSPOA_RDNODE_CUTNODE	2
#define BSPOA_RDNODE_CUTALL		3

static inline bspoanode_t* cut_rdnode_bspoa(BSPOA *g, u2i rid, int pos, int cut){
	bspoanode_t *u, *x;
	u4i nidx, header[2], nodes[2], nodecov;
	u = get_rdnode_bspoa(g, rid, pos);
	nidx = offset_bspoanodev(g->nodes, u);
	nodes[0] = nidx + 1;
	nodes[1] = nidx - 1;
	header[0] = u->header;
	header[1] = u->prev;
	nodecov   = get_nodecov_bspoa(g, u);
	if((cut & BSPOA_RDNODE_CUTNODE) && (u->next != nidx)){
		x = ref_bspoanodev(g->nodes, u->prev);
		x->next = u->next;
		x = ref_bspoanodev(g->nodes, u->next);
		x->prev = u->prev;
		u->next = nidx;
		u->prev = nidx;
		u->header = nidx;
		if(header[0] == nidx){
			x = ref_bspoanodev(g->nodes, header[1]);
			// change header of the rest
			while(1){
				x->header = header[1];
				if(x->next == header[1]) break;
				x = ref_bspoanodev(g->nodes, x->next);
			}
			x = ref_bspoanodev(g->nodes, header[1]);
			if(u->rdd){
				_mov_node_edges_bspoacore(g, u, x, nodes[0], 0, BSPOA_EMOVTYPE_KPTONE);
			} else {
				_mov_node_edges_bspoacore(g, u, x, nodes[0], 0, BSPOA_EMOVTYPE_MOVALL);
			}
			if(u->rdc){
				_mov_node_edges_bspoacore(g, u, x, nodes[1], 1, BSPOA_EMOVTYPE_KPTONE);
			} else {
				_mov_node_edges_bspoacore(g, u, x, nodes[1], 1, BSPOA_EMOVTYPE_MOVALL);
			}
		} else {
			x = ref_bspoanodev(g->nodes, header[0]);
			if(u->rdd){
				_mov_node_edges_bspoacore(g, x, u, nodes[0], 0, BSPOA_EMOVTYPE_MOVONE);
			}
			if(u->rdc){
				_mov_node_edges_bspoacore(g, x, u, nodes[1], 1, BSPOA_EMOVTYPE_MOVONE);
			}
		}
		set_nodecov_bspoa(g, x, nodecov - 1);
		set_nodecov_bspoa(g, u, 1);
	}
	if(cut & BSPOA_RDNODE_CUTEDGE){
		disconnect_rdnode_bspoa(g, rid, pos);
		disconnect_rdnode_bspoa(g, rid, pos + 1);
	}
	return u;
}

static inline bspoanode_t* merge_nodes_bspoa(BSPOA *g, bspoanode_t *n1, bspoanode_t *n2){
	bspoanode_t *ns[2], *xs[3];
	u4i nidxs[2], f;
	u2i nodecov;
	ns[0] = ref_bspoanodev(g->nodes, n1->header);
	ns[1] = ref_bspoanodev(g->nodes, n2->header);
	if(ns[0] == ns[1]) return ns[0];
	nodecov = ns[0]->cov + ns[1]->cov;
#if DEBUG
	if(ns[0]->rid == ns[1]->rid){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(ns[1]->header == g->HEAD && ns[0]->pos >= 0){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(ns[0]->header == g->HEAD && ns[1]->pos >= 0){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	if(ns[0]->cov < ns[1]->cov){
		swap_var(ns[0], ns[1]);
	} else if(ns[0]->cov > ns[1]->cov){
	} else if(ns[0]->rid > ns[1]->rid){
		swap_var(ns[0], ns[1]);
	}
	_mov_node_edges_bspoacore(g, ns[1], ns[0], MAX_U4, 0, BSPOA_EMOVTYPE_MOVALL);
	_mov_node_edges_bspoacore(g, ns[1], ns[0], MAX_U4, 1, BSPOA_EMOVTYPE_MOVALL);
	nidxs[0] = offset_bspoanodev(g->nodes, ns[0]);
	nidxs[1] = offset_bspoanodev(g->nodes, ns[1]);
	if(0){
		xs[0] = ns[0];
		xs[1] = ns[1];
		do {
			do {
				xs[2] = ref_bspoanodev(g->nodes, xs[0]->next);
				if(xs[2] == ns[0]) break;
				if(xs[2]->rid > xs[1]->rid) break;
				else if(xs[2]->rid == xs[1]->rid){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				xs[0] = xs[2];
			} while(1);
			xs[2] = xs[1];
			xs[1] = ref_bspoanodev(g->nodes, xs[1]->next);
			xs[2]->next = xs[0]->next;
			xs[2]->header  = offset_bspoanodev(g->nodes, ns[0]);
			xs[0]->next = offset_bspoanodev(g->nodes, xs[2]);
		} while(xs[1] != ns[1]);
		set_nodecov_bspoa(g, ns[0], nodecov);
	} else if(0){
		xs[0] = ns[0]; xs[1] = ns[1]; xs[2] = NULL;
		while(1){
			f = 0;
			if(xs[0]){
				if(xs[1]){
					if(xs[0]->rid == xs[1]->rid){
						fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					} else if(xs[0]->rid > xs[1]->rid){
						f = 1;
					}
				}
			} else if(xs[1]){
				f = 1;
			} else {
				break;
			}
			if(xs[2]) xs[2]->next = offset_bspoanodev(g->nodes, xs[f]);
			xs[2] = xs[f];
			if(xs[f]->next == nidxs[f]){
				xs[f] = NULL;
			} else {
				xs[f] = ref_bspoanodev(g->nodes, xs[f]->next);
			}
			xs[2]->header = nidxs[0];
			xs[2]->cov    = nodecov;
		}
		if(xs[2]) xs[2]->next = nidxs[0];
	} else {
		ns[0]->cov = nodecov;
		xs[2] = ns[1];
		do {
			xs[2]->header = nidxs[0];
			xs[2] = ref_bspoanodev(g->nodes, xs[2]->next);
		} while(xs[2] != ns[1]);
		xs[0] = ref_bspoanodev(g->nodes, ns[0]->prev);
		xs[1] = ref_bspoanodev(g->nodes, ns[1]->prev);
		ns[0]->prev = offset_bspoanodev(g->nodes, xs[1]);
		ns[1]->prev = offset_bspoanodev(g->nodes, xs[0]);
		xs[1]->next = nidxs[0];
		xs[0]->next = nidxs[1];
	}
	return ns[0];
}

static inline void cut_rdnodes_bspoa(BSPOA *g, u2i rid,  int cut){
	u4i i, len;
	len = g->seqs->rdlens->buffer[rid];
	for(i=0;i<len;i++){
		cut_rdnode_bspoa(g, rid, i, cut);
	}
}

static inline void _push_bspoa_core(BSPOA *g, u4i len, u4i *cgs, u4i ncg){
	u4i rid, off;
	rid = g->seqs->nseq - 1;
	off = g->seqs->rdoffs->buffer[rid];
	push_u4v(g->cgbs, g->cigars->size);
	append_array_u4v(g->cigars, cgs, ncg);
	push_u4v(g->cges, g->cigars->size);
	push_b4v(g->rdregs[0], 0);
	push_b4v(g->rdregs[1], len);
	encap_bitvec(g->states, g->nodes->size);
}

static inline void _add_read_bspoa_core(BSPOA *g, u2i rid){
	bspoanode_t *u, *v;
	u4i i, off, len;
	off = g->seqs->rdoffs->buffer[rid];
	len = g->seqs->rdlens->buffer[rid];
	new_node_bspoa(g, rid, -1, 4);
	push_u4v(g->ndoffs, g->nodes->size);
	for(i=0;i<len;i++){
		new_node_bspoa(g, rid, i, get_basebank(g->seqs->rdseqs, off + i));
	}
	new_node_bspoa(g, rid, len, 4);
	if(rid == 0){
		g->backbone = len;
		g->HEAD = g->ndoffs->buffer[rid] - 1;
		v = get_rdnode_bspoa(g, rid, -1);
		v->cpos = 0;
		g->TAIL = g->ndoffs->buffer[rid] + len;
		v = get_rdnode_bspoa(g, rid, len);
		v->cpos = len;
		for(i=0;i<len;i++){
			v = get_rdnode_bspoa(g, 0, i);
			v->ref   = 1;
			v->bless = 1;
			v->cpos  = i;
			connect_rdnode_bspoa(g, 0, i);
		}
		connect_rdnode_bspoa(g, 0, len);
	} else {
		u = ref_bspoanodev(g->nodes, g->HEAD);
		v = get_rdnode_bspoa(g, rid, -1);
		merge_nodes_bspoa(g, u, v);
		u = ref_bspoanodev(g->nodes, g->TAIL);
		v = get_rdnode_bspoa(g, rid, len);
		merge_nodes_bspoa(g, u, v);
	}
}

static inline void push_bspoacore(BSPOA *g, char *seq, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		push_seqbank(g->seqs, NULL, 0, seq, len);
		_push_bspoa_core(g, len, cgs, ncg);
	}
}

static inline void push_bspoa(BSPOA *g, char *seq, u4i len){ push_bspoacore(g, seq, len, NULL, 0); }

static inline void fwdbitseqpush_bspoacore(BSPOA *g, u1i *seq, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		fwdbitseqpush_seqbank(g->seqs, NULL, 0, seq, len);
		_push_bspoa_core(g, len, cgs, ncg);
	}
}

static inline void fwdbitseqpush_bspoa(BSPOA *g, u1i *seq, u4i len){ fwdbitseqpush_bspoacore(g, seq, len, NULL, 0); }

static inline void revbitseqpush_bspoacore(BSPOA *g, u1i *seq, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		revbitseqpush_seqbank(g->seqs, NULL, 0, seq, len);
		_push_bspoa_core(g, len, cgs, ncg);
	}
}

static inline void revbitseqpush_bspoa(BSPOA *g, u1i *seq, u4i len){ revbitseqpush_bspoacore(g, seq, len, NULL, 0); }

static inline void fwdbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		fwdbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		_push_bspoa_core(g, len, cgs, ncg);
	}
}

static inline void fwdbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ fwdbitpush_bspoacore(g, bits, off, len, NULL, 0); }

static inline void revbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u4i *cgs, u4i ncg){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX){
		len = num_min(len, BSPOA_RDLEN_MAX);
		revbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		_push_bspoa_core(g, len, cgs, ncg);
	}
}

static inline void revbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ revbitpush_bspoacore(g, bits, off, len, NULL, 0); }

static inline void shuffle_reads_by_kmers_bspoa(BSPOA *g){
	SeqBank *sb;
	uuhash *khash;
	u4v *kidxs;
	f4i *kords;
	//u2v *rbegs, *rends;
	uuhash_t *u;
	u8i roff;
	u4i ridx, i, ksize, kmer, kmask, rlen, khit, mincnt, *krevs;
	double logv;
	int exists;
	sb = g->seqs;
	if(sb->nseq == 0) return;
	khash = g->hash;
	kidxs = g->ords;
	ksize = g->par->ksz;
	clear_u4v(kidxs);
	clear_uuhash(khash);
	resize_b1v(g->memp, sb->nseq * sizeof(f4i));
	kords = (f4i*)g->memp->buffer;
	kmask = MAX_U4 >> ((16 - ksize) << 1);
	mincnt = g->par->refmode? 1 : num_max(2, Int(sb->nseq * 0.05));
	for(ridx=0;ridx<sb->nseq;ridx++){
		rlen = sb->rdlens->buffer[ridx];
		kmer = 0;
		roff = sb->rdoffs->buffer[ridx];
		for(i=0;i<rlen;i++){
			kmer = ((kmer << 2) | get_basebank(sb->rdseqs, roff + i)) & kmask;
			if(i + 1 < ksize) continue;
			u = prepare_uuhash(khash, kmer, &exists);
			if(exists){
				if(((u->val >> 16) & 0x7FFFU) == ridx + 1){
					u->val |= 1U << 31;
				} else {
					u->val = (u->val & 0x8000FFFFU) | ((ridx + 1) << 16);
				}
				u->val ++;
			} else {
				u->key = kmer;
				u->val = (0U << 31) | ((ridx + 1) << 16) | 1;
			}
		}
		if(g->par->refmode) break;
	}
	logv = log(1.2);
	kords[0] = 3e+38F;
	push_u4v(kidxs, 0);
	for(ridx=1;ridx<sb->nseq;ridx++){
		rlen = sb->rdlens->buffer[ridx];
		kmer = 0;
		roff = sb->rdoffs->buffer[ridx];
		khit = 0;
		for(i=0;i<rlen;i++){
			kmer = ((kmer << 2) | get_basebank(sb->rdseqs, roff + i)) & kmask;
			if(i + 1 < ksize) continue;
			u = get_uuhash(khash, kmer);
			if(u && (u->val & 0x80000000U) == 0 && (u->val & 0xFFFFU) >= mincnt){
				khit ++;
			}
		}
		if(g->par->refmode){
			kords[ridx] = ((double)khit) * logv / log(num_max(rlen, sb->rdlens->buffer[0]));
		} else {
			kords[ridx] = ((double)khit) * logv / log(rlen);
		}
		push_u4v(kidxs, ridx);
	}
	sort_array(kidxs->buffer, kidxs->size, u4i, num_cmpgt(kords[b], kords[a]));
	if(_DEBUG_LOG_ > 1){
		for(i=0;i<kidxs->size;i++){
			fprintf(_DEBUG_LOGFILE_, "SHUFFLE[%u] %u\t%u\t%0.4f\n", i, kidxs->buffer[i], sb->rdlens->buffer[kidxs->buffer[i]], kords[kidxs->buffer[i]]);
		}
	}
	encap_u4v(kidxs, kidxs->size * 2);
	krevs = kidxs->buffer + kidxs->size;
#if DEBUG
	if(1){
		u8v *rdoffs = init_u8v(32);
		u4v *shfs = init_u4v(sb->nseq * 2);
		append_u8v(rdoffs, sb->rdoffs);
		append_u4v(shfs, kidxs);
		append_u4v(shfs, kidxs);
		for(i=0;i<sb->nseq;i++) shfs->buffer[shfs->buffer[i + sb->nseq]] = i;
		shuffle_array(sb->rdoffs->buffer, kidxs->buffer, krevs, kidxs->size);
		shuffle_array(sb->rdoffs->buffer, shfs->buffer, shfs->buffer + sb->nseq, kidxs->size);
		if(memcmp(sb->rdoffs->buffer, rdoffs->buffer, rdoffs->size * sizeof(u8i))){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
#endif
	shuffle_array(sb->rdoffs->buffer, kidxs->buffer, krevs, kidxs->size);
	shuffle_array(sb->rdlens->buffer, kidxs->buffer, krevs, kidxs->size);
	// skip sb->rdatgs and sb->rdhash
	//shuffle_array(g->ndoffs->buffer, kidxs->buffer, krevs, kidxs->size);
	shuffle_array(g->cgbs->buffer, kidxs->buffer, krevs, kidxs->size);
	shuffle_array(g->cges->buffer, kidxs->buffer, krevs, kidxs->size);
}

static inline void restore_rd_orders_bspoa(BSPOA *g){
	bspoanode_t *v, *ve;
	u4i *ords, *revs, i, nseq;
	nseq = g->seqs->nseq;
	append_u4v(g->ords, g->ords);
	ords = g->ords->buffer;
	revs = ords + nseq;
	for(i=0;i<nseq;i++) ords[revs[i]] = i;
	shuffle_array(g->seqs->rdoffs->buffer, ords, revs, nseq);
	shuffle_array(g->seqs->rdlens->buffer, ords, revs, nseq);
	shuffle_array(g->ndoffs->buffer, ords, revs, nseq);
	shuffle_array(g->cgbs->buffer, ords, revs, nseq);
	shuffle_array(g->cges->buffer, ords, revs, nseq);
	for(i=0;i<nseq;i++){
		v  = get_rdnode_bspoa(g, i, -1);
		ve = get_rdnode_bspoa(g, i, g->seqs->rdlens->buffer[i]);
		while(v <= ve){
			v->rid = i;
			v ++;
		}
	}
}

static inline u4i count_node_edges_bspoa(BSPOA *g, bspoanode_t *v, int dir){
	bspoaedge_t *e;
	u4i eidx, ret;
	ret = 0;
	eidx = dir? v->erev : v->edge;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		ret ++;
	}
	return ret;
}

static inline void _del_edge_bspoa(BSPOA *g, bspoanode_t *_v, u4i eidx){
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	v = ref_bspoanodev(g->nodes, _v->header);
	e = ref_bspoaedgev(g->edges, eidx);
	w = ref_bspoanodev(g->nodes, e->node);
	_del_edge_bspoacore(g, v, eidx);
	_del_edge_bspoacore(g, w, eidx ^ 0x1);
#if DEBUG
	if(eidx == 0){ // for debug
		check_all_node_edges_bspoa(g);
	}
#endif
}

static inline void set_edges_unvisited_bspoa(BSPOA *g){
	bspoaedge_t *e;
	u4i eidx;
	for(eidx=0;eidx<g->edges->size;eidx++){
		e = ref_bspoaedgev(g->edges, eidx);
		e->vst = 0;
	}
}

static inline void print_dot_bspoa(BSPOA *g, int posbeg, int posend, u4i mincnt, FILE *out){
	bspoanode_t *n;
	bspoaedge_t *e;
	u4i nidx, eidx, nodecov;
	fprintf(out, "digraph {\n");
	fprintf(out, "rankdir=LR\n");
	fprintf(out, "N0 [label=\"BEG\"]\n");
	fprintf(out, "N1 [label=\"END\"]\n");
	for(nidx=g->TAIL+1;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->mpos < posbeg || n->mpos >= posend) continue;
#if 1
		if(n->nin == 0 && n->nou == 0) continue;
		nodecov = get_nodecov_bspoa(g, n);
		if(nodecov >= mincnt){
			fprintf(out, "N%u [label=%c%u_%d_%d_N%u color=blue]\n", nidx, bit_base_table[(n->base) & 0x03], n->mpos, n->cpos, nodecov, nidx);
		} else {
			fprintf(out, "N%u [label=%c%u_%d_%d_N%u]\n", nidx, bit_base_table[(n->base) & 0x03], n->mpos, n->cpos, nodecov, nidx);
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
		if(n->next != nidx){
			fprintf(out, "N%u -> N%u [color=magenta style=dashed]\n", nidx, n->next);
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
		if(get_bitvec(g->states, nidx) == 0) continue;
		n = ref_bspoanodev(g->nodes, nidx);
		fprintf(out, "N%u [label=\"N%u:%u:%u:%d\"]\n", nidx, nidx, n->nin, n->nct, n->vst);
		eidx = n->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(get_bitvec(g->states, e->node) == 0) continue;
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
		fprintf(out, "N%u [label=%c%u_%d_%d_N%u%s]\n", nidx, bit_base_table[(n->base) & 0x03], n->rpos, n->cpos, get_nodecov_bspoa(g, n), nidx, (nidx == ori_nidx)? " style=filled fillcolor=yellow" : "");
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		if(auxs[nidx] == 0) continue;
		n = ref_bspoanodev(g->nodes, nidx);
		eidx = n->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
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

static inline int print_node_edges_bspoa(BSPOA *g, u4i nidx, int _dir, FILE *out){
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	u4i eidx, dir, rev, ret;
	v = ref_bspoanodev(g->nodes, nidx);
	nidx = v->header;
	v = ref_bspoanodev(g->nodes, nidx);
	dir = _dir < 2? (1 << _dir) : 3;
	ret = 0;
	for(rev=0;rev<2;rev++){
		if(((dir >> rev) & 0x1) == 0) continue;
		eidx = rev? v->erev : v->edge;
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

static inline void str_cns_ruler_bspoacore(BSPOA *g, u4i mbeg, u4i mend, u4i cbeg, int colorful){
	u4i i, j, b, nseq, mrow;
	char *cp;
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
	UNUSED(colorful);
	clear_string(g->strs); encap_string(g->strs, mend - mbeg);
	cp = g->strs->string;
	j = cbeg;
	for(i=b=mbeg;i<mend;i++){
		if(g->msacols->buffer[i * mrow + nseq] < 4){
			if((j % 10) == 0){
				while(b < i){
					cp[0] = ' ';
					cp ++;
					b  ++;
				}
				if(b + 6 < mend){
					sprintf(cp, "|%05u", j);
					cp += 6;
					b  += 6;
				}
			}
			j ++;
		}
	}
	while(b < mend){ cp[0] = ' '; cp ++; b ++; }
	cp[0] = 0;
	g->strs->size = mend - mbeg;
}

static inline void str_msa_ruler_bspoacore(BSPOA *g, u4i mbeg, u4i mend, int colorful){
	bspoavar_t *x;
	u4i i, j;
	char *cp;
	UNUSED(colorful);
	clear_string(g->strs); encap_string(g->strs, mend - mbeg);
	cp = g->strs->string;
	for(i=j=mbeg;i<mend;i++){
		if((i % 10) == 0 && j + 6 <= mend){
			sprintf(cp, "|%05u", i);
			j  += 6;
			cp += 6;
		} else if(i >= j){
			cp[0] = ' ';
			cp ++;
			j ++;
		}
	}
	cp[0] = 0;
	cp = g->strs->string;
	for(i=0;i<g->var->size;i++){
		x = ref_bspoavarv(g->var, i);
		if(x->mpos >= mend) break;
		if(x->mpos >= mbeg){
			cp[x->mpos - mbeg] = '~';
		}
	}
	g->strs->size = mend - mbeg;
}

static inline u4i str_msa_seq_bspoacore(BSPOA *g, u4i mbeg, u4i mend, u2i rid, u4i rbeg, int colorful){
	bspoanode_t *v;
	bspoavar_t *x;
	u4i i, ci, nseq, mrow, roff;
	u1i *col;
	char ch;
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
	roff = rbeg;
	clear_string(g->strs);
	if(colorful){
		v = g->nodes->size? (rid < nseq? get_rdnode_bspoa(g, rid, roff) : NULL) : NULL;
		ci = 0;
		while(1){
			if(ci < g->var->size){
				x = ref_bspoavarv(g->var, ci);
				if(x->mpos >= mbeg) break;
				ci ++;
			} else {
				x = NULL;
				break;
			}
		}
		for(i=mbeg;i<mend;i++){
			col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
			if(x && x->mpos == i){
				append_string(g->strs, "\e[7m", 4);
				ci ++;
				x = (ci < g->var->size)? ref_bspoavarv(g->var, ci) : NULL;
			} else if(v && v->colorful){
				append_string(g->strs, "\e[1m", 4);
			}
			if(col[rid] <= 4 && col[rid] != col[nseq]){
				append_string(g->strs, "\e[31m", 5);
				ch = "acgt-.*"[col[rid]];
			} else {
				ch = "ACGT-.*"[col[rid]];
			}
			add_char_string(g->strs, ch);
			append_string(g->strs, "\e[0m", 4);
			if(col[rid] < 4){
				roff ++;
				v = g->nodes->size? (rid < nseq? get_rdnode_bspoa(g, rid, roff) : NULL) : NULL;
			}
		}
	} else {
		for(i=mbeg;i<mend;i++){
			col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
			if(col[rid] <= 4 && col[rid] != col[nseq]){
				ch = "acgt-.*"[col[rid]];
			} else {
				ch = "ACGT-.*"[col[rid]];
			}
			add_char_string(g->strs, ch);
			if(col[rid] < 4) roff ++;
		}
	}
	return roff;
}

static inline void str_msa_qlt_bspoacore(BSPOA *g, u4i mbeg, u4i mend, u4i row, int colorful){
	u4i i, j, nseq, mrow;
	u1i *col;
	char ch;
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
	clear_string(g->strs);
	if(colorful){
		for(j=row,i=mbeg;i<mend;i++){
			col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
			ch = '!' + col[j];
			if(colorful){
				if(col[j] < g->par->althi){
					append_string(g->strs, "\e[32m", 5);
				} else if(col[j] < g->par->qlthi){
					append_string(g->strs, "\e[31m", 5);
				}
				add_char_string(g->strs, ch);
				append_string(g->strs, "\e[0m", 4);
			} else {
				add_char_string(g->strs, ch);
			}
		}
	} else {
		for(j=row,i=mbeg;i<mend;i++){
			col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
			ch = '!' + col[j];
			add_char_string(g->strs, ch);
		}
	}
}

/*
# Attribute codes:
# 00=none 01=bold 04=underscore 05=blink 07=reverse 08=concealed
# Text color codes:
# 30=black 31=red 32=green 33=yellow 34=blue 35=magenta 36=cyan 37=white
# Background color codes:
# 40=black 41=red 42=green 43=yellow 44=blue 45=magenta 46=cyan 47=white
*/

static inline void print_msa_bspoa(BSPOA *g, const char *label, u4i mbeg, u4i mend, u4i linewidth, int colorful, FILE *out){
	u4i beg, end, nseq, mrow, i, *roffs, rend, cbeg;
	u1i *col;
	char *str;
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
	if(mend == 0 || mend > g->msaidxs->size) mend = g->msaidxs->size;
	if(linewidth == 0 || linewidth > mend - mbeg) linewidth = mend - mbeg;
	roffs = alloca((nseq + 1) * sizeof(u4i));
	memset(roffs, 0, (nseq + 1) * sizeof(u4i));
	for(beg=0;beg<mbeg;beg++){
		col = g->msacols->buffer + g->msaidxs->buffer[beg] * mrow;
		for(i=0;i<=nseq;i++){
			if(col[i] < 4) roffs[i] ++;
		}
	}
	for(;beg<mend;beg=end){
		end = num_min(mend, beg + linewidth);
		str_msa_ruler_bspoacore(g, beg, end, colorful);
		fprintf(out, "%s MSA [POS] %s\n", label, g->strs->string);
		cbeg = roffs[nseq];
		for(i=0;i<mrow;i++){
			fprintf(out, "%s MSA ", label);
			if(i <= nseq){
				rend = str_msa_seq_bspoacore(g, beg, end, i, roffs[i], colorful);
				if(i == nseq){
					fprintf(out, "[CNS] ");
				} else {
					fprintf(out, "[%03u] ", i);
				}
				fprintf(out, "%s %d\t%d\n", g->strs->string, roffs[i], rend);
				roffs[i] = rend;
			} else {
				str_msa_qlt_bspoacore(g, beg, end, i, colorful);
				if(i == nseq + 1){
					fprintf(out, "[QLT] ");
				} else {
					fprintf(out, "[ALT] ");
				}
				fprintf(out, "%s\n", g->strs->string);
			}
		}
		str_cns_ruler_bspoacore(g, beg, end, cbeg, colorful);
		fprintf(out, "%s MSA [POS] %s\n", label, g->strs->string);
		clear_string(g->strs); encap_string(g->strs, roffs[nseq] - cbeg);
		str = g->strs->string;
		for(i=cbeg;i<roffs[nseq];i++){
			str[i-cbeg] = bit_base_table[g->cns->buffer[i]];
		}
		str[i-cbeg] = 0;
		fprintf(out, "%s CNS\t%d\t%s\n", label, roffs[nseq] - cbeg, str);
		for(i=cbeg;i<roffs[nseq];i++){
			str[i-cbeg] = '!' + g->qlt->buffer[i];
		}
		str[i-cbeg] = 0;
		fprintf(out, "%s QLT\t%d\t%s\n", label, roffs[nseq] - cbeg, str);
		for(i=cbeg;i<roffs[nseq];i++){
			str[i-cbeg] = '!' + g->alt->buffer[i];
		}
		str[i-cbeg] = 0;
		fprintf(out, "%s ALT\t%d\t%s\n", label, roffs[nseq] - cbeg, str);
	}
}

static inline void dump_binary_msa_bspoa(BSPOA *g, char *metadat, u4i metalen, FILE *out){
	u4i nseq, mlen, mrow, i, dlen, dcnt;
	u1i dtag, *col;
	nseq = g->nrds;
	mlen = g->msaidxs->size;
	mrow = g->seqs->nseq + 3;
	if(metadat && metalen){
		dtag = 0x81; // META_DATA 0b10000001
		dlen = metalen;
		fwrite(&dtag, 1, 1, out);
		fwrite(&dlen, 4, 1, out);
		fwrite(metadat, 1, dlen, out);
	}
	dtag = 0x22; // MSA
	dcnt = mlen;
	dlen = nseq;
	fwrite(&dtag, 1, 1, out);
	fwrite(&dcnt, 4, 1, out);
	fwrite(&dlen, 4, 1, out);
	resize_u1v(g->qseq, mlen * 2);
	for(i=0;i<mlen;i++){
		col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
		g->qseq->buffer[i] = col[nseq + 1];
		g->qseq->buffer[mlen + i] = col[nseq + 2];
		fwrite(col, 1, nseq + 1, out); // SEQS + CNS
	}
	fwrite(g->qseq->buffer, 1, mlen * 2, out);
	dtag = 0xFF; // MSA END
	fwrite(&dtag, 1, 1, out);
}

static inline int load_binary_msa_bspoa_core(BSPOA *g, FILE *inp, String *metadat){
	u4i dlen, dcnt, i, nseq, mrow;
	u1i dtag, *col;
	clear_bspoa(g);
	if(metadat) clear_string(metadat);
	while(fread(&dtag, 1, 1, inp) == 1){
		switch(dtag){
			case 0xFF: // MSA END
				return 0;
			case 0x81: // NAME
				if(fread(&dlen, 4, 1, inp) != 1){
					fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
				}
				if(metadat){
					encap_string(metadat, dlen + 1);
					if(fread(metadat->string, 1, dlen, inp) != dlen){
						fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 2;
					}
					metadat->size = dlen;
					metadat->string[dlen] = 0;
					dlen = 0;
				}
				if(dlen && fseek(inp, dlen, SEEK_CUR) == -1){
					fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 3;
				}
				break;
			case 0x22: // MSA
				if(fread(&dcnt, 4, 1, inp) != 1){
					fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 7;
				}
				if(fread(&dlen, 4, 1, inp) != 1){
					fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 8;
				}
				g->nrds = g->nmsa = nseq = dlen;
				mrow = nseq + 3;
				resize_u4v(g->msaidxs, dcnt);
				for(i=0;i<dcnt;i++) g->msaidxs->buffer[i] = i;
				resize_u1v(g->msacols, dcnt * mrow);
				for(i=0;i<dcnt;i++){
					col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
					if(fread(col, 1, nseq + 1, inp) != nseq + 1){
						fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 9;
					}
				}
				clear_seqbank(g->seqs);
				resize_u1v(g->qseq, dcnt * 2);
				if(fread(g->qseq->buffer, 1, dcnt * 2, inp) != dcnt * 2){
					fflush(stdout); fprintf(stderr, " -- read error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 10;
				}
				for(i=0;i<dcnt;i++){
					col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
					col[nseq + 1] = g->qseq->buffer[i];
					col[nseq + 2] = g->qseq->buffer[dcnt + i];
				}
				break;
		}
	}
	return EOF;
}

static inline void post_load_binary_msa_bspoa(BSPOA *g){
	u4i dcnt, i, j, nseq, mrow;
	u1i *col;
	nseq = g->nrds;
	mrow = nseq + 3;
	dcnt = g->msaidxs->size;
	for(i=0;i<nseq;i++){
		push_seqbank(g->seqs, NULL, 0, NULL, 0);
	}
	for(i=0;i<nseq;i++){
		g->seqs->rdoffs->buffer[i] = g->seqs->rdseqs->size;
		for(j=0;j<dcnt;j++){
			col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
			if(col[i] < 4){
				bit2basebank(g->seqs->rdseqs, col[i]);
			}
		}
		g->seqs->rdlens->buffer[i] = g->seqs->rdseqs->size - g->seqs->rdoffs->buffer[i];
	}
	clear_u1v(g->cns);
	clear_u1v(g->qlt);
	clear_u1v(g->alt);
	for(i=0;i<dcnt;i++){
		col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
		if(col[nseq] < 4){
			push_u1v(g->cns, col[nseq]);
			push_u1v(g->qlt, col[nseq + 1]);
			push_u1v(g->alt, col[nseq + 2]);
		}
	}
}


// load_binary_msa_bspoa -> print_msa_bspoa -> denoising_msa_bspoa -> print_snp_bspoa
static inline int load_binary_msa_bspoa(BSPOA *g, FILE *inp, String *metadat){
	int ret;
	ret = load_binary_msa_bspoa_core(g, inp, metadat);
	if(ret == 0) post_load_binary_msa_bspoa(g);
	return ret;
}

static inline void check_aligned_nodes_bspoa(BSPOA *g){
	bspoanode_t *v, *x;
	u4i nidx, xidx;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		v = ref_bspoanodev(g->nodes, nidx);
		xidx = v->header;
		do{
			x = ref_bspoanodev(g->nodes, xidx);
			xidx = x->next;
		} while(xidx != v->header && xidx != nidx);
		if(xidx != nidx){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

static inline u4i print_aligned_nodes_bspoa(BSPOA *g, u4i nidx, FILE *out){
	bspoanode_t *v;
	u4i xidx, ret;
	ret = 0;
	xidx = nidx;
	do {
		v = ref_bspoanodev(g->nodes, xidx);
		fprintf(out, "N%d rid=%u, pos=%d, base=%d, state=%d, vst=%d, nct=%d, nin=%d, edge=%u(N%u)\n", xidx, v->rid, v->pos, v->base, (int)get_bitvec(g->states, nidx), v->vst, v->nct, v->nin, v->edge, v->edge? ref_bspoaedgev(g->edges, v->edge)->node : MAX_U4);
		ret ++;
		xidx = v->next;
	} while(xidx != nidx);
	return ret;
}

static inline u2i count_nodecov_bspoa(BSPOA *g, bspoanode_t *u){
	bspoanode_t *v;
	u2i nodecov;
	nodecov = 0;
	v = u;
	do {
		nodecov ++;
		v = ref_bspoanodev(g->nodes, v->next);
	} while(v != u);
	return nodecov;
}

static inline void check_nodecovs_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i nidx;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		v = ref_bspoanodev(g->nodes, nidx);
		if(get_nodecov_bspoa(g, v) != count_nodecov_bspoa(g, v)){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

static inline bspoanode_t* get_aligned_rdnode_bspoa(BSPOA *g, bspoanode_t *u, u2i rid){
	bspoanode_t *v;
	v = u;
	do {
		if(v->rid == rid) return v;
		v = ref_bspoanodev(g->nodes, v->next);
	} while(v != u);
	return NULL;
}

static inline void check_nodes_aligned_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i step, nidx;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		v = ref_bspoanodev(g->nodes, nidx);
		step = 0;
		if(v->header != nidx) continue;
		while(step <= g->nrds && v->next != nidx){
			v = ref_bspoanodev(g->nodes, v->next);
			step ++;
		}
		if(v->next != nidx){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
}

// keep g->seqs/cns/qlt/alt once for next call of beg_bspoa
static inline void keep_seqs_bspoa(BSPOA *g){
	g->keep_seqs = 1;
}

static inline void beg_bspoa(BSPOA *g){
	g->ncall ++;
	if((g->ncall % 32) == 0){
		renew_bspoa(g);
	} else {
		clear_bspoa(g);
	}
	if(!g->par->refmode){
		fwdbitseqpush_bspoa(g, NULL, 0);
	}
}

static inline b1i* dpalign_row_prepare_data(BSPOA *g, u4i mmidx, b1i **us, b1i **es, b1i **qs, int **ubegs){
	us[0] = g->memp->buffer + mmidx * g->mmblk;
	es[0] = (g->piecewise > 0)? us[0] + g->bandwidth : NULL;
	qs[0] = (g->piecewise > 1)? es[0] + g->bandwidth : NULL;
	ubegs[0] = (int*)(us[0] + g->bandwidth * (g->piecewise + 1));
	return ((b1i*)ubegs[0]) + (WORDSIZE + 1) * sizeof(int);
}

static inline void add_rdnodes_bspoa(BSPOA *g, u2i rid){
	bspoanode_t *u, *v;
	u4i i, seqoff, seqlen;
	fflush(stdout); fprintf(stderr, " -- TODO: revise me in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	if(g->nrds != rid){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	g->nrds ++;
	//zero2bitvec(g->rdbits);
	u = ref_bspoanodev(g->nodes, g->HEAD);
	v = new_node_bspoa(g, rid, -1, 4);
	merge_nodes_bspoa(g, u, v);
	push_u4v(g->ndoffs, g->nodes->size);
	for(i=0;i<seqlen;i++){
		v = new_node_bspoa(g, rid, i, get_basebank(g->seqs->rdseqs, seqoff + i));
	}
	u = ref_bspoanodev(g->nodes, g->TAIL);
	v = new_node_bspoa(g, rid, seqlen, 4);
	merge_nodes_bspoa(g, u, v);
	push_b4v(g->rdregs[0], 0);
	push_b4v(g->rdregs[1], seqlen);
	encap_bitvec(g->states, g->nodes->size);
}

static inline void check_rdnodes_bspoa(BSPOA *g){
	u4i rid, pos, rdlen;
	for(rid=0;rid<g->nrds;rid++){
		rdlen = g->seqs->rdlens->buffer[rid];
		if(rdlen == 0) continue;
		for(pos=0;pos<=rdlen;pos++){
			if(isconnect_rdnode_bspoa(g, rid, pos) == 0){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
}

static inline void check_unvisited_bspoacore(BSPOA *g, u4i nhead, u4i ntail, BitVec *states, int dir){
	BitVec *ndbits;
	u4v *queue;
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	u4i nidx, eidx, has;
	has = 0;
	ndbits = init_bitvec(g->nodes->size);
	queue = init_u4v(32);
	if(dir){
		push_u4v(queue, ntail);
		one_bitvec(ndbits, ntail);
	} else {
		push_u4v(queue, nhead);
		one_bitvec(ndbits, nhead);
	}
	while(queue->size){
		nidx = queue->buffer[0];
		remove_u4v(queue, 0);
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = dir? u->erev : u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(states && get_bitvec(states, e->node) == 0) continue;
			if(get_bitvec(ndbits, e->node) == 1) continue;
			v = ref_bspoanodev(g->nodes, e->node);
			if(v->vst < v->nct){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				print_node_edges_bspoa(g, e->node, !dir, stderr);
				has ++;
				one_bitvec(ndbits, e->node);
			}
			if(get_bitvec(ndbits, e->node) == 0){
				one_bitvec(ndbits, e->node);
				push_u4v(queue, e->node);
			}
		}
	}
	free_bitvec(ndbits);
	free_u4v(queue);
	if(has){
		fflush(stdout); fprintf(stderr, " -- total errors: %d in %s -- %s:%d --\n", has, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
}

static inline void check_unvisited_bspoa(BSPOA *g, int dir){
	check_unvisited_bspoacore(g, g->HEAD, g->TAIL, NULL, dir);
}

static inline u4i sel_nodes_bspoa(BSPOA *g, u4i nhead, u4i ntail, u2i ridxbeg, u2i ridxend){
	bspoanode_t *u, *v, *x;
	bspoaedge_t *e;
	u4i nidx, eidx, i, nseq;
	int j, rb, re, bonus;
	u = ref_bspoanodev(g->nodes, nhead);
	nhead = u->header;
	v = ref_bspoanodev(g->nodes, ntail);
	ntail = v->header;
	nseq = g->nmsa;
	for(i=0;i<nseq;i++){
		g->rdregs[0]->buffer[i] = MAX_B4;
		g->rdregs[1]->buffer[i] = -1;
	}
	clear_u4v(g->sels);
	reg_zeros_bitvec(g->states, 0, g->nodes->size);
	if(nhead == ntail) return 0;
	u = ref_bspoanodev(g->nodes, nhead);
	v = ref_bspoanodev(g->nodes, ntail);
	for(i=0;i<2;i++){
		nidx = i? ntail : nhead;
		u = ref_bspoanodev(g->nodes, nidx);
		x = u;
		do {
			if(x->rid >= ridxbeg && x->rid < ridxend){
				g->rdregs[i]->buffer[x->rid] = x->pos;
			}
			x = ref_bspoanodev(g->nodes, x->next);
		} while(x != u);
	}
	for(i=0;i<nseq;i++){
		rb = g->rdregs[0]->buffer[i];
		re = g->rdregs[1]->buffer[i];
		if(rb >= re) continue;
		for(j=rb;j<=re;j++){
			u = get_rdnode_bspoa(g, i, j);
			if(get_bitvec(g->states, u->header)) continue;
			push_u4v(g->sels, u->header);
			one_bitvec(g->states, u->header);
			v = header_node_bspoa(g, u);
			v->nct   = 0;
			v->vst   = 0;
		}
	}
	for(i=0;i<g->sels->size;i++){
		nidx = g->sels->buffer[i];
		if(nidx == nhead) continue;
		u = ref_bspoanodev(g->nodes, nidx);
		j = 0;
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(get_bitvec(g->states, e->node) == 0) continue;
			j |= 1;
			break;
		}
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(get_bitvec(g->states, e->node) == 0) continue;
			j |= 2;
			break;
		}
		if(j == 3){
		} else if(j == 1 || nidx == ntail){
			// add aux link to make sure all selected node can be visited from nhead
			chg_edge_bspoa(g, ref_bspoanodev(g->nodes, nhead), u, 1, NULL);
			push_u8v(g->todels, (((u8i)nhead) << 32) | nidx);
		} else if(j == 2){
			// add aux link to make sure all selected node can be visited from ntail
			chg_edge_bspoa(g, u, ref_bspoanodev(g->nodes, ntail), 1, NULL);
			push_u8v(g->todels, (((u8i)nidx) << 32) | ntail);
		}
	}
	for(i=0;i<g->sels->size;i++){
		nidx = g->sels->buffer[i];
		u = ref_bspoanodev(g->nodes, nidx);
		if(1){
			bonus = 0;
			x = u;
			do {
				bonus |= x->bless;
				x = ref_bspoanodev(g->nodes, x->next);
			} while(x != u && !bonus);
			u->bonus = bonus;
		}
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(get_bitvec(g->states, e->node) == 0) continue;
			v = ref_bspoanodev(g->nodes, e->node);
			v->nct ++;
		}
	}
#if DEBUG
	set_edges_unvisited_bspoa(g);
	clear_u4v(g->stack);
	push_u4v(g->stack, nhead);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			if(get_bitvec(g->states, e->node) == 0) continue;
			e->vst = 1;
			(e + 1)->vst = 1;
			v->vst ++;
			if(v->vst == v->nct){
				push_u4v(g->stack, e->node);
			} else if(v->vst > v->nct){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	v = ref_bspoanodev(g->nodes, ntail);
	if(v->nct == 0 || v->vst != v->nct){
		check_unvisited_bspoacore(g, nhead, ntail, g->states, 0);
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	for(i=0;i<g->sels->size;i++){
		nidx = g->sels->buffer[i];
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst = 0;
	}
#endif
	return g->sels->size;
}

static inline void prepare_rd_align_bspoa(BSPOA *g, BSPOAPar *par, u4i nhead, u4i ntail, u2i rid, int rbeg, int rend){
	seqalign_result_t rs;
	bspoanode_t *u, *v;
	u4i seqoff, reflen, *rmap, *cgs, i, ncg, x, y, j, sz, op;
	int seqlen, rpos, tb, te, *ubegs, exists;
	b1i *us, *es, *qs;
	u = ref_bspoanodev(g->nodes, nhead);
	nhead = u->header;
	v = ref_bspoanodev(g->nodes, ntail);
	ntail = v->header;
	seqoff = g->seqs->rdoffs->buffer[rid] + rbeg;
	seqlen = rend - rbeg;
	g->qlen = g->slen = seqlen;
	g->qb = 0; g->qe = g->qlen;
	clear_and_encap_u1v(g->qseq, g->qlen);
	bitseq_basebank(g->seqs->rdseqs, seqoff, g->qlen, g->qseq->buffer);
	g->qseq->size = g->qlen;
	if(g->par->refmode){
		reflen = g->backbone;
	} else {
		reflen = g->cns->size;
	}
	tb = 0; te = reflen;
	cgs = NULL; ncg = 0;
	if(par->bandwidth == 0){
		g->bandwidth = roundup_times(seqlen, WORDSIZE);
	} else {
		g->bandwidth = num_min(par->bandwidth, Int(seqlen));
		g->bandwidth = roundup_times(g->bandwidth, WORDSIZE);
	}
	if(par->bwtrigger && nhead == g->HEAD && ntail == g->TAIL){ // global
		if(g->par->refmode && g->cges->buffer[rid] > g->cgbs->buffer[rid]){
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
			x = 0; y = tb;
			tb = (tb >= Int(g->bandwidth / 2))? tb - g->bandwidth / 4 : 0;
			te = (reflen - te >= g->bandwidth / 2)? te + g->bandwidth / 4 : (reflen);
		} else if(g->cns->size && Int(roundup_times(seqlen, WORDSIZE)) > par->bandwidth){
			if(par->ksz){
				rs = kmer_striped_seqedit_pairwise(par->ksz, g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->memp, g->stack, 0);
			} else {
				rs = striped_seqedit_pairwise(g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, par->alnmode, 0, g->memp, g->stack, 0);
			}
			if(_DEBUG_LOG_){
				char *alnstr[3];
				alnstr[0] = malloc(rs.aln + 1); alnstr[1] = malloc(rs.aln + 1); alnstr[2] = malloc(rs.aln + 1);
				seqalign_cigar2alnstr(g->qseq->buffer, g->cns->buffer, &rs, g->stack, alnstr, rs.aln);
				fprintf(stderr, "#RID%d\t%d\t%d\t%d\tCNS\t%d\t%d\t%d\tmat=%d\taln=%d\n", rid, Int(g->qseq->size), rs.qb, rs.qe, Int(g->cns->size), rs.tb, rs.te, rs.mat, rs.aln);
				fprintf(stderr, "#%s\n#%s\n#%s\n", alnstr[0], alnstr[2], alnstr[1]);
				free(alnstr[0]); free(alnstr[1]); free(alnstr[2]);
			}
			g->qb = rs.qb;
			g->qe = rs.qe;
			g->slen = g->qe - g->qb;
			tb = (rs.tb >= Int(g->bandwidth / 2))? rs.tb - g->bandwidth / 4 : 0;
			te = (g->cns->size - rs.te >= g->bandwidth / 2)? Int(rs.te + g->bandwidth / 4) : Int(g->cns->size);
			cgs = g->stack->buffer;
			ncg = g->stack->size;
			x = 0; y = rs.tb;
		} else {
			g->bandwidth = roundup_times(seqlen, WORDSIZE);
		}
	} else {
		g->bandwidth = roundup_times(seqlen, WORDSIZE);
	}
	if(cgs && ncg){
		clear_and_encap_b1v(g->memp, (reflen + 1) * sizeof(u4i));
		rmap = (u4i*)(g->memp->buffer);
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
		for(i=y;i<reflen;i++){
			rmap[i] = x + (i - y + 1) * (g->slen - x) / (reflen - y + 1);
		}
#if DEBUG && 0
		for(i=1;i<reflen;i++){
			if(rmap[i] < rmap[i-1]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
#endif
		rmap[reflen] = g->slen;
		if(_DEBUG_LOG_ > 1){
			fprintf(stderr, "RMAP:");
			for(i=0;i<=reflen;i++){
				fprintf(stderr, "\t[%d]%d", i, rmap[i]);
			}
			fprintf(stderr, "\n");
		}
		// set rpos
		for(i=0;i<g->sels->size;i++){
			u = ref_bspoanodev(g->nodes, g->sels->buffer[i]);
			rpos = rmap[u->cpos] - g->bandwidth / 2;
			if(rpos < 0) rpos = 0;
			else if(g->bandwidth >= g->slen) rpos = 0;
			else if(rpos + g->bandwidth > g->slen) rpos = g->slen - g->bandwidth;
			u->rpos = rpos;
			if(u->cpos == tb && tb){
				chg_edge_bspoa(g, ref_bspoanodev(g->nodes, nhead), u, 1, &exists);
				push_u8v(g->todels, (((u8i)nhead) << 32) | offset_bspoanodev(g->nodes, u));
				tb = 0; // used once
				if(exists == 0 && get_bitvec(g->states, nhead) && get_bitvec(g->states, g->sels->buffer[i])){
					u->nct ++;
				}
			}
			if(u->cpos == te && te != Int(reflen)){
				chg_edge_bspoa(g, u, ref_bspoanodev(g->nodes, ntail), 1, &exists);
				push_u8v(g->todels, (((u8i)u->header) << 32) | ntail);
				te = reflen; // used once
				if(exists == 0 && get_bitvec(g->states, ntail) && get_bitvec(g->states, g->sels->buffer[i])){
					ref_bspoanodev(g->nodes, ntail)->nct ++;
				}
			}
		}
	} else {
		for(i=0;i<g->sels->size;i++){
			u = ref_bspoanodev(g->nodes, g->sels->buffer[i]);
			u->rpos = 0;
		}
	}
	// no bonus, no homopolymer
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[0], par->M, par->X);
	clear_and_encap_b1v(g->qprof[0], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof_hpc(g->qseq->buffer + g->qb, g->slen, g->qprof[0]->buffer, g->bandwidth, g->matrix[0], 1);
	//banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[0]->buffer, g->bandwidth, g->matrix[0]);
	// bonus, no homopolymer
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[1], par->M + par->refbonus, par->X);
	clear_and_encap_b1v(g->qprof[1], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof_hpc(g->qseq->buffer + g->qb, g->slen, g->qprof[1]->buffer, g->bandwidth, g->matrix[1], 1);
	//banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[1]->buffer, g->bandwidth, g->matrix[1]);
	// no bonus, homopolymer
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[2], par->M, par->X);
	clear_and_encap_b1v(g->qprof[2], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[2]->buffer, g->bandwidth, g->matrix[2]);
	// bonus, homopolymer
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix[3], par->M + par->refbonus, par->X);
	clear_and_encap_b1v(g->qprof[3], banded_striped_epi8_seqalign_qprof_size(g->slen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer + g->qb, g->slen, g->qprof[3]->buffer, g->bandwidth, g->matrix[3]);
	g->piecewise = banded_striped_epi8_seqalign_get_piecewise(par->O, par->E, par->Q, par->P, g->bandwidth);
	g->mmblk = roundup_times(g->bandwidth * (g->piecewise + 1) + (WORDSIZE + 1) * sizeof(int), WORDSIZE); // us, es, qs, ubegs, qoff, max_score, max_index
	g->mmcnt = 2; // reserve two blocks
	for(i=0;i<g->sels->size;i++){
		u = ref_bspoanodev(g->nodes, g->sels->buffer[i]);
		u->mmidx = g->mmcnt ++;
	}
	clear_and_encap_b1v(g->memp, g->mmcnt * g->mmblk);
	u = ref_bspoanodev(g->nodes, nhead);
	dpalign_row_prepare_data(g, u->mmidx, &us, &es, &qs, &ubegs);
	banded_striped_epi8_seqalign_piecex_row_init(us, es, qs, ubegs, NULL, par->alnmode, g->bandwidth, par->M + par->refbonus + 1, par->X, par->O, par->E, par->Q, par->P);
	g->maxscr = SEQALIGN_SCORE_MIN;
	g->maxidx = -1;
	g->maxoff = -1;
}

static inline void dpalign_row_update_bspoa(BSPOA *g, BSPOAPar *par, b1i *qp, u4i mmidx1, u4i mmidx2, u4i toff, u4i qoff1, u4i qoff2, u1i next_base){
	b1i *us[2], *es[2], *qs[2];
	int W, rh, *ubegs[2];
	W = g->bandwidth / WORDSIZE;
	dpalign_row_prepare_data(g,      0, us + 0, es + 0, qs + 0, ubegs + 0);
	__builtin_prefetch(us[0], 0);
	dpalign_row_prepare_data(g, mmidx1, us + 1, es + 1, qs + 1, ubegs + 1);
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_movx(us + 0, es + 0, qs + 0, ubegs + 0, W, qoff2 - qoff1, g->piecewise, par->M + par->refbonus + 1, par->X, par->O, par->E, par->Q, par->P);
	dpalign_row_prepare_data(g, mmidx2, us + 1, es + 1, qs + 1, ubegs + 1);
	if(qoff1 == qoff2){
		if(qoff1){
			rh = SEQALIGN_SCORE_MIN;
		} else {
			if(seqalign_mode_type(par->alnmode) == SEQALIGN_MODE_OVERLAP || toff == 0) rh = 0;
			else if(g->piecewise < 2) rh = par->O + par->E * toff;
			else rh = num_max(par->O + par->E * toff, par->Q + par->P * toff);
		}
	} else if(qoff1 + W * WORDSIZE >= qoff2){
		rh = ubegs[0][0]; // movx -> aligned
	} else {
		rh = SEQALIGN_SCORE_MIN;
	}
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_cal(qoff2, next_base, us, es, qs, ubegs, qp, par->O, par->E, par->Q, par->P, W, qoff2 - qoff1, rh, g->piecewise);
#if 0
	dpalign_row_prepare_data(g, mmidx1, us + 0, es + 0, qs + 0, ubegs + 0);
	banded_striped_epi8_seqalign_piecex_row_verify(toff, qoff1, par->alnmode, W, qoff2 - qoff1, next_base, us, es, qs, ubegs, qp, par->O, par->E, par->Q, par->P);
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

static inline seqalign_result_t alignment2graph_bspoa(BSPOA *g, BSPOAPar *par, u4i rid, u4i rbeg, u4i nhead, u4i ntail, u4i midx, int xe, String **alnstrs){
	seqalign_result_t rs;
	bspoanode_t *u, *v, *w, *n;
	bspoaedge_t *e;
	b1i *us, *es, *qs, q;
	u4i nidx, eidx, i, W, bt, ft, btc, bti;
	int x, *ubegs, Hs[3], *scrs, t, s, realn, cpos;
	nhead = ref_bspoanodev(g->nodes, nhead)->header;
	ntail = ref_bspoanodev(g->nodes, ntail)->header;
	realn = 0;
	u = get_rdnode_bspoa(g, rid, 0);
	v = get_rdnode_bspoa(g, rid, g->qlen);
	while(u < v){ u->cpos = 0; u ++; }
	W = g->bandwidth / WORDSIZE;
	if(alnstrs){
		clear_string(alnstrs[0]);
		clear_string(alnstrs[1]);
		clear_string(alnstrs[2]);
	}
	v = ref_bspoanodev(g->nodes, ntail);
	ZEROS(&rs);
	rs.qe = xe + 1;
	rs.qb = x = xe;
	nidx = midx;
	bt = MAX_U4; // to be calculated
	n = ref_bspoanodev(g->nodes, nidx);
	rs.te = n->cpos + 1;
	cpos = n->cpos;
	dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
	Hs[0] = 0;
	Hs[1] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - n->rpos);
	Hs[2] = 0;
	while(1){
		if(n->header == nhead || x < 0){
			rs.qb = x;
			rs.tb = n->cpos;
			break;
		}
		if(bt == SEQALIGN_BT_D || bt == SEQALIGN_BT2_D2){
			rs.del ++;
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
				if(get_bitvec(g->states, e->node) == 0) continue;
				w = ref_bspoanodev(g->nodes, e->node);
				if(x < Int(w->rpos) || x >= Int(w->rpos + g->bandwidth)) continue;
				dpalign_row_prepare_data(g, w->mmidx, &us, &es, &qs, &ubegs);
				Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos);
				if(bt == SEQALIGN_BT_D){
					q = g->piecewise? es[banded_striped_epi8_pos2idx(g->bandwidth, x - w->rpos)] : par->O + par->E;
				} else {
					q = qs[banded_striped_epi8_pos2idx(g->bandwidth, x - w->rpos)];
				}
				if(Hs[0] + q != Hs[1]) continue;
				n = ref_bspoanodev(g->nodes, e->node);
				if(q == ((bt == SEQALIGN_BT_D)? par->O + par->E : par->Q + par->P)){
					bt = MAX_U4;
					Hs[1] = Hs[0];
					Hs[2] = 0;
				} else {
					Hs[1] -= ((bt == SEQALIGN_BT_D)? par->E : par->P);
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
			rs.ins ++;
			if(g->piecewise == 2){
				t = num_max(par->O + par->E * Hs[2], par->Q + par->P * Hs[2]);
			} else {
				t = par->O + par->E * Hs[2];
			}
#if DEBUG
			if(x < Int(n->rpos)){ // should never happen
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			u = get_rdnode_bspoa(g, rid, rbeg + g->qb + x);
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
			u = get_rdnode_bspoa(g, rid, rbeg + g->qb + x);
			u->cpos = n->cpos;
			if(alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], offset_bspoanodev(g->nodes, n)? bit_base_table[n->base&0x03] : '^');
				add_char_string(alnstrs[2], "*|"[(u->base&0x03) == (n->base&0x03)]);
			}
			x --;
			if(offset_bspoanodev(g->nodes, n) != nhead && offset_bspoanodev(g->nodes, n) != ntail && u->base == n->base){
				merge_nodes_bspoa(g, n, u);
				rs.mat ++;
			} else {
				rs.mis ++;
			}
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
			eidx = n->erev;
			clear_u4v(g->stack);
			btc = 0;
			bti = MAX_U4;
			while(eidx){
				e = ref_bspoaedgev(g->edges, eidx);
				eidx = e->next;
				if(get_bitvec(g->states, e->node) == 0) continue;
				w = ref_bspoanodev(g->nodes, e->node);
				dpalign_row_prepare_data(g, w->mmidx, &us, &es, &qs, &ubegs);
				ft = 0;
				if(x < Int(w->rpos) || x > Int(g->bandwidth + w->rpos)){
					continue;
				} else if(x == Int(g->bandwidth + w->rpos)){
					Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos - 1);
					ft |= 1 << SEQALIGN_BT_D;
					ft |= 1 << SEQALIGN_BT2_D2;
				} else if(x == Int(w->rpos)){
					if(w->rpos == 0 && (seqalign_mode_type(par->alnmode) == SEQALIGN_MODE_OVERLAP || offset_bspoanodev(g->nodes, w) == nhead)){
						Hs[0] = ubegs[0];
						ft |= 1 << 15;
						//Hs[0] = 0;
					} else {
						Hs[0] = ubegs[0]; // useless
						ft |= 1 << SEQALIGN_BT_M; // forbid M
					}
				} else {
					Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos - 1);
				}
				//s = g->matrix[(w->base == n->base) * 2 + n->bonus][g->qseq->buffer[x + g->qb] * 4 + n->base];
				t = (x * 4 + n->base) * WORDSIZE;
				s = g->qprof[(w->base == n->base) * 2 + n->bonus]->buffer[t];
				if(ft & (1 << 15)) s -= ubegs[0];
				t = ((x - w->rpos) % W) * WORDSIZE + ((x - w->rpos) / W);
				push_u4v(g->stack, e->node);
				push_u4v(g->stack, UInt(Hs[0]));
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT_M))? UInt(SEQALIGN_SCORE_MIN) : UInt(s));
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT_D))? UInt(SEQALIGN_SCORE_MIN) : UInt(us[t] + (es? es[t] : par->E)));
				push_u4v(g->stack, (ft & (1 << SEQALIGN_BT2_D2))? UInt(SEQALIGN_SCORE_MIN) : UInt((qs? us[t] + qs[t] : SEQALIGN_SCORE_MAX)));
				scrs = (int*)(g->stack->buffer + g->stack->size - 3);
				for(i=0;i<3;i++){
					if(Hs[0] + scrs[i] == Hs[1]){
						if(e->cov > btc){
							bti = (g->stack->size << 8) | i;
							btc = e->cov;
						} else if(e->cov == btc && i == 0 && (bti & 0xFF) != 0){
							bti = (g->stack->size << 8) | i;
							btc = e->cov;
						}
					}
				}
			}
			if(bti == MAX_U4){
				bt = SEQALIGN_BT_I;
				Hs[2] = 1;
				dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
				Hs[0] = Hs[1] - us[banded_striped_epi8_pos2idx(g->bandwidth, x - n->rpos)];
			} else {
				if((bti & 0xFF) == 0){
					bt = SEQALIGN_BT_M;
					nidx = g->stack->buffer[(bti >> 8) - 5];
					Hs[1] = Int(g->stack->buffer[(bti >> 8) - 4]);
					Hs[2] = 0;
				} else if((bti & 0xFF) == 1){
					bt = SEQALIGN_BT_D;
					Hs[2] = 1;
				} else {
					bt = SEQALIGN_BT2_D2;
					Hs[2] = 1;
				}
			}
		}
	}
	rs.qb += g->qb;
	rs.qe += g->qb;
	v = get_rdnode_bspoa(g, rid, rbeg + g->qlen);
	connect_rdnode_bspoa(g, rid, rbeg + g->qlen);
	for(x=g->qlen-1;x>=0;x--){
		connect_rdnode_bspoa(g, rid, rbeg + x);
		v --;
		if(v->cpos){
			cpos = v->cpos;
		} else if(x < Int(g->qlen)){
			v->cpos = cpos;
		}
	}
	return rs;
}

static inline int align_rd_bspoacore(BSPOA *g, BSPOAPar *par, u2i rid, u4i nhead, u4i ntail){
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	u4i i, nidx, eidx, mmidx;
	b1i *us, *es, *qs;
	int *ubegs, smax, rmax, maxoff;
	UNUSED(rid);
	for(i=0;i<g->sels->size;i++){
		nidx = g->sels->buffer[i];
		u = ref_bspoanodev(g->nodes, nidx);
		u->mpos = MAX_B4 - 1;
	}
#if DEBUG
	set_edges_unvisited_bspoa(g);
#endif
	clear_u4v(g->stack);
	u = ref_bspoanodev(g->nodes, nhead);
	u->mpos = -1;
	push_u4v(g->stack, nhead);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			if(get_bitvec(g->states, e->node) == 0) continue;
#if DEBUG
			e->vst = 1;
			(e + 1)->vst = 1;
#endif
			v = ref_bspoanodev(g->nodes, e->node);
			if(u->mpos + 1 < v->mpos){
				v->mpos = u->mpos + 1;
			}
			if(e->node == ntail){
				dpalign_row_prepare_data(g, u->mmidx, &us, &es, &qs, &ubegs);
				// If reads are not high similar, the final score of graph alignment will be negative
				// but the HEAD->TAIL score is 0, we will always put new read on the graph without alignment
				// so, I force add par->E * (unaligned length) to avoid the above disaster
				{
					maxoff = num_min(g->slen, u->rpos + g->bandwidth) - 1;
					smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, maxoff - u->rpos);
					if(Int(g->slen) > maxoff + 1){
						if(g->piecewise < 2){
							smax += par->O + par->E * (g->slen - maxoff - 1);
						} else {
							smax += num_max(par->O + par->E * (g->slen - maxoff - 1), par->Q + par->P * (g->slen - maxoff - 1));
						}
					}
					smax += par->T;
					if(smax > g->maxscr){
						g->maxscr = smax;
						g->maxidx = offset_bspoanodev(g->nodes, u);
						g->maxoff = maxoff;
					}
				}
				if(seqalign_mode_type(par->alnmode) == SEQALIGN_MODE_OVERLAP){
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
				dpalign_row_update_bspoa(g, par, g->qprof[(v->base == u->base) * 2 + v->bonus]->buffer, u->mmidx, mmidx, v->mpos, u->rpos, v->rpos, v->base);
				if(v->vst){
					dpalign_row_merge_bspoa(g, mmidx, v->mmidx);
				}
				v->vst ++;
				if(v->vst == v->nct){
					if(seqalign_mode_type(par->alnmode) != SEQALIGN_MODE_GLOBAL && v->rpos + g->bandwidth >= g->slen){
						dpalign_row_prepare_data(g, v->mmidx, &us, &es, &qs, &ubegs);
						smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, g->slen - 1 - v->rpos);
						smax += par->T;
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
	v = ref_bspoanodev(g->nodes, ntail);
#if DEBUG
	if(v->vst != v->nct){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		check_unvisited_bspoacore(g, nhead, ntail, g->states, 0);
		abort();
	}
#endif
	return g->maxscr;
}

static inline seqalign_result_t align_rd_bspoa(BSPOA *g, BSPOAPar *par, int realn, u2i rid, int rbeg, int rlen){
	seqalign_result_t rs;
	String **alnstrs;
	u4i i, nhead, ntail;
	u2i ridxbeg, ridxend;
	int score;
	if(realn && rid){
		for(i=rbeg;i<UInt(rbeg+rlen);i++){
			cut_rdnode_bspoa(g, rid, i, BSPOA_RDNODE_CUTALL);
		}
	}
	clear_u8v(g->todels);
	ZEROS(&rs);
	if(rlen == 0) return rs;
	nhead = get_rdnode_bspoa(g, rid, rbeg - 1)->header;
	ntail = get_rdnode_bspoa(g, rid, rbeg + rlen)->header;
	if(realn == 0 && par->nrec){
		ridxbeg = num_max(0, Int(rid) - par->nrec - 1);
		ridxend = rid;
	} else {
		ridxbeg = 0;
		ridxend = MAX_U2;
	}
	sel_nodes_bspoa(g, nhead, ntail, ridxbeg, ridxend);
	prepare_rd_align_bspoa(g, par, nhead, ntail, rid, rbeg, rbeg + rlen);
	score = align_rd_bspoacore(g, par, rid, nhead, ntail);
	alnstrs = NULL;
	if(_DEBUG_LOG_){
		alnstrs = malloc(3 * sizeof(String*));
		alnstrs[0] = init_string(rlen * 2); alnstrs[1] = init_string(rlen * 2); alnstrs[2] = init_string(rlen * 2);
	}
	rs = alignment2graph_bspoa(g, par, rid, rbeg, nhead, ntail, g->maxidx, g->maxoff, alnstrs);
	rs.qb += g->qb;
	rs.qe += g->qb;
	rs.score = score;
	for(i=0;i<g->todels->size;i++){
		chg_edge_bspoa(g, ref_bspoanodev(g->nodes, g->todels->buffer[i] >> 32), ref_bspoanodev(g->nodes, g->todels->buffer[i] & MAX_U4), -1, NULL);
	}
	clear_u8v(g->todels);
	if(alnstrs){
		reverse_string(alnstrs[0]); reverse_string(alnstrs[1]); reverse_string(alnstrs[2]);
		fprintf(stderr, "ALIGN[%03d] len=%u band=%d aligned=%d,%d mat=%d,%0.3f tail=%d score=%d\n", rid, rlen, g->bandwidth, rs.qb + 1, rs.qe + 1, rs.mat, 1.0 * rs.mat / rlen, rs.qb + g->qlen - rs.qe, rs.score);
		fprintf(stderr, "%s\n%s\n%s\n", alnstrs[0]->string, alnstrs[2]->string, alnstrs[1]->string);
		fflush(stderr);
		free_string(alnstrs[0]); free_string(alnstrs[1]); free_string(alnstrs[2]); free(alnstrs);
	}
	return rs;
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
	u4i nseq, mrow, nidx, eidx, xidx, ready;
	int moff;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst   = 0;
		u->nct   = u->nou;
		u->inuse = 0;
		u->mpos  = 0;
	}
#if DEBUG
	set_edges_unvisited_bspoa(g);
#endif
	clear_u4v(g->stack);
	nidx = g->TAIL;
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
			if(v->vst > v->nct){
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
			if(v->inuse){
				continue; // already pushed
			}
			if(v->vst > v->nct){
				print_vstdot_bspoa(g, "1.dot");
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				check_dup_edges_bspoa(g);
				abort();
			}
			if(v->vst == v->nct){
				ready = 1;
				{
					xidx = v->next;
					moff = v->mpos;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						if(x->nct > x->vst){
							ready = 0;
							break;
						}
						if(x->mpos > moff){
							moff = x->mpos;
						}
						xidx = x->next;
					}
				}
				if(ready){
					v->mpos  = moff;
					v->inuse = 1;
					push_u4v(g->stack, e->node);
					xidx = v->next;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						x->mpos = moff;
						if(x->edge){
							push_u4v(g->stack, xidx);
							x->inuse = 1;
						}
						xidx = x->next;
					}
				}
			}
		}
	}
	if(nidx != g->HEAD){
		fprint_dot_bspoa(g, 0, MAX_U4, 1, "1.dot", NULL);
		print_seqs_bspoa(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		check_unvisited_bspoa(g, 1);
		abort();
	}
#if 0
	// tidy toword HEAD
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst   = 0;
		u->nct   = u->nin;
		u->inuse = 0;
	}
#if DEBUG
	set_edges_unvisited_bspoa(g);
#endif
	clear_u4v(g->stack);
	nidx = g->HEAD;
	push_u4v(g->stack, nidx);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		moff = MAX_B4;
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
#if DEBUG
			if(e->vst == 0){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			eidx = e->next;
			if(e->node == g->HEAD) continue;
			v = ref_bspoanodev(g->nodes, e->node);
			if(v->mpos - 1 > moff){
				moff = v->mpos - 1;
			}
		}
		if(moff != MAX_B4 && moff > u->mpos){
			u->mpos = moff;
			xidx = u->next;
			while(xidx != nidx){
				x = ref_bspoanodev(g->nodes, xidx);
				x->mpos = moff;
				xidx = x->next;
			}
		}
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
#if DEBUG
			e->vst = 1;
			(e+1)->vst = 1;
#endif
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			v->vst ++;
#if DEBUG
			if(v->vst > v->nct){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			if(v->vst == v->nct){
				v->inuse = 1;
				push_u4v(g->stack, e->node);
			}
		}
	}
	if(nidx != g->TAIL){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#else
	u = ref_bspoanodev(g->nodes, g->TAIL);
	clear_u4v(g->stack);
	eidx = u->erev;
	while(eidx){
		e = ref_bspoaedgev(g->edges, eidx);
		eidx = e->next;
		if(e->node == g->HEAD) continue;
		x = u;
		v = ref_bspoanodev(g->nodes, e->node);
		while(1){
			u4i nou = 0;
			xidx = v->edge;
			while(xidx){
				if(g->edges->buffer[xidx].node != offset_bspoanodev(g->nodes, x) && g->edges->buffer[xidx].node != g->TAIL) nou ++;
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
			xidx = v->next;
#if DEBUG
			if(v->mpos > moff){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			do {
				x = ref_bspoanodev(g->nodes, xidx);
				x->mpos = moff;
				xidx = x->next;
			} while(x != v);
			moff --;
			xidx = v->edge;
			v = NULL;
			while(xidx){
				if(g->edges->buffer[xidx].node != g->TAIL){
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
#endif
#if DEBUG
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
	v = ref_bspoanodev(g->nodes, g->HEAD);
	clear_u4v(g->msaidxs);
	for(nidx=0;Int(nidx)<v->mpos;nidx++){
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

static inline void check_msa_rdseqs_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i rid, roff, rlen, nseq, mrow, mlen, pos;
	u1i *col;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	for(rid=0;rid<nseq;rid++){
		rlen = g->seqs->rdlens->buffer[rid];
		for(roff=0;roff<rlen;roff++){
			v = get_rdnode_bspoa(g, rid, roff);
			if(v->base != get_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[rid] + roff)){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
		for(pos=roff=0;pos<mlen;pos++){
			col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(col[rid] >= 4) continue;
			if(roff >= rlen){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			if(col[rid] != get_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[rid] + roff)){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			roff ++;
		}
	}
}

static inline void del_cnsnodes_bspoa(BSPOA *g){
	int nseq, clen, i;
	nseq = g->seqs->nseq;
	if(Int(g->ndoffs->size) <= nseq){
		return;
	}
	clen = g->nodes->size - g->ndoffs->buffer[nseq] - 1;
	for(i=-1;i<=clen;i++){
		cut_rdnode_bspoa(g, nseq, i, BSPOA_RDNODE_CUTALL);
	}
	g->nodes->size = g->ndoffs->buffer[nseq] - 1;
	g->ndoffs->size = nseq;
	g->rdregs[0]->size = nseq;
	g->rdregs[1]->size = nseq;
}

static inline void add_cnsnodes_bspoa(BSPOA *g, u4i *rps){
	bspoanode_t *u, *v;
	u4i nseq, mrow, mlen, clen, pos, rid, i, noff;
	u1i *col;
	del_cnsnodes_bspoa(g);
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	clen = 0;
	noff = g->nodes->size;
	if(rps == NULL){
		encap_b1v(g->memp, nseq * sizeof(u4i));
		rps = (u4i*)(g->memp->buffer + g->memp->size);
		g->memp->size += nseq * sizeof(u4i);
	}
	memset(rps, 0, nseq * sizeof(u4i));
	v = ref_bspoanodev(g->nodes, g->HEAD);
	u = new_node_bspoa(g, nseq, -1, 4);
	merge_nodes_bspoa(g, u, v);
	push_u4v(g->ndoffs, g->nodes->size);
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		if(col[nseq] < 4){
			u = new_node_bspoa(g, nseq, clen ++, col[nseq]);
			for(rid=0;rid<nseq;rid++){
				if(col[rid] == col[nseq]){
					v = get_rdnode_bspoa(g, rid, rps[rid]);
					merge_nodes_bspoa(g, u, v);
					u->mpos = pos;
					break;
				}
			}
			if(rid == nseq){ // SHOULD never happen
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
		for(rid=0;rid<nseq;rid++){
			if(col[rid] < 4){
				rps[rid] ++;
			}
		}
	}
	v = ref_bspoanodev(g->nodes, g->TAIL);
	u = new_node_bspoa(g, nseq, clen, 4);
	merge_nodes_bspoa(g, u, v);
	for(i=0;i<=clen;i++){
		connect_rdnode_bspoa(g, nseq, i);
	}
	push_b4v(g->rdregs[0], 0);
	push_b4v(g->rdregs[1], clen);
	encap_bitvec(g->states, g->nodes->size);
}

static inline void del_msanodes_bspoa(BSPOA *g){
	int nseq, clen, i, j;
	nseq = g->seqs->nseq;
	if(Int(g->ndoffs->size) < nseq + 5){
		return;
	}
	for(j=4;j>=0;j--){
		clen = g->nodes->size - g->ndoffs->buffer[nseq + j] - 1;
		for(i=-1;i<=clen;i++){
			cut_rdnode_bspoa(g, nseq + j, i, BSPOA_RDNODE_CUTALL);
		}
		g->nodes->size = g->ndoffs->buffer[nseq + j] - 1;
		g->ndoffs->size = nseq + j;
		g->rdregs[0]->size = nseq + j;
		g->rdregs[1]->size = nseq + j;
	}
}

static inline void add_msanodes_bspoa(BSPOA *g, u4i *rps){
	bspoanode_t *u, *v;
	u4i nseq, nall, mrow, mlen, clen, pos, rid, i;
	u1i *col;
	if(g->ndoffs->size > g->seqs->nseq){
		del_msanodes_bspoa(g);
	}
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	mlen = g->msaidxs->size;
	clen = 0;
	if(rps == NULL){
		encap_b1v(g->memp, nseq * sizeof(u4i));
		rps = (u4i*)(g->memp->buffer + g->memp->size);
		g->memp->size += nseq * sizeof(u4i);
	}
	memset(rps, 0, nseq * sizeof(u4i));
	{
		u = new_node_bspoa(g, nall, -1, 4);
		v = ref_bspoanodev(g->nodes, g->HEAD);
		merge_nodes_bspoa(g, u, v);
		push_u4v(g->ndoffs, g->nodes->size);
		for(pos=0;pos<mlen;pos++){
			col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(col[nall] < 4){
				u = new_node_bspoa(g, nall, clen ++, col[nall]);
				for(rid=0;rid<nseq;rid++){
					if(col[rid] == col[nall]){
						v = get_rdnode_bspoa(g, rid, rps[rid]);
						merge_nodes_bspoa(g, u, v);
						u->mpos = pos;
						break;
					}
				}
				if(rid == nseq){ // SHOULD never happen
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
			for(rid=0;rid<nseq;rid++){
				if(col[rid] < 4){
					rps[rid] ++;
				}
			}
		}
		v = ref_bspoanodev(g->nodes, g->TAIL);
		u = new_node_bspoa(g, nall, clen, 4);
		merge_nodes_bspoa(g, u, v);
		for(i=0;i<=clen;i++){
			connect_rdnode_bspoa(g, nall, i);
		}
		push_b4v(g->rdregs[0], 0);
		push_b4v(g->rdregs[1], clen);
		encap_bitvec(g->states, g->nodes->size);
	}
	for(i=0;i<4;i++){
		u = new_node_bspoa(g, nall + 1 + i, -1, i);
		v = ref_bspoanodev(g->nodes, g->HEAD);
		merge_nodes_bspoa(g, u, v);
		push_u4v(g->ndoffs, g->nodes->size);
		for(pos=0;pos<mlen;pos++){
			u = new_node_bspoa(g, nall + 1 + i, pos, i);
			u->mpos = pos;
		}
		u = new_node_bspoa(g, nall + 1 + i, mlen, i);
		v = ref_bspoanodev(g->nodes, g->TAIL);
		merge_nodes_bspoa(g, u, v);
		push_b4v(g->rdregs[0], 0);
		push_b4v(g->rdregs[1], mlen);
		encap_bitvec(g->states, g->nodes->size);
	}
	memset(rps, 0, nseq * sizeof(u4i));
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(rid=0;rid<nseq;rid++){
			if(col[rid] < 4){
				u = get_rdnode_bspoa(g, rid, rps[rid]);
				v = get_rdnode_bspoa(g, nall + 1 + u->base, pos);
				if(u->header != v->header){
					merge_nodes_bspoa(g, u, v);
				}
				rps[rid] ++;
			}
		}
	}
}

static inline u4i msa_bspoa(BSPOA *g){
	bspoanode_t *v, *u, *x;
	bspoaedge_t *e;
	u1i *qs;
	u4i nseq, mrow, mlen, nidx, eidx, xidx, ready, rid, pos;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	sort_nodes_bspoa(g);
	mlen = g->msaidxs->size;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst = 0;
		u->nct = u->nin;
	}
#if DEBUG
	set_edges_unvisited_bspoa(g);
#endif
	clear_u4v(g->stack);
	push_u4v(g->stack, g->HEAD);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
#if DEBUG
			e->vst = 1;
			(e + 1)->vst = 1;
#endif
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			v->vst ++;
			if(v->vst == v->nct){
				ready = 1;
				xidx = v->next;
				while(xidx != e->node){
					x = ref_bspoanodev(g->nodes, xidx);
					if(x->vst < x->nct){
						ready = 0;
						break;
					}
					xidx = x->next;
				}
				if(ready){
					xidx = e->node;
					do {
						x = ref_bspoanodev(g->nodes, xidx);
						set_u1v(g->msacols, x->mpos * mrow + x->rid, x->base);
						if(x->erev){
							push_u4v(g->stack, xidx);
						}
						xidx = x->next;
					} while(xidx != e->node);
				}
			} else if(v->vst > v->nct){
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
	if(nidx != g->TAIL){
		fprint_dot_bspoa(g, 0, MAX_U4, 1, "1.dot", NULL);
		print_seqs_bspoa(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		print_node_edges_bspoa(g, nidx, 1, stderr);
		print_node_edges_bspoa(g, nidx, 1, stderr);
		check_unvisited_bspoa(g, 0);
		abort();
	}
#if DEBUG
	check_msa_rdseqs_bspoa(g);
#endif
	return g->msaidxs->size;
}

static inline u4i mask_free_ends_msa_bspoa(BSPOA *g){
	bspoanode_t *v, *u;
	bspoaedge_t *e;
	u4i dir, nseq, mrow, mlen, i, nidx, eidx, ridx, rlen, max, midx, ret;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	ret = 0;
	for(dir=0;dir<2;dir++){
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			u->vst = 0;
		}
		max = 0; midx = dir? g->TAIL : g->HEAD;
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
				if(v->next != e->node) continue;
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

// quick and dirty, major bases
static inline void simple_cns_bspoa(BSPOA *g){
	bspoanode_t *v;
	u4i nseq, nall, mrow, mlen, pos, cpos, rid, i, b;
	u2i bcnts[7], brank[7];
	u1i *qs;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
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
			} else if(bcnts[i] && bcnts[i] == bcnts[b]){
				if(brank[i] < brank[b] || b == 4){
					b = i;
				}
			}
		}
		qs[nall] = b;
		qs[nall + 1] = 0;
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
			if(qs[nall] < 4) cpos ++;
		}
	}
	ref_bspoanodev(g->nodes, g->HEAD)->cpos = 0;
	ref_bspoanodev(g->nodes, g->TAIL)->cpos = g->cns->size;
}

#define MAX_LOG_CACHE	1000
static u4i _log_caches_n = 1;
static double _log_caches[MAX_LOG_CACHE + 1];

static inline double cal_permutation_bspoa(u4i n, u4i m){
	if(n > MAX_LOG_CACHE) return 1;
	_log_caches[0] = 0;
	while(_log_caches_n <= n){
		_log_caches[_log_caches_n] = _log_caches[_log_caches_n - 1] + log(_log_caches_n);
		_log_caches_n ++;
	}
	return _log_caches[n] - _log_caches[m] - _log_caches[n - m];
}

static inline double cal_binomial_bspoa(u4i n, u4i m, double p){
	return log(p) * m + log(1 - p) * (n - m) + cal_permutation_bspoa(n, m);
}

#define BS_M_SQRT2 1.4142135623731
static inline double cal_normalCDF_bspoa(double value){
	return erfc(- value / BS_M_SQRT2) / 2;
}

static inline double sum_log_nums(int cnt, double *vals){
	double delta, sum;
	int i;
	if(cnt <= 0) return 0;
#if 1
	sum = BSPOA_MIN_LOGVAL;
	for(i=0;i<cnt;i++){
		if(vals[i] == BSPOA_MIN_LOGVAL) continue;
		if(vals[i] > sum){
			if(vals[i] >= sum + 40){
				sum = vals[i];
				continue;
			} else {
				delta = sum - vals[i];
				sum = vals[i];
			}
		} else {
			if(sum >= vals[i] + 40){
				continue;
			} else {
				delta = vals[i] - sum;
			}
		}
		delta = exp(delta);
		delta = log(1 + delta);
		sum += delta;
	}
#else
	sort_array(vals, cnt, double, num_cmpgt(b, a));
	sum = vals[0];
	for(i=1;i<cnt;i++){
		if(vals[i] == BSPOA_MIN_LOGVAL) break;
		delta = vals[i] - sum;
		if(delta <= -40) break; // overflow
		delta = exp(delta);
		delta = log(1 + delta);
		sum += delta;
	}
#endif
	return sum;
}

static int print_probs = 0;

static inline double cns_bspoa(BSPOA *g){
	typedef struct {double sc[6]; u1i bt, lb;} dp_t;
	bspoanode_t *v;
	//bspoavar_t *var;
	dp_t *dps[5], *dp[5], *lp;
	double ret, errs[10], errd, erre, p, log10, min_freq_calq;
	u1i *r, *qs, *ts, *bs[10];
	u4i a, b, c, lc, d, e, f;
	u4i nseq, nmax, nall, mrow, mlen, cpos, rid, *rps, i, mpsize, cnts[6];
	int pos;
	min_freq_calq = 0.1;
	nseq = num_min(g->nmsa, g->nrds);
	nmax = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	log10 = log(10);
	mlen = g->msaidxs->size;
	mpsize = (mlen + 1) * 5 * sizeof(dp_t);
	mpsize += nseq * sizeof(u4i); // rps
	mpsize += nseq * 10; // bs[0/1/2/3/4] is the state of the last read non-N base vs cns[A/C/G/T/-], M/I/D = 0/1/2; 5-9 is a backup
	clear_and_encap_b1v(g->memp, mpsize);
	r = (u1i*)(g->memp->buffer);
	for(i=0;i<5;i++){
		dps[i] = (dp_t*)r;
		r += (mlen + 1) * sizeof(dp_t);
		memset(dps[i][0].sc, 0, 5 * sizeof(double));
		dps[i][0].sc[5] = (i == 4)? 0 : BSPOA_MIN_LOGVAL;
		dps[i][0].bt = 4;
		dps[i][0].lb = 4;
		dps[i] ++;
	}
	rps = (u4i*)r; r += nseq * sizeof(u4i);
	for(i=0;i<10;i++){
		bs[i] = (u1i*)r; r += nseq;
		memset(bs[i], 0, nseq);
	}
	for(pos=0;pos<Int(mlen);pos++){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(a=0;a<=4;a++){
			dp[a] = dps[a] + pos;
			memset(dp[a]->sc, 0, 6 * sizeof(double));
		}
		memset(cnts, 0, 6 * sizeof(u4i));
		for(rid=0;rid<nseq;rid++){
			b = qs[rid];
			if(b > 4) continue;
			cnts[5] ++;
			cnts[b] ++;
		}
		//fprintf(stdout, "[%d]\t%d\t%d\t%d\t%d\t%d\n", pos, cnts[0], cnts[1], cnts[2], cnts[3], cnts[4]);
		for(i=0;i<5;i++){
			if(cnts[i] < (u4i)(min_freq_calq * cnts[5])){
				cnts[i] = 0;
			}
		}
		for(a=0;a<=4;a++){
			if(cnts[5] && cnts[a] == 0){
				for(i=0;i<6;i++){
					dp[a]->sc[i] = BSPOA_MIN_LOGVAL;
				}
				dp[a]->bt = 4;
				dp[a]->lb = 4;
				for(rid=0;rid<nseq;rid++){
					bs[a + 5][rid] = 0;
				}
				continue;
			}
			for(e=0;e<=4;e++){
				lp = dps[e] + pos - 1;
				c = lp->lb;
				if(cnts[5] && lp->sc[5] == BSPOA_MIN_LOGVAL){
					dp[a]->sc[e] = BSPOA_MIN_LOGVAL;
					errs[e] = BSPOA_MIN_LOGVAL;
				} else {
					for(rid=0;rid<nseq;rid++){
						b = qs[rid];
						if(b > 4) continue;
						d = bs[e][rid];
						ts = g->dptable + b * 5 + c * 25 + d * 125;
						dp[a]->sc[e] += g->dpvals[ts[a] >> 3];
					}
					errs[e] = dp[a]->sc[e] + lp->sc[5];
				}
				errs[e + 5] = errs[e];
			}
			dp[a]->sc[5] = sum_log_nums(5, errs + 5);
			dp[a]->bt = 4;
			for(e=0;e<4;e++){
				if(errs[e] > errs[dp[a]->bt]) dp[a]->bt = e;
			}
			lp = dps[dp[a]->bt] + pos - 1;
			dp[a]->lb = (a < 4)? a : lp->lb;
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
#if 0
		fprintf(stdout, "[%d]", pos);
		for(a=0;a<=4;a++){
			fprintf(stdout, "\t%f:%d", dp[a]->sc[5], dp[a]->bt);
		}
		fprintf(stdout, "\n");
#endif
		for(a=0;a<5;a++){
			memcpy(bs[a], bs[a + 5], nseq);
		}
	}
	pos = mlen - 1;
	c = 4;
	for(a=0;a<4;a++){
		if(dps[a][pos].sc[5] > dps[c][pos].sc[5]){
			c = a;
		}
	}
	ret = dps[c][pos].sc[5];
	clear_u1v(g->cns);
	clear_u1v(g->qlt);
	clear_u1v(g->alt);
	while(1){
		r = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow + nall;
		r[0] = c;
		c = dps[c][pos].bt;
		if(pos == 0) break;
		pos --;
	}
	for(i=0;i<10;i++){
		memset(bs[i], 0, nseq);
	}
	lc = 4;
	for(pos=0;pos<Int(mlen);pos++){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		c = qs[nall];
		errs[0] = dps[0][pos].sc[5];
		errs[1] = dps[1][pos].sc[5];
		errs[2] = dps[2][pos].sc[5];
		errs[3] = dps[3][pos].sc[5];
		errs[4] = dps[4][pos].sc[5];
		erre = sum_log_nums(5, errs);
		errd = dps[c][pos].sc[5];
		erre = log(1 - exp(errd - erre));
		if(c < 4) lc = c;
		erre = - (10 * (erre) / log10);
		if(print_probs){
			fprintf(stdout, "[QLT%04d]\t%c\t%0.2f\t%f\t%f\t%f\t%f\t%f\t=%f\n", pos, "ACGT-"[c], erre, dps[0][pos].sc[5], dps[1][pos].sc[5], dps[2][pos].sc[5], dps[3][pos].sc[5], dps[4][pos].sc[5], sum_log_nums(5, errs));
		}
		// cns quality
		qs[nall + 1] = (int)num_min(erre, BSPOA_QLT_MAX);
		if(1){
			//if(pos == 12061){
				//fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			//}
			memset(cnts, 0, 6 * sizeof(u4i));
			for(rid=0;rid<nmax;rid++){
				b = qs[rid];
				if(b > 4) continue;
				cnts[5] ++;
				cnts[b] ++;
			}
			a = (c + 1) % 5;
			for(e=0;e<=4;e++){
				if(e == c) continue;
				if(cnts[e] > cnts[a]) a = e;
			}
			if(0){
				e = 4;
				for(i=pos;i>0;i--){
					e = g->msacols->buffer[g->msaidxs->buffer[i - 1] * mrow + nall];
					if(e < 4) break;
				}
				d = 4;
				for(i=pos+1;i<mlen;i++){
					d = g->msacols->buffer[g->msaidxs->buffer[i] * mrow + nall];
					if(d < 4) break;
				}
				p = num_max(g->dporis[g->dptable[c + a * 5 + e * 25 + e * 125] >> 3], g->dporis[g->dptable[c + a * 5 + d * 25 + d * 125] >> 3]);
			} else {
				p = g->par->psub;
			}
			erre = 0;
			if(cnts[5] > 50 && cnts[5] * p > 5 && cnts[5] * (1 - p) > 5){
				if(0){
					for(e=0;e<cnts[a];e++){
						erre += exp(cal_binomial_bspoa(cnts[5], e, p));
					}
					double stdv = sqrt(cnts[5] * p * (1 - p));
					double perr = (cnts[a] - cnts[5] * p) / stdv;
					double errf = cal_normalCDF_bspoa(perr);
					fflush(stdout); fprintf(stdout, " --  [%d]\t%d\t%d\t%f\tNORMAL:%f\t%f in %s -- %s:%d --\n", pos, cnts[5], cnts[a], erre, perr, errf, __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
				} else {
					erre = cal_normalCDF_bspoa((cnts[a] - cnts[5] * p) / sqrt(cnts[5] * p * (1- p)));
				}
			} else {
				for(e=0;e<cnts[a];e++){
					erre += exp(cal_binomial_bspoa(cnts[5], e, p));
				}
			}
			if(erre == 0){
				errd = 0;
			} else {
				errd = - (10 * log(1 - erre) / log10);
			}
			if(print_probs){
				fprintf(stdout, "[ALT%04d]\t%c\t%c\t%0.2f\t%0.2f\t%d\t%d\n", pos, "ACGT-"[c], "ACGT-"[a], erre, errd, cnts[c], cnts[a]);
			}
			if(_DEBUG_LOG_ > 2){
				char flanks[2][3];
				for(f=2,i=pos;i>0&&f>0;i--){
					d = g->msacols->buffer[g->msaidxs->buffer[i - 1] * mrow + nall];
					if(d < 4){
						flanks[0][f - 1] = "ACGT-"[d];
						f --;
					}
				}
				for(f=0,i=pos+1;i<mlen&&f<2;i++){
					d = g->msacols->buffer[g->msaidxs->buffer[i] * mrow + nall];
					if(d < 4) flanks[1][f++] = "ACGT-"[d];
				}
				flanks[0][2] = flanks[1][2] = 0;
				fprintf(stderr, "[ALTQ] %d/%d\t%s\t%c\t%c\t%s\t%d\t%d\t%f\t%f\t%f\n", pos, mlen, flanks[0], "ACGT- "[c], "ACGT-"[a], flanks[1], cnts[1], cnts[0], p, erre, errd);
			}
			qs[nall + 2] = num_min(errd, BSPOA_QLT_MAX);
		} else {
			qs[nall + 2] = 0;
		}
		if(qs[nall + 0] < 4){
			push_u1v(g->cns, qs[nall + 0]);
			push_u1v(g->qlt, qs[nall + 1]);
			push_u1v(g->alt, qs[nall + 2]);
		}
	}
	//reverse_u1v(g->cns);
	//reverse_u1v(g->qlt);
	//reverse_u1v(g->alt);
	if(_DEBUG_LOG_ > 2){
		for(pos=0;pos<Int(mlen);pos++){
			double minp = dps[0][pos].sc[5];
			for(a=1;a<=4;a++){
				if(dps[a][pos].sc[5] < minp) minp = dps[a][pos].sc[5];
			}
			fprintf(stderr, "[CNSQ] %d\t%c\t%0.2f", pos, "ACGT- "[g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + nall]], minp);
			for(a=0;a<=4;a++){
				fprintf(stderr, "\t%c:%0.2f[", "ACGT-"[a], dps[a][pos].sc[5] - minp);
				for(e=0;e<=4;e++){
					if(e == dps[a][pos].bt) fprintf(stderr, "*");
					fprintf(stderr, "%c:%0.2f", "acgt-"[e], dps[a][pos].sc[e] + dps[e][pos+1].sc[5] - minp);
					if(e < 4) fprintf(stderr, ",");
				}
				fprintf(stderr, "]");
			}
			fprintf(stderr, "\n");
		}
	}
	// set node_t->cpos to the position on cns, useful in bandwidth
	if(g->nodes->size){
		for(rid=0;rid<nall;rid++){
			cpos = 0;
			//v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + g->seqs->rdlens->buffer[rid] - 1);
			v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid]);
			for(pos=0;pos<Int(mlen);pos++){
				if(g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + rid] < 4){
					v->cpos = cpos;
					v ++;
				}
				if(g->msacols->buffer[g->msaidxs->buffer[pos] * mrow + nall] < 4) cpos ++;
			}
		}
		ref_bspoanodev(g->nodes, g->HEAD)->cpos = 0;
		ref_bspoanodev(g->nodes, g->TAIL)->cpos = g->cns->size;
	}
	return ret;
}

static inline int diag2idx(int W, int xi, int yi){
	int rowlen = W * WORDSIZE + 2;
	int rowbeg = 1;
	int half = W * WORDSIZE / 2;
	int moff = (xi + yi);
	int mdir = (moff & 0x1);
	int midx = (xi - yi - mdir) / 2 + half;
	if(moff < 0 || midx < 0 || midx >= W * WORDSIZE){
		return -1;
	}
	return rowlen * (moff + 1) + rowbeg + midx;
}

static inline void maxmat_dp_diag_rowcal_init(int W, int mbeg, u1i *matrix[2]){
	int rowlen, rowbeg, half, idx;
	rowlen = W * WORDSIZE + 2;
	rowbeg = 1;
	half   = W * WORDSIZE / 2;
	idx = (mbeg + mbeg) * rowlen;
	memset(matrix[0] + idx, 0, rowlen);
	memset(matrix[1] + idx, 0, rowlen);
	matrix[0][idx + rowbeg + half - 1] = MAX_U1;
	matrix[1][idx + rowbeg + half    ] = MAX_U1;
}

static inline int maxmat_dp_diag_rowcal_prepare(int x, int y, int mlen, int W, u1i *matrix[2], u1i *_seqs[2], u1i *_mats[2][4], u1i *rows[2][2], u1i *seqs[2], u1i *mats[2][4]){
	int rowlen, rowbeg, half, xb, yb, moff, mdir, midx;
	rowlen = W * WORDSIZE + 2;
	rowbeg = 1;
	half   = W * WORDSIZE / 2;
	moff = x + y;
	mdir = (moff & 0x1);
	midx = (x - y - mdir) / 2 + half;
	xb = x - midx;
	yb = mlen - 1 - (y + midx);
	rows[0][0] = matrix[0] + rowlen * moff + rowbeg;
	rows[0][1] = matrix[1] + rowlen * moff + rowbeg;
	rows[1][0] = matrix[0] + rowlen * (moff + 1) + rowbeg;
	rows[1][1] = matrix[1] + rowlen * (moff + 1) + rowbeg;
	seqs[0]    = _seqs[0] + xb;
	seqs[1]    = _seqs[1] + yb;
	mats[0][0] = _mats[0][0] + xb;
	mats[0][1] = _mats[0][1] + xb;
	mats[0][2] = _mats[0][2] + xb;
	mats[0][3] = _mats[0][3] + xb;
	mats[1][0] = _mats[1][0] + yb;
	mats[1][1] = _mats[1][1] + yb;
	mats[1][2] = _mats[1][2] + yb;
	mats[1][3] = _mats[1][3] + yb;
	return midx;
}

static inline void maxmat_dp_diag_rowcal_verify(int x, int y, int W, int dir, u1i *matrix[2], u1i *rows[2][2], u1i *seqs[2], u1i *mats[2][4], int display){
	int i, len, h, s, u, v, b;
	len = W * WORDSIZE;
	if(display){
		fprintf(stderr, "[%d:%d]", x, y);
	}
	for(i=0;i<len;i++){
		s = (seqs[1][i] < 4? mats[0][seqs[1][i]][i] : 0) + (seqs[0][i] < 4? mats[1][seqs[0][i]][i] : 0);
		if(s > Int(MAX_U1)) s = MAX_U1;
		if(dir){
			u = rows[0][0][i + 1];
			v = rows[0][1][i];
		} else {
			u = rows[0][0][i];
			v = rows[0][1][i - 1];
		}
		{
			int idx;
			idx = diag2idx(W, x + i - W * WORDSIZE / 2, y - i + W * WORDSIZE / 2);
			if(idx != -1 && matrix[0][idx] != rows[1][0][i]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			if(idx != -1 && matrix[1][idx] != rows[1][1][i]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			idx = diag2idx(W, x + i - W * WORDSIZE / 2 - 1, y - i + W * WORDSIZE / 2);
			if(idx != -1 && matrix[1][idx] != v){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			idx = diag2idx(W, x + i - W * WORDSIZE / 2, y - i + W * WORDSIZE / 2 - 1);
			if(idx != -1 && matrix[0][idx] != u){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
		if(s < u){
			b = SEQALIGN_BT_D;
			h = u;
		} else {
			h = s;
			b = SEQALIGN_BT_M;
		}
		if(h < v){
			h = v;
			b = SEQALIGN_BT_I;
		}
		if(rows[1][0][i] != h - v){
			fflush(stdout); fprintf(stderr, " -- [%d:%d] in %s -- %s:%d --\n", x, y, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		if(rows[1][1][i] != h - u){
			fflush(stdout); fprintf(stderr, " -- [%d:%d] in %s -- %s:%d --\n", x, y, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		if(display){
			fprintf(stderr, "\t[%d]%d:%d:%d:%c", i, s, u, v, "MID"[b]);
		}
	}
	if(display) fprintf(stderr, "\n");
}

static xint MM_EPI8_ALL0;
static xint MM_EPI8_ALL1;
static xint MM_EPI8_ALL2;
static xint MM_EPI8_ALL3;

static inline void maxmat_dp_diag_rowcal(int W, int dir, u1i *rows[2][2], u1i *seqs[2], u1i *mats[2][4]){
	xint h, u, v, x, y;
	int i;
	for(i=0;i<W;i++){
		h = mm_loadu(((xint*)seqs[1]) + i);
		x = mm_setzero();
		x = mm_blendv(x, mm_loadu(((xint*)mats[0][0]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL0));
		x = mm_blendv(x, mm_loadu(((xint*)mats[0][1]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL1));
		x = mm_blendv(x, mm_loadu(((xint*)mats[0][2]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL2));
		x = mm_blendv(x, mm_loadu(((xint*)mats[0][3]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL3));
		h = mm_loadu(((xint*)seqs[0]) + i);
		y = mm_setzero();
		y = mm_blendv(y, mm_loadu(((xint*)mats[1][0]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL0));
		y = mm_blendv(y, mm_loadu(((xint*)mats[1][1]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL1));
		y = mm_blendv(y, mm_loadu(((xint*)mats[1][2]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL2));
		y = mm_blendv(y, mm_loadu(((xint*)mats[1][3]) + i), mm_cmpeq_epi8(h, MM_EPI8_ALL3));
		h = mm_adds_epu8(x, y);
		if(dir){ // left
			u = mm_loadu((xint*)(rows[0][0] + i * WORDSIZE + 1));
			v = mm_loadu((xint*)(rows[0][1] + i * WORDSIZE));
		} else { // down
			u = mm_loadu((xint*)(rows[0][0] + i * WORDSIZE));
			v = mm_loadu((xint*)(rows[0][1] + i * WORDSIZE - 1));
		}
		h = mm_max_epu8(h, u);
		h = mm_max_epu8(h, v);
		mm_storeu(((xint*)rows[1][0]) + i, mm_subs_epu8(h, v));
		mm_storeu(((xint*)rows[1][1]) + i, mm_subs_epu8(h, u));
	}
	if(dir){
		rows[1][0][-1] = MAX_U1;
		rows[1][1][-1] = 0;
		rows[1][0][W * WORDSIZE] = 0;
		rows[1][1][W * WORDSIZE] = 0;
	} else {
		rows[1][0][-1] = 0;
		rows[1][1][-1] = 0;
		rows[1][0][W * WORDSIZE] = 0;
		rows[1][1][W * WORDSIZE] = MAX_U1;
	}
}

static inline void print_diag_region(int W, int mbeg, int mend, u1i *matrix[2]){
	int x, y, idx, u, v;
	for(x=mbeg;x<mend;x++){
		printf("%d\t", x);
		for(y=mbeg;y<mend;y++){
			idx = diag2idx(W, x, y);
			if(idx == -1){
				u = v = -1;
			} else {
				u = matrix[0][idx];
				v = matrix[1][idx];
			}
			printf("%d:%d\t", u, v);
		}
		printf("\n");
	}
}

static inline int remsa_pedit_rd_bspoacore(BSPOA *g, u2i rid, u4i rbeg, u4i rend, u1i *matrix[2], u1i *_seqs[2], u1i *_mats[2][4], int mlen, int mbeg, int mend, int W){
	bspoanode_t *u, *v;
	int i, x, y, xi, yi, xx, dir, roff, rowlen, half, rowbeg, scr;
	u1i *rows[2][2], *seqs[2], *mats[2][4], *tmp;
	UNUSED(rbeg);
	rowlen = W * WORDSIZE + 2;
	rowbeg = 1;
	half   = W * WORDSIZE / 2;
	tmp = alloca(W * WORDSIZE);
	maxmat_dp_diag_rowcal_init(W, mbeg, matrix);
	x = y = mbeg;
	for(i=x+y;;i++){
		dir = (i & 0x1);
		maxmat_dp_diag_rowcal_prepare(x, y, mlen, W, matrix, _seqs, _mats, rows, seqs, mats);
		maxmat_dp_diag_rowcal(W, dir, rows, seqs, mats);
		//maxmat_dp_diag_rowcal_verify(x, y, W, dir, matrix, rows, seqs, mats, 0);
		if(dir) y ++; // NOTE: next row
		else x ++;
		if(x >= mend) break;
	}
#if 0
	void simple_eval_score(int beg, int end){
		int xi, yi, scr;
		scr = 0;
		xi = end - 1;
		yi = end - 1;
		while(1){
			i = xi + yi;
			if(i < beg + beg) break;
			dir = (i & 0x1);
			xx = maxmat_dp_diag_rowcal_prepare(xi, yi, mlen, W, matrix, _seqs, _mats, rows, seqs, mats);
			x = xi - xx + half;
			y = yi + xx - half;
			maxmat_dp_diag_rowcal_verify(x, y, W, dir, matrix, rows, seqs, mats, 0);
			{
				int h;
				h = (seqs[1][xx] < 4? mats[0][seqs[1][xx]][xx] : 0) + (seqs[0][xx] < 4? mats[1][seqs[0][xx]][xx] : 0);
				scr += h;
				xi --;
				yi --;
			}
		}
		fflush(stdout); fprintf(stderr, " -- R%d %d-%d -> %d in %s -- %s:%d --\n", rid, beg, end, scr, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	simple_eval_score(mbeg, mend);
#endif
	scr = 0;
	xi = mend - 1;
	yi = mend - 1;
	roff = rend;
	//print_diag_region(W, mend - 10, mend, matrix);
	while(xi >= 0 && yi >= 0){
		//simple_eval_score(xi + 1, mend);
		i = xi + yi;
		if(i < mbeg + mbeg) break;
		dir = (i & 0x1);
		xx = maxmat_dp_diag_rowcal_prepare(xi, yi, mlen, W, matrix, _seqs, _mats, rows, seqs, mats);
		x = xi - xx + half;
		y = yi + xx - half;
		//maxmat_dp_diag_rowcal_verify(x, y, W, dir, matrix, rows, seqs, mats, 0);
		if(xx < 0 || xx >= W * WORDSIZE){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		{
			int s, h, e, f, bt;
			h = (seqs[1][xx] < 4? mats[0][seqs[1][xx]][xx] : 0) + (seqs[0][xx] < 4? mats[1][seqs[0][xx]][xx] : 0);
			if(h > Int(MAX_U1)) h = MAX_U1;
			if(dir){
				e = rows[0][0][xx + 1];
				f = rows[0][1][xx];
			} else {
				e = rows[0][0][xx];
				f = rows[0][1][xx - 1];
			}
#if DEBUG
			if(f + rows[1][0][xx] != e + rows[1][1][xx]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			if(rows[1][0][xx] != matrix[0][diag2idx(W, xi, yi)]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			if(rows[1][1][xx] != matrix[1][diag2idx(W, xi, yi)]){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			s = f + rows[1][0][xx];
			if(s == f && !(xx == 0 && dir == 0)){
				if(seqs[0][xx] < 4) roff --;
				xi --;
				bt = 1;
			} else if(s == e){
				yi --;
				bt = 2;
			} else if(s == h){
				if(seqs[0][xx] < 4){
					roff --;
					u = get_rdnode_bspoa(g, g->seqs->nseq + 1 + seqs[0][xx], yi);
					v = get_rdnode_bspoa(g, rid, roff);
					//v->colorful = 1;
#ifdef DEBUG
					if(v->base != seqs[0][xx]){
						fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					merge_nodes_bspoa(g, u, v);
				}
				scr += s;
				bt = 0;
				xi --;
				yi --;
			} else {
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			//printf("RID[%d] %c\t%d\t%d\tscr=%d\ts=%d\th=%d\te=%d\tf=%d\n", rid, "MID"[bt], xi, yi, scr, s, h, e, f);
		}
	}
#if 0
	if(roff > Int(rbeg) + 2){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	return scr;
}

static inline u4i calc_profile_mismatches_bspoa_core(u4i nseq, u1i *v1, u1i *v2){
	xint a, b, c, d;
	u4i i, xx, yy, rc, rm, mask;
	rc = nseq / WORDSIZE;
	rm = nseq % WORDSIZE;
	xx = 0;
	if(v1 == v2){
		if(1){
			for(yy=i=0;i<nseq;i++){
				if(v1[i] >= 4) yy ++;
			}
		}
		d = mm_set1_epi8(3);
		for(i=0;i<rc;i++){
			a = mm_loadu(((xint*)v1) + i);
			c = mm_cmpgt_epi8(a, d);
			mask = mm_movemask_epi8(c);
			xx += __builtin_popcount(mask);
		}
		if(rm){
			a = mm_loadu(((xint*)v1) + i);
			c = mm_cmpgt_epi8(a, d);
			mask = mm_movemask_epi8(c);
			xx += __builtin_popcount(mask << (32 - rm));
		}
		if(yy != xx){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	} else {
		if(1){
			for(yy=i=0;i<nseq;i++){
				if(((v1[i] ^ v2[i]) & 4) == 0) yy ++;
			}
		}
		d = mm_set1_epi8(4);
		for(i=0;i<rc;i++){
			a = mm_loadu(((xint*)v1) + i);
			b = mm_loadu(((xint*)v2) + i);
			c = mm_xor(a, b);
			c = mm_cmpgt_epi8(d, c);
			mask = mm_movemask_epi8(c);
			xx += __builtin_popcount(mask);
		}
		if(rm){
			a = mm_loadu(((xint*)v1) + i);
			b = mm_loadu(((xint*)v2) + i);
			c = mm_xor(a, b);
			c = mm_cmpgt_epi8(d, c);
			mask = mm_movemask_epi8(c);
			xx += __builtin_popcount(mask << (32 - rm));
		}
		if(yy != xx){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
	return xx;
}

//TODO: unfinished
static inline void optimise_profile_bspoa(BSPOA *g, u1i *prof[2][4], u4i *matrix, int bandwidth){
	u4i nseq, mrow, h, e, f, *rows[2], SMAX;
	int mlen, i, j, x, y, jb, je, W, HW, bt;
	u1i *p1, *p2;
	UNUSED(prof);
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	W = bandwidth + 1;
	HW = bandwidth / 2;
	SMAX = MAX_U4 >> 2;
	rows[1] = matrix;
	for(i=0;i<HW;i++) rows[1][i] = SMAX;
	for(;i<W;i++) rows[1][i] = SMAX;
	for(i=0;i<mlen;i++){
		rows[0] = rows[1];
		rows[1] += W;
		jb = i > HW? 0 : HW - i;
		je = i + HW < mlen? bandwidth : mlen + HW - i;
		f = rows[1][jb - 1] = SMAX;
		p2 = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
		for(j=jb;j<je;j++){
			p1 = g->msacols->buffer + g->msaidxs->buffer[j + i - HW] * mrow;
			h = calc_profile_mismatches_bspoa_core(nseq, p1, p2);
			h = h + rows[0][j - 1];
			e = rows[0][j] + 1;
			if(h > e) h = e;
			if(h > f) h = f;
			rows[1][j] = h;
			f = h + 1;
		}
		rows[1][je] = SMAX;
	}
	x = y = mlen - 1;
	while(y >= 0){
		rows[0] = matrix + W * y;
		rows[1] = rows[0] + W;
		j = y - (x - HW);
		if(j < 0 || j >= bandwidth){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		p1 = g->msacols->buffer + g->msaidxs->buffer[x] * mrow;
		p2 = g->msacols->buffer + g->msaidxs->buffer[y] * mrow;
		h = calc_profile_mismatches_bspoa_core(nseq, p1, p2);
		h += rows[0][j - 1];
		e = rows[0][j] + 1;
		f = rows[1][j - 1] + 1;
		if(h > e){
			h = e;
			bt = 1;
		} else {
			bt = 0;
		}
		if(h > f){
			h = f;
			bt = 2;
		}
		//TODO: revise prof
		if(bt == 0){
			x --;
			y --;
		} else if(bt == 1){
			y --;
		} else {
			x --;
		}
	}
}

static inline void remsa_pedits_bspoa(BSPOA *g, u4i bandwidth, int hpadjust, int all){
	typedef union { struct { u4i base:3, cns:1, off:12, bcnt:16; }; u4i val; } hp_base_t;
	hp_base_t PB;
	bspoanode_t *v, *u;
	b1v *mmp;
	u4i rid, nseq, nall, mrow, mlen, W, HW, pos, lpos, p, i, j, rdlen, mbeg, mend, mpsize, cc, bc, mc, cnts[4], scr;
	u4i qb, qe;
	u1i *col, *seqs[2], *matrix[5], *mats[2][4], lc;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	mlen = g->msaidxs->size;
	if(nseq >= MAX_U1) return;
	MM_EPI8_ALL0 = mm_set1_epi8(0);
	MM_EPI8_ALL1 = mm_set1_epi8(1);
	MM_EPI8_ALL2 = mm_set1_epi8(2);
	MM_EPI8_ALL3 = mm_set1_epi8(3);
	//print_msa_bspoa(g, "RWALIGN", 0, 0, 0, 1, stdout);
	bandwidth = roundup_times(bandwidth, WORDSIZE);
	W = bandwidth / WORDSIZE;
	HW = bandwidth / 2;
	//mpsize = nseq * sizeof(u4i);
	//clear_and_encap_b1v(g->memp, mpsize);
	//rps = (u4i*)g->memp->buffer;
	add_msanodes_bspoa(g, NULL);
	if(all && nseq < g->seqs->nseq){
		mmp = adv_init_b1v(1024, 0, WORDSIZE, 0);
	} else {
		mmp = NULL;
	}
	mpsize = 2 * roundup_times(mlen + bandwidth, WORDSIZE); // seqs[]
	mpsize += 5 * roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE); // matrix[]
	mpsize += 2 * 4 * roundup_times(mlen + bandwidth, WORDSIZE); // mats[][]
	clear_and_encap_b1v(g->memp, mpsize);
	memset(g->memp->buffer, 0, mpsize);
	seqs[0] = (u1i*)g->memp->buffer;
	seqs[1] = seqs[0] + roundup_times(mlen + bandwidth, WORDSIZE);
	matrix[0] = seqs[1] + roundup_times(mlen + bandwidth, WORDSIZE);
	matrix[1] = matrix[0] + roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE);
	matrix[2] = matrix[1] + roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE);
	matrix[3] = matrix[2] + roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE);
	matrix[4] = matrix[3] + roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE);
	mats[0][0] = matrix[4] + roundup_times((2 * mlen + 1) * (bandwidth + 2), WORDSIZE) + HW;
	mats[0][1] = mats[0][0] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[0][2] = mats[0][1] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[0][3] = mats[0][2] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[1][0] = mats[0][3] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[1][1] = mats[1][0] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[1][2] = mats[1][1] + roundup_times(mlen + bandwidth, WORDSIZE);
	mats[1][3] = mats[1][2] + roundup_times(mlen + bandwidth, WORDSIZE);
	seqs[0] += HW;
	seqs[1] += HW;
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		seqs[1][mlen - 1 - pos] = col[nall];
		for(rid=0;rid<nseq;rid++){
			if(col[rid] < 4){
				mats[1][col[rid]][pos] ++;
			}
		}
	}
	if(hpadjust){
		//adjust [cns=4] minor bases to right-side identitical cns
		for(pos=0;pos<mlen;pos++){
			if((lc = seqs[1][mlen - 1 - pos]) < 4){
				for(i=pos;i;i--){
					if(seqs[1][mlen - 1 - (i - 1)] < 4) break;
					if(mats[1][lc][i - 1] && UInt(mats[1][lc][i - 1]) + UInt(mats[1][lc][pos]) <= MAX_U1){
						mats[1][lc][pos] = mats[1][lc][pos] + mats[1][lc][i - 1];
						mats[1][lc][i - 1] = 0;
					}
				}
			}
		}
		// adjust bases and cnts within a cns-homopolymer
		lc = 4; mc = 0; memset(cnts, 0, 4 * sizeof(u4i));
		clear_u4v(g->stack);
		col = NULL;
		for(lpos=pos=0;pos<=mlen;pos++){
			if(pos == mlen || ((col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow)[nall] < 4 && col[nall] != lc)){
				sort_array(g->stack->buffer, g->stack->size, hp_base_t, num_cmpgt(a.base, b.base));
				for(i=p=0;i<=g->stack->size;i++){
					if(i < g->stack->size && ((hp_base_t)g->stack->buffer[i]).base == ((hp_base_t)g->stack->buffer[p]).base){
						continue;
					}
					cc = cnts[TOTYPE(g->stack->buffer[p], hp_base_t).base];
					if(TOTYPE(g->stack->buffer[p], hp_base_t).base == lc){
						sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgtx(b.cns, a.cns, a.off, b.off));
						for(j=p;cc&&j<i;j++){
							PB = TOTYPE(g->stack->buffer[j], hp_base_t);
							if(PB.cns == 0) break;
							bc = num_min(cc, mc);
							mats[1][PB.base][lpos + PB.off] = bc;
							cc -= bc;
						}
						while(p < j){
							PB = TOTYPE(g->stack->buffer[p], hp_base_t);
							if(j - p + mats[1][PB.base][lpos + PB.off] <= MAX_U1){
								mats[1][PB.base][lpos + PB.off] += (j - p);
							}
							p ++;
						}
						p = j;
						sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgt(b.bcnt, a.bcnt));
						for(j=p;cc&&j<i;j++){
							PB = TOTYPE(g->stack->buffer[j], hp_base_t);
							bc = num_min(cc, mc);
							mats[1][PB.base][lpos + PB.off] = bc;
							cc -= bc;
						}
					} else {
						sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgt(b.bcnt, a.bcnt));
						for(j=p;cc&&j<i;j++){
							PB = TOTYPE(g->stack->buffer[j], hp_base_t);
							bc = num_min(cc, mc);
							mats[1][PB.base][lpos + PB.off] = bc;
							cc -= bc;
						}
					}
					p = i;
				}
			}
			if(pos == mlen) break;
			if(col[nall] < 4 && col[nall] != lc){
				lc = col[nall]; mc = 0;
				memset(cnts, 0, 4 * sizeof(u4i));
				lpos = pos;
				clear_u4v(g->stack);
			}
			for(i=0;i<4;i++){
				if(mats[1][i][pos]){
					if(mats[1][i][pos] > mc) mc = mats[1][i][pos];
					cnts[i] += mats[1][i][pos];
					PB.base = i;
					PB.cns  = (i == col[nall]);
					PB.off  = pos - lpos;
					PB.bcnt = mats[1][i][pos];
					push_u4v(g->stack, PB.val);
					mats[1][i][pos] = 0;
				}
			}
		}
	}
#if 0
	for(lc=0;lc<4;lc++){
		if(1){
			printf("[%c]              ", "ACGT"[lc]);
			for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i] / 100) % 10);
			printf("\n");
			printf("[%c]              ", "ACGT"[lc]);
			for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i] / 10) % 10);
			printf("\n");
			printf("[%c]              ", "ACGT"[lc]);
			for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i]) % 10);
			printf("\n");
			printf("\n");
		}
	}
	fflush(stdout);
	exit(1);
#endif
	for(lc=0;lc<4;lc++){
		i = 0; j = mlen - 1;
		while(i < j){
			swap_var(mats[1][lc][i], mats[1][lc][j]);
			i ++; j --;
		}
	}
	for(rid=0;rid<g->seqs->nseq;rid++){
		rdlen = g->seqs->rdlens->buffer[rid];
		if(rdlen == 0) continue;
		memset(seqs[0] - HW, 4, (mlen + bandwidth)); // TODO: when query read is very short, should only memset the span of read
		memset(mats[0][0] - HW, 0, 4 * (mlen + bandwidth));
		if(rid >= nseq){ // aligning read aginst consensus seq to assign mpos for each bases in aligned region
			seqalign_result_t rs;
			u8i iter;
			int op, x, y, z;
			if(!all) break;
			g->qlen = g->slen = rdlen;
			g->qb = 0; g->qe = g->qlen;
			resize_u1v(g->qseq, rdlen);
			bitseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[rid], g->qlen, g->qseq->buffer);
			if(g->par->ksz){
				rs = kmer_striped_seqedit_pairwise(g->par->ksz, g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, mmp, g->stack, 0);
			} else {
				rs = striped_seqedit_pairwise(g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->par->alnmode, 0, mmp, g->stack, 0);
			}
			g->qb = rs.qb;
			g->qe = rs.qe;
			g->slen = g->qe - g->qb;
			for(i=0;i<rdlen;i++){
				v = get_rdnode_bspoa(g, rid, i);
				v->mpos = 0;
			}
			iter = 0;
			x = rs.qb;
			y = rs.tb;
			z = 0;
			while((op = seqalign_iter_cigars(g->stack->buffer, g->stack->size, &iter)) != -1){
				switch(op){
					case 0:
					case 7:
					case 8: // match/mistmatch
					u = get_rdnode_bspoa(g, g->seqs->nseq, y); // cns node
					v = get_rdnode_bspoa(g, rid, x); // read node
					v->mpos = u->mpos;
					x ++; y ++;
					break;
					case 1:
					case 4: // insertion
					x ++;
					break;
					case 2:
					case 3: // deletion
					y ++;
					break;
					case 5: // H, used to locate the begin and end
					y ++;
					break;
				}
			}
			// revise leading
			y = -1;
			for(x=rs.qb;x<rs.qe;x++){
				v = get_rdnode_bspoa(g, rid, x);
				if(v->mpos){
					y = v->mpos;
					break;
				}
			}
			for(z=x;z>=rs.qb;z--){
				v = get_rdnode_bspoa(g, rid, z);
				if(y + z - x < 0) break;
				v->mpos = y + z - x;
			}
			qb = ++ z;
			// infering mpos
			for(;x<rs.qe;x++){
				v = get_rdnode_bspoa(g, rid, x);
				if(v->mpos == 0){
					v->mpos = y ++;
				} else if(v->mpos < y){
					v->mpos = y ++;
				} else {
					y = v->mpos + 1;
				}
				if(y >= Int(mlen)) break;
			}
			// revise tailing
			qe = x;
		} else {
			qb = 0; qe = rdlen;
		}
		{
			lc = 4; cc = 0;
			scr = 0;
			for(i=qe;i>qb;i--){
				v = cut_rdnode_bspoa(g, rid, i - 1, BSPOA_RDNODE_CUTALL);
				seqs[0][v->mpos] = v->base;
				scr += mats[1][v->base][mlen - 1 - v->mpos];
				if(v->base == lc){
					if(cc < MAX_U1) cc ++;
					mats[0][v->base][v->mpos] = cc;
				} else {
					lc = v->base;
					cc = 0;
				}
			}
			v = get_rdnode_bspoa(g, rid, qb);
			mbeg = v->mpos;
			v = get_rdnode_bspoa(g, rid, qe - 1);
			mend = v->mpos + 1;
		}
		scr = remsa_pedit_rd_bspoacore(g, rid, qb, qe, matrix, seqs, mats, mlen, mbeg, mend, W);
		connect_rdnodes_bspoa(g, rid);
	}
	del_msanodes_bspoa(g);
	if(all) g->nrds = g->seqs->nseq;
	if(mmp) free_b1v(mmp);
}

static inline u4i remsa_edit_rd_bspoacore(BSPOA *g, u2i rid, u4i rend, u2i *bcnts, u4i *breps, u1i *seqs[2], int *rows[2], u1i *matrix, int mbeg, int mend, int W){
	bspoanode_t *u, *v;
	int *cur, *lst, h, e, s, SMIN, scr;
	int x, y, xi, bt, roff, HW;
	u2i *bs;
	u1i *seq, *cns, *mtx, lb, cb;
	SMIN = - (MAX_U4 >> 2);
	HW = W / 2; // W MUST be even
	lst = rows[0];
	cur = rows[1];
	seq = seqs[0];
	cns = seqs[1];
	for(x=0;x<HW;x++) lst[x] = SMIN;
	for(x=HW;x<W;x++) lst[x] = 0;
	lst[W] = cur[W] = SMIN;
	for(y=mbeg;y<mend;y++){
		mtx = matrix + y * W;
		bs = bcnts + y * 5;
		cb = cns[y] << 2;
		s = SMIN;
		for(x=0;x<W;x++){
			xi = y + x - HW;
			bt = SEQALIGN_BT_I;
			lb = seq[xi];
			h = lst[x] + bs[lb] + ((breps[xi] >> cb) & 0xf);
			e = lst[x + 1];
			if(s < h){
				s = h;
				bt = SEQALIGN_BT_M;
			}
			if(s < e){
				s = e;
				bt = SEQALIGN_BT_D;
			}
			cur[x] = s;
			mtx[x] = bt;
		}
		swap_var(lst, cur);
	}
	scr = lst[0];
	bt = SEQALIGN_BT_M;
	y = mend - 1;
	x = mend - 1;
	roff = rend;
	while(1){
		xi = x - y + HW;
		if(xi < 0 || xi >= W){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		bt = matrix[y * W + xi];
		if(bt == SEQALIGN_BT_M){
			if(seq[x] < 4){
				roff --;
				u = get_rdnode_bspoa(g, g->seqs->nseq + seq[x], y);
				v = get_rdnode_bspoa(g, rid, roff);
				//v->colorful = 1;
				if(v->base != seq[x]){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				merge_nodes_bspoa(g, u, v);
			}
			x --;
			y --;
		} else if(bt == SEQALIGN_BT_I){
			if(seq[x] < 4){
				roff --;
			}
			x --;
		} else {
			y --;
		}
		if(x < mbeg) break;
		if(y < mbeg) break;
	}
	return scr;
}

static inline void remsa_edits_bspoa(BSPOA *g, u4i W){
	typedef union { struct { u4i base:3, cns:1, off:12, bcnt:16; }; u4i val; } hp_base_t;
	hp_base_t PB;
	bspoanode_t *v;
	u4i nseq, nall, mrow, mlen, HW, pos, lpos, p, i, j, rdlen, mbeg, mend, mpsize, *rps, *breps, cc, bc, mc, cnts[4];
	int *rows[2], scr;
	u2i *bcnts, *bs, rid;
	u1i *col, *seqs[2], *matrix, lc;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	mlen = g->msaidxs->size;
	W = (W + 1) & (~0x1U); // make sure it is even
	HW = W >> 1;
	if(mlen < W) return;
	//print_msa_bspoa(g, "RWALIGN", 0, 0, 0, 1, stdout);
	mpsize = nseq * sizeof(u4i);
	clear_and_encap_b1v(g->memp, mpsize);
	rps = (u4i*)g->memp->buffer;
	add_msanodes_bspoa(g, rps);
	mpsize = 2 * (W + 1) * sizeof(u4i);
	mpsize += mlen * sizeof(u4i);
	mpsize += 5 * mlen * sizeof(u2i);
	mpsize += mlen * W * sizeof(u1i);
	mpsize += 2 * (mlen + W) * sizeof(u1i);
	clear_and_encap_b1v(g->memp, mpsize);
	memset(g->memp->buffer, 0, mpsize);
	rows[0] = (int*)g->memp->buffer;
	rows[1] = rows[0] + W + 1;
	breps   = (u4i*)(rows[1] + W + 1);
	bcnts   = (u2i*)(breps + mlen);
	matrix   = (u1i*)(bcnts + 5 * mlen);
	seqs[0]  = matrix + mlen * W;
	seqs[1]  = seqs[0] + mlen + W;
	memset(seqs[0], 4, HW * sizeof(u1i));
	seqs[0] += HW;
	memset(seqs[0] + mlen, 4, HW * sizeof(u1i));
	memset(seqs[1], 4, HW * sizeof(u1i));
	seqs[1] += HW;
	memset(seqs[1] + mlen, 4, HW * sizeof(u1i));
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		seqs[1][pos] = col[nall];
		bs = bcnts + pos * 5;
		for(rid=0;rid<nall;rid++){
			if(col[rid] < 4){
				bs[col[rid]] ++;
			}
		}
	}
	//adjust [cns=4] minor bases to right-side identitical cns
	for(pos=0;pos<mlen;pos++){
		if((lc = seqs[1][pos]) < 4){
			for(i=pos;i;i--){
				if(seqs[1][i-1] < 4) break;
				if(bcnts[(i - 1) * 5 + lc]){
					bcnts[pos * 5 + lc] += bcnts[(i - 1) * 5 + lc];
					bcnts[(i - 1) * 5 + lc] = 0;
				}
			}
		}
	}
	// adjust bases and cnts within a cns-homopolymer
	lc = 4; mc = 0; memset(cnts, 0, 4 * sizeof(u4i));
	clear_u4v(g->stack);
	col = NULL;
	for(lpos=pos=0;pos<=mlen;pos++){
		if(pos == mlen || ((col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow)[nall] < 4 && col[nall] != lc)){
			sort_array(g->stack->buffer, g->stack->size, hp_base_t, num_cmpgt(a.base, b.base));
			for(i=p=0;i<=g->stack->size;i++){
				if(i < g->stack->size && ((hp_base_t)g->stack->buffer[i]).base == ((hp_base_t)g->stack->buffer[p]).base){
					continue;
				}
				cc = cnts[TOTYPE(g->stack->buffer[p], hp_base_t).base];
				if(TOTYPE(g->stack->buffer[p], hp_base_t).base == lc){
					sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgtx(b.cns, a.cns, a.off, b.off));
					for(j=p;cc&&j<i;j++){
						PB = TOTYPE(g->stack->buffer[j], hp_base_t);
						if(PB.cns == 0) break;
						bc = num_min(cc, mc);
						bcnts[(lpos + PB.off) * 5 + PB.base] = bc;
						cc -= bc;
					}
					while(p < j){
						PB = TOTYPE(g->stack->buffer[p], hp_base_t);
						bcnts[(lpos + PB.off) * 5 + PB.base] += (j - p);
						p ++;
					}
					p = j;
					sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgt(b.bcnt, a.bcnt));
					for(j=p;cc&&j<i;j++){
						PB = TOTYPE(g->stack->buffer[j], hp_base_t);
						bc = num_min(cc, mc);
						bcnts[(lpos + PB.off) * 5 + PB.base] = bc;
						cc -= bc;
					}
				} else {
					sort_array(g->stack->buffer + p, i - p, hp_base_t, num_cmpgt(b.bcnt, a.bcnt));
					for(j=p;cc&&j<i;j++){
						PB = TOTYPE(g->stack->buffer[j], hp_base_t);
						bc = num_min(cc, mc);
						bcnts[(lpos + PB.off) * 5 + PB.base] = bc;
						cc -= bc;
					}
				}
				p = i;
			}
		}
		if(pos == mlen) break;
		if(col[nall] < 4 && col[nall] != lc){
			lc = col[nall]; mc = 0;
			memset(cnts, 0, 4 * sizeof(u4i));
			lpos = pos;
			clear_u4v(g->stack);
		}
		bs = bcnts + pos * 5;
		for(i=0;i<4;i++){
			if(bs[i]){
				if(bs[i] > mc) mc = bs[i];
				if(i == lc){
					cnts[i] += bs[i];
				} else {
					//cnts[i] += bs[i] - 1; // descrease non-cns bases cnt
					cnts[i] += bs[i];
				}
				PB.base = i;
				PB.cns  = (i == col[nall]);
				PB.off  = pos - lpos;
				PB.bcnt = bs[i];
				push_u4v(g->stack, PB.val);
				bs[i] = 0;
			}
		}
	}
#if 0
	for(lc=0;lc<4;lc++){
		if(1){
			printf("[%c]              ", "ACGT"[lc]);
			//for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i] / 100) % 10);
			for(i=0;i<mlen;i++) printf("%d", (bcnts[i * 5 + lc] / 100) % 10);
			printf("\n");
			printf("[%c]              ", "ACGT"[lc]);
			//for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i] / 10) % 10);
			for(i=0;i<mlen;i++) printf("%d", (bcnts[i * 5 + lc] / 10) % 10);
			printf("\n");
			printf("[%c]              ", "ACGT"[lc]);
			//for(i=0;i<mlen;i++) printf("%d", (mats[1][lc][i]) % 10);
			for(i=0;i<mlen;i++) printf("%d", (bcnts[i * 5 + lc]) % 10);
			printf("\n");
			printf("\n");
		}
	}
	fflush(stdout);
	exit(1);
#endif
	for(rid=0;rid<nseq;rid++){
		rdlen = g->seqs->rdlens->buffer[rid];
		if(rdlen == 0) continue;
		//if(rid != 1) continue;
		memset(seqs[0], 4, mlen * sizeof(u1i));
		memset(breps, 0, mlen * sizeof(u4i));
		lc = 4; cc = 0;
		for(i=rdlen;i;i--){
			v = cut_rdnode_bspoa(g, rid, i - 1, BSPOA_RDNODE_CUTALL);
			seqs[0][v->mpos] = v->base;
			if(v->base == lc){
				if(cc < 0xf) cc ++;
				breps[v->mpos] = cc << (v->base << 2);
			} else {
				lc = v->base;
				cc = 0;
			}
		}
		v = get_rdnode_bspoa(g, rid, 0);
		mbeg = v->mpos;
		v = get_rdnode_bspoa(g, rid, rdlen - 1);
		mend = v->mpos + 1;
		scr = remsa_edit_rd_bspoacore(g, rid, rdlen, bcnts, breps, seqs, rows, matrix, mbeg, mend, W);
		//fflush(stdout); fprintf(stderr, " -- R%d = %d in %s -- %s:%d --\n", rid, scr, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		connect_rdnodes_bspoa(g, rid);
	}
	del_msanodes_bspoa(g);
}

static inline void end_bspoa(BSPOA *g){
	seqalign_result_t rs;
	int i;
	u2i rid;
	clear_u1v(g->cns);
	clear_u1v(g->qlt);
	clear_u1v(g->alt);
	if(g->par->refmode){
		resize_u1v(g->cns, g->seqs->rdlens->buffer[0]);
		resize_u1v(g->qlt, g->seqs->rdlens->buffer[0]);
		resize_u1v(g->alt, g->seqs->rdlens->buffer[0]);
		bitseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[0], g->seqs->rdlens->buffer[0], g->cns->buffer);
		memset(g->qlt->buffer, 0, g->seqs->rdlens->buffer[0]);
		memset(g->alt->buffer, 0, g->seqs->rdlens->buffer[0]);
	}
	if(g->seqs->nseq <= 1){
		return;
	}
	if(g->par->shuffle) shuffle_reads_by_kmers_bspoa(g);
	if(g->par->seqcore){
		g->nmsa = num_min(g->seqs->nseq, g->par->seqcore);
	} else {
		g->nmsa = g->seqs->nseq;
	}
	for(rid=0;rid<g->seqs->nseq;rid++){
		_add_read_bspoa_core(g, rid);
	}
	g->nrds = 1;
	for(rid=1;rid<g->nmsa;rid++){
		if(!g->par->refmode && g->par->bwtrigger){
			msa_bspoa(g);
			simple_cns_bspoa(g);
			if(_DEBUG_LOG_ > 1){
				print_msa_bspoa(g, __FUNCTION__, 0, 0, 0, 1, _DEBUG_LOGFILE_);
			}
		}
		rs = align_rd_bspoa(g, g->par, 0, rid, 0, g->seqs->rdlens->buffer[rid]);
		g->nrds ++;
	}
	// generate MSA and CNS
	for(i=0;i<g->par->realn;i++){
		msa_bspoa(g);
		//simple_cns_bspoa(g);
		cns_bspoa(g);
		//print_msa_bspoa(g, "BSALIGN", 0, 0, 0, 1, stdout);
		if(g->par->editbw < 0){ // no parallel mode
			remsa_edits_bspoa(g, - g->par->editbw);
		} else {
			remsa_pedits_bspoa(g, g->par->editbw / 2, 1, (i + 1 == g->par->realn));
		}
	}
	if(g->par->shuffle) restore_rd_orders_bspoa(g);
	//print_probs = 1;
	msa_bspoa(g);
	cns_bspoa(g);
	//print_probs = 0;
}

static inline void fix_tenon_mortise_msa_bspoa2(BSPOA *g){
	bspoanode_t *v;
	u4i rid, nseq, mrow, mlen, i, pos, *rps;
	u2i pairs[2][4], cut;
	u1i *cols[2], pmaxs[2];
	nseq = g->seqs->nseq;
	mrow = g->seqs->nseq + 3;
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

static inline void tidy_msa_bspoa(BSPOA *g){
	u4i nseq, nall, realnseq, mrow, mlen, i, pos, lpos, e, f;
	u2i bcnts[6], rid, gap, lst;
	u1i *col, *bss, qlt, alt, lc, lq, m1, m2;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	realnseq = (nseq && g->seqs->rdlens->size && g->seqs->rdlens->buffer[0])? nseq : nseq - 1;
	mrow = nall + 3;
	mlen = g->msaidxs->size;
	void find_top2_bases(int calc_n){
		memset(bcnts, 0, 6 * sizeof(u2i));
		for(rid=0;rid<nseq;rid++){
			if(col[rid] <= 4){
				bcnts[col[rid]] ++;
				bcnts[5] ++;
			}
		}
		if(bcnts[0] >= bcnts[1]){
			m1 = 0; m2 = 1;
		} else {
			m1 = 1; m2 = 0;
		}
		for(i=2;i<(calc_n? 5 : 4);i++){
			if(bcnts[i] > bcnts[m1]){
				m2 = m1;
				m1 = i;
			} else if(bcnts[i] > bcnts[m2]){
				m2 = i;
			}
		}
	}
	bss = NULL; lst = 0; lc = 4; lq = 0; lpos = MAX_U4;
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		qlt = col[nall + 2];
		//if(pos == 1494){
			//fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		//}
		if(qlt < g->par->althi) continue;
		find_top2_bases(1);
		alt = 4; gap = 0;
		if(m1 == 4 && bcnts[m2]){
			gap = bcnts[m1];
			alt = m2;
		} else if(m2 == 4 && bcnts[m1]){
			alt = m1;
			gap = bcnts[m2];
		}
		if(alt == 4 || gap < Int(0.1 * bcnts[5]) || bcnts[alt] < Int(0.1 * bcnts[5])) continue;
		if(lpos == MAX_U4){
		} else if(alt == lc){
			if(qlt < lq) continue;
		} else if(lst >= Int(0.75 * bcnts[alt])){
			e = lpos;
			while(e < pos){
				bss = g->msacols->buffer + g->msaidxs->buffer[e + 1] * mrow;
				if(bss[nall] < 4 && bss[nall] != lc) break;
				e ++;
			}
			f = pos;
			while(f > e){
				bss = g->msacols->buffer + g->msaidxs->buffer[f - 1] * mrow;
				if(bss[nall] < 4 && bss[nall] != alt) break;
				f --;
			}
			if(e + 1 >= f){
				bss = g->msacols->buffer + g->msaidxs->buffer[lpos] * mrow;
				for(e=f=rid=0;rid<nseq;rid++){
					if(col[rid] < 4 && bss[rid] < 4 && col[rid] != bss[rid]) e ++; // if col[rid] == bss[rid], maybe homopolymer insertion
					else if(col[rid] < 4 || bss[rid] < 4) f ++;
				}
				//fprintf(stdout, "TENON_MORTISE\t%d\t%d\tfixed=%d\tunfix=%d\n", lpos, pos, f, e);
				if(f >= Int(0.75 * bcnts[5]) && e <= Int(0.40 * num_min(gap, lst))){
					//FIX[0].pos  = lpos;
					//FIX[0].fix  = f;
					//FIX[0].base = lc;
					//FIX[1].pos  = pos;
					//FIX[1].fix  = f;
					//FIX[1].base = alt;
					//push_u8v(g->heap, ((u8i*)FIX)[0]);
					//push_u8v(g->heap, ((u8i*)FIX)[1]);
					for(rid=0;rid<nseq;rid++){
						if(col[rid] == alt && bss[rid] == 4){
							bss[rid] = alt;
							col[rid] = 4;
						}
					}
					lpos = MAX_U4; lst = 0; lc = 4; lq = 0;
					continue;
				}
			}
		}
		lpos = pos; lst = gap; lc = alt; lq = qlt;
	}
	cns_bspoa(g);
}

static inline void call_snvs_bspoa(BSPOA *g){
	bspoavar_t *var;
	u4i nseq, nall, realnseq, mrow, mlen, i, pos, lpos;
	u2i bcnts[6], rid, mincov, *acnts;
	u1i *col, m1, m2;
	u4i pcnt, j, k, covn, altn, acnt, qual;
	float pmax, pmin, pinc, perr, pexp, prob, *psums;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	realnseq = (nseq && g->seqs->rdlens->size && g->seqs->rdlens->buffer[0])? nseq : nseq - 1;
	mrow = nall + 3;
	mlen = g->msaidxs->size;
	void find_top2_bases(int calc_n){
		memset(bcnts, 0, 6 * sizeof(u2i));
		for(rid=0;rid<nseq;rid++){
			if(col[rid] <= 4){
				bcnts[col[rid]] ++;
				bcnts[5] ++;
			}
		}
		if(bcnts[0] >= bcnts[1]){
			m1 = 0; m2 = 1;
		} else {
			m1 = 1; m2 = 0;
		}
		for(i=2;i<(calc_n? 5 : 4);i++){
			if(bcnts[i] > bcnts[m1]){
				m2 = m1;
				m1 = i;
			} else if(bcnts[i] > bcnts[m2]){
				m2 = i;
			}
		}
	}
	// calc probs
	pinc = 0.0005;
	pcnt = 100;
	pmax = pinc * pcnt;
	pmin = 0.01; // min prob to be considered in calculation
	acnt = realnseq * (realnseq / 2 + 1);
	resize_b1v(g->memp, pcnt * sizeof(float) + acnt * sizeof(u2i));
	zeros_b1v(g->memp);
	psums = (float*)g->memp->buffer;
	acnts = (u2i*)(psums + pcnt);
	mincov = num_max(2, realnseq * g->par->min_covfrq);
	for(lpos=pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		find_top2_bases(0);
		if(bcnts[m1] + bcnts[m2] >= mincov){
			covn = bcnts[5];
			altn = bcnts[m2];
			acnts[altn * realnseq + covn - 1] ++;
		}
	}
	for(i=1;i<acnt;i++){
		if(acnts[i] == 0) continue;
		altn = i / realnseq;
		covn = (i % realnseq) + 1;
		pexp = 1.0 * altn / covn;
		j = Int(pexp / pinc);
		if(j > 0 && j < pcnt){
			for(k=0;k<j;k++){
				perr = pexp - pinc * k;
				prob = exp(cal_binomial_bspoa(covn, altn, perr));
				psums[j - k] += acnts[i] * prob;
				if(prob <= pmin) break;
			}
			for(k=1;j+k<pcnt;k++){
				perr = pexp + pinc * k;
				prob = exp(cal_binomial_bspoa(covn, altn, perr));
				psums[j + k] += acnts[i] * prob;
				if(prob <= pmin) break;
			}
		}
	}
	perr = 1; pexp = 0.01;
	for(i=0;i<pcnt;i++){
		if(perr < psums[i]){
			pexp = i * pinc;
			perr = psums[i];
		}
	}
	if(_DEBUG_LOG_){
		fprintf(stdout, "[MSA_SNV_ERR_PROB:%.4f]", pexp);
		for(i=0;i<pcnt;i++){
			fprintf(stdout, " [%.4f]:%.4f", pinc * i, psums[i]);
		}
		fprintf(stdout, "\n");
	}
	// call SNVs
	clear_bspoavarv(g->var);
	for(lpos=pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		find_top2_bases(0);
		//if(m1 < 4 && m2 < 4 && bcnts[m1] + bcnts[m2] >= UInt(realnseq * g->par->min_covfrq) && bcnts[m2] >= g->par->min_varcnt && bcnts[m2] >= (u4i)(bcnts[5] * g->par->min_varfrq)){
		if(m1 < 4 && m2 < 4 && bcnts[m2] >= g->par->min_varcnt && bcnts[m1] + bcnts[m2] >= mincov){
			qual = - ((prob = cal_binomial_bspoa(bcnts[5], bcnts[m2], pexp)) / log(10));
			if(qual > 1000) qual = 1000;
			if(qual >= g->par->min_snvqlt){
				var = next_ref_bspoavarv(g->var);
				var->cpos = lpos;
				var->mpos = pos;
				var->covn = bcnts[5];
				var->refn = bcnts[m1];
				var->refb = m1;
				var->altn = bcnts[m2];
				var->altb = m2;
				var->qual = qual;
			}
		}
		if(col[nall] < 4) lpos ++;
	}
	if(0){
		for(i=0;i<g->var->size;i++){
			var = ref_bspoavarv(g->var, i);
			fprintf(stdout, " -- pos=%d,%d cov=%d,%d m1=%c:%d m2=%c:%d qual=%d --\n", var->mpos, var->cpos, realnseq, var->covn, "ACGT-"[var->refb], var->refn, "ACGT-"[var->altb], var->altn, var->qual);
		}
	}
}

#define call_snv_bspoa(g) call_snvs_bspoa(g)

static inline void print_snvs_bspoa(BSPOA *g, char *label, FILE *out){
	bspoavar_t *var;
	u4i i, j, fsz, fct, rid, nseq, mrow, mlen;
	u1i *col;
	char flanks[4][6], *genotypes;
	fsz = 5;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	clear_string(g->strs);
	encap_string(g->strs, nseq);
	genotypes = g->strs->string;
	for(i=0;i<g->var->size;i++){
		var = ref_bspoavarv(g->var, i);
		fct = num_min(var->cpos, fsz);
		memcpy(flanks[0], g->cns->buffer + var->cpos - fct, fct);
		for(j=0;j<fct;j++) flanks[0][j] = bit_base_table[(int)flanks[0][j]];
		flanks[0][fct] = 0;
		memcpy(flanks[2], g->qlt->buffer + var->cpos - fct, fct);
		for(j=0;j<fct;j++) flanks[2][j] += '!';
		flanks[2][fct] = 0;
		fct = num_min(g->cns->size - var->cpos - 1, fsz);
		memcpy(flanks[1], g->cns->buffer + var->cpos + 1, fct);
		for(j=0;j<fct;j++) flanks[1][j] = bit_base_table[(int)flanks[1][j]];
		flanks[1][fct] = 0;
		memcpy(flanks[3], g->qlt->buffer + var->cpos + 1, fct);
		for(j=0;j<fct;j++) flanks[3][j] += '!';
		flanks[3][fct] = 0;
		fprintf(out, "%s SNP\t", label);
		col = g->msacols->buffer + g->msaidxs->buffer[var->mpos] * mrow;
		for(rid=0;rid<nseq;rid++){
			genotypes[rid] = "ACGT-.*"[(int)col[rid]];
		}
		genotypes[nseq] = '\0';
		fprintf(out, "%d\t%d\t%s\t%s\t%c\t%d\t%c\t%d\t%s\t%s\t%d\t%d\t%s\n", var->cpos, var->mpos, flanks[0], flanks[2], bit_base_table[var->refb], var->refn, bit_base_table[var->altb], var->altn, flanks[1], flanks[3], var->covn, var->qual, genotypes);
	}
}

static inline void check_msa_bspoa(BSPOA *g){
	u4i nseq, mrow, mlen, rid, pos, len;
	u1i c;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
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
		xidx = v->next;
		while(xidx != node){
			x = ref_bspoanodev(g->nodes, xidx);
			ncov ++;
			xidx = x->next;
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

#define BSPOA_HSP_MINLEN	3 // MUST >= 3
// 1, base not match cns
// 2, node->cov > 1, means it indeed affact the POA process
// 3, OR nonbase but cns base
// 4, edge->cov > 1, means it indeed affact the POA process
// 5, three nonplolymer bases as boundaray
static inline u4i gen_rd_lsp_bspoa(BSPOA *g, BSPOAPar *par, u2i rid){
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	bspoalsp_t *lsp;
	int hsp[6], tot;
	float score;
	u4i nseq, mrow, mlen, rdlen, rbeg, rpos, pos, ret;
	u1i *qs, q, c, ld, lc, f, state;
	UNUSED(par);
	nseq = g->nrds;
	mrow = nseq + 3;
	mlen = g->msaidxs->size;
	rdlen = g->seqs->rdlens->buffer[rid];
	ret = 0;
	// find all low scoring pairs (LSP)
	score = 0;
	rbeg = rpos = 0;
	memset(hsp, 0, 6 * sizeof(int)); // 0:lst hsp rpos, 1:lst hsp mpos, 2:tmp hsp rpos, 3:tmp hsp mpos, 4:hsp len, 5:lsp cnt
	lc = 4; ld = 0;
	state = 0;
	tot = 0;
	for(pos=0;pos<mlen;pos++){
		qs = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		q = qs[rid]  < 4? qs[rid]  : 4;
		c = qs[nseq] < 4? qs[nseq] : 4;
		f = g->dptable[c + q * 5 + lc * 25 + ld * 125]; // Maybe I need to reverse pos to get more accurate probs
		score += g->dpvals[f >> 3] * qs[nseq + 1]; // probs * cns_phred_quality
		ld = f & 0x7;
		if(q == c){
			if(q == 4){
				continue;
			}
			if(state == 0){
				state = 2;
			}
			if(q != lc){
				if(state == 2){
					hsp[4] ++;
					if(hsp[4] == 2){
						hsp[2] = rpos;
						hsp[3] = pos;
					} else if(hsp[4] == BSPOA_HSP_MINLEN){
						state = 1;
						if(hsp[0] && hsp[5]){
							if((lsp = tail_bspoalspv(g->lsp)) && lsp->rid == rid && hsp[0] <= Int(lsp->rbeg + lsp->rlen + 2)){
								lsp->scr += score;
								lsp->rlen = hsp[2] - lsp->rbeg;
								lsp->mlen = hsp[3] - lsp->mbeg;
								tot += hsp[2] - hsp[0];
								if(_DEBUG_LOG_ > 1){
									fprintf(_DEBUG_LOGFILE_, "LSP: MERGED RID=%d\tscr=%0.4f\trbeg=%d\trlen=%d\tmbeg=%d\tmlen=%d\t", lsp->rid, lsp->scr, lsp->rbeg, lsp->rlen, lsp->mbeg, lsp->mlen);
									println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[lsp->rid] + lsp->rbeg, lsp->rlen, _DEBUG_LOGFILE_);
								}
							} else {
								ret ++;
								lsp = next_ref_bspoalspv(g->lsp);
								lsp->rid   = rid;
								lsp->scr   = score;
								lsp->rbeg  = hsp[0];
								lsp->rlen  = hsp[2] - hsp[0];
								lsp->mbeg  = hsp[1];
								lsp->mlen  = hsp[3] - hsp[1];
								tot += lsp->rlen;
								if(_DEBUG_LOG_ > 1){
									fprintf(_DEBUG_LOGFILE_, "LSP: RID=%d\tscr=%0.4f\trbeg=%d\trlen=%d\tmbeg=%d\tmlen=%d\t", lsp->rid, lsp->scr, lsp->rbeg, lsp->rlen, lsp->mbeg, lsp->mlen);
									println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[lsp->rid] + lsp->rbeg, lsp->rlen, _DEBUG_LOGFILE_);
								}
							}
						}
					}
				}
				if(state == 1){
					hsp[2] = rpos;
					hsp[3] = pos;
					hsp[5] = 0;
					score = 0;
				}
			}
		} else {
			hsp[4] = 0;
			if(state == 1){
				hsp[0] = hsp[2];
				hsp[1] = hsp[3];
			}
			state = 0;
			v = get_rdnode_bspoa(g, rid, rpos);
			if(q < 4){
				if(get_nodecov_bspoa(g, v) > 1){
					hsp[5] ++;
				}
			} else {
				w = get_rdnode_bspoa(g, rid, rpos - 1);
				e = get_edge_bspoa(g, w, v);
				if(e && e->cov > 1){
					hsp[5] ++;
				}
			}
		}
		if(c < 4) lc = c;
		if(q < 4){
			rpos ++;
			if(rpos >= rdlen) break;
		}
	}
	if(_DEBUG_LOG_){
		fprintf(_DEBUG_LOGFILE_, "LSP RID[%d] cnt=%d\ttot=%d\trdlen=%d\n", rid, ret, tot, rdlen);
	}
	return ret;
}

static inline float cal_rd_lsp_score_bspoa(BSPOA *g, bspoalsp_t *lsp){
	float scr;
	u4i nseq, nall, mrow, p, pos, a, b, c, d, f;
	u1i *col;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	c = 4;
	d  = 0;
	scr = 0;
	for(p=0;p<lsp->mlen;p++){
		pos = p + lsp->mbeg;
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		a = col[nall];
		b = col[lsp->rid];
		if(b > 4) continue;
		if(a >= 4 && b >= 4) continue;
		f = g->dptable[a + b * 5 + c * 25 + d * 125];
		scr += g->dpvals[f >> 3];
		d = f & 0x7;
		if(col[nall] < 4) c = col[nall];
	}
	return - scr;
}

static inline u4i gen_lsps_bspoa(BSPOA *g, BSPOAPar *par){
	bspoalsp_t *lsp;
	bspoanode_t *v, *w;
	bspoaedge_t *e;
	u8i mpsize;
	u4i nseq, nall, mrow, ret, wsz;
	u4i i, pos, rid, *roffs, *hsps[8], LSP, cnts[6], idxs[6], revs[6];
	u1i *mem, *col, *states, lc, q, a, b, c, x;
	wsz = 5;
	if(g->msaidxs->size < wsz) return 0;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
	mpsize = 0;
	mpsize += roundup_times((1 + 8) * nseq * sizeof(u4i) + nseq * sizeof(u1i), WORDSIZE);
	clear_and_encap_b1v(g->memp, mpsize);
	mem = (u1i*)g->memp->buffer;
	memset(mem, 0, mpsize);
	roffs = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[0] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[1] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[2] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[3] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[4] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[5] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[6] = (u4i*)mem; mem += nseq * sizeof(u4i);
	hsps[7] = (u4i*)mem; mem += nseq * sizeof(u4i);
	states  = (u1i*)mem; mem += nseq * sizeof(u1i);
	memset(states, 1, nseq * sizeof(u1i));
	ret = 0;
	lc = 4;
	for(pos=0;pos<g->msaidxs->size;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		q   = col[nall + 1];
		a   = col[nall + 2];
		c   = num_min(col[nall], 4);
		if(q < par->qlthi || a > par->althi){
			LSP = 1;
			memset(cnts, 0, 6 * sizeof(u4i));
			for(rid=0;rid<nseq;rid++){
				cnts[col[rid]] ++;
			}
			if(0){
				for(i=0;i<6;i++) idxs[i] = i;
				sort_array(idxs, 6, u4i, num_cmpgt(cnts[b], cnts[a]));
				for(i=0;i<6;i++) revs[idxs[i]] = (i + 1);
				for(rid=0;rid<nseq;rid++){
					//if(col[rid] == c) continue;
					if(hsps[6][rid] > 1){
						hsps[6][rid] += revs[col[rid]];
					} else {
						hsps[6][rid] = 1 + revs[col[rid]];
					}
				}
			} else {
				for(rid=0;rid<nseq;rid++){
					if(col[rid] == c){
						if(hsps[6][rid] == 0) hsps[6][rid] = 1;
					} else if(hsps[6][rid] > 1){
						hsps[6][rid] += cnts[col[rid]];
					} else {
						hsps[6][rid] = 1 + cnts[col[rid]];
					}
				}
			}
		} else {
			LSP = 0;
		}
		for(rid=0;rid<nseq;rid++){
			b = num_min(col[rid], 4);
			if(b == 4 && c == 4) continue;
			if(LSP || b != c){
				if(!LSP){
					x = 0;
					v = get_rdnode_bspoa(g, rid, roffs[rid]);
					if(b < 4){
						if(get_nodecov_bspoa(g, v) > 1){
							if(c < 4){
								x = 1;
							} else {
								do {
									w = get_rdnode_bspoa(g, rid, roffs[rid] - 1);
									if(w->base == b) break;
									w = get_rdnode_bspoa(g, rid, roffs[rid] + 1);
									if(w->base == b) break;
									x = 1;
								} while(0);
							}
						}
					} else {
						w = get_rdnode_bspoa(g, rid, roffs[rid] - 1);
						e = get_edge_bspoa(g, w, v);
						if(e && e->cov > 1){
							x = 1;
						}
					}
					if(x == 1){
						hsps[5][rid] ++;
						hsps[4][rid] = 0;
					}
				}
				if(LSP || (b != c && b != 4 && c != 4)){
					if(b != c && b != 4 && c != 4 && hsps[6][rid] == 0){
						hsps[6][rid] = 1;
					}
					hsps[4][rid] = 0;
					if(states[rid] == 1){
						hsps[0][rid] = hsps[2][rid];
						hsps[1][rid] = hsps[3][rid];
					}
					states[rid] = 0;
				}
				hsps[7][rid] = 1;
			} else {
				if(states[rid] == 0) states[rid] = 2;
				if(c != lc){
					hsps[4][rid] ++;
					if(states[rid] == 2){
						if(hsps[4][rid] >= 2 && hsps[7][rid] == 0){
							hsps[2][rid] = roffs[rid];
							hsps[3][rid] = pos;
							hsps[5][rid] = 0;
							states[rid] = 3;
						}
					}
					if(states[rid] == 3){
						if(hsps[4][rid] >= wsz){
							states[rid] = 1;
							if(hsps[0][rid] && hsps[6][rid]){
								ret ++;
								lsp = next_ref_bspoalspv(g->lsp);
								lsp->rid   = rid;
								lsp->rbeg  = hsps[0][rid];
								lsp->rlen  = hsps[2][rid] - hsps[0][rid];
								lsp->mbeg  = hsps[1][rid];
								lsp->mlen  = hsps[3][rid] - hsps[1][rid];
								lsp->scr   = hsps[6][rid] * 100 + hsps[5][rid];
								lsp->scr   += cal_rd_lsp_score_bspoa(g, lsp);
								if(_DEBUG_LOG_ > 1){
									fprintf(_DEBUG_LOGFILE_, "LSP: RID=%d\tscr=%0.4f\trbeg=%d\trlen=%d\tmbeg=%d\tmlen=%d\t", lsp->rid, lsp->scr, lsp->rbeg, lsp->rlen, lsp->mbeg, lsp->mlen);
									println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[lsp->rid] + lsp->rbeg, lsp->rlen, _DEBUG_LOGFILE_);
								}
							}
						}
					}
					if(states[rid] == 1){
						if(hsps[4][rid] >= wsz && hsps[7][rid] == 0){
							hsps[2][rid] = roffs[rid];
							hsps[3][rid] = pos;
						}
						hsps[5][rid] = 0;
						hsps[6][rid] = 0;
					}
				}
				hsps[7][rid] = 0;
			}
			if(b < 4) roffs[rid] ++;
		}
		if(c < 4) lc = c;
	}
	return ret;
}

static inline void remsa_lsps_bspoa(BSPOA *g, BSPOAPar *par){
	bspoalsp_t *lsp;
	bspoanode_t *v;
	u4i i, pos, mrow, *roffs;
	u2i rid, nseq, nall;
	u1i *col, lc;
	nseq = g->nrds;
	nall = (g->seqs->nseq == 0)? nseq : g->seqs->nseq;
	mrow = nall + 3;
#if DEBUG
	check_rdnodes_bspoa(g);
	check_msa_rdseqs_bspoa(g);
	check_aligned_nodes_bspoa(g);
	check_nodecovs_bspoa(g);
#endif
	for(i=0;i<g->nodes->size;i++){
		v = ref_bspoanodev(g->nodes, i);
		v->colorful = 0;
	}
	clear_and_encap_b1v(g->memp, nseq * sizeof(u4i));
	roffs = (u4i*)g->memp->buffer;
	add_cnsnodes_bspoa(g, roffs);
	// set blessed nodes
	if(1){
		for(i=0;i<g->cns->size;i++){
			v = get_rdnode_bspoa(g, nseq, i);
			v->bless = 1;
		}
	} else {
		memset(roffs, 0, nseq * sizeof(u4i));
		lc = 4;
		for(pos=0;pos<g->msaidxs->size;pos++){
			col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			for(rid=0;rid<nseq;rid++){
				if(col[rid] >= 4) continue;
				if(col[rid] == col[nall] && col[nall] != lc){
					v = get_rdnode_bspoa(g, rid, roffs[rid]);
					v->bless = 1;
				}
				roffs[rid] ++;
			}
			if(col[nall] < 4) lc = col[nall];
		}
	}
	clear_bspoalspv(g->lsp);
	if(1){
		gen_lsps_bspoa(g, par);
	} else {
		for(rid=0;rid<nseq;rid++){
			gen_rd_lsp_bspoa(g, par, rid);
		}
	}
	u4i mkey = MAX_U4;
	if(mkey != MAX_U4){
		u4i j;
		for(i=j=0;i<g->lsp->size;i++){
			lsp = ref_bspoalspv(g->lsp, i);
			if(lsp->mbeg < mkey && lsp->mbeg + lsp->mlen > mkey){
				g->lsp->buffer[j ++] = *lsp;
			}
		}
		g->lsp->size = j;
	}
	sort_list(g->lsp, num_cmpgt(a.scr, b.scr));
	for(i=0;i<g->lsp->size;i++){
		lsp = ref_bspoalspv(g->lsp, i);
		for(pos=lsp->rbeg;pos<lsp->rbeg+lsp->rlen;pos++){
			v = cut_rdnode_bspoa(g, lsp->rid, pos, BSPOA_RDNODE_CUTEDGE);
			// no bless to the lsp node
			v->bless = 0;
			v->colorful = 1;
		}
	}
	if(_DEBUG_LOG_){
		print_msa_bspoa(g, __FUNCTION__, 0, 0, 0, 1, _DEBUG_LOGFILE_);
	}
	for(i=0;i<g->lsp->size;i++){
		lsp = ref_bspoalspv(g->lsp, i);
		if(_DEBUG_LOG_){
			fprintf(_DEBUG_LOGFILE_, "LSP: RID=%d\tscr=%0.4f\trbeg=%d\trlen=%d\tmbeg=%d\tmlen=%d\t", lsp->rid, lsp->scr, lsp->rbeg, lsp->rlen, lsp->mbeg, lsp->mlen);
			println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[lsp->rid] + lsp->rbeg, lsp->rlen, _DEBUG_LOGFILE_);
		}
		if(lsp->rid == 4 && lsp->rbeg == 1380){
			fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		for(pos=lsp->rbeg;pos<lsp->rbeg+lsp->rlen;pos++){
			v = cut_rdnode_bspoa(g, lsp->rid, pos, BSPOA_RDNODE_CUTNODE);
		}
		align_rd_bspoa(g, par, 1, lsp->rid, lsp->rbeg, lsp->rlen);
	}
	del_cnsnodes_bspoa(g);
	msa_bspoa(g);
	cns_bspoa(g);
	if(_DEBUG_LOG_){
		print_msa_bspoa(g, __FUNCTION__, 0, 0, 0, 1, _DEBUG_LOGFILE_);
	}
#if DEBUG
	check_rdnodes_bspoa(g);
	check_msa_rdseqs_bspoa(g);
	check_aligned_nodes_bspoa(g);
	check_nodecovs_bspoa(g);
#endif
}

static inline int revise_seq_joint_point(u4v *cigars, int *qe, int *te){
	u4i i, op, ln, max;
	int qq, q, tt, t;
	q = t = 0;
	qq = tt = 0;
	max = 0;
	for(i=1;i<=cigars->size;i++){
		op = cigars->buffer[cigars->size - i] & 0xF;
		ln = cigars->buffer[cigars->size - i] >> 4;
		if(op == 0){
			if(ln > max){
				qq = q; tt = t;
				max = ln;
			}
			q += ln;
			t += ln;
		} else if(op == 1){
			q += ln;
		} else {
			t += ln;
		}
	}
	*qe -= qq;
	*te -= tt;
	return 1;
}

// ret[0]: revised end of seq1
// ret[1]  revised beg of seq2
static inline seqalign_result_t cat_cns_seqs(int ret[2], u1v *seq1, u1v *seq2, int overlap, b1v *mempool, u4v *cigars, int M, int X, int O, int E){
	seqalign_result_t RS;
	b1i matrix[16];
	int qb, qe, tb, te, ol, maxl;
	if(seq1->size == 0 || seq2->size == 0){
		ZEROS(&RS);
		ret[0] = seq1->size;
		ret[1] = 0;
		return RS;
	}
	qb = 0; qe = seq1->size;
	tb = 0; te = seq2->size;
	if(qe > overlap) qb = qe - overlap;
	if(te > overlap) te = overlap;
	ol = num_min(qe - qb, te - tb);
	banded_striped_epi8_seqalign_set_score_matrix(matrix, M, X);
	RS = banded_striped_epi8_seqalign_pairwise(seq1->buffer + qb, qe - qb, seq2->buffer + tb, te - tb, mempool, cigars, SEQALIGN_MODE_OVERLAP, 0, matrix, O, E, 0, 0, 0);
	if(RS.aln < Int(0.5 * overlap) || RS.mat < Int(RS.aln * 0.9)){
		// full length alignment
		maxl = num_min(seq1->size, seq2->size);
		maxl = num_min(maxl, overlap * 4);
		qb = 0; qe = seq1->size;
		tb = 0; te = seq2->size;
		if(qe > maxl) qb = qe - maxl;
		if(te > maxl) te = maxl;
		ol = num_min(qe - qb, te - tb);
		RS = banded_striped_epi8_seqalign_pairwise(seq1->buffer + qb, qe - qb, seq2->buffer + tb, te - tb, mempool, cigars, SEQALIGN_MODE_OVERLAP, 0, matrix, O, E, 0, 0, 0);
	}
	RS.qb += qb;
	RS.qe += qb;
	RS.tb += tb;
	RS.te += tb;
	ret[0] = RS.qe;
	ret[1] = RS.te;
	revise_seq_joint_point(cigars, ret + 0, ret + 1);
	return RS;
}

// TODO: debug MERGE MSA

static inline void beg_merge_bspoa(BSPOA *dg, BSPOA *mg){
	beg_bspoa(dg);
	beg_bspoa(mg);
	clear_u4v(dg->sels); // cns2rdnode
	push_u4v(dg->sels, 0);
	push_u4v(dg->sels, 0);
}

static inline void push_merge_bspoa(BSPOA *dg, BSPOA *mg, u1i *msa, u4i nseq, u4i mlen){
	u4i roff, noff, *rpos, *xoffs, mpos, mrow;
	u4i ridx, nidx, nidxs[4];
	u1i *col;
	mrow = nseq + 3;
	roff = dg->seqs->nseq;
	noff = dg->nodes->size;
	clear_and_encap_b1v(dg->memp, 2 * nseq * sizeof(u4i));
	rpos  = (u4i*)dg->memp->buffer;
	xoffs = rpos + nseq;
	memset(rpos, 0, 2 * nseq * sizeof(u4i));
	for(ridx=0;ridx<=nseq;ridx++){
		clear_u1v(mg->cns);
		for(mpos=0;mpos<mlen;mpos++){
			col = msa + mpos * mrow;
			if(col[ridx] < 4){
				push_u1v(mg->cns, col[ridx]);
			}
		}
		if(ridx == nseq){
			fwdbitseqpush_bspoa(mg, mg->cns->buffer, mg->cns->size);
		} else {
			fwdbitseqpush_bspoa(dg, mg->cns->buffer, mg->cns->size);
			xoffs[ridx] = mg->cns->size;
		}
	}
	for(ridx=0;ridx<nseq;ridx++){
		xoffs[ridx] = dg->nodes->size + 1;
		add_rdnodes_bspoa(dg, roff + ridx);
	}
	push_u4v(dg->sels, 0);
	for(mpos=0;mpos<mlen;mpos++){ // merge inner nodes
		col = msa + mpos * mrow;
		memset(nidxs, 0, 4 * sizeof(u4i));
		for(ridx=0;ridx<nseq;ridx++){
			if(col[ridx] >= 4) continue;
			nidx = xoffs[ridx] + rpos[ridx];
			if(nidxs[col[ridx]]){
				merge_nodes_bspoa(dg, ref_bspoanodev(dg->nodes, nidxs[col[ridx]]), ref_bspoanodev(dg->nodes, nidx));
			} else {
				nidxs[col[ridx]] = nidx;
			}
			rpos[ridx] ++;
		}
		if(col[nseq] < 4){
			if(nidxs[col[nseq]]){
				push_u4v(dg->sels, nidxs[col[nseq]]);
			} else {
				push_u4v(dg->sels, 0);
			}
		}
	}
	push_u4v(dg->sels, 0);
	for(ridx=0;ridx<nseq;ridx++){
		connect_rdnodes_bspoa(dg, roff + ridx);
	}
}

static inline void end_merge_bspoa(BSPOA *dg, BSPOA *mg){
	bspoanode_t *v;
	u4i *rpos, mpos, mlen, nseq, mrow;
	u4i i, ridx, ridb, nidx, nidxs[4];
	u1i *col;
	end_bspoa(mg);
	nseq = mg->nrds;
	mlen = mg->msaidxs->size;
	mrow = nseq + 3;
	clear_and_encap_b1v(dg->memp, nseq * sizeof(u4i));
	rpos = (u4i*)dg->memp->buffer;
	memset(rpos, 0, nseq * sizeof(u4i));
	ridb = 1;
	for(mpos=0;mpos<mlen;mpos++){ // merge inter nodes
		col = mg->msacols->buffer + mg->msaidxs->buffer[mpos] * mrow;
		memset(nidxs, 0, 4 * sizeof(u4i));
		for(ridx=ridb;ridx<nseq;ridx++){
			if(col[ridx] >= 4) continue;
			v = get_rdnode_bspoa(mg, ridx, rpos[ridx]);
			nidx = offset_bspoanodev(mg->nodes, v);
			nidx = dg->sels->buffer[nidx];
			if(nidxs[col[ridx]]){
				merge_nodes_bspoa(dg, ref_bspoanodev(dg->nodes, nidxs[col[ridx]]), ref_bspoanodev(dg->nodes, nidx));
			} else {
				nidxs[col[ridx]] = nidx;
			}
			rpos[ridx] ++;
		}
	}
	dg->nrds = dg->seqs->nseq;
	msa_bspoa(dg);
	for(i=0;i<(unsigned)dg->par->realn;i++){
		cns_bspoa(dg);
		remsa_edits_bspoa(dg, dg->par->editbw);
		msa_bspoa(dg);
	}
	cns_bspoa(dg);
}

/**
 * bspoa subset
 */

static inline u4i subset_bspoa(BSPOA *g, BitVec *rdbits, BSPOA *d){
	bspoanode_t *u, *v;
	u4i i, j, rid, ridx, rdoff, rdlen;
	u4i mrow, mlen, pos;
	u1i *col, *bss;
	if(g->nrds != g->seqs->nseq){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	d->keep_seqs = 0;
	d->par->refmode = 1;
	one_bitvec(rdbits, 0);
	clear_bspoa(d);
	for(rid=0;rid<g->seqs->nseq;rid++){
		if(get_bitvec(rdbits, rid) == 0) continue;
		d->nrds ++;
	}
	d->nmsa = d->nrds;
	// copy seqs
	for(rid=0;rid<g->seqs->nseq;rid++){
		if(get_bitvec(rdbits, rid) == 0) continue;
		rdoff = g->seqs->rdoffs->buffer[rid];
		rdlen = g->seqs->rdlens->buffer[rid];
		fwdbitpush_bspoa(d, g->seqs->rdseqs->bits, rdoff, rdlen);
	}
	if(0){
		// copy nodes
		for(rid=0;rid<d->seqs->nseq;rid++){
			_add_read_bspoa_core(d, rid);
		}
		for(rid=ridx=0;rid<g->seqs->nseq;rid++){
			if(ridx && get_bitvec(rdbits, rid) == 0) continue;
			rdlen = g->seqs->rdlens->buffer[rid];
			u = get_rdnode_bspoa(g, rid, 0);
			v = get_rdnode_bspoa(d, ridx, 0);
			for(i=0;i<rdlen;i++){
				v->mpos = u->mpos;
			}
			ridx ++;
		}
	}
	append_u4v(d->msaidxs, g->msaidxs);
	mlen = d->msaidxs->size;
	mrow = d->nrds + 3;
	resize_u1v(d->msacols, mrow * mlen);
	zeros_u1v(d->msacols);
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + (g->seqs->nseq + 3) * pos;
		bss = d->msacols->buffer + mrow * pos;
		for(i=j=0;i<g->seqs->nseq;i++){
			if(get_bitvec(rdbits, i) == 0) continue;
			bss[j++] = col[i];
		}
	}
	cns_bspoa(d);
	return d->nrds;
}

static inline void calc_msa_rdregs_bspoa(BSPOA *g, u4i *rbegs, u4i *rends){
	u4i pos, i, nseq, mlen, mrow;
	u1i *col;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	for(i=0;i<nseq;i++){
		if(g->seqs->rdlens->buffer[i] == 0){
			rbegs[i] = rends[i] = 0;
			continue;
		}
		for(pos=0;pos<mlen;pos++){
			col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(col[i] < 4) break;
		}
		rbegs[i] = pos;
		for(pos=mlen;pos;pos--){
			col = g->msacols->buffer + g->msaidxs->buffer[pos - 1] * mrow;
			if(col[i] < 4) break;
		}
		rends[i] = pos;
	}
}

static inline void calc_cns_rdregs_bspoa(BSPOA *g, u4i *rbegs, u4i *rends){
	u4i pos, i, j, nseq, mlen, mrow;
	u1i *col;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	for(i=0;i<nseq;i++){
		if(g->seqs->rdlens->buffer[i] == 0){
			rbegs[i] = rends[i] = 0;
			continue;
		}
		for(pos=j=0;pos<mlen;pos++){
			col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
			if(col[i] < 4) break;
			if(col[nseq] < 4) j ++;
		}
		rbegs[i] = j;
		for(j=0,pos=mlen;pos;pos--){
			col = g->msacols->buffer + g->msaidxs->buffer[pos - 1] * mrow;
			if(col[i] < 4) break;
			if(col[nseq] < 4) j ++;
		}
		rends[i] = g->cns->size - j;
	}
}

static inline void clip_reads_msa_bspoa_core(BSPOA *g, u4i *roffs, u4i *rbegs, u4i *rends){
	u4i pos, i, nseq, mlen, mrow;
	u1i *col;
	nseq = g->nrds;
	mrow = g->seqs->nseq + 3;
	mlen = g->msaidxs->size;
	for(pos=0;pos<mlen;pos++){
		col = g->msacols->buffer + g->msaidxs->buffer[pos] * mrow;
		for(i=0;i<nseq;i++){
			if(col[i] < 4) roffs[i] ++;
			if(roffs[i] <= rbegs[i] || roffs[i] > rends[i]){
				col[i] = 5;
			}
		}
	}
}

static inline void clip_reads_msa_bspoa(BSPOA *g, u4i clip_head, u4i clip_tail){
	u4i *roffs, *rbegs, *rends, nseq, i;
	if(clip_head + clip_tail == 0) return;
	nseq = g->nrds;
	resize_u4v(g->stack, nseq);
	zeros_u4v(g->stack);
	roffs = g->stack->buffer;
	encap_b4v(g->rdregs[0], nseq);
	rbegs = (u4i*)g->rdregs[0]->buffer;
	encap_b4v(g->rdregs[1], nseq);
	rends = (u4i*)g->rdregs[1]->buffer;
	for(i=0;i<nseq;i++) rbegs[i] = clip_head;
	for(i=0;i<nseq;i++) rends[i] = num_max(g->seqs->rdlens->buffer[i], clip_tail) - clip_tail;
	clip_reads_msa_bspoa_core(g, roffs, rbegs, rends);
}

#if DEBUG
static inline void include_debug_funs_bspoa(BSPOA *g){
	if((u8i)g == MAX_U8){
		print_vstdot_bspoa(g, "1.dot");
		print_aligned_nodes_bspoa(g, 0, stdout);
	}
}
#endif

#endif
