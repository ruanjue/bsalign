#ifndef BANDED_STRIPED_SIMD_RECURRENT_ALIGNMENT_GRAPH_MSA_CNS_RJ_H
#define BANDED_STRIPED_SIMD_RECURRENT_ALIGNMENT_GRAPH_MSA_CNS_RJ_H

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
static int bspoa_cns_debug = 0;

//#define BSPOACNS_NO_RECUR	1

#define SEQALIGN_MODE_POA	16

#define BSPOA_RDLEN_MAX	0x7FF8
#define BSPOA_RDCNT_MAX	0x3FFF
#define BSPOA_HEAD_NODE	0
#define BSPOA_TAIL_NODE	1

#define BSPOA_VST_MAX	MAX_U2

typedef struct {
	u2i rid, pos;
	u1i base, state;
	u2i vst;
	u2i nin, nou;
	u2i nct, cov;
	u2i cpos, rpos;
	u4i edge, erev;
	u4i aligned, header;
	u4i mmidx;
} bspoanode_t;
define_list(bspoanodev, bspoanode_t);

// put pair of edges together, idx % 2 = 0, 1
typedef struct {
	u4i node;
	u4i cov:30, is_aux:1, mask:1;
	u4i next;
} bspoaedge_t;
define_list(bspoaedgev, bspoaedge_t);

typedef struct {
	int refmode; // 0: no
	int alnmode; // SEQALIGN_MODE_GLOBAL
	int nrec; // 20, new read will align against previous nrec reads on POA
	int bwtigger; // 10, min aligned reads to tigger bandwidth for next read
	int bandwidth; // 128, when tiggering bandwidth, first generate cns, then align read against cns to get the offset of each node' band
	int M, X, O, E, Q, P; // 2, -6, -3, -2, -8, -1
	float psub, pins, pdel, pext, hins, hdel; // 0.05, 0.05, 0.05, 0.50, 0.15, 0.25
} BSPOAPar;

static const BSPOAPar DEFAULT_BSPOA_PAR = (BSPOAPar){0, SEQALIGN_MODE_GLOBAL, 20, 10, 128, 2, -6, -3, -2, -8, -1, 0.05, 0.05, 0.05, 0.50, 0.15, 0.25};

typedef struct {
	u4i coff:29, bt:3;
	float max;
} bspoacns_t;

typedef struct {
	SeqBank  *seqs;
	u4v      *ndoffs;
	u2v      *sbegs, *sends; // suggested [beg, end) on ref(1st seq) in reference-based mode
	bspoanodev *nodes;
	bspoaedgev *edges;
	u4v      *ecycs;
	BSPOAPar *par;
	int piecewise, bwtigger;
	u4i  bandwidth, qlen, nrds; // real bandwidth used in alignment
	u1v *qseq;
	b1i matrix[16];
	b1v *qprof;
	b1v *memp; // memory pool
	u8i mmblk, mmcnt;
	int maxscr, maxidx, maxoff;
	u4v *stack;
	u2i backbone;
	u4i msa_len;
	u1v *msa;
	u2v *hcovs;
	u1v *cns;
	u8i ncall, naln, nnode, nsse;
} BSPOA;

static inline BSPOA* init_bspoa(BSPOAPar par){
	BSPOA *g;
	g = malloc(sizeof(BSPOA));
	g->seqs = init_seqbank();
	g->ndoffs = init_u4v(1024);
	g->sbegs = init_u2v(32);
	g->sends = init_u2v(32);
	g->nodes = init_bspoanodev(16 * 1024);
	g->edges = init_bspoaedgev(16 * 1024);
	g->ecycs = init_u4v(32);
	g->par   = malloc(sizeof(BSPOAPar));
	memcpy(g->par, &par, sizeof(BSPOAPar));
	g->par->bandwidth = roundup_times(g->par->bandwidth, WORDSIZE);
	g->piecewise = 1;
	g->nrds      = 0;
	g->bwtigger  = 0;
	g->bandwidth = 0;
	g->qseq  = init_u1v(1024);
	g->qprof = adv_init_b1v(4 * 1024, 0, WORDSIZE, WORDSIZE);
	g->memp  = init_b1v(1024);
	g->mmblk = 0;
	g->mmcnt = 0;
	g->stack = init_u4v(32);
	g->backbone = 0;
	g->msa_len = 0;
	g->msa = init_u1v(16 * 1024);
	g->hcovs = init_u2v(4 * 1024);
	g->cns = init_u1v(1024);
	g->ncall = 0;
	g->naln = 0;
	g->nnode = 0;
	g->nsse = 0;
	return g;
}

static inline void renew_bspoa(BSPOA *g){
	free_seqbank(g->seqs); g->seqs = init_seqbank();
	renew_u4v(g->ndoffs, 1024);
	renew_u2v(g->sbegs, 32);
	renew_u2v(g->sends, 32);
	renew_bspoanodev(g->nodes, 16 * 1024);
	renew_bspoaedgev(g->edges, 16 * 1024);
	renew_u4v(g->ecycs, 32);
	renew_u1v(g->qseq, 1024);
	renew_b1v(g->qprof, 4 * 1024);
	renew_b1v(g->memp, 4 * 1024);
	g->nrds  = 0;
	g->mmblk = 0;
	g->mmcnt = 0;
	renew_u4v(g->stack, 32);
	renew_u1v(g->msa, 16 * 1024);
	renew_u2v(g->hcovs, 4 * 1024);
	renew_u1v(g->cns, 1024);
}

static inline void free_bspoa(BSPOA *g){
	free_seqbank(g->seqs);
	free_u4v(g->ndoffs);
	free_u2v(g->sbegs);
	free_u2v(g->sends);
	free_bspoanodev(g->nodes);
	free_bspoaedgev(g->edges);
	free_u4v(g->ecycs);
	free_u1v(g->qseq);
	free_b1v(g->qprof);
	free_b1v(g->memp);
	free_u4v(g->stack);
	free_u1v(g->msa);
	free_u2v(g->hcovs);
	free_u1v(g->cns);
	free(g->par);
	free(g);
}

static inline void push_bspoacore(BSPOA *g, char *seq, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX && len){
		len = num_min(len, BSPOA_RDLEN_MAX);
		push_seqbank(g->seqs, NULL, 0, seq, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void push_bspoa(BSPOA *g, char *seq, u4i len){ push_bspoacore(g, seq, len, 0, 0); }

static inline void fwdbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX && len){
		len = num_min(len, BSPOA_RDLEN_MAX);
		fwdbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void fwdbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ fwdbitpush_bspoacore(g, bits, off, len, 0, 0); }

static inline void revbitpush_bspoacore(BSPOA *g, u8i *bits, u8i off, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < BSPOA_RDCNT_MAX && len){
		len = num_min(len, BSPOA_RDLEN_MAX);
		revbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void revbitpush_bspoa(BSPOA *g, u8i *bits, u8i off, u4i len){ revbitpush_bspoacore(g, bits, off, len, 0, 0); }

static inline void print_dot_bspoa(BSPOA *g, u4i posbeg, u4i posend, FILE *out){
	bspoanode_t *n;
	bspoaedge_t *e;
	u4i nidx, eidx;
	fprintf(out, "digraph {\n");
	fprintf(out, "rankdir=LR\n");
	fprintf(out, "N0 [label=\"BEG\"]\n");
	fprintf(out, "N1 [label=\"END\"]\n");
	for(nidx=BSPOA_TAIL_NODE+1;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->pos < posbeg || n->pos >= posend) continue;
		fprintf(out, "N%u [label=R%u_%u_%c]\n", nidx, n->rid, n->pos, bit_base_table[(n->base) & 0x03]);
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_bspoanodev(g->nodes, nidx);
		if(n->pos < posbeg || n->pos >= posend) continue;
		if(n->aligned != nidx){
			fprintf(out, "N%u -> N%u [color=magenta style=dashed]\n", nidx, n->aligned);
		}
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
			fprintf(out, "N%u -> N%u [label=%u]\n", nidx, e->node, e->cov);
		}
	}
	fprintf(out, "}\n");
}

static inline void fprint_dot_bspoa(BSPOA *g, u4i posbeg, u4i posend, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_dot_bspoa(g, posbeg, posend, out);
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
		fprintf(out, "E%u N%u -> N%u(R%u_%u_%u)\n", eidx, nidx, e->node, w->rid, w->pos, w->base);
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

static inline void print_msa_bspoa(BSPOA *g, FILE *out){
	char *str;
	u1i *cns;
	u4i i, j, b, e, c, n, nseq;
	nseq = g->nrds;
	str = malloc(g->msa_len + 1);
	fprintf(out, "[POS] ");
	for(i=j=0;i<g->msa_len;i++){
		if((i % 10) == 0){
			fprintf(out, "|%04u", i + 1);
			j += 5;
		} else if(i >= j){
			putc(' ', out);
			j ++;
		}
	}
	fputc('\n', out);
	for(i=0;i<nseq+1;i++){
		if(i == nseq){
			fprintf(out, "[CNS] ");
		} else {
			fprintf(out, "[%03u] ", i);
		}
		b = i * g->msa_len;
		e = b + g->msa_len;
		n = 0;
		for(j=b;j<e;j++){
			c = g->msa->buffer[j];
			if(c < 4) n ++;
			str[j-b] = "ACGT-acgt*"[c];
		}
		str[e-b] = 0;
		fputs(str, out);
		if(i == nseq){
			fprintf(out, "\t%u\t%u\n", (u4i)g->cns->size, n);
		} else {
			fprintf(out, "\t%u\t%u\n", g->seqs->rdlens->buffer[i], n);
		}
	}
	fprintf(out, "[POS] ");
	cns = ref_u1v(g->msa, g->msa_len * nseq);
	for(i=j=b=0;i<g->msa_len;i++){
		if(cns[i] < 4){
			j ++;
			if((j % 10) == 1){
				while(b < i){
					fputc(' ', out);
					b ++;
				}
				fprintf(out, "|%04u", j);
				b += 5;
			}
		}
	}
	fprintf(out, "\n");
	if(1){
		u4i divs[5];
		u2i *hcovs;
		hcovs = g->hcovs->buffer;
		divs[0] = 10000;
		divs[1] = 1000;
		divs[2] = 100;
		divs[3] = 10;
		divs[4] = 1;
		for(i=0;i<4;i++){
			for(j=0;j<g->msa_len;j++){
				if(hcovs[j] < divs[i + 1]){
					str[j] = ' ';
				} else {
					str[j] = '0' + ((hcovs[j] % divs[i]) / divs[i + 1]);
				}
			}
			str[j] = 0;
			fprintf(out, "[%c  ] %s\n", "HCOV"[i], str);
		}
	}
	free(str);
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

static inline bspoanode_t* add_node_bspoa(BSPOA *g, u2i rid, u2i pos, u1i base, int realn){
	bspoanode_t *u;
	if(realn){
		u = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + g->seqs->rdlens->buffer[rid] - 1 - pos);
	} else {
		u = next_ref_bspoanodev(g->nodes);
	}
	ZEROS(u);
	u->rid  = rid;
	u->pos  = pos;
	u->base = base;
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
	u4i eidx, *eptr;
	int exs;
#if DEBUG
	if(u == v){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	//if(offset_bspoanodev(g->nodes, u) == 2186 && offset_bspoanodev(g->nodes, v) == 1104){
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
	encap_bspoaedgev(g->edges, 2);
	e = f = NULL;
	eptr = &u->edge;
	eidx = u->edge;
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
		e->cov ++;
		r->cov ++;
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
	del_edge_bspoacore(g, v, eidx + 0);
	del_edge_bspoacore(g, w, eidx + 1);
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
	bspoanode_t *v, *w, *y;
	UNUSED(rid);
	v = x;
	y = u;
	do {
		if(v->base == u->base){
			y = ref_bspoanodev(g->nodes, v->header);
			break;
		}
		v = ref_bspoanodev(g->nodes, v->aligned);
	} while(v != x);
	if(y != u){
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
	u->aligned = v->aligned;
	v->aligned = offset_bspoanodev(g->nodes, u);
	u->header  = offset_bspoanodev(g->nodes, y);
	return y;
}

static inline void del_read_nodes_bspoa(BSPOA *g, u2i rid){
	bspoanode_t *u, *v, *x;
	bspoaedge_t *e;
	u4i nidx, nide, xidx, rdlen, eidx, ret;
	ret = 0;
	u = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
	nidx = g->ndoffs->buffer[rid];
	rdlen = g->seqs->rdlens->buffer[rid];
	nide = nidx + rdlen;
	for(;;nidx++){
		if(nidx == nide) nidx = BSPOA_HEAD_NODE;
		v = ref_bspoanodev(g->nodes, nidx);
		xidx = v->header;
		x = ref_bspoanodev(g->nodes, xidx);
		// find edge from u -> x
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
		//if(nidx == 2186){
			//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		//}
		e->cov --;
		(e + 1)->cov --;
		if(e->cov == 0 && e->is_aux == 0){ // keep aux edge
			del_edge_bspoa(g, x, eidx);
		}
		u = x;
		if(nidx == BSPOA_HEAD_NODE) break;
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

static inline void beg_bspoacore(BSPOA *g, BaseBank *cns, u8i off, u2i len, int clear_all){
	bspoanode_t *head, *tail, *u, *v;
	bspoaedge_t *e;
	u4i i, bb;
	g->bwtigger = 0;
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
	head->aligned = BSPOA_HEAD_NODE;
	tail = next_ref_bspoanodev(g->nodes);
	ZEROS(tail);
	tail->base = 4;
	tail->aligned = BSPOA_TAIL_NODE;
	u = head;
	// add backbone nodes
	g->backbone = len;
	for(i=0;i<len;i++){
		bb = get_basebank(cns, off + i);
		v = add_node_bspoa(g, g->seqs->nseq, i, bb, 0);
		e = add_edge_bspoa(g, u, v, 0, 1, NULL);
		u = v;
	}
	e = add_edge_bspoa(g, u, tail, 0, 1, NULL);
	if(clear_all){
		clear_seqbank(g->seqs);
		clear_u2v(g->sbegs);
		clear_u2v(g->sends);
		clear_u1v(g->cns);
	}
}

static inline void beg_bspoa(BSPOA *g){
	if((g->ncall % 16) == 0){
		renew_bspoa(g);
	}
	beg_bspoacore(g, NULL, 0, 0, 1);
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
	u4i seqoff, seqlen;
	u4i nidx, eidx, i, j, x, y, op, sz;
	u2i *rmap;
	int rpos;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	g->qlen = seqlen;
	clear_and_encap_u1v(g->qseq, g->qlen);
	bitseq_basebank(g->seqs->rdseqs, seqoff, g->qlen, g->qseq->buffer);
	g->qseq->size = g->qlen;
	if(g->bwtigger && Int(roundup_times(seqlen, WORDSIZE)) > g->par->bandwidth){
		g->bandwidth = g->par->bandwidth;
		rs = striped_epi2_seqedit_pairwise(g->qseq->buffer, g->qseq->size, g->cns->buffer, g->cns->size, g->memp, g->stack, g->par->alnmode, 0);
		// cigar2rpos
		clear_and_encap_b1v(g->memp, (g->cns->size + 1) * sizeof(u2i));
		rmap = (u2i*)g->memp->buffer;
		rmap[0] = 0;
		for(i=1;Int(i)<rs.tb;i++) rmap[i] = (rs.qb - 1) / (rs.tb - i);
		x = rs.qb;
		y = rs.tb;
		for(i=0;i<g->stack->size;i++){
			op = g->stack->buffer[i] & 0xf;
			sz = g->stack->buffer[i] >> 4;
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
			}
		}
		for(i=y;i<g->cns->size;i++){
			rmap[i] = x + (i - y + 1) * (g->qseq->size - x) / (g->cns->size - y + 1);
		}
		rmap[g->cns->size] = g->qseq->size;
		// set rpos
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			rpos = rmap[u->cpos] - g->bandwidth / 2;
			if(rpos < 0) rpos = 0;
			else if(rpos + g->bandwidth > g->qlen){
				rpos = g->qlen - g->bandwidth;
			}
			u->rpos = rpos;
		}
	} else {
		g->bandwidth = roundup_times(seqlen, WORDSIZE);
		for(nidx=0;nidx<g->nodes->size;nidx++){
			u = ref_bspoanodev(g->nodes, nidx);
			u->rpos = 0;
		}
	}
	banded_striped_epi8_seqalign_set_score_matrix(g->matrix, g->par->M, g->par->X);
	clear_and_encap_b1v(g->qprof, banded_striped_epi8_seqalign_qprof_size(g->qlen, g->bandwidth));
	banded_striped_epi8_seqalign_set_query_prof(g->qseq->buffer, g->qlen, g->qprof->buffer, g->bandwidth, g->matrix);
	// prepare space for nodes and edges
	encap_bspoanodev(g->nodes, seqlen + 2);
	encap_bspoaedgev(g->edges, seqlen + 2);
	// calculate nct
	for(i=0;i<g->nodes->size;i++){
		v = ref_bspoanodev(g->nodes, i);
		v->mmidx = 0;
		v->vst  = 0;
		v->cov  = 0;
		if(g->par->nrec == 0){
			v->state = 1;
			v->nct = v->nin;
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

static inline void dpalign_row_update_bspoa(BSPOA *g, u4i mmidx1, u4i mmidx2, u4i qoff1, u4i qoff2, u1i next_base){
	b1i *us[2], *es[2], *qs[2];
	int W, rh, *ubegs[2];
	if(next_base > 3){
		fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	W = g->bandwidth / WORDSIZE;
	dpalign_row_prepare_data(g,      0, us + 0, es + 0, qs + 0, ubegs + 0);
	__builtin_prefetch(us[0], 0);
	dpalign_row_prepare_data(g, mmidx1, us + 1, es + 1, qs + 1, ubegs + 1);
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_movx(us, es, qs, ubegs, W, qoff2 - qoff1, g->piecewise, g->par->X, g->par->O, g->par->E, g->par->Q, g->par->P);
	dpalign_row_prepare_data(g, mmidx2, us + 1, es + 1, qs + 1, ubegs + 1);
	if(qoff1 == qoff2){
		rh = (qoff1 == 0 && (g->par->alnmode == SEQALIGN_MODE_OVERLAP || mmidx1 == g->nodes->buffer[BSPOA_HEAD_NODE].mmidx))? 0 : SEQALIGN_SCORE_MIN;
	} else if(qoff1 + W * WORDSIZE >= qoff2){
		//rh = banded_striped_epi8_seqalign_getscore(us[0], ubegs[0], W, qoff2 - qoff1 - 1);
		rh = ubegs[0][0];
	} else {
		rh = SEQALIGN_SCORE_MIN;
	}
	__builtin_prefetch(us[1], 1);
	banded_striped_epi8_seqalign_piecex_row_cal(qoff2, next_base, us, es, qs, ubegs, g->qprof->buffer, g->par->O, g->par->E, g->par->Q, g->par->P, W, qoff2 - qoff1, rh, g->piecewise);
	if(0){
		dpalign_row_prepare_data(g, mmidx1, us + 0, es + 0, qs + 0, ubegs + 0);
		banded_striped_epi8_seqalign_piecex_row_verify(mmidx1 != g->nodes->buffer[BSPOA_HEAD_NODE].mmidx, qoff1, g->par->alnmode, W, qoff2 - qoff1, next_base, us, es, qs, ubegs, g->qprof->buffer, g->par->O, g->par->E, g->par->Q, g->par->P);
	}
}

static inline void dpalign_row_merge_bspoa(BSPOA *g, u4i mmidx1, u4i mmidx2){
	b1i *us[3], *es[3], *qs[3];
	int *ubegs[3];
	dpalign_row_prepare_data(g, mmidx1, us + 0, es + 0, qs + 0, ubegs + 0);
	__builtin_prefetch(us[0], 0);
	dpalign_row_prepare_data(g, mmidx2, us + 1, es + 1, qs + 1, ubegs + 1);
	__builtin_prefetch(us[1], 1);
	dpalign_row_prepare_data(g, mmidx2, us + 2, es + 2, qs + 2, ubegs + 2);
	//dpalign_row_prepare_data(g, 0, us + 2, es + 2, qs + 2, ubegs + 2);
	//if(g->bandwidth == 64 && mmidx1 == 1 && mmidx2 == 1026){
		//fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
	banded_striped_epi8_seqalign_piecex_row_merge(us, es, qs, ubegs, g->bandwidth / WORDSIZE, g->piecewise);
}

static inline u2i get_rdbase_bspoa(BSPOA *g, u4i rid, u4i pos){
	u4i seqoff;
	seqoff = g->seqs->rdoffs->buffer[rid];
	return get_basebank(g->seqs->rdseqs, seqoff + pos);
}

static inline int _alignment2graph_bspoa(BSPOA *g, u4i rid, u4i midx, int xe, int *_mat, int *_tail, String **alnstrs){
	bspoanode_t *u, *v, *w, *n;
	bspoaedge_t *e;
	b1i *us, *es, *qs;
	u4i nidx, eidx, i, seqoff, cpos, W, bt;
	int x, xb, *ubegs, Hs[2], t, s, seqlen, realn, mat, tail;
	realn = 0;
	seqlen = g->seqs->rdlens->buffer[rid];
	seqoff = g->seqs->rdoffs->buffer[rid];
	encap_bspoanodev(g->nodes, seqlen);
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
	tail = seqlen - x - 1;
	if(x + 1 < (int)seqlen){
		for(i=seqlen-1;(int)i>x;i--){
			u = add_node_bspoa(g, rid, i, get_rdbase_bspoa(g, rid, i), realn);
			u->cpos = cpos;
			add_edge_bspoa(g, u, v, 1, 0, NULL);
			if(alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], '-');
				add_char_string(alnstrs[2], '-');
			}
			v = u;
		}
	}
	xb = x;
	nidx = midx;
	bt = 0;
	while(1){
		if(nidx == BSPOA_HEAD_NODE || x < 0){
			xb = x;
			while(x >= 0){
				w = add_node_bspoa(g, rid, x, get_rdbase_bspoa(g, rid, x), realn);
				w->cpos = cpos;
				e = add_edge_bspoa(g, w, v, 1, 0, NULL);
				if(alnstrs){
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
		n = ref_bspoanodev(g->nodes, nidx);
		dpalign_row_prepare_data(g, n->mmidx, &us, &es, &qs, &ubegs);
		Hs[1] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - n->rpos);
		eidx = n->erev;
		clear_u4v(g->stack);
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			w = ref_bspoanodev(g->nodes, e->node);
			if(w->state == 0) continue;
			dpalign_row_prepare_data(g, w->mmidx, &us, &es, &qs, &ubegs);
			s = g->matrix[g->qseq->buffer[x] * 4 + n->base];
			if(x < w->rpos || x > Int(g->bandwidth + w->rpos)){
				push_u4v(g->stack, e->node);
				push_u4v(g->stack, SEQALIGN_SCORE_MIN);
				push_u4v(g->stack, SEQALIGN_SCORE_MIN);
				push_u4v(g->stack, SEQALIGN_SCORE_MIN);
				continue;
			} else if(x == w->rpos){
				Hs[0] = ubegs[0];
				if(w->rpos == 0 && (g->par->alnmode == SEQALIGN_MODE_OVERLAP || offset_bspoanodev(g->nodes, w) == BSPOA_HEAD_NODE)){
				} else {
					s = SEQALIGN_SCORE_MIN;
				}
			} else {
				Hs[0] = banded_striped_epi8_seqalign_getscore(us, ubegs, W, x - w->rpos - 1);
			}
			t = ((x - w->rpos) % W) * WORDSIZE + ((x - w->rpos) / W);
			push_u4v(g->stack, e->node);
			push_u4v(g->stack, UInt(Hs[0] + s));
			push_u4v(g->stack, UInt(Hs[0] + us[t] + (es? es[t] : g->par->E)));
			push_u4v(g->stack, UInt(Hs[0] + (qs? us[t] + qs[t] : SEQALIGN_SCORE_MAX)));
		}
		bt = SEQALIGN_BT_I;
		for(i=0;i<g->stack->size;i+=4){
			if(Hs[1] == Int(g->stack->buffer[i + 1])){
				bt = SEQALIGN_BT_M;
				nidx = g->stack->buffer[i + 0];
				break;
			}
		}
		if(bt == SEQALIGN_BT_I){
			for(i=0;i<g->stack->size;i+=4){
				if(Hs[1] == Int(g->stack->buffer[i + 2]) || Hs[1] == Int(g->stack->buffer[i + 3])){
					bt = SEQALIGN_BT_D;
					nidx = g->stack->buffer[i + 0];
					break;
				}
			}
		}
		if(bt == SEQALIGN_BT_M){
			u = add_node_bspoa(g, rid, x, g->qseq->buffer[x], realn);
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
		} else if(bt == SEQALIGN_BT_D){
			if(alnstrs){
				add_char_string(alnstrs[0], '-');
				add_char_string(alnstrs[1], bit_base_table[n->base&0x03]);
				add_char_string(alnstrs[2], '-');
			}
		} else {
			u = add_node_bspoa(g, rid, x, g->qseq->buffer[x], realn);
			u->cpos = cpos;
			e = add_edge_bspoa(g, u, v, 1, 0, NULL);
			if(alnstrs){
				add_char_string(alnstrs[0], bit_base_table[u->base&0x03]);
				add_char_string(alnstrs[1], '-');
				add_char_string(alnstrs[2], '-');
			}
			if(x == 27 && rid == 10){
				fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			x --;
			v = u;
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
	int *ubegs, smax, rmax;
	seqlen = g->seqs->rdlens->buffer[rid];
	prepare_rd_align_bspoa(g, rid);
	clear_u4v(g->stack);
	push_u4v(g->stack, BSPOA_HEAD_NODE);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			f = e + 1; // v -> u
			eidx = e->next;
			//if(e->node == BSPOA_TAIL_NODE) tail = 1;
			v = ref_bspoanodev(g->nodes, e->node);
			if(v->state == 0) continue;
			if(e->node == BSPOA_TAIL_NODE){
				dpalign_row_prepare_data(g, u->mmidx, &us, &es, &qs, &ubegs);
				if(g->par->alnmode == SEQALIGN_MODE_GLOBAL){
					smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, g->qlen - 1 - v->rpos);
					if(smax > g->maxscr){
						g->maxscr = smax;
						g->maxidx = offset_bspoanodev(g->nodes, u);
						g->maxoff = g->qlen - 1;
					}
				} else {
					rmax = banded_striped_epi8_seqalign_row_max(us, ubegs, g->bandwidth / WORDSIZE, &smax);
					if(smax > g->maxscr){
						g->maxscr = smax;
						g->maxidx = offset_bspoanodev(g->nodes, u);
						g->maxoff = rmax + v->rpos;
					}
				}
				v->vst ++;
			} else {
				mmidx = v->vst? 1 : v->mmidx;
				if(v->rpos < u->rpos){
					fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				//if(rid == 1826 && u->mmidx == 2 && mmidx == 2246){
					//fflush(stdout); fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				//}
				dpalign_row_update_bspoa(g, u->mmidx, mmidx, u->rpos, v->rpos, v->base);
				if(v->vst){
					dpalign_row_merge_bspoa(g, mmidx, v->mmidx);
				}
				v->vst ++;
				if(v->vst == v->nct){
					if((g->par->alnmode & 0x7) != SEQALIGN_MODE_GLOBAL && v->rpos + g->bandwidth >= g->qlen){
						dpalign_row_prepare_data(g, v->mmidx, &us, &es, &qs, &ubegs);
						smax = banded_striped_epi8_seqalign_getscore(us, ubegs, g->bandwidth / WORDSIZE, g->qlen - 1 - v->rpos);
						if(smax > g->maxscr){
							g->maxscr = smax;
							g->maxidx = offset_bspoanodev(g->nodes, v);
							g->maxoff = g->qlen - 1;
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
	g->naln ++;
	return g->maxscr;
}

static inline int align_rd_bspoa(BSPOA *g, u2i rid){
	String **alnstrs;
	int xb, xe, rlen, score, mat, bad;
	score = align_rd_bspoacore(g, rid);
	rlen = g->seqs->rdlens->buffer[rid];
	xe = g->maxoff;
	xb = 0;
	alnstrs = NULL;
	if(bspoa_cns_debug > 1){
		alnstrs = malloc(3 * sizeof(String*));
		alnstrs[0] = init_string(rlen * 2);
		alnstrs[1] = init_string(rlen * 2);
		alnstrs[2] = init_string(rlen * 2);
	}
	xb = _alignment2graph_bspoa(g, rid, g->maxidx, xe, &mat, &bad, alnstrs);
	if(alnstrs){
		reverse_string(alnstrs[0]);
		reverse_string(alnstrs[1]);
		reverse_string(alnstrs[2]);
		fprintf(stderr, "ALIGN[%03d] len=%u ref=%d,%d band=%d aligned=%d,%d mat=%d,%0.3f tail=%d score=%d\n", rid, rlen, g->sbegs->buffer[rid], g->sends->buffer[rid], g->bandwidth, xb + 1, xe + 1, mat, 1.0 * mat / rlen, bad, score);
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

static inline u4i msa_bspoa(BSPOA *g){
	bspoanode_t *v, *u, *x;
	bspoaedge_t *e;
	u1i *r;
	u4i i, ridx, nseq, nidx, eidx, xidx, moff, ready, beg, end;
	nseq = g->nrds;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst   = 0;
		u->state = 0;
		u->cpos  = 0;
	}
	clear_u4v(g->stack);
	nidx = BSPOA_TAIL_NODE;
	push_u4v(g->stack, nidx);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_bspoanodev(g->nodes, nidx);
		eidx = u->erev;
		while(eidx){
			e = ref_bspoaedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_bspoanodev(g->nodes, e->node);
			if(u->cpos + 1 > v->cpos){
				v->cpos = u->cpos + 1;
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
			if(v->state) continue; // already pushed
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
					moff = v->cpos;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						if(x->nou > x->vst){
							ready = 0;
							break;
						}
						if(x->cpos > moff){
							moff = x->cpos;
						}
						xidx = x->aligned;
					}
				}
				if(ready){
					v->cpos  = moff;
					v->state = 1;
					push_u4v(g->stack, e->node);
					xidx = v->aligned;
					while(xidx != e->node){
						x = ref_bspoanodev(g->nodes, xidx);
						x->cpos = moff;
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
		fprint_dot_bspoa(g, 0, MAX_U4, "1.dot", NULL);
		print_seqs_bspoa(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	v = ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE);
	g->msa_len = v->cpos;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_bspoanodev(g->nodes, nidx);
		u->vst = 0;
		u->cpos = g->msa_len - 1 - u->cpos;
	}
	clear_and_encap_u1v(g->msa, g->msa_len * (nseq + 2));
	g->msa->size = g->msa_len * (nseq + 1);
	memset(g->msa->buffer, 4, g->msa->size);
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
			if(v->vst >= v->nin){
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
						if(xidx != BSPOA_TAIL_NODE){
							set_u1v(g->msa, x->rid * g->msa_len + x->cpos, (x->base) & 0x03);
						}
						if(x->edge){
							push_u4v(g->stack, xidx);
						}
						xidx = x->aligned;
					} while(xidx != e->node);
				}
			}
		}
	}
	clear_and_encap_u2v(g->hcovs, g->msa_len + 8 + nseq);
	memset(g->hcovs->buffer, 0, (g->msa_len + 8 + nseq) * sizeof(u2i));
	for(ridx=0;ridx<nseq;ridx++){
		r = ref_u1v(g->msa, g->msa_len * ridx);
		beg = 0;
		while(beg < g->msa_len && r[beg] == 4) beg ++;
		end = g->msa_len;
		while(end && r[end - 1] == 4) end --;
		for(i=beg;i<end;i++) g->hcovs->buffer[i] ++;
	}
	return g->msa_len;
}

static inline float cns_bspoa(BSPOA *g){
	typedef struct {float sc; int bt;} dp_t;
	bspoanode_t *v;
	dp_t *dps[5];
	float M, X, I, D, E, HI, HD, *table, s[5], logp;
	u1i a, b, c, d, *r, **qs;
	u4i nseq, mlen, cpos, rid, i, mpsize;
	int pos;
	M  = 0;
	X  = log(g->par->psub);
	I  = log(g->par->pins);
	D  = log(g->par->pdel);
	E  = log(g->par->pext);
	HI = log(g->par->hins);
	HD = log(g->par->hdel);
	nseq = g->nrds;
	mlen = g->msa_len;
	mpsize = roundup_times(5 * 5 * 5 * 5 * sizeof(float), 8); // scoring table
	mpsize += roundup_times(nseq * sizeof(u1i*), 8);
	mpsize += (mlen + 1) * 5 * sizeof(dp_t);
	clear_and_encap_b1v(g->memp, mpsize);
	r = (u1i*)g->memp->buffer;
	table = (float*)r; r += roundup_times(5 * 5 * 5 * 5 * sizeof(float), 8);
	qs = (u1i**)r; r += roundup_times(nseq * sizeof(u1i*), 8);
	for(i=0;i<nseq;i++) qs[i] = g->msa->buffer + mlen * i;
	for(i=0;i<5*5*5*5;i++){
		a = i / 125;
		b = (i % 125) / 25;
		c = (i % 25) / 5;
		d = i % 5;
		if(d == b){
			table[i] = M;
		} else if(d < 4){
			if(b == 4){
				if(d == a){
					table[i] = HI;
				} else if(c < 4){
					table[i] = I;
				} else {
					table[i] = E;
				}
			} else {
				table[i] = X;
			}
		} else {
			if(c == 4){
				table[i] = E;
			} else if(a == b){
				table[i] = HD;
			} else {
				table[i] = D;
			}
		}
	}
	for(i=0;i<5;i++){
		dps[i] = (dp_t*)r;
		r += (mlen + 1) * sizeof(dp_t);
		memset(dps[i], 0, (mlen + 1) * sizeof(dp_t));
	}
	for(pos=mlen-1;pos>=0;pos--){
		for(b=0;b<=4;b++){
			for(a=0;a<=4;a++){
				s[a] = dps[a][pos + 1].sc;
				for(rid=0;rid<nseq;rid++){
					s[a] += table[a * 125 + b * 25 + (pos + 1 == Int(mlen)? 4 : qs[rid][pos + 1]) * 5 + qs[rid][pos]];
				}
			}
			c = 0;
			for(a=1;a<=4;a++){
				if(s[a] > s[c]) c = a;
			}
			dps[b][pos].sc = s[c];
			dps[b][pos].bt = c;
		}
	}
	c = 0;
	for(a=1;a<=4;a++){
		if(dps[a][0].sc > dps[c][0].sc){
			c = a;
		}
	}
	logp = dps[c][0].sc;
	r = g->msa->buffer + mlen * nseq;
	clear_u1v(g->cns);
	for(i=0;i<mlen;i++){
		r[i] = c;
		if(r[i] < 4){
			push_u1v(g->cns, r[i]);
		}
		c = dps[c][i].bt;
	}
	// set node_t->cpos to the position on cns, useful in bandwidth
	for(rid=0;rid<nseq;rid++){
		cpos = 0;
		v = ref_bspoanodev(g->nodes, g->ndoffs->buffer[rid] + g->seqs->rdlens->buffer[rid] - 1);
		for(pos=0;pos<Int(mlen);pos++){
			if(qs[rid][pos] < 4){
				v->cpos = cpos;
				v --;
			}
			if(r[pos] < 4) cpos ++;
		}
	}
	ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE)->cpos = 0;
	ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE)->cpos = g->cns->size;
	return logp;
}

static inline void end_bspoa(BSPOA *g){
	bspoanode_t *u, *v;
	bspoaedge_t *e;
	u4i i, nidx, reflen, margin, tigger;
	int score;
	clear_u1v(g->cns);
	g->nrds = 0;
	if(g->seqs->nseq == 0) return;
	clear_u4v(g->ndoffs);
	nidx = 2;
	for(i=0;i<g->seqs->nseq;i++){
		push_u4v(g->ndoffs, nidx);
		nidx += g->seqs->rdlens->buffer[i];
	}
	score = align_rd_bspoa(g, 0);
	g->nrds ++;
	if(g->par->refmode){
		reflen = g->seqs->rdlens->buffer[0];
		if(g->par->bwtigger > 0){
			g->bwtigger = 1;
			clear_and_encap_u1v(g->cns, reflen);
			bitseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[0], reflen, g->cns->buffer);
			g->cns->size = reflen;
			for(i=0;i<reflen;i++){
				u = ref_bspoanodev(g->nodes, 2 + reflen - 1 - i);
				u->cpos = i;
			}
		}
		// add edges to refbeg and refend, to make sure reads can be aligned quickly and correctly within small bandwidth
		u = ref_bspoanodev(g->nodes, BSPOA_HEAD_NODE);
		margin = 5;
		for(i=0;i<g->sbegs->size;i++){
			if(g->sbegs->buffer[i] < margin || g->sbegs->buffer[i] + margin >= reflen) continue;
			nidx = BSPOA_TAIL_NODE + 1 + (reflen - 1 - g->sbegs->buffer[i]);
			v = ref_bspoanodev(g->nodes, nidx);
			e = add_edge_bspoa(g, u, v, 0, 1, NULL);
		}
		v = ref_bspoanodev(g->nodes, BSPOA_TAIL_NODE);
		for(i=0;i<g->sends->size;i++){
			if(g->sends->buffer[i] < margin || g->sends->buffer[i] + margin >= reflen) continue;
			nidx = BSPOA_TAIL_NODE + 1 + (reflen - g->sends->buffer[i]);
			u = ref_bspoanodev(g->nodes, nidx);
			e = add_edge_bspoa(g, u, v, 0, 1, NULL);
		}
	}
	tigger = g->par->bwtigger;
	for(;g->nrds<g->seqs->nseq;g->nrds++){
		if(g->par->refmode == 0 && tigger > 0 && (g->nrds % tigger) == 0){
			// generate MSA+CNS per bwtigger reads, the cpos of node is needed to perform banded alignment
			g->bwtigger = 1;
			tigger = 2 * tigger;
			msa_bspoa(g);
			cns_bspoa(g);
		}
		score = align_rd_bspoa(g, g->nrds);
	}
	// generate MSA and CNS
	msa_bspoa(g);
	cns_bspoa(g);
}

#endif
