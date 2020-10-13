#include "dna.h"
#include "filereader.h"
#include "bsalign.h"
#include "bspoa.h"
#include <stdlib.h>
#include <regex.h>

int usage(){
	fprintf(stdout,
	"Program: bsalign\n"
	"Version: %s\n"
	"Author : Jue Ruan <ruanjue@caas.cn>\n"
	"Usage  : bsalign <cmd> [options]\n"
	"\n"
	"commands:\n"
	" align       Pairwise alignment implemented by 8-bit encoded Banded Striped SIMD\n"
	" edit        Pairwise alignment using edit distance implemented by 2-bit encoded banded Striped algorithm\n"
	" poa         Multiple alignment implemented by 8-bit encoded Banded Striped SIMD Partial Order Alignment\n"
	"Tips:\n"
	"# To invoke affine gap cost pairwise/multiple alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 0 -P 0\n"
	"# To invoke 2-piece gap cost pairwise/multiple alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 8 -P 1\n"
	"\n", TOSTR(VERSION)
	);
	return 1;
}

int usage_align(){
	fprintf(stdout,
	"Usage: bsalign align [options] <input fasta/q.gz file>\n"
	" -m <string> Align mode: global/extend/overlap, [overlap]\n"
	" -W <int>    Bandwidth, 0: full length, [0]\n"
	" -M <int>    Score for match, [2]\n"
	" -X <int>    Penalty for mismatch, [6]\n"
	" -O <int>    Penalty for gap open, [3]\n"
	" -E <int>    Penalty for gap extension, [2]\n"
	" -Q <int>    Penalty for gap2 open, [0]\n"
	" -P <int>    Penalty for gap2 extension, [0]\n"
	" -R <int>    Repeat times (for benchmarking) [1]\n"
	" -v          Verbose\n"
	"# To invoke linear gap cost pairwise alignment\n"
	"  -M 2 -X 6 -O 0 -E 3 -Q 0 -P 0\n"
	"# To invoke affine gap cost pairwise alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 0 -P 0\n"
	"# To invoke 2-piece gap cost pairwise alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 8 -P 1\n"
	);
	return 1;
}

int usage_edit(){
	fprintf(stdout,
	"Usage: bsalign edit [option] <input fasta/q.gz file>\n"
	" -m <string> Align mode: global/extend/overlap/kmer, [global]\n"
	"             in overlap and extend mode, disable bandwidth\n"
	" -W <int>    Bandwidth, 0: full length, [0]\n"
	" -k <int>    Kmer size (<=15), [13]\n"
	" -v          Verbose\n"
	" -R <int>    Repeat times (for benchmarking) [1]\n"
	);
	return 1;
}

int usage_poa(){
	fprintf(stdout,
	"Usage: bsalign poa [option] <input fasta/q.gz file>\n"
	" -o <string> Consensus fasta file, [NULL]\n"
	" -m <string> Align mode: global/extend/overlap, [global]\n"
	" -W <int>    Bandwidth, 0: full length, [128]\n"
	" -M <string> Score for match for poa and local realignment [2,1]\n"
	" -X <string> Penalty for mismatch, [6,2]\n"
	" -O <string> Penalty for gap open, [3,0]\n"
	" -E <string> Penalty for gap extension, [2,1]\n"
	" -Q <string> Penalty for gap2 open, [0,0]\n"
	" -P <string> Penalty for gap2 extension, [0,0]\n"
	" -G <string> misc parameters for POA, <tag>=<val>\n"
	"             Defaults: refmode=0,refbonus=1,nrec=20,kmer=15,trigger=1,remsa=1,rma_win=5,qltlo=30,qlthi=35,\n"
	"                       psub=0.05,pins=0.05,pdel=0.10,piex=0.25,pdex=0.30,hins=0.10,hdel=0.20\n"
	"              refmode: whether the first sequences is reference sequence, useful in polishing\n"
	"              refbonus: base match score on reference will be M + refbonus\n"
	"              nrec: every query read is aligning against previous <nrec> reads on graph, 0 to all the previous\n"
	"              trigger: when <trigger> > 0 and <-W> < query length, genrates CNS per after <trigger> reads, and trigger banded alignment\n"
	"              remsa: based on consensus sequence, invoke local realignment\n"
	"              rma_win: min length of flinking high quality cns bases\n"
	"              qltlo: trigger local remsa when cns quality <= qltlo\n"
	"              qlthi: high cns quality\n"
	"              psub/pins/pdel/piex/pdex: for consensus, probs. of mis/ins/del/ins_ext/del_ext\n"
	"              hins/hdel: probs of ins/del in homopolymer region\n"
	" -L          Print MSA in 'one seq one line'\n"
	" -R <int>    Repeat times (for benchmarking) [1]\n"
	" -v          Verbose\n"
	"# different sequencing error pattern\n"
	" tunes psub,pins,pdel,piex,pdex, hins, hdel\n"
	"# Treats with reads having large offsets with each other in multiple alignment\n"
	" set '-G trigger=0' to disable banded alignment, will be slower\n"
	);
	return 1;
}

int main_edit(int argc, char **argv){
	FileReader *fr;
	SeqBank *seqs;
	BioSequence *seq;
	seqalign_result_t rs;
	u4v *cigars;
	u1v *qseq, *tseq;
	b1v *mempool;
	char *alnstr[3], *str, *tok;
	int c, repm, repn, verbose, strn, mode, W, ksz;
	repm = 1;
	verbose = 0;
	W = 0;
	ksz = 13;
	mode = SEQALIGN_MODE_GLOBAL;
	while((c = getopt(argc, argv, "hm:k:W:R:v")) != -1){
		switch(c){
			case 'm': 
			str = optarg;
			while(str && *str){
				tok = index(str, ','); if(tok) *tok = '\0';
				if(strcasecmp(str, "GLOBAL") == 0) mode = SEQALIGN_MODE_GLOBAL;
				else if(strcasecmp(str, "OVERLAP") == 0) mode = SEQALIGN_MODE_OVERLAP;
				else if(strcasecmp(str, "EXTEND") == 0) mode = SEQALIGN_MODE_EXTEND;
				else if(strcasecmp(str, "KMER") == 0) mode = SEQALIGN_MODE_KMER;
				else return usage_align();
				if(tok) str = tok + 1;
				else break;
			}
			break;
			case 'W': W = atoi(optarg); break;
			case 'k': ksz = atoi(optarg); break;
			case 'R': repm = atoi(optarg); break;
			case 'v': verbose ++; break;
			default: return usage_edit();
		}
	}
	if(optind < argc){
		fr = open_all_filereader(argc - optind, argv + optind, 0);
	} else {
		return usage_edit();
	}
	if(mode == SEQALIGN_MODE_OVERLAP && W){
		fprintf(stderr, " ** disable band in bsalign-edit's overlap mode ** \n");
		W = 0;
	}
	seqs = init_seqbank();
	seq  = init_biosequence();
	mempool = adv_init_b1v(1024, 0, WORDSIZE, 0);
	qseq   = init_u1v(1024);
	tseq   = init_u1v(1024);
	cigars = init_u4v(64);
	alnstr[0] = NULL;
	alnstr[1] = NULL;
	alnstr[2] = NULL;
	strn = 0;
	while(readseq_filereader(fr, seq)){
		if(seq->seq->size == 0) continue;
		push_seqbank(seqs, seq->tag->string, seq->tag->size, seq->seq->string, seq->seq->size);
		if(seqs->nseq == 2){
			clear_and_encap_u1v(qseq,  seqs->rdlens->buffer[0]);
			bitseq_basebank(seqs->rdseqs, seqs->rdoffs->buffer[0], seqs->rdlens->buffer[0], qseq->buffer);
			qseq->size = seqs->rdlens->buffer[0];
			clear_and_encap_u1v(tseq,  seqs->rdlens->buffer[1]);
			bitseq_basebank(seqs->rdseqs, seqs->rdoffs->buffer[1], seqs->rdlens->buffer[1], tseq->buffer);
			tseq->size = seqs->rdlens->buffer[1];
			for(repn=1;repn<repm;repn++){ // for benchmarking
				if(mode == SEQALIGN_MODE_KMER){
					rs = kmer_striped_seqedit_pairwise(ksz, qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, verbose);
				} else {
					rs = striped_seqedit_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mode, W, mempool, cigars, verbose);
				}
			}
			if(mode == SEQALIGN_MODE_KMER){
				rs = kmer_striped_seqedit_pairwise(ksz, qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, verbose);
			} else {
				rs = striped_seqedit_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mode, W, mempool, cigars, verbose);
			}
			if(rs.mat){
				if(strn < rs.aln){
					strn = rs.aln;
					alnstr[0] = realloc(alnstr[0], strn + 1);
					alnstr[1] = realloc(alnstr[1], strn + 1);
					alnstr[2] = realloc(alnstr[2], strn + 1);
				}
				if(verbose){
					u4i ci;
					fflush(stdout);
					fprintf(stderr, "CIGAR: %d\t", rs.aln);
					for(ci=0;ci<cigars->size;ci++){
						if((cigars->buffer[ci] >> 4) == 1){
							fprintf(stderr, "%c", "MIDNSHP=X*"[cigars->buffer[ci] & 0xf]);
						} else {
							fprintf(stderr, "%d%c", cigars->buffer[ci] >> 4, "MIDNSHP=X*"[cigars->buffer[ci] & 0xf]);
						}
					}
					fprintf(stderr, "\n");
				}
				seqalign_cigar2alnstr(qseq->buffer, tseq->buffer, &rs, cigars, alnstr, strn);
				fprintf(stdout, "%s\t%d\t+\t%d\t%d\t%s\t%d\t+\t%d\t%d\t", seqs->rdtags->buffer[0], Int(qseq->size), rs.qb, rs.qe, seqs->rdtags->buffer[1], Int(tseq->size), rs.tb, rs.te);
				fprintf(stdout, "%d\t%.3f\t%d\t%d\t%d\t%d\n", rs.score, 1.0 * rs.mat / rs.aln, rs.mat, rs.mis, rs.ins, rs.del);
				fprintf(stdout, "%s\n%s\n%s\n", alnstr[0], alnstr[2], alnstr[1]);
			}
			clear_seqbank(seqs);
		}
	}
	free_u1v(qseq);
	free_u1v(tseq);
	free_b1v(mempool);
	free_u4v(cigars);
	if(strn){
		free(alnstr[0]);
		free(alnstr[1]);
		free(alnstr[2]);
	}
	free_biosequence(seq);
	close_filereader(fr);
	free_seqbank(seqs);
	return 0;
}

int main_align(int argc, char **argv){
	FileReader *fr;
	SeqBank *seqs;
	BioSequence *seq;
	BSPOAPar par;
	seqalign_result_t rs;
	u4v *cigars;
	u1v *qseq, *tseq;
	b1v *mempool;
	b1i mtx[16];
	char *alnstr[3], *str, *tok;
	int c, mode, _W, W, repm, repn, verbose, strn;
	par  = DEFAULT_BSPOA_PAR;
	mode = SEQALIGN_MODE_OVERLAP;
	_W = 0;
	par.M = 2; par.X = -6; par.O = -3; par.E = -2; par.Q = 0; par.P = 0;
	repm = 1;
	verbose = 0;
	while((c = getopt(argc, argv, "hm:W:M:X:O:E:Q:P:R:v")) != -1){
		switch(c){
			case 'm':
			str = optarg;
			while(str && *str){
				tok = index(str, ','); if(tok) *tok = '\0';
				if(strcasecmp(str, "GLOBAL") == 0) mode = SEQALIGN_MODE_GLOBAL;
				else if(strcasecmp(str, "EXTEND") == 0) mode = SEQALIGN_MODE_EXTEND;
				else if(strcasecmp(str, "OVERLAP") == 0) mode = SEQALIGN_MODE_OVERLAP;
				else return usage_align();
				if(tok) str = tok + 1;
				else break;
			}
			break;
			case 'W': par.bandwidth = _W = atoi(optarg); break;
			case 'M': par.M = atoi(optarg); break;
			case 'X': par.X = - atoi(optarg); break;
			case 'O': par.O = - atoi(optarg); break;
			case 'E': par.E = - atoi(optarg); break;
			case 'Q': par.Q = - atoi(optarg); break;
			case 'P': par.P = - atoi(optarg); break;
			case 'R': repm = atoi(optarg); break;
			case 'v': verbose ++; break;
			default: return usage_align();
		}
	}
	if(optind < argc){
		fr = open_all_filereader(argc - optind, argv + optind, 0);
	} else {
		return usage_align();
	}
	banded_striped_epi8_seqalign_set_score_matrix(mtx, par.M, par.X);
	seqs = init_seqbank();
	seq  = init_biosequence();
	mempool = adv_init_b1v(1024, 0, WORDSIZE, 0);
	qseq   = init_u1v(1024);
	tseq   = init_u1v(1024);
	cigars = init_u4v(64);
	alnstr[0] = NULL;
	alnstr[1] = NULL;
	alnstr[2] = NULL;
	strn = 0;
	while(readseq_filereader(fr, seq)){
		if(seq->seq->size == 0) continue;
		push_seqbank(seqs, seq->tag->string, seq->tag->size, seq->seq->string, seq->seq->size);
		if(seqs->nseq == 2){
			if(_W <= 0) W = roundup_times(seqs->rdlens->buffer[0], 16);
			else W = _W;
			clear_and_encap_u1v(qseq,  seqs->rdlens->buffer[0]);
			bitseq_basebank(seqs->rdseqs, seqs->rdoffs->buffer[0], seqs->rdlens->buffer[0], qseq->buffer);
			qseq->size = seqs->rdlens->buffer[0];
			clear_and_encap_u1v(tseq,  seqs->rdlens->buffer[1]);
			bitseq_basebank(seqs->rdseqs, seqs->rdoffs->buffer[1], seqs->rdlens->buffer[1], tseq->buffer);
			tseq->size = seqs->rdlens->buffer[1];
			for(repn=1;repn<repm;repn++){ // for benchmarking
				rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode, W, mtx, par.O, par.E, par.Q, par.P, verbose);
			}
			rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode, W, mtx, par.O, par.E, par.Q, par.P, verbose);
			if(rs.mat){
				if(strn < rs.aln){
					strn = rs.aln;
					alnstr[0] = realloc(alnstr[0], strn + 1);
					alnstr[1] = realloc(alnstr[1], strn + 1);
					alnstr[2] = realloc(alnstr[2], strn + 1);
				}
				if(verbose){
					u4i ci;
					fprintf(stderr, "CIGAR: %d\t", rs.aln);
					for(ci=0;ci<cigars->size;ci++){
						if((cigars->buffer[ci] >> 4) == 1){
							fprintf(stderr, "%c", "MIDNSHP=X*"[cigars->buffer[ci] & 0xf]);
						} else {
							fprintf(stderr, "%d%c", cigars->buffer[ci] >> 4, "MIDNSHP=X*"[cigars->buffer[ci] & 0xf]);
						}
					}
					fprintf(stderr, "\n");
				}
				seqalign_cigar2alnstr(qseq->buffer, tseq->buffer, &rs, cigars, alnstr, strn);
				fprintf(stdout, "%s\t%d\t+\t%d\t%d\t%s\t%d\t+\t%d\t%d\t", seqs->rdtags->buffer[0], Int(qseq->size), rs.qb, rs.qe, seqs->rdtags->buffer[1], Int(tseq->size), rs.tb, rs.te);
				fprintf(stdout, "%d\t%.3f\t%d\t%d\t%d\t%d\n", rs.score, 1.0 * rs.mat / rs.aln, rs.mat, rs.mis, rs.ins, rs.del);
				fprintf(stdout, "%s\n%s\n%s\n", alnstr[0], alnstr[2], alnstr[1]);
				fflush(stdout);
			}
			clear_seqbank(seqs);
		}
	}
	free_u1v(qseq);
	free_u1v(tseq);
	free_b1v(mempool);
	free_u4v(cigars);
	if(strn){
		free(alnstr[0]);
		free(alnstr[1]);
		free(alnstr[2]);
	}
	free_biosequence(seq);
	close_filereader(fr);
	free_seqbank(seqs);
	return 0;
}

int main_poa(int argc, char **argv){
	FileReader *fr;
	FILE *out;
	SeqBank *seqs;
	BioSequence *seq;
	BSPOAPar par, rpar;
	BSPOA *g, *lg;
	regex_t reg;
	regmatch_t mats[3];
	char *str, *tok, *cnsfn;
	int c, repm, repn, mline, verbose;
	int msabeg, msaend, msacnt, rmabeg, rmaend;
	par  = DEFAULT_BSPOA_PAR;
	rpar = DEFAULT_BSPOA_PAR;
	 par.ksz = 13; par.alnmode = SEQALIGN_MODE_OVERLAP;  par.M = 2;  par.X = -6; par.O = -3;  par.E = -2;  par.Q = 0;  par.P = 0;
	rpar.ksz = 0; rpar.alnmode =  SEQALIGN_MODE_GLOBAL; rpar.M = 1; rpar.X = -2; rpar.O = 0; rpar.E = -1; rpar.Q = 0; rpar.P = 0;
	mline = 1;
	repm = 1;
	verbose = 0;
	msabeg = 0;
	msaend = -1;
	msacnt = 3;
	rmabeg = 0;
	rmaend = -1;
	cnsfn = NULL;
	c = regcomp(&reg, "([a-zA-Z_]+?)=([.0-9]+?)", REG_EXTENDED);
	if(c){
		char regtag[14];
		regerror(c, &reg, regtag, 13);
		fprintf(stderr, " -- REGCOMP: %s --\n", regtag); fflush(stderr);
		exit(1);
	}
	while((c = getopt(argc, argv, "hvo:m:W:M:X:O:E:Q:P:G:LT:R:")) != -1){
		switch(c){
			case 'h': return usage_poa();
			case 'v': verbose ++; break;
			case 'o': cnsfn = optarg; break;
			case 'm':
			str = optarg;
			while(str && *str){
				tok = index(str, ','); if(tok) *tok = '\0';
				if(strcasecmp(str, "GLOBAL") == 0) par.alnmode = SEQALIGN_MODE_GLOBAL;
				else if(strcasecmp(str, "EXTEND") == 0) par.alnmode = SEQALIGN_MODE_EXTEND;
				else if(strcasecmp(str, "OVERLAP") == 0) par.alnmode = SEQALIGN_MODE_OVERLAP;
				else return usage_poa();
				if(tok) str = tok + 1;
				else break;
			}
			break;
			case 'W': par.bandwidth = atoi(optarg); break;
			case 'M': par.M = strtol(optarg, &tok, 10);   if(tok && tok[1]) rpar.M = strtol(tok + 1, &tok, 10); break;
			case 'X': par.X = - strtol(optarg, &tok, 10); if(tok && tok[1]) rpar.X = - strtol(tok + 1, &tok, 10); break;
			case 'O': par.O = - strtol(optarg, &tok, 10); if(tok && tok[1]) rpar.O = - strtol(tok + 1, &tok, 10); break;
			case 'E': par.E = - strtol(optarg, &tok, 10); if(tok && tok[1]) rpar.E = - strtol(tok + 1, &tok, 10); break;
			case 'Q': par.Q = - strtol(optarg, &tok, 10); if(tok && tok[1]) rpar.Q = - strtol(tok + 1, &tok, 10); break;
			case 'P': par.P = - strtol(optarg, &tok, 10); if(tok && tok[1]) rpar.P = - strtol(tok + 1, &tok, 10); break;
			case 'G':
				str = optarg;
				while(1){
					if(regexec(&reg, str, 3, mats, 0)) break;
					if(strncasecmp("psub", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.psub = atof(str + mats[2].rm_so);
					else if(strncasecmp("pins", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pins = atof(str + mats[2].rm_so);
					else if(strncasecmp("pdel", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pdel = atof(str + mats[2].rm_so);
					else if(strncasecmp("piex", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.piex = atof(str + mats[2].rm_so);
					else if(strncasecmp("pdex", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pdex = atof(str + mats[2].rm_so);
					else if(strncasecmp("hins", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.hins = atof(str + mats[2].rm_so);
					else if(strncasecmp("hdel", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.hdel = atof(str + mats[2].rm_so);
					else if(strncasecmp("nrec", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.nrec = atof(str + mats[2].rm_so);
					else if(strncasecmp("kmer", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.ksz = atoi(str + mats[2].rm_so);
					else if(strncasecmp("trigger", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.bwtrigger = atof(str + mats[2].rm_so);
					else if(strncasecmp("refmode", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.refmode = atoi(str + mats[2].rm_so);
					else if(strncasecmp("refbonus", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.refbonus = atoi(str + mats[2].rm_so);
					else if(strncasecmp("remsa", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.remsa = atoi(str + mats[2].rm_so);
					else if(strncasecmp("rma_win", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.rma_win = atoi(str + mats[2].rm_so);
					else if(strncasecmp("qltlo", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.qltlo = atoi(str + mats[2].rm_so);
					else if(strncasecmp("qlthi", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.qlthi = atoi(str + mats[2].rm_so);
					else {
						fprintf(stderr, "Unknown parameter: %s\n", str);
						return 1;
					}
					str += mats[0].rm_eo;
				}
				break;
			case 'T':
				str = optarg;
				while(1){
					if(regexec(&reg, str, 3, mats, 0)) break;
					if(strncasecmp("msabeg", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) msabeg = atoi(str + mats[2].rm_so);
					else if(strncasecmp("msaend", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) msaend = atoi(str + mats[2].rm_so);
					else if(strncasecmp("msacnt", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) msacnt = atoi(str + mats[2].rm_so);
					else if(strncasecmp("rmabeg", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rmabeg = atoi(str + mats[2].rm_so);
					else if(strncasecmp("rmaend", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rmaend = atoi(str + mats[2].rm_so);
					else {
						fprintf(stderr, "Unknown parameter: %s\n", str);
						return 1;
					}
					str += mats[0].rm_eo;
				}
			case 'L': mline = 0; break;
			case 'R': repm = atoi(optarg); break;
			default: return usage_poa();
		}
	}
	regfree(&reg);
	if(optind < argc){
		fr = open_all_filereader(argc - optind, argv + optind, 0);
	} else {
		return usage_poa();
	}
	if(cnsfn){
		out = open_file_for_write(cnsfn, NULL, 1);
	} else {
		out = NULL;
	}
	seqs = init_seqbank();
	seq = init_biosequence();
	_DEBUG_LOG_ = verbose;
	g  = init_bspoa(par);
	lg = init_bspoa(rpar);
	beg_bspoa(g);
	while(readseq_filereader(fr, seq)){
		if(seq->seq->size == 0) continue;
		push_bspoa(g, seq->seq->string, seq->seq->size);
	}
	end_bspoa(g);
	for(repn=1;repn<repm;repn++){ // for benchmarking
		beg_bspoacore(g, NULL, 0, 0);
		end_bspoa(g);
	}
	if(par.remsa){
		int recnt;
		for(recnt=0;recnt<par.remsa;recnt++){
			if(_DEBUG_LOG_){
				fprintf(stdout, "## TOBE MSA\n");
				print_msa_sline_bspoa(g, stdout);
			}
			remsa_bspoa(g, lg);
		}
	}
	if(rmaend >= rmabeg){
		local_remsa_bspoa(g, rmabeg, rmaend, NULL, NULL, NULL, NULL, NULL, lg);
		print_msa_sline_bspoa(lg, stderr);
	}
	if(mline){
		print_msa_mline_bspoa(g, stdout);
	} else {
		print_msa_sline_bspoa(g, stdout);
	}
	print_snp_bspoa(g, stdout);
	if(out){
		fprintf(out, ">cns_seq\n%s\n", g->strs->string);
		close_file(out);
	}
	if(msaend >= msabeg){
		FILE *out;
		out = open_file_for_write("1.dot", NULL, 1);
		print_dot_bspoa(g, msabeg, msaend, msacnt, out);
		fclose(out);
	}
	free_bspoa(g);
	free_bspoa(lg);
	free_biosequence(seq);
	close_filereader(fr);
	free_seqbank(seqs);
	return 0;
}

int main(int argc, char **argv){
	if(argc < 2){
		return usage();
	}
	if(strcasecmp("EDIT", argv[1]) == 0) return main_edit(argc - 1, argv + 1);
	if(strcasecmp("ALIGN", argv[1]) == 0) return main_align(argc - 1, argv + 1);
	if(strcasecmp("POA", argv[1]) == 0) return main_poa(argc - 1, argv + 1);
	fprintf(stderr, " -- Unknown command '%s' -- \n", argv[0]);
	return 1;
}
