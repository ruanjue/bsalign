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
	" bsalign implements banded-striped-SIMD alignment, supports both pairwise and multiple sequence alignment\n"
	"Usage: bsalign [options] <input files>\n"
	"options:\n"
	" -h          Show this document\n"
	" -m <string> align mode: {global/extend/*overlap*},{edit/*align*/poa} [overlap,align]\n"
	"             edit mode is to find minimal DNA edit distance, ultra-fast, but only supports global alignmment\n"
	"             align mode is normal pairwise sequence alignment\n"
	"             poa mode is multiple sequences alignment\n"
	"             accepts 'edit,global', 'poa,extend' and others\n"
	" -W <int>    Bandwidth, 0: full length of query [0]\n"
	"             in POA mode, when set trigger bigger than zero to invoke banded alignment, the default option -W is 128\n"
	" -M <int>    Score for match, [2]\n"
	" -X <int>    Penalty for mismatch, [6]\n"
	" -O <int>    Penalty for gap open, [3]\n"
	" -E <int>    Penalty for gap extension, [2]\n"
	" -Q <int>    Penalty for gap2 open, [0]\n"
	" -P <int>    Penalty for gap2 extension, [0]\n"
	" -G <sting>  parameters for POA, <tag>=<val>\n"
	"             Defaults: refmode=0,refbonus=1,nrec=20,trigger=10,remsa=1,rma_win=5,qltlo=30,qlthi=35,\\"
	"                       M=1,X=2,O=0,E=1,Q=0,P=0,psub=0.05,pins=0.05,pdel=0.10,piex=0.25,pdex=0.30,hins=0.10,hdel=0.20\n"
	"              refmode: whether the first sequences is reference sequence, useful in polishing\n"
	"              refbonus: base match score on reference will be M + refbonus\n"
	"              nrec: every query read is aligning against previous <nrec> reads on graph, 0 to all the previous\n"
	"              trigger: when <trigger> > 0 and <-W> < query length, genrates CNS per after <trigger> reads, and trigger banded alignment\n"
	"              remsa: based on consensus sequence, do MSA again in <-W> band\n"
	"              rma_win: min length of flinking high quality cns bases\n"
	"              qltlo: trigger local remsa when cns quality <= qltlo\n"
	"              qlthi: high cns quality\n"
	"              M/X/O/E/Q/P: score in local realignment, can be different with whole poa alignment\n"
	"              psub/pins/pdel/piex/pdex: probs. of mis/ins/del/ins_ext/del_ext\n"
	"              hins/hdel: probs of ins/del in homopolymer region\n"
	"             To polish long reads' consensus with short reads, you might set\n"
	"              -G refmode=1,remsa=0,trigger=0,nrec=0,psub=0.02,pins=0.005,pdel=0.005,piex=0.002,pdex=0.002,hins=0.005,hdel=0.005\n"
	" -R <int>    repeat times (for benchmarking) [1]\n"
	" -v          Verbose\n"
	"Tips:\n"
	"# To invoke affine gap cost pairwise/multiple alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 0 -P 0\n"
	"# To invoke 2-piece gap cost pairwise/multiple alignment\n"
	"  -M 2 -X 6 -O 3 -E 2 -Q 8 -P 1\n"
	"# different sequencing error pattern\n"
	" tunes psub,pins,pdel,piex,pdex, hins, hdel\n"
	"# Treats with reads having large offsets with each other in multiple alignment\n"
	" set '-G trigger=0' to disable banded alignment, will be slower\n"
	"\n", TOSTR(VERSION)
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	SeqBank *seqs;
	BioSequence *seq;
	BSPOAPar par, rpar;
	regex_t reg;
	regmatch_t mats[3];
	char *str, *tok;
	//int c, mode, _W, W, mat, mis, gapo1, gape1, gapo2, gape2, repm, repn, verbose;
	int c, mode, _W, W, repm, repn, verbose;
	int msabeg, msaend, msacnt, rmabeg, rmaend;
	par  = DEFAULT_BSPOA_PAR;
	rpar = DEFAULT_BSPOA_PAR;
	mode = SEQALIGN_MODE_OVERLAP;
	_W = 0;
	par.M = 2; par.X = -6; par.O = -2; par.E = -2; par.Q = 0; par.P = 0;
	rpar.M = 1; rpar.X = -2; rpar.O = 0; rpar.E = -1; rpar.Q = 0; rpar.P = 0;
	repm = 1;
	verbose = 0;
	msabeg = 0;
	msaend = -1;
	msacnt = 3;
	rmabeg = 0;
	rmaend = -1;
	c = regcomp(&reg, "([a-zA-Z_]+?)=([.0-9]+?)", REG_EXTENDED);
	if(c){
		char regtag[14];
		regerror(c, &reg, regtag, 13);
		fprintf(stderr, " -- REGCOMP: %s --\n", regtag); fflush(stderr);
		exit(1);
	}
	while((c = getopt(argc, argv, "hvm:W:M:X:O:E:Q:P:G:T:R:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'v': verbose ++; break;
			case 'm':
			str = optarg;
			while(str && *str){
				tok = index(str, ','); if(tok) *tok = '\0';
				if(strcasecmp(str, "GLOBAL") == 0) mode = (mode & 0xF8) | SEQALIGN_MODE_GLOBAL;
				else if(strcasecmp(str, "EXTEND") == 0) mode = (mode & 0xF8) | SEQALIGN_MODE_EXTEND;
				else if(strcasecmp(str, "OVERLAP") == 0) mode = (mode & 0xF8) | SEQALIGN_MODE_OVERLAP;
				else if(strcasecmp(str, "EDIT") == 0) mode = (mode & 0x7) | SEQALIGN_MODE_EDIT;
				else if(strcasecmp(str, "ALIGN") == 0) mode = mode & 0x7;
				else if(strcasecmp(str, "POA") == 0) mode = (mode & 0x7) | SEQALIGN_MODE_POA;
				else return usage();
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
					else if(strncasecmp("trigger", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.bwtrigger = atof(str + mats[2].rm_so);
					else if(strncasecmp("refmode", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.refmode = atoi(str + mats[2].rm_so);
					else if(strncasecmp("refbonus", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.refbonus = atoi(str + mats[2].rm_so);
					else if(strncasecmp("remsa", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.remsa = atoi(str + mats[2].rm_so);
					else if(strncasecmp("rma_win", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.rma_win = atoi(str + mats[2].rm_so);
					else if(strncasecmp("qltlo", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.qltlo = atoi(str + mats[2].rm_so);
					else if(strncasecmp("qlthi", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.qlthi = atoi(str + mats[2].rm_so);
					else if(strncasecmp("M", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.M = atoi(str + mats[2].rm_so);
					else if(strncasecmp("X", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.X = - atoi(str + mats[2].rm_so);
					else if(strncasecmp("O", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.O = - atoi(str + mats[2].rm_so);
					else if(strncasecmp("E", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.E = - atoi(str + mats[2].rm_so);
					else if(strncasecmp("Q", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.Q = - atoi(str + mats[2].rm_so);
					else if(strncasecmp("P", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) rpar.P = - atoi(str + mats[2].rm_so);
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
			case 'R': repm = atoi(optarg); break;
			default: return usage();
		}
	}
	regfree(&reg);
	if(optind < argc){
		fr = open_all_filereader(argc - optind, argv + optind, 0);
	} else {
		return usage();
		//fr = open_filereader(NULL, 0);
	}
	seqs = init_seqbank();
	seq = init_biosequence();
	if((mode & SEQALIGN_MODE_POA)){
		BSPOA *g;
		bspoa_cns_debug = verbose;
		g = init_bspoa(par);
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
			BSPOA *lg;
			int recnt;
			lg = init_bspoa(rpar);
			for(recnt=0;recnt<par.remsa;recnt++){
				if(bspoa_cns_debug){
					fprintf(stdout, "## TOBE MSA\n");
					print_msa_sline_bspoa(g, stdout);
				}
				remsa_bspoa(g, lg);
			}
			free_bspoa(lg);
		}
		if(rmaend >= rmabeg){
			BSPOA *lg;
			lg = init_bspoa(rpar);
			local_remsa_bspoa(g, rmabeg, rmaend, NULL, NULL, NULL, NULL, NULL, lg);
			//print_msa_sline_bspoa(lg, stdout);
			free_bspoa(lg);
		}
		print_msa_mline_bspoa(g, stdout);
		if(msaend >= msabeg){
			FILE *out;
			out = open_file_for_write("1.dot", NULL, 1);
			print_dot_bspoa(g, msabeg, msaend, msacnt, out);
			fclose(out);
		}
		free_bspoa(g);
	} else {
		seqalign_result_t rs;
		char *alnstr[3]; int strn;
		u4v *cigars;
		u1v *qseq, *tseq;
		b1v *mempool;
		b1i mtx[16];
		banded_striped_epi8_seqalign_set_score_matrix(mtx, par.M, par.X);
		mempool = init_b1v(1024);
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
					if(mode & SEQALIGN_MODE_EDIT){
						rs = striped_epi2_seqedit_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, verbose);
					} else {
						rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, W, mtx, par.O, par.E, par.Q, par.P, verbose);
					}
				}
				if(mode & SEQALIGN_MODE_EDIT){
					rs = striped_epi2_seqedit_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, verbose);
				} else {
					rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, W, mtx, par.O, par.E, par.Q, par.P, verbose);
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
	}
	free_biosequence(seq);
	close_filereader(fr);
	free_seqbank(seqs);
	return 0;
}
