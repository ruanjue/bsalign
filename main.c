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
	"Author : Jue Ruan <ruanjue@gmail.com>\n"
	" bsalign reads sequences from STDIN and perform banded-striped-SIMD alignment\n"
	"Usage: bsalign [options]\n"
	"options:\n"
	" -h          Show this document\n"
	" -m <string> align mode: {global/extend/*overlap*},{edit/*align*/poa} [overlap,align]\n"
	"             edit mode is to find minimal DNA edit distance, no other option\n"
	"             align mode is normal pairwise sequence alignment, supposes options -W/M/X/O/E/Q/P\n"
	"             poa mode is multiple sequences alignment, supposes options -W/M/X/O/E/Q/P and -G\n"
	"             accepts 'edit,global', 'poa,extend' and others\n"
	" -W <int>    Bandwidth, 0: full length of query [0]\n"
	"             in POA mode, default option -W 128, see -G tigger=10\n"
	" -M <int>    Score for match, [2]\n"
	" -X <int>    Penalty for mismatch, [6]\n"
	" -O <int>    Penalty for gap open, [3]\n"
	" -E <int>    Penalty for gap extension, [2]\n"
	" -Q <int>    Penalty for gap2 open, [8]\n"
	" -P <int>    Penalty for gap2 extension, [1]\n"
	" -G <sting>  parameters for POA, <tag>=<val>\n"
	"             Defaults: tigger=10,nrec=20,refmode=0,psub=0.05,pins=0.05,pdel=0.05,pext=0.50,hins=0.15,hdel=0.25\n"
	"              tigger: when <tigger> > 0 and <-W> < query length, genrates CNS per after <tigger> reads, and tigger banded alignment\n"
	"              nrec: every query read is aligning against previous <nrec> reads on graph, 0 to all the previous\n"
	"              ref: whether the first sequences is reference sequence, useful in polishing\n"
	"              psub/pins/pdel/pext: probs. of mis/ins/del/ext\n"
	"              hins/hdel: probs of ins/del in homopolymer region\n"
	"             To polish long reads' consensus with short reads, you might set\n"
	"              -G refmode=1,tigger=0,nrec=0,psub=0.02,pins=0.005,pdel=0.005,pext=0.002,hins=0.005,hdel=0.005\n"
	" -R <int>    repeat times (for benchmarking) [1]\n"
	" -v          Verbose\n"
	"\n", TOSTR(VERSION)
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	SeqBank *seqs;
	BioSequence *seq;
	BSPOAPar par;
	regex_t reg;
	regmatch_t mats[3];
	char *str, *tok;
	int c, mode, _W, W, mat, mis, gapo1, gape1, gapo2, gape2, repm, repn, verbose;
	par = DEFAULT_BSPOA_PAR;
	mode = SEQALIGN_MODE_OVERLAP;
	_W = 0;
	mat = 2; mis = -6; gapo1 = -3; gape1 = -2; gapo2 = -8; gape2 = -1;
	repm = 1;
	verbose = 0;
	c = regcomp(&reg, "([a-zA-Z_]+?)=([.0-9]+?)", REG_EXTENDED);
	if(c){
		char regtag[14];
		regerror(c, &reg, regtag, 13);
		fprintf(stderr, " -- REGCOMP: %s --\n", regtag); fflush(stderr);
		exit(1);
	}
	while((c = getopt(argc, argv, "hvm:W:M:X:O:E:Q:P:G:R:")) != -1){
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
			case 'M': par.M = mat = atoi(optarg); break;
			case 'X': par.X = mis = - atoi(optarg); break;
			case 'O': par.O = gapo1 = - atoi(optarg); break;
			case 'E': par.E = gape1 = - atoi(optarg); break;
			case 'Q': par.Q = gapo2 = - atoi(optarg); break;
			case 'P': par.P = gape2 = - atoi(optarg); break;
			case 'G':
				str = optarg;
				while(1){
					if(regexec(&reg, str, 3, mats, 0)) break;
					if(strncasecmp("psub", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.psub = atof(str + mats[2].rm_so);
					else if(strncasecmp("pins", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pins = atof(str + mats[2].rm_so);
					else if(strncasecmp("pdel", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pdel = atof(str + mats[2].rm_so);
					else if(strncasecmp("pext", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.pext = atof(str + mats[2].rm_so);
					else if(strncasecmp("hins", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.hins = atof(str + mats[2].rm_so);
					else if(strncasecmp("hdel", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.hdel = atof(str + mats[2].rm_so);
					else if(strncasecmp("nrec", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.nrec = atof(str + mats[2].rm_so);
					else if(strncasecmp("tigger", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.bwtigger = atof(str + mats[2].rm_so);
					else if(strncasecmp("refmode", str + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so) == 0) par.refmode = atoi(str + mats[2].rm_so);
					str += mats[0].rm_eo;
				}
				break;
			case 'R': repm = atoi(optarg); break;
			default: return usage();
		}
	}
	regfree(&reg);
	fr = open_filereader(NULL, 0);
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
			beg_bspoacore(g, NULL, 0, 0, 0);
			end_bspoa(g);
		}
		print_msa_bspoa(g, stdout);
		free_bspoa(g);
	} else {
		seqalign_result_t rs;
		char *alnstr[3]; int strn;
		u4v *cigars;
		u1v *qseq, *tseq;
		b1v *mempool;
		b1i mtx[16];
		banded_striped_epi8_seqalign_set_score_matrix(mtx, mat, mis);
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
						rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, W, mtx, gapo1, gape1, gapo2, gape2, verbose);
					}
				}
				if(mode & SEQALIGN_MODE_EDIT){
					rs = striped_epi2_seqedit_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, verbose);
				} else {
					rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, mempool, cigars, mode & 0x7, W, mtx, gapo1, gape1, gapo2, gape2, verbose);
				}
				if(rs.mat){
					if(strn < rs.aln){
						strn = rs.aln;
						alnstr[0] = realloc(alnstr[0], strn + 1);
						alnstr[1] = realloc(alnstr[1], strn + 1);
						alnstr[2] = realloc(alnstr[2], strn + 1);
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
