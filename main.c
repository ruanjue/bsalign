#include "bsalign.h"
#include "filereader.h"
#include <stdlib.h>

int usage(){
	fprintf(stdout,
	"Program: bsalign\n"
	"Version: %s\n"
	"Author : Jue Ruan <ruanjue@gmail.com>\n"
	" bsalign read every pairs of two sequences from STDIN and perform banded-striped-SIMD-overlap alignment\n"
	"Usage: bsalign [options]\n"
	"options:\n"
	" -h          Show this document\n"
	" -W <int>    Bandwidth, 0: full band, [0]\n"
	" -M <int>    Score for match, [2]\n"
	" -X <int>    Penalty for mismatch, [6]\n"
	" -O <int>    Penalty for gap open, [3]\n"
	" -E <int>    Penalty for gap extension, [2]\n"
	" -Q <int>    Penalty for gap2 open, [0]\n"
	" -P <int>    Penalty for gap2 extension, [0]\n"
	" -v          Verbose\n"
	" Gap-weighting functions\n"
	" If P < E and Q + P > O + E, 2-piecewise affine gap cost\n"
	"  e.g. -M 2 -X 6 -O 3 -E 2 -Q 18 -P 1\n"
	" Else if O > 0, 1-piecewise affine gap cost\n"
	"  e.g. -M 2 -X 6 -O 3 -E 2 -Q 0 -P 0\n"
	" Else, linear gap cost\n"
	"  e.g. -M 2 -X 6 -O 0 -E 5 -Q 0 -P 0\n", TOSTR(VERSION)
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	SeqBank *seqs;
	BioSequence *seq;
	String *alnstr[3];
	seqalign_result_t rs;
	u4v *begs;
	b1i mtx[16];
	int c, _W, W, mat, mis, gapo1, gape1, gapo2, gape2, verbose;
	_W = 0;
	mat = 2; mis = -6; gapo1 = -3; gape1 = -2; gapo2 = 0; gape2 = 0;
	verbose = 0;
	while((c = getopt(argc, argv, "hvW:M:X:O:E:Q:P:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'v': verbose ++; break;
			case 'W': _W = atoi(optarg); break;
			case 'M': mat = atoi(optarg); break;
			case 'X': mis = - atoi(optarg); break;
			case 'O': gapo1 = - atoi(optarg); break;
			case 'E': gape1 = - atoi(optarg); break;
			case 'Q': gapo2 = - atoi(optarg); break;
			case 'P': gape2 = - atoi(optarg); break;
			default: return usage();
		}
	}
	fr = open_filereader(NULL, 0);
	seqs = init_seqbank();
	seq = init_biosequence();
	alnstr[0] = init_string(64);
	alnstr[1] = init_string(64);
	alnstr[2] = init_string(64);
	begs  = init_u4v(64);
	{
		b1v *qprof, *rows, *btds;
		banded_striped_epi8_seqalign_set_score_matrix(mtx, mat, mis);
		qprof = init_b1v(1024);
		rows  = init_b1v(1024);
		btds  = init_b1v(1024);
		while(readseq_filereader(fr, seq)){
			push_seqbank(seqs, seq->tag->string, seq->tag->size, seq->seq->string, seq->seq->size);
			if(seqs->nseq == 2){
				if(_W <= 0) W = roundup_times(seqs->rdlens->buffer[0], 16);
				else W = _W;
				rs = banded_striped_epi8_seqalign_pairwise_overlap(seqs->rdseqs, seqs->rdoffs->buffer[0], seqs->rdlens->buffer[0], 
					seqs->rdoffs->buffer[1], seqs->rdlens->buffer[1], qprof, rows, btds, begs, W, mtx, gapo1, gape1, gapo2, gape2, alnstr, verbose);
				fprintf(stdout, "%s\t%d\t+\t%d\t%d\t%s\t%d\t+\t%d\t%d\t", seqs->rdtags->buffer[0], seqs->rdlens->buffer[0], rs.qb, rs.qe,
					seqs->rdtags->buffer[1], seqs->rdlens->buffer[1], rs.tb, rs.te);
				fprintf(stdout, "%d\t%.3f\t%d\t%d\t%d\t%d\n", rs.score, 1.0 * rs.mat / (rs.mat + rs.mis + rs.ins + rs.del) ,rs.mat, rs.mis, rs.ins, rs.del);
				fprintf(stdout, "%s\n%s\n%s\n", alnstr[0]->string, alnstr[2]->string, alnstr[1]->string);
				clear_seqbank(seqs);
			}
		}
		free_b1v(qprof);
		free_b1v(rows);
		free_b1v(btds);
	}
	free_u4v(begs);
	free_biosequence(seq);
	close_filereader(fr);
	free_string(alnstr[0]);
	free_string(alnstr[1]);
	free_string(alnstr[2]);
	free_seqbank(seqs);
	return 0;
}
