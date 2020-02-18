#include "bsalign.h"
#include "dna.h"
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
	" -m <string> align mode: global, extend, overlap [overlap]\n"
	" -W <int>    Bandwidth, 0: full length of query [0]\n"
	"             set it to a small value if the begins of two reads are nearly aligned\n"
	" -M <int>    Score for match, [2]\n"
	" -X <int>    Penalty for mismatch, [6]\n"
	" -O <int>    Penalty for gap open, [3]\n"
	" -E <int>    Penalty for gap extension, [2]\n"
	" -Q <int>    Penalty for gap2 open, [0]\n"
	" -P <int>    Penalty for gap2 extension, [0]\n"
	" -R <int>    repeat times (for benchmarking) [1]\n"
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
	int c, mode, _W, W, mat, mis, gapo1, gape1, gapo2, gape2, repm, repn, verbose;
	mode = SEQALIGN_MODE_OVERLAP;
	_W = 0;
	mat = 2; mis = -6; gapo1 = -3; gape1 = -2; gapo2 = 0; gape2 = 0;
	repm = 1;
	verbose = 0;
	while((c = getopt(argc, argv, "hvm:W:M:X:O:E:Q:P:R:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'v': verbose ++; break;
			case 'm':
			if(strcasecmp(optarg, "GLOBAL") == 0) mode = SEQALIGN_MODE_GLOBAL;
			else if(strcasecmp(optarg, "EXTEND") == 0) mode = SEQALIGN_MODE_EXTEND;
			else if(strcasecmp(optarg, "OVERLAP") == 0) mode = SEQALIGN_MODE_OVERLAP;
			else return usage();
			break;
			case 'W': _W = atoi(optarg); break;
			case 'M': mat = atoi(optarg); break;
			case 'X': mis = - atoi(optarg); break;
			case 'O': gapo1 = - atoi(optarg); break;
			case 'E': gape1 = - atoi(optarg); break;
			case 'Q': gapo2 = - atoi(optarg); break;
			case 'P': gape2 = - atoi(optarg); break;
			case 'R': repm = atoi(optarg); break;
			default: return usage();
		}
	}
	fr = open_filereader(NULL, 0);
	seqs = init_seqbank();
	seq = init_biosequence();
	{
		seqalign_result_t rs;
		char *alnstr[3]; int strn;
		u4v *begs, *cigars;
		u1v *qseq, *tseq;
		b1v *qprof, *rows, *btds;
		b1i mtx[16];
		banded_striped_epi8_seqalign_set_score_matrix(mtx, mat, mis);
		qseq   = init_u1v(1024);
		tseq   = init_u1v(1024);
		qprof  = adv_init_b1v(1024, 0, 32, 32); // 32-bytes aligned memory, there is also 32 allocated bytes before qprof->buffer to support accessing buffer[-1]
		rows   = adv_init_b1v(1024, 0, 32, 0);
		btds   = adv_init_b1v(1024, 0, 32, 0);
		begs   = init_u4v(64);
		cigars = init_u4v(64);
		alnstr[0] = NULL;
		alnstr[1] = NULL;
		alnstr[2] = NULL;
		strn = 0;
		while(readseq_filereader(fr, seq)){
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
					rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, qprof, rows, btds, begs, cigars, mode, W, mtx, gapo1, gape1, gapo2, gape2, verbose);
				}
				rs = banded_striped_epi8_seqalign_pairwise(qseq->buffer, qseq->size, tseq->buffer, tseq->size, qprof, rows, btds, begs, cigars, mode, W, mtx, gapo1, gape1, gapo2, gape2, verbose);
				if(strn < rs.aln){
					strn = rs.aln;
					alnstr[0] = realloc(alnstr[0], strn + 1);
					alnstr[1] = realloc(alnstr[1], strn + 1);
					alnstr[2] = realloc(alnstr[2], strn + 1);
				}
				banded_striped_epi8_seqalign_cigar2alnstr(qseq->buffer, tseq->buffer, &rs, cigars, alnstr, strn);
				fprintf(stdout, "%s\t%d\t+\t%d\t%d\t%s\t%d\t+\t%d\t%d\t", seqs->rdtags->buffer[0], Int(qseq->size), rs.qb, rs.qe, seqs->rdtags->buffer[1], Int(tseq->size), rs.tb, rs.te);
				fprintf(stdout, "%d\t%.3f\t%d\t%d\t%d\t%d\n", rs.score, 1.0 * rs.mat / rs.aln, rs.mat, rs.mis, rs.ins, rs.del);
				fprintf(stdout, "%s\n%s\n%s\n", alnstr[0], alnstr[2], alnstr[1]);
				clear_seqbank(seqs);
			}
		}
		free_u1v(qseq);
		free_u1v(tseq);
		free_b1v(qprof);
		free_b1v(rows);
		free_b1v(btds);
		free_u4v(begs);
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
