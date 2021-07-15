#ifndef __GET_OPTIONS_RJ
#define __GET_OPTIONS_RJ
#include "chararray.h"
#include "list.h"
#include "hashset.h"
#include <getopt.h>
#include <regex.h>

#define PROG_OPT_TYPE_SET	0   // b
#define PROG_OPT_TYPE_INT	1   // i
#define PROG_OPT_TYPE_INC	2   // a
#define PROG_OPT_TYPE_FLT	3   // f
#define PROG_OPT_TYPE_STR	4   // c
#define PROG_OPT_TYPE_INTS	5   // I
#define PROG_OPT_TYPE_FLTS	6   // F
#define PROG_OPT_TYPE_STRS	7   // C

typedef struct {
	char   *opt_tag;
	u4i     parent;
	char    opt_alt;
	int     opt_typ, opt_req, opt_cnt;
	char    opt_sep;
	char   *opt_dsc;
	b8i     val_int;
	b8i     val_inc;
	f8i     val_flt;
	char   *val_str;
	b8v    *val_ints;
	f8v    *val_flts;
	cplist *val_strs;
} prog_opt_t;
define_list(progoptv, prog_opt_t);

typedef struct {
	progoptv *opts;
	cuhash   *hash;
} PROGOPT;

static inline int parse_val_progopt(PROGOPT *pg, u4i optidx, char *valstr, int init){
	prog_opt_t *opt;
	char *ptr, *tok;
	u4i i;
	opt = ref_progoptv(pg->opts, optidx);
	if(opt->opt_cnt == 0){
		switch(opt->opt_typ){
			case PROG_OPT_TYPE_INTS: clear_b8v(opt->val_ints); break;
			case PROG_OPT_TYPE_FLTS: clear_f8v(opt->val_flts); break;
			case PROG_OPT_TYPE_STRS:
				for(i=0;i<opt->val_strs->size;i++){
					if(opt->val_strs->buffer[i]) free(opt->val_strs->buffer[i]);
				}
				clear_cplist(opt->val_strs);
				break;
		}
	}
	if(!init){
		opt->opt_cnt ++;
	}
	if((valstr == NULL || valstr[0] == '\0') && (opt->opt_typ != PROG_OPT_TYPE_SET && opt->opt_typ != PROG_OPT_TYPE_INC)){
		return 0;
	}
	switch(opt->opt_typ){
		case PROG_OPT_TYPE_SET: opt->val_int = (init == 0); break;
		case PROG_OPT_TYPE_INT: opt->val_int = atoll(valstr); break;
		case PROG_OPT_TYPE_INC: opt->val_int += (init == 0); break;
		case PROG_OPT_TYPE_FLT: opt->val_flt = strtod(valstr, NULL); break;
		case PROG_OPT_TYPE_STR: if(opt->val_str) free(opt->val_str); opt->val_str = strdup(valstr); break;
		case PROG_OPT_TYPE_INTS:
			ptr = valstr;
			while(1){
				b8i v;
				tok = index(ptr, opt->opt_sep);
				v = strtoll(ptr, NULL, 10);
				push_b8v(opt->val_ints, v);
				if(tok == NULL || *tok == '\0') break;
				ptr = tok + 1;
			}
			break;
		case PROG_OPT_TYPE_FLTS:
			ptr = valstr;
			while(1){
				f8i v;
				tok = index(ptr, opt->opt_sep);
				v = strtold(ptr, NULL);
				push_f8v(opt->val_flts, v);
				if(tok == NULL || *tok == '\0') break;
				ptr = tok + 1;
			}
			break;
		case PROG_OPT_TYPE_STRS:
			ptr = valstr;
			while(1){
				tok = index(ptr, opt->opt_sep);
				push_cplist(opt->val_strs, strndup(ptr, tok? (size_t)(tok - ptr) : strlen(ptr)));
				if(tok == NULL || *tok == '\0') break;
				ptr = tok + 1;
			}
			break;
		default:
			fflush(stdout); fprintf(stderr, " ** PROG_OPT * bad type [%d] \n", opt->opt_typ); fflush(stderr);
			return 1;
	}
	return 0;
}

static inline u4i add_progopt(PROGOPT *pg, char *tag, u4i parent, int type, char alt, char *val, char sep, char *dsc){
	prog_opt_t *opt, *prt;
	u4i optidx;
	if(exists_cuhash(pg->hash, tag)){
		fflush(stdout); fprintf(stderr, " ** PROG_OPT * duplicated option '%s'\n", tag); fflush(stderr);
	}
	optidx = pg->opts->size;
	opt = next_ref_progoptv(pg->opts);
	if(parent){
		ZEROS(opt);
		prt = ref_progoptv(pg->opts, parent);
		opt->opt_tag = strdup(tag);
		opt->opt_typ = type;
		opt->parent  = parent;
		opt->opt_req = prt->opt_req;
		put_cuhash(pg->hash, (cuhash_t){opt->opt_tag, optidx});
	} else {
		opt->opt_cnt = 0;
		opt->opt_tag = strdup(tag);
		opt->parent  = 0;
		opt->opt_alt = alt;
		opt->opt_typ = type;
		opt->opt_req = (type != PROG_OPT_TYPE_SET && type != PROG_OPT_TYPE_INC);
		opt->opt_sep = sep;
		opt->opt_dsc = strdup(dsc);
		opt->val_int = 0;
		opt->val_flt = 0;
		opt->val_str = NULL;
		opt->val_ints = init_b8v(4);
		opt->val_flts = init_f8v(4);
		opt->val_strs = init_cplist(4);
		put_cuhash(pg->hash, (cuhash_t){opt->opt_tag, optidx});
		parse_val_progopt(pg, optidx, val, 1);
	}
	return parent? parent : optidx;
}

static inline PROGOPT* init_progopt(){
	PROGOPT *pg;
	//prog_opt_t *opt;
	pg = malloc(sizeof(PROGOPT));
	pg->opts = init_progoptv(4);
	pg->hash = init_cuhash(11);
	add_progopt(pg, "thereisabeautifulworld", 0, PROG_OPT_TYPE_SET, 0, NULL, 0, "howtoyourheart");
	//ZEROS(opt = next_ref_progoptv(pg->opts));
	//opt->opt_tag = "thereisabeautifulworld";
	//opt->opt_req = 0;
	return pg;
}

static inline int adds_progopt(PROGOPT *pg, const char *txt){
	String *tag, *val, *dsc;
	char alt, sep;
	int type, ret, line, optidx;
	char *lb, *le, *ptr;
	tag = init_string(16);
	val = init_string(16);
	dsc = init_string(16);
	lb = ptr = (char*)txt;
	ret = 0;
	alt = '\0'; sep = ','; type = PROG_OPT_TYPE_SET; clear_string(tag); clear_string(val); clear_string(dsc);
	optidx = 0;
	line = 0;
	while(1){
		if(ptr[0] != '\0' && ptr[0] != '\n'){
			ptr ++;
			continue;
		}
		line ++;
		le = ptr;
		ptr = lb;
		if(ptr[0] == ' '){
			append_string(dsc, lb + 1, le - lb - 1);
		} else {
			if(tag->size){
				optidx = add_progopt(pg, tag->string, optidx, type, alt, val->size? val->string : NULL, sep, dsc->string);
				ret ++;
			}
			alt = '\0'; sep = ','; type = PROG_OPT_TYPE_SET; clear_string(tag); clear_string(val); clear_string(dsc);
			while(1){
				{
					if(ptr < le){
						if(ptr[0] == '='){
							// current option is a alias of the previous
							ptr ++;
							lb ++;
						} else {
							optidx = 0;
						}
					}
					while(ptr < le && ptr[0] != ' ' && ptr[0] != '\t') ptr ++;
					if(ptr == lb) break;
					append_string(tag, lb, ptr - lb);
					while(ptr < le && (ptr[0] == ' ' || ptr[0] == '\t')) ptr ++;
					lb = ptr;
				}
				{
					while(ptr < le && ptr[0] != ' ' && ptr[0] != '\t') ptr ++;
					if(ptr == lb) break;
					if(ptr > lb + 1){
						fflush(stdout); fprintf(stderr, " ** PROG_OPT * bad option type '%c%c': line %d\n", lb[0], lb[1], line); fflush(stderr);
						break;
					}
					switch(lb[0]){
						case 'b': type = PROG_OPT_TYPE_SET; break;
						case 'i': type = PROG_OPT_TYPE_INT; break;
						case 'a': type = PROG_OPT_TYPE_INC; break;
						case 'f': type = PROG_OPT_TYPE_FLT; break;
						case 'c': type = PROG_OPT_TYPE_STR; break;
						case 'I': type = PROG_OPT_TYPE_INTS; break;
						case 'F': type = PROG_OPT_TYPE_FLTS; break;
						case 'C': type = PROG_OPT_TYPE_STRS; break;
						default:
							fflush(stdout); fprintf(stderr, " ** PROG_OPT * bad option type '%c': line %d\n", lb[0], line); fflush(stderr);
					}
					while(ptr < le && (ptr[0] == ' ' || ptr[0] == '\t')) ptr ++;
					lb = ptr;
				}
				{
					while(ptr < le && ptr[0] != ' ' && ptr[0] != '\t') ptr ++;
					if(ptr == lb) break;
					if(ptr > lb + 1){
						fflush(stdout); fprintf(stderr, " ** PROG_OPT * bad shortcut '%c%c': line %d\n", lb[0], lb[1], line); fflush(stderr);
						break;
					}
					alt = (lb[0] == '*')? '\0' : lb[0];
					while(ptr < le && (ptr[0] == ' ' || ptr[0] == '\t')) ptr ++;
					lb = ptr;
				}
				{
					while(ptr < le && ptr[0] != ' ' && ptr[0] != '\t') ptr ++;
					if(ptr == lb) break;
					if(ptr - lb == 1 && lb[0] == '*'){
					} else {
						append_string(val, lb, ptr - lb);
					}
					while(ptr < le && (ptr[0] == ' ' || ptr[0] == '\t')) ptr ++;
					lb = ptr;
				}
				{
					while(ptr < le && ptr[0] != ' ' && ptr[0] != '\t') ptr ++;
					if(ptr == lb) break;
					if(ptr > lb + 1){
						fflush(stdout); fprintf(stderr, " ** PROG_OPT * bad SEP '%c%c': line %d\n", lb[0], lb[1], line); fflush(stderr);
						break;
					}
					sep = (lb[0] == '*')? '\0' : lb[0];
					while(ptr < le && (ptr[0] == ' ' || ptr[0] == '\t')) ptr ++;
					lb = ptr;
				}
				break;
			}
		}
		if(le[0] == '\0') break;
		lb = ptr = le + 1;
	}
	free_string(tag);
	free_string(val);
	free_string(dsc);
	return ret;
}

static inline int parse_progopt(PROGOPT *pg, int argc, char **argv){
	prog_opt_t *opt;
	struct option *popts;
	char *alts, *ptr;
	u4i i;
	int c, optidx;
	if(pg->opts->size == 0) return 0;
	popts = alloca((pg->opts->size + 1) * sizeof(struct option));
	alts  = alloca(pg->opts->size * 2 + 1);
	ptr = alts;
	for(i=0;i<pg->opts->size;i++){
		opt = ref_progoptv(pg->opts, i);
		popts[i].name    = opt->opt_tag;
		popts[i].has_arg = opt->opt_req;
		popts[i].flag    = NULL;
		popts[i].val     = opt->opt_alt? opt->opt_alt : Int(256 + i);
		if(opt->opt_alt){
			*ptr = opt->opt_alt;
			if(opt->opt_req){
				ptr ++;
				*ptr = ':';
			}
			ptr ++;
		}
	}
	popts[i].name    = NULL;
	popts[i].has_arg = 0;
	popts[i].flag    = NULL;
	popts[i].val     = -1;
	ptr = '\0';
	optidx = -1;
	while((c = getopt_long(argc, argv, alts, popts, &optidx)) != -1){
		if(optidx == -1){
			for(i=1;i<pg->opts->size;i++){
				opt = ref_progoptv(pg->opts, i);
				if(opt->opt_alt == c){
					optidx = i;
					break;
				}
			}
			if(optidx == -1){
				continue;
			}
		}
		opt = ref_progoptv(pg->opts, optidx);
		if(opt->parent){
			optidx = opt->parent;
		}
		parse_val_progopt(pg, optidx, optarg, 0);
		optidx = -1;
	}
	return optind;
}

static inline prog_opt_t* get_progopt(PROGOPT *pg, char *tag, const char *file, const char *func, const int line){
	prog_opt_t *opt;
	u4i optidx;
	optidx = getval_cuhash(pg->hash, tag);
	if(optidx == MAX_U4){
		fflush(stdout); fprintf(stderr, " -- \e[7mUnknown option '%s' in %s -- %s:%d\e[0m --\n", tag, func, file, line); fflush(stderr);
		optidx = 0;
	}
	opt = ref_progoptv(pg->opts, optidx);
	if(opt->parent){
		optidx = opt->parent;
	}
	return ref_progoptv(pg->opts, optidx);
}

#define getopt_int(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_int)
#define OPTINT(tag) getopt_int(opts, tag)
#define getopt_ints(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_ints)
#define OPTINTS(tag) getopt_ints(opts, tag)
#define getopt_flt(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_flt)
#define OPTFLT(tag) getopt_flt(opts, tag)
#define getopt_flts(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_flts)
#define OPTFLTS(tag) getopt_flts(opts, tag)
#define getopt_str(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_str)
#define OPTSTR(tag) getopt_str(opts, tag)
#define getopt_strs(pg, tag) (get_progopt(pg, tag, __FILE__, __FUNCTION__, __LINE__)->val_strs)
#define OPTSTRS(tag) getopt_strs(opts, tag)

static inline void free_progopt(PROGOPT *pg){
	prog_opt_t *opt;
	u4i i, j;
	for(i=1;i<pg->opts->size;i++){
		opt = ref_progoptv(pg->opts, i);
		if(opt->opt_tag) free(opt->opt_tag);
		if(opt->opt_dsc) free(opt->opt_dsc);
		if(opt->val_str) free(opt->val_str);
		if(opt->val_ints) free_b8v(opt->val_ints);
		if(opt->val_flts) free_f8v(opt->val_flts);
		if(opt->val_strs){
			for(j=0;j<opt->val_strs->size;j++){
				if(opt->val_strs->buffer[j]) free(opt->val_strs->buffer[j]);
			}
			free_cplist(opt->val_strs);
		}
	}
	free_progoptv(pg->opts);
	free_cuhash(pg->hash);
	free(pg);
}

#endif
