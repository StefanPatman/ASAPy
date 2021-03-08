#ifndef OLDFNS_H
#define OLDFNS_H

void print_distmat(  struct DistanceMatrix distmat  );
void html_error(int n);
void f_html_error(int nb, char *ledir,FILE *f);

void distancesimple(struct FastaSeq *mesSeqs,int l,struct  DistanceMatrix  my_mat, char *ledir, FILE *fres);
void distanceJC69 (struct FastaSeq *mesSeqs, int l, struct  DistanceMatrix  mymat, char *ledir, FILE *fres);
void distanceK80 (struct FastaSeq *mesSeqs,int l,struct  DistanceMatrix  my_mat, char *ledir, FILE *fres);
void distanceTN93(struct FastaSeq *mesSeqs,int l,struct  DistanceMatrix  my_mat, char *ledir, FILE *fres);
void clean_str(char *ch);
void print_groups_newick( Composante my_comp, DistMat mat  , char *lastring, FILE *f2, char *ledir,FILE *fres);
void free_distmat(  struct DistanceMatrix mat );
void transition_transversion_sequences(char *seq1, char *seq2, long L, long *tsi, long *tsv);
void readMatrixMega_string(char *data,struct DistanceMatrix *my_mat,char *ledir,FILE *fres);
void readMatrixMegaCVS_string(char *data,struct DistanceMatrix *my_mat,char *ledir,FILE *fres);
void remplace(char *name,char c,char newc);

int check_names(struct FastaSeq *mesSeq, int nbseq,char *ledir,FILE *fres);
int ReadFastaSequence( FILE *f, struct FastaSeq *laseq,int *len_seq);
int check_compat(char *s1,char *s2,int l);
int search_delim(char *cgiinput, char *delim,char *ledir,FILE *fres);
int check_valid_name(char *name);

long del_sequences(char *seq1, char *seq2, long L);

short compare_DNA( char s1, char s2 );

struct DistanceMatrix compute_dis(FILE *f,int method,float ts_tv, int *len_seq, char *ledir, FILE *fres);
struct DistanceMatrix GetDistMat (int nseq, struct FastaSeq *mesSeqs, int method,float ts_tv, char *ledir, FILE *fres);
struct DistanceMatrix read_distmat(  FILE *file ,float ts_tv,char *ledir,FILE *fres);
struct DistanceMatrix read_distmat_string( char *data ,int fmega,char *ledir,FILE *fres);
struct DistanceMatrix  read_fasta_and_compute_dis(char *input,int method,float ts_tv,char *ledir,FILE *fres, int *ll);

char IsTransition( char nt1, char nt2 );
char IsTransversion( char nt1, char nt2 );
char **getMultPartData(int *nb,char *nf,char *ledir,FILE *fres);

double P_given_t_R( double t, double R );
double Q_given_t_R( double t, double R );
double compute_logL_given_t_R( long nsites, long n_tsv, long n_tsi, double t, double R );
double compute_k80( long nsites, long n_tsv, long n_tsi );
double find_ML_t_given_R( double R, long nsites, long n_tsv, long n_tsi );



#endif