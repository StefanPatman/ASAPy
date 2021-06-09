/*ASAP:Agglomerate specimens by automatic process*/
/*
	Copyright (C) 2015-2016 G Achaz/ S Brouillet

	This progam is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 	for more information, please contact guillaume achaz <guillaume.achaz@mnhnfr>/<sophie.brouillet@mnhn.fr>

*/
/******
        file     : asap
        function : Agglomerate specimens by automatic process
                   command line version

        created  : Feb 2016


        author   : madamesophie


*****/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <float.h>
#include "asap.h"
#include "asap_core.h"
#include "oldfns.h"
//#include "drawMat.h"
#include "gdtosvg.h"
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <unistd.h>
#include <dirent.h>
#define strcasestr strstr
#endif


#define ASAP_CL


#ifdef MACOSX
#include <float.h>
#elif defined(SGI)
#include <limits.h>
#elif defined(LINUX)
#include <values.h>
#endif

#ifdef _PROTOTYPES
#undef _PROTOTYPES
#include <float.h>
#define _PROTOTYPES
#else
#include <float.h>
#endif
#define NBCHARMALLOC 256

void usage(char *arg)
{
	fprintf(stderr, "/*\n\tAgglomerate Specimens by Automatic Processing\n*/\n");
	fprintf(stderr, "syntax is '%s [-h] [options] distance_matrix or fasta file'\n", arg);
	fprintf(stderr, "\tfile is EITHER a distance matrix in phylip format OR aligned sequences in fasta format\n");

	fprintf(stderr,"Options are:\n\
	\t-h    : this help\n\
	\t-r #  : nbr of replicates for statistical tests (default is 10^4)\n\
	\t-b #  : nbr of low-pvalues to be reported (0.001 default)\n\
	\t-m    : if present the distance Matrix is supposed to be MEGA CVS (other formats than mega are guessed)\n\
	\t-a    : output all files: all probabilities, tree and graph files [Better with -o option]\n\
	\t-d #  : distance (0: Kimura-2P, 1: Jukes-Cantor --default--, 2: Tamura-Nei 3:simple distance)\n\
	\t-o #  : directory where results files are written (default is where the script is run)\n\
	\t-l #	: original length of seqs if a distance matrix was provided (default value 600)\n\
	\t-t #  : transition/transversion (for Kimura) default:2\n\
	\t-n #  : nbr of best scores to be kept (default 10))\n\
	\t-x #  : seed value\n");

/*	\t-p #  : slope ponderation for asap score calculation: default is 0.5\n\*/

	exit(1);


}

//returns the position of c in string l 1 to length(l) return 0 if not
int myIndex(char *l, char c)
{
int i,lo=strlen(l);

for (i=1;i<=lo;i++)
	if (l[i-1]==c)
	return(i);
return(0);
}


/*my own fgets which reallocs sizeof line if needed*/
char *my_get_line(char *ligne,FILE *f_in,int *nbcharmax)
{
char c;
int nbc=0;

	while (1)
		{
		 c=fgetc(f_in);
		// if (feof(f_in))
		 //	printf("EOF (1) detected wrong format\n"),exit(1);

		 if (c=='\n'|| c==10 || c=='\r' || feof(f_in)){
 			ligne[nbc]='\0';

 			break;
 			}
 		ligne[nbc++]=c;
 		if (nbc== *nbcharmax)
 			{
 			*nbcharmax= *(nbcharmax)+NBCHARMALLOC;
 			ligne=realloc(ligne, sizeof(char)*(*nbcharmax));
 			}
 		}


return(ligne);
}



int check_nbr(char *nb)
{
	int i, l=strlen(nb);
	//printf("%s-->%d\n",nb,l);
	for (i=0;i<l;i++)
		{
		//printf("car:%c %d -->%d\n",nb[i],nb[i], isdigit(nb[i]));
		if (isdigit(nb[i])) continue;
		if (nb[i]=='.')continue;
		if (nb[i]=='+')continue;
		if (nb[i]=='-')continue;
		if (nb[i]=='e')continue;
		if (nb[i]=='E')continue;
		if (nb[i]==' ')continue;
		return (0);
		}
return (1);
}





void read_mega10(FILE *f_in,struct DistanceMatrix *my_mat)
{
int nb=0,a,b,c;
int nbcharmax=NBCHARMALLOC,to_alloc=0;
char *ligne,letter,nombre [128];
int i,l,w=0;;
long ppos;
//float ff;
//long posit;

	printf("CVS MEGA X FILE\n");fflush(stdout);

	ligne=(char *)malloc(sizeof(char)*nbcharmax);
	*ligne='\0';


	ligne=my_get_line(ligne,f_in,&nbcharmax);
	l=strlen(ligne);
	//printf("ligne:%d char\n",l);
	for (i=0;i<l;i++)
		if (ligne[i]==',')
			nb++;

	nb=nb+1;
	printf("%ld seq\n",nb);fflush(stdout);

	my_mat->n = nb;

	my_mat->names = (char **)malloc( (size_t) sizeof(char *)*my_mat->n );
	if( ! my_mat->names )fprintf(stderr, "read_distmat: cannot allocate my_mat.names, bye"), exit(4);

/*	for(a=0;a<my_mat->n; a++){
		my_mat->names[a] = (char *)malloc( (size_t) sizeof(char)*(SIZE_NAME_DIST +1));
		if( ! my_mat->names[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat.names[%d], bye",a), exit(4);
	}
*/
	my_mat->dist = (double **)malloc( (size_t) sizeof(double *)*my_mat->n );
	if( ! my_mat->dist)fprintf(stderr, "read_distmat: cannot allocate my_mat.dist, bye"), exit(4);
	for(a=0;a<my_mat->n; a++){
		my_mat->dist[a] = (double *)malloc( (size_t) sizeof(double)*my_mat->n );
		if( ! my_mat->dist[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat.dist[%d], bye",a), exit(4);
		}


for (a=0;a<my_mat->n;a++){
		c=0;
		ppos=ftell(f_in);
		to_alloc=0;
		while( (letter=fgetc(f_in)) != ','){ //count length of title
		to_alloc++;
		}
		fseek(f_in,ppos,SEEK_SET);
		my_mat->names[a]=(char *)malloc(sizeof(char)*(to_alloc+2));
		while( (letter=fgetc(f_in)) != ','){
				my_mat->names[a][c] = (char)letter;
				c++;
			}
		//printf("%d %d\n",c,to_alloc);
		my_mat->names[a][c]='\0';
	//printf("--->%s\n",my_mat->names[a]);
		for (b=0;b<a;b++)
			{

			c=0;
			while( (letter=fgetc(f_in)) != ',' && !feof(f_in)){
				if (letter=='?'){
				fprintf(stderr,"**Warning distance between %s and %s is unknown,exiting<BR>\n",my_mat->names[a],my_mat->names[b]);exit(1);
				}

				nombre[c]=(char) letter;
				c++;

				}
			if (feof(f_in)) printf("your MEGA format can't be read. asap only reads MEGA 6 or MEGA X csv format\n"),exit(1);
	    	nombre[c]='\0';
	    	//
	    //printf("%s %d %d %s\n",my_mat->names[a],a,b,nombre);
	    	if (c==0)
	    		my_mat->dist[b][a]=my_mat->dist[a][b]=0;
	    	else
	    		{
	    		if(check_nbr(nombre)==0) {printf("your MEGA format can't be read. asap only reads MEGA 6 or MEGA X csv format\n");exit(1);}
				my_mat->dist[b][a]=my_mat->dist[a][b]=strtod(nombre,NULL);
				if (my_mat->dist[b][a]<0)
					{
					if (w==0){printf("negative distances have been found. Asap can't deal with that, your distance will be assumed to be 0\n");w=1;}
					my_mat->dist[b][a]=my_mat->dist[a][b]=0;
					}
				}
			}

		 my_mat->dist[a][a]=0;

		while (letter != 10  && letter!=13 && letter !='\n'&& !feof(f_in))/* go to end of line*/
			{letter=fgetc(f_in);}
		if (feof(f_in) && b!=a)
			printf("%d %d pb reading matrix CVS, asap only reads MEGA 6 or MEGA X csv format\n",a,b),exit(1);

		}
/*for (a=0;a<my_mat->n;a++)
{
printf("%s ",my_mat->names[a]);
for (b=0;b<my_mat->n;b++)
	printf("%5.4f ",my_mat->dist[a][b]);
printf("\n");
}	*/

free(ligne);

	}




/*Read CVS mega matrix which is the default for MEGA 5*/
void readMatrixMegaCVS(FILE *f_in,struct DistanceMatrix *my_mat)
{
int nb=0,a,b,c;
int nbcharmax=NBCHARMALLOC,to_alloc=0;
char *ligne,letter,nombre [128];
long ppos;
//float ff;
//long posit;

	printf("CVS MEGA FILE\n");fflush(stdout);



	ligne=(char *)malloc(sizeof(char)*nbcharmax);
		*ligne='\0';
	while (1)
		{
		ligne=my_get_line(ligne,f_in,&nbcharmax);
		//printf("%d ->%s\n",nb,ligne);
		if(strncmp(ligne,"Table",5)==0 || feof(f_in)) break;
		if (strlen(ligne)>2)
			{nb++;}

		//if(strncmp(ligne,"Table",5)==0 || feof(f_in)) break;
		}

	rewind(f_in);
	my_mat->n = nb;
		printf("%ld seq\n",my_mat->n);fflush(stdout);


	my_mat->names = (char **)malloc( (size_t) sizeof(char *)*my_mat->n );
	if( ! my_mat->names )fprintf(stderr, "read_distmat: cannot allocate my_mat.names, bye"), exit(4);

/*	for(a=0;a<my_mat->n; a++){
		my_mat->names[a] = (char *)malloc( (size_t) sizeof(char)*(SIZE_NAME_DIST +1));
		if( ! my_mat->names[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat.names[%d], bye",a), exit(4);
	}
*/
	my_mat->dist = (double **)malloc( (size_t) sizeof(double *)*my_mat->n );
	if( ! my_mat->dist)fprintf(stderr, "read_distmat: cannot allocate my_mat.dist, bye"), exit(4);
	for(a=0;a<my_mat->n; a++){
		my_mat->dist[a] = (double *)malloc( (size_t) sizeof(double)*my_mat->n );
		if( ! my_mat->dist[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat.dist[%d], bye",a), exit(4);
		}

/*now read */

for (a=0;a<my_mat->n;a++){
		c=0;
		ppos=ftell(f_in);
		to_alloc=0;
		while( (letter=fgetc(f_in)) != ','){ //count length of title
		to_alloc++;
		}
		fseek(f_in,ppos,SEEK_SET);
		my_mat->names[a]=(char *)malloc(sizeof(char)*(to_alloc+1));
		while( (letter=fgetc(f_in)) != ','){
				my_mat->names[a][c] = (char)letter;
				c++;
			}

		my_mat->names[a][c]='\0';
		//printf("--->%s\n",my_mat->names[a]);
		for (b=0;b<a;b++)
			{
			c=0;
			while( (letter=fgetc(f_in)) != ',' && !feof(f_in)){
				if (letter=='?'){
				fprintf(stderr,"**Warning distance between %s and %s is unknown,exiting<BR>\n",my_mat->names[a],my_mat->names[b]);exit(1);
				}

				nombre[c]=(char) letter;
				c++;
			}

	    	nombre[c]='\0';
	    	//
	    //printf("%s %d %d %s\n",my_mat->names[a],a,b,nombre);
	    	if (c==0)
	    		my_mat->dist[b][a]=my_mat->dist[a][b]=0;
	    	else
	    	{
	    	if(check_nbr(nombre)==0) {printf("your MEGA format can't be read. asap only reads MEGA 6 or MEGA X csv format\n");exit(1);}
			my_mat->dist[b][a]=my_mat->dist[a][b]=strtod(nombre,NULL);
			}
			}
		 my_mat->dist[a][a]=0;

	while (letter != 10  && letter!=13 && letter !='\n'&& !feof(f_in))/* go to end of line*/
		{letter=fgetc(f_in);}
	if (feof(f_in) && b!=a)
		printf("%d %d pb reading matrix CVS, asap only reads MEGA 6 or MEGA X csv format\n",a,b),exit(1);

	}
//for (a=0;a<my_mat->n;a++)

//printf("ok\n");
free(ligne);
//printf("all done\nRETURN");
}


/*MEGA matrix is a plague because output can be customize a lot..  */
void readMatrixMega(FILE *f_in,struct DistanceMatrix *my_mat)
{

	int a,b,nbc=0,c,n;

	char *ligne,letter,nombre[16];

//	int nbcol=0;;
	int lower=-1;
	int nbcharmax=NBCHARMALLOC;
	int lindex=0;


	ligne=(char *)malloc(sizeof(char)*nbcharmax);

	my_mat->n=0;
	my_mat->names=NULL;
	my_mat->dist=NULL;

	printf("Read Mega Format\n");

	//read the header
	while (1)
		 {
			fscanf(f_in,"%[^\n]\n",ligne);

#ifdef _WIN32
			// Naive solution to replace strcasestr: capitalize everything
			    char *s = ligne;
			    while (*s) {
			      *s = toupper((unsigned char) *s);
			      s++;
			    }
#endif

			if (feof(f_in)) printf("pb reading file...\n"),exit(1);

		 	if (strcasestr(ligne," OF TAXA :") != NULL)
				my_mat->n=atoi(strchr(ligne,':')+1);

			if (strcasestr(ligne,"NTAXA=") !=NULL)
				my_mat->n=atoi(strchr(strcasestr(ligne,"NTAXA="),'=')+1);

			if (strcasestr(ligne,"DATAFORMAT=")!=NULL)
				{
				if (strcasestr(ligne,"LOWERLEFT")!=NULL)
					lower=1;
				else
					if (strcasestr(ligne,"UPPERRIGHT")!=NULL)
						lower=0;
					else
					printf("Unknown data format\n"),exit(1);
				}
			if (*ligne!='!' && strchr(ligne,';'))// we have reach the species desc line
				break;

			}


	printf("%ld data\n",my_mat->n);

	if (my_mat->n ==0) printf("asap was not able to read your MEGA file: [TAXA] number not in the header\n"),exit(1);


	nbc=0;


//do some memory initialisation
	my_mat->names = (char **)malloc( sizeof(char *)* my_mat->n );
	if( ! my_mat->names )fprintf(stderr, "read_distmat: cannot allocate my_mat->names, bye"), exit(4);

/*	for(a=0;a<my_mat->n; a++){
		my_mat->names[a] = (char *)malloc( sizeof(char)*SIZE_NAME_DIST +1);
		if( ! my_mat->names[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat->names[%d], bye",a), exit(4);
	}*/

	my_mat->dist = (double **)malloc( sizeof(double *)* my_mat->n );
	if( ! my_mat->dist )fprintf(stderr, "read_distmat: cannot allocate my_mat->dist, bye"), exit(4);
	for(a=0;a<my_mat->n; a++){
		my_mat->dist[a] = (double *)malloc( sizeof(double)* my_mat->n );
		if( ! my_mat->dist[a] )
			fprintf(stderr, "read_distmat: cannot allocate my_mat->dist[%d], bye",a), exit(4);
		}


	a=0;


//read species name
	while (1)
		{
			lindex=0;
			do
				fscanf(f_in,"%[^\n]\n",ligne);
			while (strlen(ligne)<=1); //skip white lines if needed

			if (strlen(ligne)<=1) break;

 			if (strchr(ligne,'#')!=0)
 					lindex=myIndex(ligne,'#');
 				else
 					{
 					if (strchr(ligne,']'))
 				 		lindex=myIndex(ligne,']');
					else
 						lindex=0;//printf("cant read species \n"),exit(1);
 					}
 			n=strlen(ligne+lindex);
 			my_mat->names[a]= (char *)malloc( sizeof(char)*(n+1));
 			strncpy(my_mat->names[a],ligne+lindex,n);
 			my_mat->names[a][n]='\0';

 									/*names with ( stink */
 			if (strchr(my_mat->names[a],'('))
 				remplace(my_mat->names[a],'(','_');
 			if (strchr(my_mat->names[a],')'))
 				remplace(my_mat->names[a],')','_');


 			a++;

 			if (a==my_mat->n)
 				break;

		}



	do {
		letter=fgetc(f_in);
		if (feof(f_in)) printf("error reading values\n"),exit(1);
		}
	while (letter!=']');	//last line read should be very long but some empty lines occur ....


letter=fgetc(f_in); //be sure we areon  line 1 of matrix
for (a=0;a<my_mat->n;a++){
		c=0;
		while( letter != ']' && !feof(f_in)) //reading after the name.
			letter=fgetc(f_in);

		if (feof(f_in))printf("problem reading your file\n"),exit(1);

		for (b=0;b<=a;b++)
			{
			c=0;
			while( (letter=fgetc(f_in)) == ' ');
			if (feof(f_in) ) break;
			while ( (letter != ' ') && (letter!='\n') && (letter != 10 ) && (letter!=13) && (letter!='[')){
				if (letter==',') letter='.';
				if (letter=='?')
				{
				fprintf(stderr,"**Warning distance between %s and %s is unknown,exiting<BR>\n",my_mat->names[a],my_mat->names[b]);exit(1);
				}


				nombre[c]=(char) letter;
//				printf("%d %c ",letter,letter);
				c++;
				if (c>15) {printf("too much char %d \n",letter);break;}

				letter=fgetc(f_in);
				if (feof(f_in)) break;
				}
	    	nombre[c]='\0';
	    	if (c==0)
	    		my_mat->dist[b][a]=my_mat->dist[a][b]=0;
	    	else
				my_mat->dist[b][a]=my_mat->dist[a][b]=strtod(nombre,NULL);

			}

		while (letter != 10  && letter != ']'  && letter!=13 && letter !='\n'&& !feof(f_in))/* go to end of line*/
			{letter=fgetc(f_in);}
		if (a!=my_mat->n -1 && feof(f_in))
			printf("pb reading matrix CVS\n"),exit(1);

	}

	free(ligne);


}




/*--------------------------------------------------*/
#ifndef ismodule

int main(int argc, char**argv)
{

	/* a bunch of beautiful structures usefull for the code */

	DistMat mat;               /* The matrix distance */
	DistPair *ListDistance;    /* distance for all sequence pairs */
	Composante comp;           /* The whole graph that describes the species */
	Results *scores;
	Tabcompo *strucompo;       /* each elemnt store how many groups and how many sequences in each group */
	Node *zenodes;             /* Nodes of the hierarchical clusterting tree */
	Parameter asap_param;  		/*stuff for asap*/
Spart *myspar;
	int i,
//		*grp,
//	    nb_pairs,
	    nbresults = 0,
	    firstpart,
//	    n_best=0,
	    color_ori=5,

	    seed_asap=-1,
	    *no_node;       // report for each sequence, its current node


	FILE *f_in,

	     *fgroups,
//	 	 *ffgroups,
	     *svgout;

	char *fout,
	     *dirfiles,
	     *file_data,
//	      *nametree,
	     *fname,
	     *namegroups,
	     *simple_name;
	char *meth[5]={ "K80_Kimura","JC69_Jukes-Cantor","N93_Tamura-Nei" , "Simple_Dist"};

//	char nametree[512];

	time_t t1, t2, t3, t5;

	char c;
		char thedate[80];

	short int imethode = 1,fmeg = 0, withallfiles = 0;//imethode1 for Jukes
	int last_node;
	//int fmeg2=0;
	float maxDist,
	      min,
	      ts_tv = 2.0;     /* default value for the trans/transv rates for Kimura 2-p */
	int nbBestAsap=10;
	double best_score, echx, echy,max_score,min_score;

	int widthKlado;
	//float seuil_pvalue=0.05;

       extern char *optarg;           /* for options parisng */
        extern  int optind;

		float minAsapDist=0.005,maxAsapDist=0.05;
struct stat     statbuf;
	struct tm      *tm;

	struct stat st = {0};

	stat(argv[0], &statbuf);
    tm = localtime(&statbuf.st_mtime);
	 strftime(thedate,80,"%FT%T", tm); // 11/19/20 - 05:34PM



	/*
		init
	 */
	dirfiles = NULL;
	t1 = time(NULL);

	asap_param.pond_pente=0.1;
	asap_param.pond_score=0.5;
	asap_param.replicates=1000;
	asap_param.seuil_pvalue=0.001;
	//asap_param.ledir="";
	asap_param.fres=stderr;
	asap_param.lenSeq=600;
	/*
		Header
	*/
	fprintf(stderr, "/*\n");
	fprintf(stderr, "\tASAP (Agglomerate Specimens by Automatic Processing)\n");
	fprintf(stderr, "\twill delineate species in your dataset in a few moments.\n");
	fprintf(stderr, "\tRemember that the final cut remains yours!\n");
	fprintf(stderr, "*/\n");


	/*
		parse options
	*/
	while ( (c = getopt(argc, argv, "o:l:n:p:d:amhr:b:")) != -1 ) {

		switch (c) {
			case 'a':
				withallfiles = 1;                     /* all files are output  default is just graphic files */
				break;

			case 'd':
				imethode = atoi(optarg);              /* nbr choosing dist method */
				break;

			case 'o':								/*dir where results files are written*/
				dirfiles = malloc((strlen(optarg) + 2) * sizeof(char));
				strcpy(dirfiles, optarg);
				if (dirfiles[strlen(dirfiles)-1]!='/')
					strcat(dirfiles,"/");

				if (stat(dirfiles, &st) == -1) {
					mkdir(dirfiles, 0700);
				}

				break;

			case 'h':
				usage(argv[0]);
				break;

			case 'l':
				asap_param.lenSeq=atoi(optarg);
				break;

			case 't':
				ts_tv = atof(optarg);		/* trans/trav rate */
				break;

			case 'n':
				nbBestAsap = atoi(optarg);		/* nb scores */
				break;


			case 'm':
				fmeg = 1;			/*if present format mega CVS*/
				break;

			case 'M':
				fmeg = 2;			/*if present format mega 5*/
				break;


			case 'r':
				asap_param.replicates = atoi(optarg);			/* for statistical testing */
				break;

			case 'b':
				asap_param.seuil_pvalue = atof(optarg);			/* limit for results to be reported */
				break;

			case 'p':
				asap_param.pond_pente = atof(optarg);			/* limit for results to be reported */
				break;

			case 'x':
				seed_asap=atoi(optarg); /* give a seed */

			default:
				usage(argv[0]);
		}

	}

	if (argc - optind != 1)usage(argv[0]), exit(1);

	if (seed_asap== -1)
	srand( time(NULL) );
	else
	srand(seed_asap);
	file_data = argv[optind];

	if (strrchr(file_data,'/')!=NULL)
	{
		int len_sim_nam= strrchr(file_data,'/') - file_data;
		simple_name=malloc(sizeof(char)*(len_sim_nam+1));
		sprintf(simple_name,"%.*s",len_sim_nam,strrchr(file_data,'/')+1);
//		nametree=malloc(sizeof(char)*(len_sim_nam+1));
	}
	else
	{simple_name=malloc(sizeof(char)*(strlen(file_data)+1));sprintf(simple_name,"%s",file_data);
//	nametree=malloc(sizeof(char)*(strlen(file_data)+1));
	}


	if (dirfiles == NULL)
	{
		dirfiles = (char *) malloc( (size_t) sizeof(char) * 3);
		if(!dirfiles)fprintf(stderr, "main: cannot allocate dirfiles bye\n"), exit(2);
		dirfiles[0] = '.'; dirfiles[1] = '/';dirfiles[2] = '\0';


	}

	/*

	*/
	f_in = fopen(file_data, "r");
	if (f_in == NULL)fprintf(stderr,"cannot open the file_data, bye\n"), exit(1);

	fout = (char * )malloc( (size_t) sizeof(char) * (strlen(simple_name)+ strlen(dirfiles) + 5));
	if (!fout)fprintf(stderr, "main: cannot allocate fout bye\n"), exit(2);

	namegroups=malloc(sizeof(char)*( (strlen (dirfiles) + strlen (simple_name) +20)));
	sprintf(namegroups,"%s%s.groups.svg",dirfiles, simple_name);
	sprintf(fout, "%s%s.all", dirfiles, simple_name);

	asap_param.f_out = fopen(fout, "w+");
	asap_param.fres=stdout;
	asap_param.web=0;
	if (asap_param.f_out == NULL)fprintf(stderr,"cannot open the output file %s, bye\n", fout), exit(1);

	fname = (char *) malloc( (size_t) sizeof(char) * (strlen (dirfiles) + strlen (simple_name) +5) );
	sprintf(fname, "%s%s.svg", dirfiles, simple_name);

	svgout = fopen(fname, "w");
	if (svgout == NULL)fprintf(stderr, "cannot open the graphic output file %s, bye\n", fname), exit(1);


	/*
		Read or build the distance matrix
	*/
	c = fgetc(f_in);
	rewind(f_in);
	if ( c == '>'){
		fprintf(stderr, "> asap is reading the fasta file and computing the distance matrix\n");


	//mat = compute_dis(f_in, imethode, ts_tv, &(asap_param.lenSeq),"",stdout);
		mat = compute_dis(f_in, imethode, ts_tv, &(asap_param.lenSeq),asap_param);
	}
	else
	{
		fprintf(stderr, "assuming a distance matrix file. Reading the matrix\n");
		if (fmeg==0)
			mat = read_distmat(f_in, ts_tv, NULL, NULL);
		else
			if (fmeg==1)
				{
					c=fgetc(f_in);
					fputc(c,f_in);

					if (c==',')
					{
						read_mega10(f_in,&mat);printf("done 10\n");
					}
					else
					readMatrixMegaCVS(f_in,&mat);
				}
			//else
				//readMatrixMega(f_in,&mat);


	}
	fclose(	f_in);
	fprintf(stderr,"End of matrix distance\n");
	myspar=malloc(sizeof(Spart)*mat.n);

	/*for (i=0;i<mat.n;i++)
		{
			//myspar[i].name=malloc(sizeof(char)*strlen( mat.names[i])+1);
			//strcpy(myspar[i].name,mat.names[i]);

			myspar[i].specie=malloc(sizeof(int)* nbBestAsap);

		}*/

	if (mat.n<MAXSPECIESGRAPH)
	widthKlado=WIDTHCLADO/3;
	else
	widthKlado=WIDTHCLADO;
	fprintf(stderr, "  %ld input sequences\n", mat.n);

	t2 = time(NULL);

	/*
		Get memory for needed struct
	*/

	asap_param.nbpairs = (mat.n * (mat.n - 1)) / 2;

	ListDistance = (DistPair *) malloc( (size_t) sizeof(DistPair) *  asap_param.nbpairs);
	if (!ListDistance)fprintf(stderr, "main: cannot allocate  ListDistance bye\n"), exit(2);
//
	no_node = (int *)malloc( (size_t) sizeof(int) * mat.n);               //indice qui me donne pour une seuqnec i a quel noeud elle est liée à chaque etape de lagglutnage
	if (!no_node)fprintf(stderr, "main: MEMORY ERROR error can allocate nonode bye\n"), exit(2);

	zenodes = (Node *) malloc( (size_t) sizeof(Node) * ((mat.n*2)-1));
	if (!zenodes)fprintf(stderr, "main: MEMORY ERROR error can allocate  zenodes bye\n"), exit(2);


	strucompo = (Tabcompo *) malloc( (size_t) sizeof(Tabcompo) * mat.n);
	if (!strucompo)fprintf(stderr, "main: cannot allocate  strucompo bye\n"), exit(2);


	scores = (Results *) malloc(  (size_t) sizeof(Results) * mat.n); //not enough juste a first dim reajusted when nbresulkts greater than mat.n
	if (!scores)fprintf(stderr, "main: cannot allocate  scores bye\n"), exit(2);
	//build the sorted struct from smallest dist to higher dist
for (i=0;i<mat.n;i++)
			scores[i].listNodes=malloc(sizeof(int)*mat.n);


	initcomp(&comp, mat.n, stderr, "");
	inittabcompo(strucompo, mat.n, stderr, "");       /* the structures are oversized currently */
	initNodes(stderr, zenodes, mat, "");
//FILE *fdeb;
//	fdeb=fopen("/Applications/MAMP/htdocs/temp/debug.txt","w");
	/*
		Set the first n nodes to their id --the leaves--
	*/
	for (i = 0; i < mat.n; i++)
		no_node[i] = i;

//int nb_pairs=(mat.n*(mat.n-1))/2;
	/*
		from the distance matrix, build a sorted list of pairwise_distance, min and max
	*/
	mattolist(ListDistance , &mat , &maxDist, &min);


		//for (i=0;i<nb_pairs;i++)
		//	fprintf(fdeb,"%d %f %d %d\n",i,ListDistance[i].d,ListDistance[i].a,ListDistance[i].b);


	nbresults = 0;

	last_node = mat.n - 1;

	/*
		Run ASAP core
	*/

	fprintf(stderr,"> asap is building and testing all partitions\n  ");

// int do_agglutine(DistMat mat, Composante *comp, DistPair *ListDist, Results *scores, Tabcompo *strucompo, int nb_pairs, FILE *f_out,double *best, int *fi, FILE *ff, Node *zenodes, int *list_node, int *lastnode, char *ledir, int lenSeq, int replicates,float seuil_pvalue,float pond_pente)


	//nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo, nb_pairs, f_out, &best_score, &firstpart, stderr, zenodes, no_node, &last_node, "", len_seq, replicates,seuil_pvalue,pond_pente);
	nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo,  &best_score, &firstpart,  zenodes, no_node, &last_node,asap_param);


		/*if (fdeb!=NULL)
		{
		fprintf(fdeb,"%d res\n",nbresults);
		for (i=0;i<nbresults;i++)
			fprintf(fdeb,"%d %d %d\n",i,scores[i].nbspec,scores[i].nbspecRec);
		fclose (fdeb);
	}*/

	qsort(scores,nbresults,sizeof (Results ),compareProba);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_proba=i+1;

	qsort(scores,nbresults,sizeof (Results ),compareParameter);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_pente=i+1;


	max_score=0;
	scores[0].score=(scores[0].rank_pente*(1.0-asap_param.pond_score))+(scores[0].rank_proba*asap_param.pond_score);
	min_score=scores[0].score;
	for (i = 0; i < nbresults+1; i++)
	{
		scores[i].score=(scores[i].rank_pente*(1.0-asap_param.pond_score))+(scores[i].rank_proba*asap_param.pond_score);
		if (max_score <scores[i].score)
			max_score =scores[i].score;
		if (min_score >scores[i].score)
			min_score =scores[i].score;
	}
	qsort(scores,nbresults,sizeof (Results ),compareRang);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_general=i+1;


	fprintf(stderr, "\n> %d Best scores (probabilities evaluated with seq length:%d)\n",nbBestAsap,asap_param.lenSeq);
	fprintf(stderr, "  distance  #species   #spec w/rec  p-value pente score\n");
	int nb_B=(nbresults<nbBestAsap)?nbresults:nbBestAsap;
	for (i = 0; i < nb_B; i++)

			{
						char toStar=' ';
					if (scores[i].d_jump>=minAsapDist && scores[i].d_jump<=maxAsapDist)
						toStar='*';

			 fprintf(stderr, "%c%8.4f %8d  %12d  %.3e %e \t%f \n",
			    	toStar,
			 		scores[i].d_jump,
			 		scores[i].nbspec,
			      	scores[i].nbspecRec,
			       	scores[i].proba,

			       	 scores[i].other_parameter *100,
			       	 scores[i].score);
			}

	/*if (withallfiles)
		//ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue);
	ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue,myspar,mat.n);*/
	printf("creating histo in %s\n",dirfiles);
	//createSVGhisto(dirfiles,mat,20,scores, nbresults,WORKDIR_CL);
	createSVGhisto(dirfiles,mat,20,scores, nbresults,"");

	/*
		That
	*/

	fprintf(stderr, "> asap is creating text and graphical output\n");
	qsort(scores,nbresults,sizeof (Results ),compareSpecies);


	fprintf(svgout, "<svg xmlns=\"http://www.w3.org/2000/svg\" onload=\"init(evt)\" ");
	fprintf(svgout, "width=\"%d\" height=\"%ld\" >\n", widthKlado + MARGECLADO + 20, HAUTEURCOURBE + MARGECLADO + ( mat.n * SIZEOFTEXT));


	CreateCurve2(scores, nbresults, dirfiles, simple_name, NULL, maxDist,  svgout,mat.n,max_score,min_score,widthKlado,minAsapDist,maxAsapDist);
	clearalltab(strucompo, &comp, mat.n);



	echy = mat.n * SIZEOFTEXT;
	echx = widthKlado / (float)maxDist;

	print_clado(zenodes, last_node, NULL, echx, echy, (widthKlado - 100) / zenodes[last_node].round, 0,0);

	draw_clado(zenodes, svgout, last_node, mat.n,widthKlado);

	color_clado(zenodes, last_node,&color_ori);


	fprintf(svgout, "</svg>\n");
	fclose(svgout);
	fgroups=fopen(namegroups,"w");
	if (fgroups==NULL)
	printf("Cant write %s\n",namegroups);
	else
	draw_nico(zenodes, fgroups, mat.n,scores,nbresults,asap_param.seuil_pvalue,10,last_node,widthKlado);
//printf("write \n");
	for (i=0;i<mat.n;i++)
		{
			int n=zenodes[i].first_to_draw; //assign ed in print_clado
			if (n>=mat.n)printf("error***** %d \n",n);
			myspar[n].name=malloc(sizeof(char)*strlen( mat.names[i])+1);
			strcpy(myspar[n].name,mat.names[i]);
			//printf("i:%d %d %s %s\n",i,n,myspar[n].name,mat.names[i]);
			myspar[i].specie=malloc(sizeof(int)* nb_B+1);

		}
		/*for (i=0;i<mat.n;i++)
		printf("i:%d %s\n",i,myspar[i].name);*/
	//qsort(scores,nbresults,sizeof (Results ),compareRang);
if (withallfiles)
{
	int **o_sp;
		//ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue);
	qsort(scores,nbresults,sizeof (Results ),compareRang);


	o_sp=malloc(sizeof(int*)*nb_B);
	for (i=0;i<nb_B;i++)
		o_sp[i]=malloc(sizeof(int)*2);
	//fprintf(stderr,"go ecrit\n");
	ecrit_fichier_texte( dirfiles,nb_B,nbresults, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue,myspar,mat.n,last_node);
	//fprintf(stderr,"go order\n");
	order_spart(o_sp,nb_B,myspar,mat.n);
	//fprintf(stderr,"go create\n");
	CreateSpartFile(myspar,dirfiles,nb_B,simple_name,stdout,scores,mat.n,thedate,"",meth[imethode],o_sp);
	for (i=0;i<nb_B;i++)
		free(o_sp[i]);
	free(o_sp);
}
	fprintf(stderr, "> results were write \n");

	t5 = time(NULL);

	resetcomp(&comp, mat.n);
	fprintf(stderr, "  partition results are logged in: \n");
	fprintf(stderr, "\tFull results: %s\n",fout);
	fprintf(stderr, "\tThe csv file of rank x: %sgroupe_x\n",dirfiles);
	fprintf(stderr, "\tThe result file of rank x: %sgroup_x\n",dirfiles);
	fprintf(stderr, "\tThe graphic outputs are: %s*.svg\n",dirfiles);
	fprintf(stderr, "\tSpart file is: %s%s.spart\n", dirfiles,simple_name);

	free (fout);
//	free(fileNex);
//	free(newickStringOriginal);
//	free(newickString);
	free (ListDistance);
	freecomp(&comp, mat.n);
	free_distmat(mat);



	t3 = time(NULL);

	fprintf(stderr, "> asap computation times were:\n");

	fprintf(stderr,"  %2ldm %2lds to read file and compute distance\n", (t2 - t1) / 60, (t2 - t1) % 60);
	fprintf(stderr,"  %2ldm %2lds to compute and test all partitions\n", (t3 - t2) / 60, (t3 - t2) % 60);
//	fprintf(stderr,"  %2ldm %2lds to build an nj tree\n", (t5 - t4) / 60, (t5 - t4) % 60);
	fprintf(stderr,"  --------------\n");
	fprintf(stderr,"  %2ldm %2lds total\n", (t5 - t1) / 60, (t5 - t1) % 60);

//	fclose(	f_out);

	return 0;
}

#endif
