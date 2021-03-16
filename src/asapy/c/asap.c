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
#ifndef _WIN32
#include <unistd.h>
#include <dirent.h>
#endif
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>





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


void usage(char *arg)
{
	fprintf(stderr, "/*\n\tAgglomerate Specimens by Automatic Processing\n*/\n");
	fprintf(stderr, "syntax is '%s [-h] [options] distance_matrix or fasta file'\n", arg);
	fprintf(stderr, "\tfile is EITHER a distance matrix in phylip format OR aligned sequences in fasta format\n");

	fprintf(stderr,
	        "Options are:\n\
	\t-h    : this help\n\
	\t-r #  : nbr of replicates for statistical tests (default is 10^4)\n\
	\t-b #  : nbr of low-pvalues to be reported (clicable)\n\
	\t-m    : if present the distance Matrix is supposed to be MEGA CVS (other formats are guessed)\n\
	\t-a    : output all files: all probabilities, tree and graph files //  Better with -o option\n\
	\t-d #  : distance (0: Kimura-2P, 1: Jukes-Cantor --default--, 2: Tamura-Nei 3:simple distance)\n\
	\t-o #  : directory where results files are written (default is where the script is run)\n\
	\t-l #	: original length of seqs if a distance matrix was provided (default value 600)\n\
	\t-t #  : transition/transversion (for Kimura) default:2\n\
	\t-p #  : slope ponderation for asap score calculation: default is 0.5\n\
	\t-x #  : seed value\n");



	exit(1);


}

/*--------------------------------------------------*/
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


//	char nametree[512];

	time_t t1, t2, t3, t5;

	char c;


	short int imethode = 1,fmeg = 0, withallfiles = 0;//imethode1 for Jukes
	int last_node;

	float maxDist,
	      min,
	      ts_tv = 2.0;     /* default value for the trans/transv rates for Kimura 2-p */

	double best_score, echx, echy,max_score,min_score;

	int widthKlado;
	//float seuil_pvalue=0.05;

       extern char *optarg;           /* for options parisng */
        extern  int optind;

		float minAsapDist=0.005,maxAsapDist=0.05;


	struct stat st = {0};


	/*
		init
	 */
	dirfiles = NULL;
	t1 = time(NULL);

	asap_param.pond_pente=0.1;
	asap_param.pond_score=0.5;
	asap_param.replicates=1000;
	asap_param.seuil_pvalue=0.05;
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
	while ( (c = getopt(argc, argv, "o:l:p:d:amhr:b:")) != -1 ) {

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

			case 'm':
				fmeg = 1;			/*if present format mega CVS*/
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
		mat = compute_dis(f_in, imethode, ts_tv, &(asap_param.lenSeq),"",stdout);
	}
	else
	{
		fprintf(stderr, "assuming a distance matrix file. Reading the matrix\n");
		if (fmeg==0)
			mat = read_distmat(f_in, ts_tv, NULL, NULL);
		else
			fprintf(stderr, "MEGA format not yet implemented\n");
	}
	fclose(	f_in);
	fprintf(stderr,"End of matrix distance\n");
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

	/*
		Set the first n nodes to their id --the leaves--
	*/
	for (i = 0; i < mat.n; i++)
		no_node[i] = i;


	/*
		from the distance matrix, build a sorted list of pairwise_distance, min and max
	*/
	mattolist(ListDistance , &mat , &maxDist, &min);


	//	for (i=0;i<nb_pairs;i++)
	//		fprintf(stderr,"%d %f %d %d\n",i,ListDistance[i].d,ListDistance[i].a,ListDistance[i].b);


	nbresults = 0;

	last_node = mat.n - 1;

	/*
		Run ASAP core
	*/

	fprintf(stderr,"> asap is building and testing all partitions\n  ");

// int do_agglutine(DistMat mat, Composante *comp, DistPair *ListDist, Results *scores, Tabcompo *strucompo, int nb_pairs, FILE *f_out,double *best, int *fi, FILE *ff, Node *zenodes, int *list_node, int *lastnode, char *ledir, int lenSeq, int replicates,float seuil_pvalue,float pond_pente)


	//nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo, nb_pairs, f_out, &best_score, &firstpart, stderr, zenodes, no_node, &last_node, "", len_seq, replicates,seuil_pvalue,pond_pente);
	nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo,  &best_score, &firstpart,  zenodes, no_node, &last_node,asap_param);

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


	fprintf(stderr, "\n> 10 Best scores (probabilities evaluated with seq length:%d)\n",asap_param.lenSeq);
	fprintf(stderr, "  distance  #species   #spec w/rec  p-value pente score\n");
	int nb_B=(nbresults<10)?nbresults:10;
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

	if (withallfiles)
		ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue);

	printf("creating histo in %s\n",dirfiles);
	createSVGhisto(dirfiles,mat,20,scores, nbresults,WORKDIR_CL);

	/*
		That
	*/

	fprintf(stderr, "> asap is creating text and graphical output\n");
	qsort(scores,nbresults,sizeof (Results ),compareSpecies);


	fprintf(svgout, "<svg xmlns=\"http://www.w3.org/2000/svg\" onload=\"init(evt)\" ");
	fprintf(svgout, "width=\"%d\" height=\"%ld\" >\n", widthKlado + MARGECLADO + 20, HAUTEURCOURBE + MARGECLADO + ( mat.n * SIZEOFTEXT));


	CreateCurve2(scores, nbresults, dirfiles, simple_name, NULL, maxDist,  svgout,mat.n,max_score,min_score,widthKlado,minAsapDist,maxAsapDist);
	clearalltab(strucompo, &comp, mat.n);

	resetcomp(&comp, mat.n);

	echy = mat.n * SIZEOFTEXT;
	echx = widthKlado / (float)maxDist;

	print_clado(zenodes, last_node, NULL, echx, echy, (widthKlado - 100) / zenodes[last_node].round, 0,0);

	draw_clado(zenodes, svgout, last_node, mat.n,widthKlado);

	color_clado(zenodes, last_node,&color_ori);


//	draw_bestlines(zenodes, last_node,svgout, scores, echx, MARGECLADO + ( mat.n * SIZEOFTEXT) + HAUTEURCOURBE, nbresults,min_pvalue);

	fprintf(svgout, "</svg>\n");
	fclose(svgout);
	fgroups=fopen(namegroups,"w");
	if (fgroups==NULL)
	printf("Cant write %s\n",namegroups);
	else
	draw_nico(zenodes, fgroups, mat.n,scores,nbresults,asap_param.seuil_pvalue,10,last_node,widthKlado);
	fprintf(stderr, "> results were write\n");



	t5 = time(NULL);

if (withallfiles)
{
	// all files after some work
}

	fprintf(stderr, "  partition results are logged in %s and groups.txt and the graphic output is in %s and %s\n", fout,fname,namegroups);
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
