
/*

Extension Module for ASAP

Consider using multi-phase extension module initialization instead:
https://www.python.org/dev/peps/pep-0489/

*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include "asap.h"
#include "oldfns.h"

#define SIGN( a ) ( ( (a) > 0 )?1: ( ((a)==0)?0:-1)  )
static int Increase(const void *v1, const void *v2){  	return (int)SIGN( *((double *)v1) - *((double *)v2));  };
#undef SIGN

// Set var = dict[str], do nothing if key does not exist.
// On failure, sets error indicator and returns -1.
// Return 0 on success.
// Beware: parsing strings allocates memory!
int parseItem(PyObject *dict, const char *str, const char t, void *var) {

	PyObject *item;
	PyObject *value;
	switch(t){
		case 'b':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(int *)var = PyObject_IsTrue(item);
			if (*(int *)var == -1) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected boolean value for key '%s'", str);
				return -1;
			}
			break;
		case 'i':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(int *)var = (int) PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected integer value for key '%s'", str);
				return -1;
			}
			break;
		case 'd':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(double *)var = (double) PyFloat_AsDouble(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected double value for key '%s'", str);
				return -1;
			}
			break;
		case 'f':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(float *)var = (float) PyFloat_AsDouble(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected float value for key '%s'", str);
				return -1;
			}
			break;
		case 's':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			value = PyUnicode_AsEncodedString(item, "utf-8", "~E~");
			if (value == NULL) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected string value for key '%s'", str);
				return -1;
			}
			const char *bytes = PyBytes_AS_STRING(value);
			*(const char **)var = malloc(strlen (bytes) + 1);
			strcpy(*(const char **)var, bytes);
			Py_XDECREF(value);
			break;
		default:
			PyErr_Format(PyExc_TypeError, "parseItem: Unexpected type: %c", t);
			return -1;
		}
	return 0;
}

static PyObject *
asap_main(PyObject *self, PyObject *args) {

	PyObject *dict;
	PyObject *item;

	int withlogfile=0;
	int stdout_bak = -1;
	int stderr_bak = -1;
	fpos_t stdout_pos;
	fpos_t stderr_pos;

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
	     *dirfiles_default = ".",
	     *file_data,
//	      *nametree,
	     *fname,
	     *namegroups,
	     *simple_name;


//	char nametree[512];

	time_t t1, t2, t3, t5;

	char c;


	int imethode = 1, fmeg = 0, withallfiles = 0;//imethode1 for Jukes
	int last_node;

	float maxDist,
	      min,
	      ts_tv = 2.0;     /* default value for the trans/transv rates for Kimura 2-p */

	double best_score, echx, echy,max_score,min_score;

	int widthKlado;
	//float seuil_pvalue=0.05;

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

	// Accept a dictionary-like python object
	if (!PyArg_ParseTuple(args, "O", &dict))
		return NULL;
	if (!PyDict_Check(dict)) {
		PyErr_SetString(PyExc_TypeError, "asap_main: Argument must be a dictionary");
		return NULL;
	}

	if (parseItem(dict, "file", 's', &file_data)) return NULL;

	if (!file_data) {
		PyErr_SetString(PyExc_KeyError, "asap_main: Mandatory key: 'file'");
		return NULL;
	}

	f_in = fopen(file_data, "r");
	if (f_in==NULL) {
		PyErr_Format(PyExc_FileNotFoundError, "asap_main: Input file not found: '%s'", file_data);
		return NULL;
	}

	if (parseItem(dict, "out", 's', &dirfiles)) return NULL;
	if (!dirfiles) dirfiles = dirfiles_default;

	if (dirfiles[strlen(dirfiles)-1]!='/')
		strcat(dirfiles,"/");

	if (stat(dirfiles, &st) == -1) {
		mkdir(dirfiles, 0700);
	}

	if (parseItem(dict, "logfile", 'b', &withlogfile)) return NULL;

	if (withlogfile) {
		sprintf(file_data,"%s/asap.log",dirfiles);
		printf("> Redirecting stdout/stderr to file: %s\n", file_data);
		printf("BEFORE REDIRECT: \n");
		printf("stdout %d\n", stdout);
		printf("stderr %d\n", stderr);
		fflush(stdout);
		fflush(stderr);
		fgetpos(stdout, &stdout_pos);
		fgetpos(stderr, &stderr_pos);
		stdout_bak = dup(fileno(stdout));
		stderr_bak = dup(fileno(stderr));
		FILE *dout = freopen(file_data,"w",stdout);
		#ifdef _WIN32
		// Open stderr or else dup2 will fail for windowed app
		FILE *derr = freopen("NUL:","w",stderr);
		#endif
		int ddup = dup2(fileno(stdout), fileno(stderr));
		printf("AFTER REDIRECT: \n");
		printf("stdout %d\n", stdout);
		printf("stderr %d\n", stderr);
		printf("ddup %d\n", ddup);
		printf("stdout_bak %p\n", stdout_bak);
		printf("stderr_bak %p\n", stderr_bak);
		if ((dout == NULL) || (ddup < 0)) {
			PyErr_SetString(PyExc_SystemError, "asap_main: Failed to redirect output, aborting.");
			return NULL;
		}
	}

	/*
		Header
	*/
	fprintf(stderr, "/*\n");
	fprintf(stderr, "\tASAP (Agglomerate Specimens by Automatic Processing)\n");
	fprintf(stderr, "\twill delineate species in your dataset in a few moments.\n");
	fprintf(stderr, "\tRemember that the final cut remains yours!\n");
	fprintf(stderr, "*/\n");

	// Print these here so they are redirected if needed
	printf("> file = %s\n", file_data);
	printf("> dirfiles = %s\n", dirfiles);
	printf("> withlogfile = %i\n", withlogfile);

	if (parseItem(dict, "method", 'i', &imethode)) return NULL;
	printf("> imethode = %i\n", imethode);

	if (parseItem(dict, "sequence_length", 'i', &(asap_param.lenSeq))) return NULL;
	printf("> asap_param.lenSeq = %i\n", asap_param.lenSeq);

	if (parseItem(dict, "replicates", 'i', &(asap_param.replicates))) return NULL;
	printf("> asap_param.replicates = %i\n", asap_param.replicates);

	if (parseItem(dict, "seed", 'i', &seed_asap)) return NULL;
	printf("> seed_asap = %i\n", seed_asap);

	if (parseItem(dict, "rate", 'f', &ts_tv)) return NULL;
	printf("> ts_tv = %f\n", ts_tv);

	if (parseItem(dict, "all", 'b', &withallfiles)) return NULL;
	printf("> withallfiles = %i\n", withallfiles);

	if (parseItem(dict, "mega", 'b', &fmeg)) return NULL;
	printf("> fmeg = %i\n", fmeg);


	printf("\n> Begin ASAP core:\n\n");

	if (seed_asap== -1)
		srand( time(NULL) );
	else
		srand(seed_asap);

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

	Py_INCREF(Py_None);
	return Py_None;
}


static PyMethodDef AsapMethods[] = {
  {"main",  asap_main, METH_VARARGS,
   "Run ASAP for given parameters."},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyDoc_STRVAR(asap_doc,
"Assemble species by automatic partitioning.");

static struct PyModuleDef asapmodule = {
  PyModuleDef_HEAD_INIT,
  "asap",   /* name of module */
  asap_doc, /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module,
               or -1 if the module keeps state in global variables. */
  AsapMethods
};

PyMODINIT_FUNC
PyInit_asapc(void)
{
	PyObject *m = NULL;
  m = PyModule_Create(&asapmodule);
	// if (m != NULL) {
	// 	if (PyModule_AddStringConstant(m, "separator", "/")) {
	// 		Py_XDECREF(m);
	// 		m = NULL;
	// 	}
	// }
	return m;
}
