/*
	Copyright (C) 2015-2016 G Achaz/ S Brouillet

	This program is free software; you can redistribute it and/or
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

 	for more information, please contact guillaume achaz <guillaume.achaz@mnhn.fr>/<sophie.brouillet@mnhn.fr>

*/
/******
        file     : oldfns
        function : All fns for asap , common with ABGD slightly changed because need to write errors on the results html file

        created  : November 2015


        author   : madamesophie


*****/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>  /* errno */
#include "asap.h"
#include "oldfns.h"




/*check if we have at least one common symbol beetween the 2 seqs*/
/*--------------------------------------------*/
void print_distmat(  struct DistanceMatrix distmat  ) {

	int a, b;
	printf("Distance Matrix:\n");
	for (a = 0; a < distmat.n; a++) {

		printf("[%d]%.10s", a + 1, distmat.names[a]);

		for (b = 0; b < distmat.n; b++) {
			printf("  %f", distmat.dist[a][b]);
		}

		printf("<BR>\n");

	}


}

char *ltrim(char *s)
{
    while(isspace(*s)) s++;
    return s;
}

char *rtrim(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *trim(char *s)
{
    return rtrim(ltrim(s)); 
}

void f_html_error(int nb,char *ledir,FILE *f )
{
fprintf(f,"******ASAP Error %d ****<BR>",nb);
switch (nb)
	{
	case 1:
	case 2:
	case 3:
	case 4:
	
	case 5:
	
	case 10:	
	fprintf(f,"Pb while reading CGI data. If this problem persists please contact <A HREF=\"mailto:%s\"> web site administrator</A>\n",WEB_ADMIN) ;
	break;
	
	case 8:
	fprintf(f,"Please try to rename your file with \".txt\" extension. If the problem persists or if the file extension is already \".txt\", please contact <A HREF=\"mailto:%s\"> web site administrator</A>\n",WEB_ADMIN) ;

	case 11:
	fprintf(f,"Your data appears to be a matrix distance but some characters have been found instead of numerical value\n") ;
	fprintf(f,"(inf or sup matrix are not supported _ yet_)\n") ;	break;
	
	case 20:
	case 56:
	case 144:
	fprintf(f,"can't open file. please contact <A HREF=\"mailto:%s\"> web site administrator</A>\n",WEB_ADMIN) ;
	break;

	case 60:
	fprintf(f,"Your data contains at least one other symbol than WSYKVHDB<BR>Please correct it\n");
	break;
	
    case 50:
    case 66:
	fprintf(f,"Weird FASTA format . Please check your data\n");
    break;       

    case 95:    
	fprintf(f,"Please check the P value you entered\n");
    break;   
    
    case 99:    
	fprintf(f,"FASTA sequences not all same size . Please check your data\n");
    break;   
    
    case 100:
	fprintf(f,"Length of seq appears to be 0\n");
    break; 
    
  	case 105:
	fprintf(f,"Either you have a CVS MEGA matrix and didn't check the option in previous page, <BR>either your matrix is not well formated\n");
    break;   	
    
	case 110:
	fprintf(f,"ASAP can't deal with sequences which names are only composed of digits, or with \": , ( )\" present in names. <BR>Please add at least one alphabetic character in sequences names or remove others\n");
	break;
	
    case 111:
	case 112:
    case 155:
    case 55:
    case 44:
   fprintf(f,"Memory pb. can't malloc. Try with smaller data file\n");
    break; 
   
   	case 177:
     fprintf(f,"Minimum number is 2. Bye\n");
    break; 
 
   
    case 200:
    case 300:
	fprintf(f,"pb while choosing distance model\n");
    break;
    case 255:
   fprintf(f,"newick file has some pb.. This should never arrives so.. try again\n");
    break;
    
    case 257:
   fprintf(f,"MEGA distance file not well formated. Please resave your file with MEGA\n");
    break;
    
    case 344:
  
   fprintf(f,"MEGA CSV distance file not well formated. Please resave your file with MEGA\n");
    break;
     
    case 343: 
    fprintf(f,"Pb reading MEGA CSV distance file . Keyword \"Table\" not found at bottom of the file\n");
    break;
    
    case 888: //error is already written
    break; 
    
    default:
  fprintf(f,"Undocumented error \n");
    }
fclose (f), exit_properly(ledir);
 
 }



/*--------------------------------------------*/
void html_error(int n)
{
	printf("error %d\n", n), exit(1);
	
	
	
}
/*check that names are uniques in alignment*/

/*--------------------------------------------*/
/*check that names are uniques in alignment*/
int check_names(struct FastaSeq *mesSeq, int nbseq, char *ledir, FILE *fres)
{
	int i, j;
	for (i = 0; i < nbseq - 1; i++)
		for (j = i + 1; j < nbseq; j++)
			if (strcmp(mesSeq[i].name, mesSeq[j].name) == 0)
			{
				fprintf(fres,"seq %d: %s and seq %d: %s have same name\n", i + 1, mesSeq[i].name, j + 1, mesSeq[j].name);
				fclose (fres), exit_properly(ledir);

			}
	return (1);
}

/*--------------------------------------------*/
long del_sequences(char *seq1, char *seq2, long L) {

	long i, del = 0;

	for (i = 0; i < L; i++)
		if ( *(seq1 + i) == '-' || *(seq2 + i) == '-')
			del++;

	return del;
}
/********************

SBGD FNS
*********************/

short compare_DNA( char s1, char s2 ) {

	/*
		if we are in standard DNA
	*/
	if (  (s1 == 'A' || s1 == 'C' || s1 == 'G' || s1 == 'T') &&  (s2 == 'A' || s2 == 'C' || s2 == 'G' || s2 == 'T') ) {
		if	(s1 == s2 )return 1;
		else			return 0;
	}

	/*
		What about deletion
	*/

	if (s1 == '-' || s2 == '-' || s1 == '+' || s2 == '+') {                  /* consider deletion as a diff only if the alignemt has been recoded */

		if (  (s1 == '-' && s2 == '+') || (s2 == '-' && s1 == '+') )
			return 0;
		else
			return 1;
	}

	/*
		if there is an N
	*/
	if ( s1 == 'N' || s2 == 'N')
		return 1;


	/*
		if only one is an A
	*/
	if (s1 == 'A') {
		if (  s2 == 'M' || s2 == 'R' || s2 == 'W' || s2 == 'V' || s2 == 'H' || s2 == 'D' )	return 1;
		else 							return 0;
	}
	if (s2 == 'A') {
		if (s1 == 'M' || s1 == 'R' || s1 == 'W' || s1 == 'V' || s1 == 'H' || s1 == 'D')	return 1;
		else 							return 0;
	}

	/*
		if only one is a C
	*/
	if (s1 == 'C') {
		if ( s2 == 'M' || s2 == 'S' || s2 == 'Y' || s2 == 'V' || s2 == 'H' || s2 == 'B' )	return 1;
		else 							return 0;
	}
	if (s2 == 'C') {
		if (  s1 == 'M' || s1 == 'S' || s1 == 'Y' || s1 == 'V' || s1 == 'H' || s1 == 'B' )	return 1;
		else 							return 0;
	}

	/*
		if only one is a G
	*/
	if (s1 == 'G') {
		if ( s2 == 'R' || s2 == 'S' || s2 == 'K' || s2 == 'V' || s2 == 'D' || s2 == 'B' )	return 1;
		else 							return 0;
	}
	if (s2 == 'G') {
		if ( s1 == 'R' || s1 == 'S' || s1 == 'K' || s1 == 'V' || s1 == 'D' || s1 == 'B' )	return 1;
		else 							return 0;
	}

	/*
		if only one is a T
	*/
	if (s1 == 'T') {
		if ( s2 == 'W' || s2 == 'Y' || s2 == 'K' || s2 == 'H' || s2 == 'D' || s2 == 'B' )	return 1;
		else 							return 0;
	}
	if (s2 == 'T') {
		if ( s1 == 'W' || s1 == 'Y' || s1 == 'K' || s1 == 'H' || s1 == 'D' || s1 == 'B' )	return 1;
		else 							return 0;
	}


	/*
		More complex case
	*/
	if ( (s1 == 'M' && s2 == 'K') || (s1 == 'K' && s2 == 'M')  )
		return 0;

	if ( (s1 == 'R' && s2 == 'Y') || (s1 == 'Y' && s2 == 'R')  )
		return 0;

	if ( (s1 == 'W' && s2 == 'S') || (s1 == 'S' && s2 == 'W')  )
		return 0;


	/*
		Else there is at least one overlap
	*/
	return 1;

}
/********************

	Distance Matrix

*********************/

/*Read one fasta sequence in a file pointer store it in a fastaseq struct
returns 0 if some pbs or some pbs and 1 if everything ok*/
int ReadFastaSequence( FILE *f, struct FastaSeq *laseq, int *lenseq)
{
	char *name;
	char *seq;
	int   n, c;
	int   nalloc;
	char *nucs = "ATGC-+NMRWSYKVHDB";
	nalloc = 128;
	c = fgetc(f);
	if (c != '>')
	{return 0;}
	n = 0;
	name = malloc(sizeof(char) * 128);
	while (1) {
		c = fgetc(f);
		if (c == '\n' || c == '\r' || c == 10 || c == 13  ) //do not store weird chars in names
			break;
		name[n++] = c;
		if (n > 127)
		{
			nalloc += 128;
			name = realloc(name, sizeof(char) * nalloc);
		}

	}

	name[n] = '\0';

	laseq->name = malloc(sizeof(char) * (n + 1));
	strcpy(laseq->name, name);

	seq = malloc(sizeof(char) * 128);      /* allocate seq in blocks of 128 residues */
	nalloc = 128;
	n = 0;

	while (1)
	{
		c = fgetc(f);
		if (c == EOF )
			break;
		if (c == '>' )
		{ungetc(c, f); break;} //put back in the stream the next new seq indicator
		if ( c != '\n' && c != '\r' && c != '\t' && c != ' ')
		{
			if (strchr(nucs, toupper(c)) == NULL) html_error(60); /*weird symbol found*/

			seq[n++] = toupper(c);
			if (nalloc == n)
			{
				nalloc += 128;
				seq = realloc(seq, sizeof(char) * nalloc);
			}
		}
	}

	seq[n] = '\0';
	*lenseq = n;
    
	laseq->seq = malloc(sizeof(char) * n + 1);
	strcpy(laseq->seq, trim(seq));

	free(seq);
	free(name);
	if (c == EOF)
		return (0);
	else
		return (1);
}
//substitute c by newc in string name
void remplace(char *name, char c, char newc)
{
	int i = 0;
//printf("remplacer %ss %c %c<BR>\n",name,c,newc);fflush(stdout);
	while (name[i] != '\0')
	{
		if (name[i] == c)
			name[i] = newc;
		i++;
	}
//printf("remplacé %s <BR>\n",name);fflush(stdout);
}


/*------------------------------------------------------*/
/*only numbers are not good!!*/
int check_valid_name(char *name)
{
	int tt = 0, i = 0;
//clean_str(name);
	tt = strlen(name);
	if (tt == 0)
		return 0;
    if (strchr(name,'>')!=NULL || strchr(name,'<')!=NULL)
            return 0;
	while (*(name + i) != '\0')
	{
		if (isdigit(*(name + i)))
			tt--;
		i++;
	}


	return tt;
}

/*------------------------------------------------------*/
/*as the form is multipart we have to find the string which is going to separate the
differents fields. delimiter are different among the browsers _for fun_
return number of fields*/
int search_delim(char *cgiinput, char *delim, char *ledir, FILE *fres)
{
	int i = 0, lequel = 0;
	char *begin, *end;
	char *achercher[2] = {"WebKitFormBoundary", "-----------------------------"}; // delimitors are parts of either first either 2nd string

	for (i = 0; i < 2; i++)
	{
		begin = strstr(cgiinput, achercher[i]); // find begining of delimitor
		if (begin != NULL)
			break;
	}
	lequel = i;
	if (i == 2)
	{
		printf("<BR>%s<BR>", cgiinput);
		html_error( 4 );
	}

	end = begin + strlen(achercher[lequel]); //look for delimitor 's end
	i = 0;

	while ((*(end + i) != ' ') && (*(end + i) != 13) && (*(end + i) != 10) && (*(end + i) != '\0'))
		i++;
	if (*(end + i) == '\0' || i >= 256) {fprintf(fres, "PB With your browser delimiter not found %d %.128s <BR>\n", i, end), fclose (fres), exit_properly(ledir);} //too much char in the delim

	// copy the clean delimitor
	strncpy(delim, begin, strlen(achercher[lequel]) + i);
	delim[strlen(achercher[lequel]) + i] = '\0';
	// count how many fields
	for (begin = cgiinput, end = begin, i = 0; end != NULL; i++, end = strstr(begin, delim), begin = end + 1)
		;

	i--;
	return (i - 1); //last delim is for ending the form
}


/*------------------------------------------------------*/
/** Read the CGI input and place all name/val pairs into list.        **/
/** Returns list containing name1, value1, name2, value2, ... , NULL  **/
char **getMultPartData(int *nb, char *nf, char *ledir, FILE *fres) {
	int i, j, k, ii ;

	int content_length;
	char *cgiinput ;
	char **cgivars ;
	int paircount ;
	char delim[256];
	char *nextdelim, *newName;
	char *begin, *end;
	char *bb;
	char *cont_type[3] = {"application/octet-stream", "text/html", "text/plain"};
//    FILE *ftemp;

	// perform some validity tests on cgi input
	if ( !(content_length = atoi(getenv("CONTENT_LENGTH"))) ) f_html_error( 1 ,ledir,fres);
	if ( !(cgiinput = (char *) malloc(content_length + 1)) )f_html_error( 2,ledir,fres);
	if (!fread(cgiinput, content_length, 1, stdin)) f_html_error( 3,ledir,fres );
//    if (strstr(getenv("HTTP_USER_AGENT"),"MSIE")!=NULL)
//        	{printf("<HR>EXPLORER found. It was not tested but you may be able to have results<HR><BR>");}

	cgiinput[content_length] = '\0' ;

//	DEBUG CGI WITH THIS
	/*	FILE *ftemp=fopen("/temp/cgivars","w");
		if (ftemp==NULL) exit(1);
		fprintf(ftemp,"%s\n",cgiinput);
		fclose(ftemp);
		ftemp=fopen("/temp/command","w");
		if (ftemp==NULL) exit(1);
		fprintf(ftemp,"export CONTENT_LENGTH=%d\n<BR> progname<tmp/cgivars<BR>",content_length);fflush(stdout);
		fclose(ftemp);
		exit(1);*/
//END of DEBUG

	paircount = search_delim(cgiinput, delim, ledir, fres); /*search the delimiter and how many of them*/
	cgivars = (char * *) malloc(2 * paircount * sizeof(char *)) ;
	end = strstr(cgiinput, delim);
	begin = cgiinput;
//	f=fopen("/temp/ttt","w");if (f!=NULL){fprintf(f,"%d\n%s",content_length,cgiinput);fclose(f);}
	bb = strstr(cgiinput, "filename="); //get the name of the file user has given
	i = 0;
	if (bb != NULL)
	{
		bb = bb + strlen("filename=") + 1;
		while (*(bb + i) != '"' && i < 256)
		{
			*(nf + i) = *(bb + i);
			i++;
		}
		*(nf + i) = '\0'; //file name is written into nf (initialised before call)
	}
	else
		sprintf(nf, "Pasted data"); //otherwise it's pasted data

	/*this part read and splits the string into an array cgivars organized as followed:
	  cgivars[i] contains the name of the field cgiVars[i+1] contains the value(s)*/
	for (i = 0;;)
	{
		end = strstr(begin, delim); //looking for next delim
		if (end == NULL || strlen(end) - 2 <= strlen(delim) )
			break;
		newName = strchr(end, '='); // looking for begining of name
		if (newName == NULL)
			fprintf(fres, "PB parsing form <BR>\n"), fclose (fres), exit_properly(ledir);



		newName += 2; // skip the = and the "
		k = 0;
		while (*(newName + k) != '"' && (newName + k) != NULL) k++; // look for end of name
		if ((newName + k) == NULL)		fprintf(fres, "PB parsing form <BR>\n"), fclose (fres), exit_properly(ledir);
		if (k == 0) printf("arrrg %d\n", i); //should never arrives exept if someone put a field with no name in the html form...
		cgivars[i] = malloc (sizeof(char) * (k + 2));
		strncpy(cgivars[i], newName, k);
		cgivars[i][k] = '\0';


		if (strstr(cgivars[i], "ubmit") != NULL)
		{ paircount = i;  break;} //last field is button submit and we don't care

		// do something special if file has been transmitted by the form because mime type is not the same for all browers _so fun_
		if (strstr(cgivars[i], "fileabgd") != NULL)
		{
			for (ii = 0; ii < 3; ii++)
			{
				newName = strstr(end, cont_type[ii]);
				if (newName != NULL)
					break;
			}
			if (newName == NULL) 		fprintf(fres, "PB parsing form <BR>\n"), fclose (fres), exit_properly(ledir);
			newName += strlen(cont_type[ii]); // so newname points on the real begining of the file...
			k = 0;
		}
		nextdelim = strstr(end + 1, delim); // look for end of data
		j = (nextdelim - (newName + k)); // how many chars are going to be written
		if (j < 0) j = 1;
		cgivars[i + 1] = malloc (sizeof(char) * (j + 1));

		while (*(newName + k) == '"' && (newName + k) != NULL) k++; //skip trailing "


		strncpy(cgivars[i + 1], newName + k, j);

		cgivars[i + 1][j] = '\0';


		begin = end + 1;
		i = i + 2;

	}
	/*clean the input i've parsed everything...*/
	free(cgiinput) ;

	*nb = i / 2;
	return cgivars ;

}

/*--------------------------------------------*/
void clean_str(char *ch)
{
	int i, j;
	char *chp;
	chp = malloc (sizeof (char) * strlen(ch) + 1);
	for (i = 0, j = 0; i < strlen(ch) - 1; i++)
		if (ch[i] == ' ' && ch[i + 1] == ' ')
			i++;
		else
			chp[j++] = ch[i];
	chp[j] = '\0';
	strcpy(ch, chp);
	free(chp);
}

/*--------------------------------------------*/
void print_groups_newick( Composante my_comp, DistMat mat  , char *lastring, FILE *f2, char *ledir,FILE *fres) {

	int i, j, k = 0, ng = 1;
	char nom[100], *bou;
	char chiffre[10];
	if (fres == NULL || fres==stdout) fres = stderr;
	for (i = 0; i < mat.n; i++) {
		for (j = 0; j < my_comp.n_in_comp[i]; j++)
		{
			k = my_comp.comp[i][j];


			if (strlen(mat.names[k]) > 100) exit(1);
			sprintf(nom, "%s", mat.names[k]);

			bou = strcasestr(lastring, nom);

			if (bou == NULL) //should never arrives
			{fprintf(fres, "ERROR: debug is %s cant be find in \n%s \n", nom, lastring); f_html_error(888,ledir,fres);}
			bou += strlen(nom) + 1;
			sprintf(chiffre, "%d", ng);
			*bou++ = '|'; *bou++ = 'g'; *bou++ = 'r'; *bou++ = 'o'; *bou++ = 'u'; *bou++ = 'p'; *bou++ = ' ';
			for (k = 0; k < strlen(chiffre); k++)*bou++ = chiffre[k];

		}
		if (my_comp.n_in_comp[i] != 0)
			ng++;

	}
	clean_str(lastring);
	fprintf(f2, "%s\n", lastring);


}

struct DistanceMatrix compute_dis(FILE *f, int method, float ts_tv, int *len_seq,char *ledir,FILE *fres)
{
	struct FastaSeq *mesSeq;

	int i = 1;;
	int nalloc = 256;
	int nseq = 0;
	struct DistanceMatrix my_mat;   /* store distance matrix, names and matrix size */


	mesSeq = (struct FastaSeq *)malloc (sizeof (struct FastaSeq ) * nalloc);

	while (i)
	{
//	fprintf(stderr,"%d seq\n",i);
		i = ReadFastaSequence(f, &mesSeq[nseq], len_seq);
		nseq++;
		if (nseq == nalloc)
		{
			nalloc += 256;
			mesSeq = realloc(mesSeq, sizeof (struct FastaSeq ) * nalloc);
			if (mesSeq == NULL) {fprintf(fres,"not enough memory\n"); f_html_error(888,ledir,fres);}
		}
	}

	fprintf(stderr,"done read\n");

	if (check_names(mesSeq, nseq,ledir,fres) == 0)
		{fprintf(fres,"Two seqs found with same name. Exit\n"); f_html_error(888,ledir,fres);}

//	fprintf(stderr, "get dist mat\n");
	my_mat = GetDistMat(nseq, mesSeq, method, ts_tv,ledir,fres);
	fprintf(stderr, "done mat\n");


	for (i = 0; i < nseq; i++)
	{free(mesSeq[i].seq); free(mesSeq[i].name);}
	free(mesSeq);
	return my_mat;
}



/*compute a very tricky distance for 2 sequences*/
void distancesimple(struct FastaSeq *mesSeqs, int l, struct  DistanceMatrix  my_mat,char *ledir,FILE *fres)

{
	char *s1, *s2, c1, c2;;
	double v = 0;
	int i, a, b, ncor = 0;
	int nseq = my_mat.n;

	if (l == 0)
		f_html_error(100,ledir,fres);

	for (a = 0; a < nseq - 1; a++)
	{
		s1 = mesSeqs[a].seq;
		my_mat.dist[a][a] = 0;
		for (b = a + 1; b < nseq; b++)
		{
			v = 0; ncor = 0;
			s2 = mesSeqs[b].seq;
			if (check_compat(s1, s2, l) == 0)
			{fprintf(fres,"<H2>Sequence %s and %s have no common site. Distance can't be computed. Bye <BR>", my_mat.names[a], my_mat.names[b]); fflush(stdout); exit(1);}

			for (i = 0; i < l; i++)
			{
				c1 = toupper(*(s1 + i));
				c2 = toupper(*(s2 + i));
				if (compare_DNA(c1, c2) == 0)
					v = v + 1;
				if ( ( (*(s1 + i) ) == '-') || ((*(s2 + i) ) == '-') ||	((*(s1 + i) ) == 'N') || ((*(s2 + i) ) == 'N' ))
					ncor++;



			}
			v = ((v) / (double)(l - ncor));
			if (isnan(v))v=1.0;
			my_mat.dist[a][b] = my_mat.dist[b][a] = v;
		}
	}
}

/*compute distance according to Jukes Cantor method*/
/*do not take in consideration gaps or N*/


void distanceJC69 (struct FastaSeq *mesSeqs, int l, struct  DistanceMatrix  mymat,char *ledir,FILE *fres)
//double distanceJC69 (char *s1,char *s2, int l)
{
	double v = 0, h;
	int i, newl = 0;
	char c1, c2;
	int a, b;
	char *s1, *s2;
	int nseq = mymat.n;


	if (l == 0)
		f_html_error(100,ledir,fres);

	for (a = 0; a < nseq - 1; a++)
	{
		s1 = mesSeqs[a].seq;
		mymat.dist[a][a] = 0;
		for (b = a + 1; b < nseq; b++)
		{
			s2 = mesSeqs[b].seq;
			newl = 0; v = 0;
			if (check_compat(s1, s2, l) == 0)
			{fprintf(fres,"<H2>Sequence %s and %s have no common site. Distance can't be computed. Bye <BR>", mymat.names[a], mymat.names[b]); fflush(stdout); exit(1);}

			for (i = 0; i < l; i++)
			{
				c1 = toupper(*(s1 + i));
				c2 = toupper(*(s2 + i));
				if (c1 != '-' && c1 != 'N' && c2 != '-' && c2 != 'N')
					newl++;
				if (compare_DNA(c1, c2) == 0)
					v = v + 1;
			}


			v = (v) / (double)(newl);
			if (v>=0.75)v=0.74;//SOFIZ
				h = (-3.0 / 4.0) * log(1.0 - ((4.0 / 3.0) * v));
		//	if (isnan(h)) h=1.0; 
			if (h == -0)
				h = 0;
			mymat.dist[a][b] = mymat.dist[b][a] = h;

		}
	}
}
char IsTransition( char nt1, char nt2 ) {

	if (nt1 == 'T') {
		if ( nt2 == 'C' || nt2 == 'S' ||  nt2 == 'M' || nt2 == 'V')
			return 1;
		else
			return 0;
	}
	if (nt1 == 'C') {
		if ( nt2 == 'T' || nt2 == 'W' ||  nt2 == 'K' || nt2 == 'D')
			return 1;
		else
			return 0;
	}
	if (nt1 == 'A') {
		if ( nt2 == 'G' || nt2 == 'S' ||  nt2 == 'K' || nt2 == 'B')
			return 1;
		else
			return 0;
	}
	if (nt1 == 'G') {
		if ( nt2 == 'A' || nt2 == 'W' ||  nt2 == 'M' || nt2 == 'H')
			return 1;
		else
			return 0;
	}
	return (0); //should never go there but compiler complains
}

/*
	Test if it is a transition
	return 1 if yes, 0 if not.
*/
char IsTransversion( char nt1, char nt2 ) {

	if ( compare_DNA(nt1, nt2) == 0 && IsTransition(nt1, nt2 ) == 0 )
		return 1;
	else
		return 0;

}
void transition_transversion_sequences(char *seq1, char *seq2, long L, long *tsi, long *tsv) {

	long i;

	*tsi = *tsv = 0;

	for ( i = 0 ; i < L ; i++ )
		if ( compare_DNA( *(seq1 + i), *(seq2 + i) ) == 0 ) {
			if ( IsTransition( *(seq1 + i), *(seq2 + i) ) == 1 )
				(*tsi)++;
			else
				(*tsv)++;
		}

	return;
}

double P_given_t_R( double t, double R ) {
	return 0.25 - 0.5 * exp( - t * (2 * R + 1) / (R + 1) ) + 0.25 * exp(- 2 * t / (R + 1));
}
double Q_given_t_R( double t, double R ) {
	return 0.5 - 0.5 * exp(- 2 * t / (R + 1));
}
double compute_logL_given_t_R( long nsites, long n_tsv, long n_tsi, double t, double R ) {

	return nsites * log(0.25) + \
	       (nsites - n_tsv - n_tsi) * log( 1.00 - P_given_t_R(t, R)  - Q_given_t_R(t, R) ) + \
	       n_tsi * log( P_given_t_R(t, R) ) + \
	       n_tsv * log( Q_given_t_R(t, R) );
}

double compute_k80( long nsites, long n_tsv, long n_tsi ) {

	double Q = n_tsv / (double)nsites,
	       P = n_tsi / (double)nsites;

	return -0.25 * log( (1 - 2 * Q) * (1 - 2 * P - Q) * (1 - 2 * P - Q) );

}

double find_ML_t_given_R( double R, long nsites, long n_tsv, long n_tsi ) {

	double t = compute_k80( nsites, n_tsv, n_tsi );     /* seed it with an empirical value from the data -- Kimura 80 */


	double epsilon = 1e-7;
	double eps = 1e-3;

	/*	printf("dist is %.10f\n", t );*/

	while ( eps >= epsilon ) {

		while ( compute_logL_given_t_R( nsites, n_tsv, n_tsi, t + eps, R ) > compute_logL_given_t_R( nsites, n_tsv, n_tsi, t, R )  )
			t += eps;
		while ( compute_logL_given_t_R( nsites, n_tsv, n_tsi, t - eps, R ) > compute_logL_given_t_R( nsites, n_tsv, n_tsi, t, R )  )
			t -= eps;

		eps *= 0.1;

	}

	return t;
}

void distanceK80 (struct FastaSeq *mesSeqs, int l, struct  DistanceMatrix  my_mat,char *ledir,FILE *fres) {
	int i, j;
	long tsi, tsv;
	long del;
	int nseq = my_mat.n;
	double h;

	for (i = 0; i < nseq; i++) {

		my_mat.dist[i][i] = 0;

		for (j = i + 1; j < nseq; j++) {
			if (check_compat(mesSeqs[i].seq, mesSeqs[j].seq, l) == 0)
			{fprintf(fres,"<H2>Sequence %s and %s have no common site. Distance can't be computed. Bye <BR>", my_mat.names[i], my_mat.names[j]); f_html_error(888,ledir,fres);}
			transition_transversion_sequences(mesSeqs[i].seq, mesSeqs[j].seq, l, &tsi, &tsv);
			del = del_sequences(mesSeqs[i].seq, mesSeqs[j].seq, l);
			h=find_ML_t_given_R( my_mat.ratio_ts_tv, l - del, tsv, tsi );
			//if (isnan(h)) h=1.0; //SOFIZ PABON
			my_mat.dist[i][j] = my_mat.dist[j][i] = h;
			if (my_mat.dist[i][j] == -0) //happens sometimes
				my_mat.dist[i][j] = my_mat.dist[j][i] = 0;


		}
	}


}



/*check if we have at least one common symbol beetween the 2 seqs*/
int check_compat(char *s1, char *s2, int l)
{
	int i, newl = 0;
	char c1, c2;
	for (i = 0; i < l; i++)
	{
		c1 = toupper(*(s1 + i));
		c2 = toupper(*(s2 + i));
		if ((c1 == '-' && c2 != '-') || (c2 == '-' && c1 != '-') || (c2 == '-' && c1 == '-'))
			newl++;
	}



	return (l - newl);
}




/*compute distance according to Tmura Nei method*/
/*do not take in consideration gaps or N*/
void distanceTN93(struct FastaSeq *mesSeqs, int l, struct  DistanceMatrix  my_mat,char *ledir,FILE *fres)
{
	double v = 0;
	double transitions = 0, transversions = 0,  q, p1, p2, ga, gg, gc, gt, gr, gy;
	//double p;
	double transitionsag = 0, transitionsct = 0;
	char c1, c2;
	double f[5];
	char nuc[5] = "ACGT-";
	int i, newl = 0;
	char *s1, *s2;
	int a, b;
	int nseq = my_mat.n;

	if (l == 0)
		f_html_error(100,ledir,fres);

	for (a = 0; a < nseq - 1; a++)
	{
		s1 = mesSeqs[a].seq;
		my_mat.dist[a][a] = 0;
		for (b = a + 1; b < nseq; b++)
		{
			s2 = mesSeqs[b].seq;
			if (check_compat(s1, s2, l) == 0)
			{fprintf(fres,"<H2>Sequence %s and %s have no common site. Distance can't be computed. Bye <BR>", my_mat.names[a], my_mat.names[b]); f_html_error(888,ledir,fres);}
			newl = 0; v = 0;
			for (i = 0; i < 5; i++)
				f[i] = 0;
			for (i = 0; i < l; i++)
			{
				c1 = toupper(*(s1 + i));
				c2 = toupper(*(s2 + i));
				if (strchr(nuc, c1) && strchr(nuc, c2) )
				{
					f[(int)(strchr(nuc, c1) - nuc)]++;
					f[(int)(strchr(nuc, c2) - nuc)]++;
				}

				if (compare_DNA(c1, c2) == 0)
				{
					v++;
					if ((c1 == 'A' && c2 == 'G') || (c2 == 'A' && c1 == 'G') )transitionsag++;
					else if ((c1 == 'C' && c2 == 'T') || (c2 == 'C' && c1 == 'T')) transitionsct++;

				}
			}


			transversions = v - transitions;
			q = transversions / (double)newl;
			//p = transitions / (double)newl; not used?????????
			p1 = transitionsag / (double)newl;
			p2 = transitionsct / (double)newl;
			ga = f[0] / (2.0 * newl);
			gc = f[1] / (2.0 * newl);
			gg = f[2] / (2.0 * newl);
			gt = f[3] / (2.0 * newl);
			gr = ga + gg;
			gy = gc + gt;
			v = ((-2 * ga * gg / gr) * log (1.0 - ((gr / (2.0 * ga * gg)) * p1) - ((1 / (2 * gr)) * q))) -
			    (( (2 * gt * gc) / gy) * log(1.0 - ((gy / (2.0 * gt * gc)) * p2) - ((1 / (2 * gy)) * q))) -
			    ((2.0 * ((gr * gy) - ( (ga * gg * gy) / gr) - ((gt * gc * gr) / gy))) * log (1.0 - ((1.0 / (2.0 * gr * gy)) * q))) ;
			if (v == -0)
				v = 0;
			my_mat.dist[a][b] = my_mat.dist[b][a] = v;

		}
	}
}
struct DistanceMatrix GetDistMat(int nseq, struct FastaSeq *mesSeqs, int method, float ts_tv, char *ledir, FILE *fres)
{

	struct DistanceMatrix my_mat;                  /* store distance matrix, names and matrix size */
	void (*distance) (struct FastaSeq *, int , struct DistanceMatrix ,char *ledir, FILE *fres) = NULL;      /* pointeur de fonction */;
	int a;
	int length;

	length = strlen(mesSeqs[0].seq);

	switch (method) {

	case 0:
		distance = distanceK80;
		break;

	case 1:
		distance = distanceJC69;
		break;

	case 2:
		distance = distanceTN93;
		break;

	case 3:
		distance = distancesimple;
		break;
	}

	my_mat.names = (char **)malloc( (size_t) sizeof(char *)*nseq + 1);
	my_mat.n = nseq;
	my_mat.ratio_ts_tv = ts_tv;


	if ( ! my_mat.names )fprintf(fres,"read_distmat: cannot allocate my_mat.names, bye<BR>"), fclose (fres), exit_properly(ledir);

	for (a = 0; a < my_mat.n; a++) {
//		my_mat.names[a] = (char *)malloc( (size_t) sizeof(char)*SIZE_NAME_DIST +1);
		my_mat.names[a] = (char *)malloc( (size_t) sizeof(char) * (strlen(mesSeqs[a].name) + 1));
		if ( ! my_mat.names[a] )
			fprintf(fres,"read_distmat: cannot allocate my_mat.names[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);

		strcpy(my_mat.names[a], mesSeqs[a].name);

//	my_mat.names[a][SIZE_NAME_DIST]='\0';
	}

	my_mat.dist = (double **)malloc( (size_t) sizeof(double *)*my_mat.n );
	if ( ! my_mat.dist)fprintf(fres,"read_distmat: cannot allocate my_mat.dist, bye<BR>"), fclose (fres), exit_properly(ledir);
	for (a = 0; a < my_mat.n; a++) {
		my_mat.dist[a] = (double *)malloc( (size_t) sizeof(double) * my_mat.n );
		if ( ! my_mat.dist[a] )
			fprintf(fres,"read_distmat: cannot allocate my_mat.dist[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);
	}


	distance(mesSeqs, length, my_mat, ledir, fres);

	return my_mat;

}

/*--------------------------------------------------*/
void readMatrixMega_string(char *data, struct DistanceMatrix *my_mat, char *ledir, FILE *fres) {

	int a;
	int nbc = 0;
	int nbl = 0;
	int l = 0;
	double v;
	char *ptr, *nptr;;
	int lower = -1;
	char name [256];


	my_mat->n = 0;
	my_mat->names = NULL;
	my_mat->dist = NULL;
	printf("Format Mega detected<BR>\n");	fflush (stdout);
	ptr = (char *)strcasestr((const char *)data, "#mega") + 5; /*mega format 1st line begins by keyword mega*/
	if (ptr == NULL)fprintf(fres, "ReadMatrixMega_string: your matrice is not recognized "), fclose (fres), exit_properly(ledir);

	if (strcasestr ((const char *) ptr, "dataformat") != NULL)
	{
		if (strcasestr( ptr, "lowerleft") != NULL)
			lower = 1;
		else if (strcasestr( ptr, "upperight") != NULL)
			lower = 0;
		else
			fprintf(fres, "ReadMatrixMega_string: your matrice is not recognized "), fclose (fres), exit_properly(ledir);
	}
	else
		fprintf(fres, "ReadMatrixMega_string: your matrice is not recognized"), fclose (fres), exit_properly(ledir);


	ptr = strcasestr((const char *)data, "of Taxa :");
	if (ptr != NULL)
		my_mat->n = atol(strchr(ptr, ':') + 1);
	else
	{
		ptr = strcasestr(data, "ntaxa=");
		if (ptr != NULL)
			my_mat->n = atol(strchr(strcasestr(data, "ntaxa="), '=') + 1);
		else
			fprintf(fres, "ReadMatrixMega_string: Nbr of taxa is not found"), fclose (fres), exit_properly(ledir);
	}

	my_mat->names = (char **)malloc( (size_t) sizeof(char *)* my_mat->n + 1 );
	if ( ! my_mat->names )fprintf(fres, "ReadMatrixMega_string: cannot allocate my_mat.names, bye<BR>"), fclose (fres), exit_properly(ledir);
	/*	for(a=0;a<my_mat->n; a++){
			my_mat->names[a] = (char *)malloc( (size_t) sizeof(char)*SIZE_NAME_DIST+1 );
			if( ! my_mat->names[a] )
				printf( "read_distmat: cannot allocate my_mat.names[%d], bye<BR>",a), exit(4);
		}*/

	my_mat->dist = (double **)malloc( (size_t) sizeof(double *)*my_mat->n + 1 );
	if ( ! my_mat->dist )fprintf(fres, "ReadMatrixMega_string: cannot allocate my_mat.dist, bye<BR>"), fclose (fres), exit_properly(ledir);

	for (a = 0; a < my_mat->n; a++) {
		my_mat->dist[a] = (double *)malloc( (size_t) sizeof(double) * my_mat->n + 1 );
		if ( ! my_mat->dist[a] )
			fprintf(fres, "ReadMatrixMega_string: cannot allocate my_mat.dist[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);
	}
//	printf("%ld taxa <BR>\n",my_mat->n);fflush(stdout);
	nptr = strchr(ptr, ';') + 1; //next ; will be names
	ptr = nptr;

	for (a = 0; a < my_mat->n; a++)
	{
		ptr = strchr(nptr, '#');
		if (ptr == NULL)   fprintf(fres, "ReadMatrixMega_string:symbol '#' was not found bye<BR>"), fclose (fres), exit_properly(ledir);

		sscanf(ptr + 1, "%s\n", name);
		l = strlen(name) + 1;
		my_mat->names[a] = (char *)malloc( (size_t) sizeof(char) * l );
		if (!my_mat->names[a])
			fprintf(fres, "ReadMatrixMega_string: cannot allocate my_mat.names[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);
		strcpy(my_mat->names[a], name);
//		my_mat->names[a][strlen(name)]='.';
//		my_mat->names[a][strlen(name)+1]='\0';

		//here change the name because some names are fully included in others and this causes pb withbionj

		if (strchr(my_mat->names[a], '(') != NULL)
			remplace(my_mat->names[a], '(', '_');

		if (strchr(my_mat->names[a], ')') != NULL)
			remplace(my_mat->names[a], ')', '_');
        
        if (strchr(my_mat->names[a], '<') != NULL)
			remplace(my_mat->names[a], '<', '_');
        if (strchr(my_mat->names[a], '>') != NULL)
			remplace(my_mat->names[a], '>', '_');
        
		nptr = strstr(ptr, name) + 1;
		ptr = nptr;

	}
//
	ptr = strchr(nptr, ']') + 1;
	nptr = ptr;
//at this point we point on 1st line of matrix; we can read it

	if (lower)
		for (nbl = 0; nbl < my_mat->n; nbl++)
		{
			nptr = strchr(ptr, ']') + 1;
			ptr = nptr;
			for (nbc = 0; nbc < nbl; nbc++)
			{
				while (*ptr == ' ' ) ptr++;
				if (*ptr == '?') 	//strangely mega puts ? when distance is null
				{

					fprintf(fres, "ReadMatrixMega_string:Warning distance between %s and %s is unknown,exiting<BR>\n", my_mat->names[nbl], my_mat->names[nbc]), fclose (fres), exit_properly(ledir);
				}
				sscanf(ptr, "%lf", &v);
//				printf("-->%d %d%f (%.10s)<BR>",nbl,nbc,v,ptr);fflush(stdout);
				my_mat->dist[nbl][nbc] = my_mat->dist[nbc][nbl] = v; //fill the matrix symetrically
				while (*ptr != ' ' && *ptr != '\n' && *ptr != '\r') ptr++;

			}
			fflush(stdout);
			my_mat->dist[nbl][nbl] = 0;
		}
	else //uper matrix
	{
		for (nbl = 0; nbl < my_mat->n - 1; nbl++)
		{
			nptr = strchr(ptr, ']') + 1;
			ptr = nptr;
			for (nbc = nbl + 1; nbc < my_mat->n - nbl; nbc++)
			{
				while (*ptr == ' ' ) ptr++;
				if (*ptr == '?')
				{
					fprintf(fres, "ReadMatrixMega_string:Warning distance between %s and %s is unknown,exiting<BR>\n", my_mat->names[nbl], my_mat->names[nbc]), fclose (fres), exit_properly(ledir);
				}

				else
					sscanf(ptr, "%lf", &v);
				my_mat->dist[nbl][nbc] = my_mat->dist[nbc][nbl] == v; //fill the matrix symetrically
				while (*ptr != ' ' && *ptr != '\n' && *ptr != '\r') ptr++;
			}


			my_mat->dist[nbl][nbl] = 0;
		}
		my_mat->dist[my_mat->n - 1][my_mat->n - 1] = 0;
	}
//printf("read mega done<BR>");fflush (stdout);
//printf("%d<BR>\n",my_mat->n);
//printf("il est temps que ca se termine\n<BR>");	fflush (stdout);

}



/*--------------------------------------------------*/
void  readMatrixMegaCVS_string(char *data, struct DistanceMatrix *my_mat, char *ledir, FILE *fres)
{
	int nb = 0, a = 0, b = 0, n = 0;
	char *pt = data;
	char *p, *ll;

	float d;
	//fprintf(fres," MEGA 5 format CSV detected. <BR>\n");
	fflush(stdout);
	while (strncmp(pt, "Table", 5) != 0 && *pt != '\0') //hoping that keyword Table is the end of the matrix
	{

		if (*pt == '\n' && ( (*(pt + 1) == '\r') || (*(pt + 1) == 10) || (*(pt + 1) == 13) ) ) if (nb != 0) break;
		//if (*pt == '\n' && isalnum(*(pt + 1))) {nb++;  }
		if (*pt == '[') {nb++;  }
//	if (*pt=='\n') printf("***RC %d %d<BR>",*pt,*(pt+1));
		pt++;
	}
	if (*pt == '\0')
		f_html_error(343,ledir,fres);

	if (nb <=2)
		f_html_error(344,ledir,fres);

	my_mat->n = nb;
	

	my_mat->names = (char **)malloc( (size_t) sizeof(char *) * (my_mat->n + 1) );
	if ( ! my_mat->names )

		fprintf(fres, "readMatrixMegaCVS_string: cannot allocate my_mat.names, bye<BR>"), fclose (fres), exit_properly(ledir);

	my_mat->dist = (double **)malloc( (size_t) sizeof(double *) * (my_mat->n + 1) );
	if ( ! my_mat->dist )
		fprintf(fres, "readMatrixMegaCVS_string:read_distmat: cannot allocate my_mat.dist, bye<BR>"), fclose (fres), exit_properly(ledir);
	for (a = 0; a < my_mat->n; a++) {
		my_mat->dist[a] = (double *)malloc( (size_t) sizeof(double) * (my_mat->n + 1) );
		if ( ! my_mat->dist[a] )
			fprintf(fres, "readMatrixMegaCVS_string: cannot allocate my_mat.dist[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);
	}
//Now read everything... rewind the string before...
//printf("%d read<BR>\n",nb);fflush(stdout);
	pt = data;
	for (a = 0; a < my_mat->n; a++)
	{
		if (*pt == ',' || *pt == 10 || *pt == 13 || *pt == ' ') pt++; //skip any no need chars
		ll = strchr(pt, ',');
		if (ll != NULL)
		{
			n = strchr(pt, ',') - pt;
//			n=MINI(strchr(pt,',')-pt,SIZE_NAME_DIST);
			my_mat->names[a] = (char *)malloc( (size_t) sizeof(char) * (n + 1) );
			strncpy(my_mat->names[a], pt, n);
			my_mat->names[a][n] = '\0';
//			strncpy(my_mat->names[a],pt, n);
			n = ll - pt;
		}
		else
			f_html_error(344,ledir,fres);

		pt += n + 1; /*+1 should advance on 1st val*/
		if (strchr(my_mat->names[a], '(') != NULL)
			remplace(my_mat->names[a], '(', '_');

		if (strchr(my_mat->names[a], ')') != NULL)
			remplace(my_mat->names[a], ')', '_');
        if (strchr(my_mat->names[a], '<') != NULL)
			remplace(my_mat->names[a], '<', '_');
		if (strchr(my_mat->names[a], '>') != NULL)
			remplace(my_mat->names[a], '>', '_');
		if (strchr(my_mat->names[a], '[') != NULL)
			remplace(my_mat->names[a], '[', '_');
		if (strchr(my_mat->names[a], ']') != NULL)
			remplace(my_mat->names[a], ']', '_');
//		printf("%s<BR>\n",my_mat->names[a]);fflush(stdout);
		for (b = 0; b <= a; b++)
		{
			if (*pt == ',')
			{my_mat->dist[a][b] = 0; pt++;}
			else
			{
				if (*pt == '?')
				{
					fprintf(fres, "readMatrixMegaCVS_string:**Warning distance between %s and %s is unknown,exiting<BR>\n", my_mat->names[a], my_mat->names[b]), fclose (fres), exit_properly(ledir);
				}
				else
					sscanf(pt, "%f", &d);
				my_mat->dist[a][b] = my_mat->dist[b][a] = d;
				p = strchr(pt, ',');
				if (p == NULL) f_html_error(344,ledir,fres);
				pt = p + 1;
			}
			//fprintf(fres,"%d %d %f<BR>\n",a,b,my_mat->dist[a][b]);fflush(stdout);

		}
	}
	

}

/*------------------------------------------------------*/
/*
	Takes a distance file as an input (phylip format)
	Return a struct with a distance matrix
*/

struct DistanceMatrix read_distmat_string( char *data , int fmega, char *ledir, FILE *fres) {

	int a, n;
	int nbc = 0;
	int nbl = 0;
	struct DistanceMatrix my_mat;   /* store distance matrix, names and matrix size */
	char *ptr;
	char *end;

	my_mat.n = 0;
	my_mat.names = NULL;
	my_mat.dist = NULL;
	//fprintf(fres,"read_distmat_string<BR> %d %s",fmega,data),fclose (fres), exit_properly(ledir);
	fflush(stdout);
	if (fmega == 5)
		readMatrixMegaCVS_string(data, &my_mat, ledir, fres);
	else if (strcasestr((const char *)data, "#mega") != NULL )
		readMatrixMega_string(data, &my_mat, ledir, fres);
	else /* classic phylip distance matrix*/
	{

		while (!isdigit(*data) && *data != '\0') data++;	//skip any  trailing spaces
		if (*data == '\0') f_html_error(10,ledir,fres);

		if (*data == '>') //fasta detetected but tested before
		{
			fprintf(fres, "read_distmat_string: FASTA file pb "), fclose (fres), exit_properly(ledir);
		}
		// read the number of species in the distance matrix

		ptr = strtok(data, " \n\r");
		my_mat.n = atol(ptr);
		
		/*
			Get memory
		*/
		my_mat.names = (char **)malloc( (size_t) sizeof(char *)* my_mat.n + 1 );
		if ( ! my_mat.names )fprintf(fres, "read_distmat_string:  cannot allocate my_mat.names, bye<BR>"), fclose (fres), exit_properly(ledir);
		/*		for(a=0;a<my_mat.n; a++){
					my_mat.names[a] = (char *)malloc( (size_t) sizeof(char)*SIZE_NAME_DIST+1 );
					if( ! my_mat.names[a] )
						printf( "read_distmat: cannot allocate my_mat.names[%d], bye<BR>",a), exit(4);
				}
		*/
		my_mat.dist = (double **)malloc( (size_t) sizeof(double *)*my_mat.n + 1 );
		if ( ! my_mat.dist )fprintf(fres, "read_distmat_string: cannot allocate my_mat.dist, bye<BR>"), fclose (fres), exit_properly(ledir);
		for (a = 0; a < my_mat.n; a++) {
			my_mat.dist[a] = (double *)malloc( (size_t) sizeof(double) * my_mat.n + 1 );
			if ( ! my_mat.dist[a] )
				fprintf(fres, "read_distmat_string: cannot allocate my_mat.dist[%d], bye<BR>", a), fclose (fres), exit_properly(ledir);
		}

		//read the first taxon name
		ptr = strtok(NULL,  " \n\r\t");
//		n=MINI(strlen(ptr),SIZE_NAME_DIST);
		n = strlen(ptr);
		my_mat.names[0] = (char *)malloc( (size_t) sizeof(char) * (n + 1) );
		strncpy(my_mat.names[0], ptr, n);  //st si initialised out of the for
		my_mat.names[0][n] = '\0';
		// read all values and other taxon names
		//fprintf(fres," (%d) %d %d [%s][%d] %d %d %d %d %d %d %d %d %d %d <BR>\n",n,nbl,nbc,ptr,ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],
		//ptr[6],ptr[7],ptr[8],ptr[9],ptr[10]);
		while (1)
		{
			/*extract the next token*/
			ptr = strtok(NULL,  " \n\r\t\009");
			//fprintf(fres,"%d %d %s<BR>\n",nbl,nbc,ptr);
			if (ptr == NULL || nbl> my_mat.n) break;
			
			if (nbc == my_mat.n) // we have found a new name so store it
			{
				nbl++;
				if (nbl == my_mat.n) break;
				n = strlen(ptr);
				//			n=strlen(ptr),SIZE_NAME_DIST);
				my_mat.names[nbl] = (char *)malloc( (size_t) sizeof(char) * (n + 1) );
				if (my_mat.names[nbl]==NULL) fprintf(fres, "Not enough memory for %d names<BR>",nbl), fclose (fres), exit_properly(ledir);
				strncpy(my_mat.names[nbl], ptr, n);  //st si initialised out of the for
			
				
				my_mat.names[nbl][n] = '\0';
				nbc = 0;
			}
			else
			{
				double tempo;
				tempo=strtod(ptr, &end);
				my_mat.dist[nbl][nbc] = tempo; //store the value found and
				if (strlen(end) >0) fprintf(fres, " Pb in your matrix file only symetrical matrix are supported %d %d (%f) %s %s<BR>",nbl,nbc,tempo,ptr,end), fclose (fres), exit_properly(ledir);
				nbc++;
					
			}
		}
        if (ptr==NULL && nbl != my_mat.n -1)
		if (nbc != my_mat.n || nbl != my_mat.n)
			fprintf(fres, "read_distmat_string:  Pb in your matrix file number of species (%ld) dont match line (%d) or col(%d)bye<BR>",my_mat.n,nbl,nbc), fclose (fres), exit_properly(ledir);


	}

//	printf("------------------<BR>\n");fflush(stdout);


	return my_mat;

}

/*------------------------------------------------------*/
/*  take a fasta file as input and compute distance as method
	Returns a struct with a distance matrix */
struct DistanceMatrix  read_fasta_and_compute_dis(char *input, int method, float ts_tv, char *ledir, FILE *fres, int *ll)
{
	struct FastaSeq *mesSeq; /*store all sequences in a struct with name and seq*/
	char *myseq;
	char *nucs = "ATGC-+NMRWSYKVHDB";


	int i;
	int n = 0, nseq = -1, k = 0, redo = 1, length = 0, nbr = 1, nk;
	struct DistanceMatrix my_mat;   /* store distance matrix, names and matrix size */
	int *agarder;
	/*some usefull mallocs*/

	mesSeq = malloc(sizeof(struct FastaSeq) * 128);
	myseq = malloc(sizeof(char) * 100);
	agarder = malloc(sizeof(int) * 100);

	while ((*input != '>') && *input != '\0') input++;	//skip any  trailing spaces

	if (*input == '\0')
		fprintf(fres, "READFASTA:File sequence is empty\n"), fclose (fres), exit_properly(ledir);
	//fprintf(fres,">>> %.10s\n<BR>",input);
	// going to parse sequences from input to mesSeqs
	while ( *input != 0 ) {
		if (*input == '>' ) {          /* if any, the previous seq ends */
			if (nseq == 0)length = n;
			if (nseq != -1) // not the fisrt seq so sizeseq is ok
			{
				if (n != length) fprintf(fres, "READFASTA:All sequences must be same size (check sequence %d)\n", nseq + 1), fclose (fres), exit_properly(ledir);; //at least one seq with different length is found
				myseq[n] = '\0';
				mesSeq[nseq].seq = malloc(sizeof(char) * length + 1);
				if (!mesSeq[nseq].seq)	fprintf(fres, "READFASTA:MEMORY ERROR error can allocate mesSeq\n"), fclose (fres), exit_properly(ledir);
				;
				strcpy(mesSeq[nseq].seq, myseq);
				n = 0;
			}

			nseq++;
			if (nseq >= 128 * redo) //need to add some more memory
			{
				mesSeq = realloc(mesSeq, sizeof(struct FastaSeq) * (++redo) * 128);
				if (!mesSeq) fprintf(fres, "READFASTA:MEMORY ERROR error can reallocate mesSeq\n"), fclose (fres), exit_properly(ledir);;
			}
			k = 0;
			input++; //skip the trailing ">"
			//	fprintf(fres,">> %.10s\n<BR>",input);
			while ( *(input + k) != 10 && *(input + k) != 13 && *(input + k) != 0) 	k++;
			if (*(input + k) == 0) fprintf(fres, "READFASTA:Pb with File, Please check your data\n"), fclose (fres), exit_properly(ledir);;
			/*			    if (k>SIZE_NAME_DIST)
						    	nk=SIZE_NAME_DIST;
						    else*/
			nk = k;

			mesSeq[nseq].name = malloc(sizeof(char) * nk + 2);
			if (!mesSeq[nseq].name) fprintf(fres, "READFASTA:MEMORY ERROR error can allocate mesSeq names\n"), fclose (fres), exit_properly(ledir);
            
            char *seqnametocheck = malloc(sizeof(char) * nk + 2);
            strncpy(seqnametocheck, input, nk);
            seqnametocheck[nk] = '\0';
            char *seqname_clean=trim(seqnametocheck);
            mesSeq[nseq].name= malloc(sizeof(char) * nk + 2);
			strcpy(mesSeq[nseq].name, seqname_clean);
            free(seqnametocheck);
			//mesSeq[nseq].name[nk] = '\0';
			if (check_valid_name(mesSeq[nseq].name) == 0)
			{
                if (strchr(mesSeq[nseq].name,'>')!=NULL || strchr(mesSeq[nseq].name,'<')!=NULL)
                    fprintf(fres, "READFASTA:Pb reading  Name Seq %d contains > or < check your entries <BR>", nseq + 1); 
                else
				fprintf(fres, "READFASTA:Pb reading  Name Seq %d : %s <BR>", nseq + 1, mesSeq[nseq].name);
                fclose (fres), exit_properly(ledir);

			//	fhtml_error(110);
			}
			while ( *(input + k) == 10 || *(input + k) == 13 ) k++;
			//fprintf(fres,">> %.10s\n<BR>",input);
			input = input + k-1; //because +1
			//fprintf(fres,"+k %.10s\n<BR>",input);
			
			k = 0;

		}//end of >

		else if (*input != 10 && *input != 13 && *input != 32 && *input != 9 ) //store the sequence
		{
			//fprintf(fres,"in else %.10s\n<BR>",input);
			if (strchr(nucs, toupper(*input)) == NULL)  fprintf(fres, "READFASTA:Some unexpected symbol was found:(%c)\n", *input), fclose (fres), exit_properly(ledir); /*weird symbol found*/
			myseq[n++] = toupper(*input);
			if (nseq == 0) /*only for the 1st seq*/
				if (n % 100 == 0) /*give another 100 extras char if not enough room for seqs*/
				{
					myseq = realloc(myseq, sizeof(char) * (++nbr) * 100);
					agarder = realloc(agarder, sizeof(int) * nbr * 100);
					if (!myseq) fprintf(fres, "READFASTA:MEMORY ERROR cant reallocate myseq\n"), fclose (fres), exit_properly(ledir);
				}
		}
		//fprintf(fres,"----> %.10s\n<BR>",input);
			input++;
		//fprintf(fres,"----> %.10s\n<BR>",input);
	}

	mesSeq[nseq].seq = malloc(sizeof(char) * length + 1);
	if (!mesSeq[nseq].seq)fprintf(fres, "READFASTA:MEMORY ERROR error can allocate mesSeq[%d].seq\n", nseq), fclose (fres), exit_properly(ledir);
	*ll = length;
	//copy last seq
	strncpy(mesSeq[nseq].seq, myseq, length);
	mesSeq[nseq].seq[length] = '\0';
	nseq++;

	//This option is de activated in online command version
	if (nseq > MAXSEQTOT)
    {
		fprintf(fres, "Your fasta file contains more than %d sequences, it can't be run on this server. You must download the command line or reduce the number of sequences in your file<BR>\n", MAXSEQTOT);
        fclose (fres);
        exit_properly(ledir);
    }
	// compute distance fowllowing method
    if (nseq <=2)
    {
        fprintf(fres, "You need at least 3 sequences to run  asap<BR>\n");
        fclose (fres); 
        exit_properly(ledir);
    }
	
    if (check_names(mesSeq, nseq,ledir, fres) == 0)
	   {
        fprintf(fres, "Two seqs found with same name. Please correct your alignment<BR>\n");
        fclose (fres); 
        exit_properly(ledir);
        }

	//for (i=0;i<nseq;i++)
	//	fprintf(fres,"%s\n<BR>",mesSeq[i].seq);
	
	my_mat = GetDistMat(nseq, mesSeq, method, ts_tv,ledir, fres);

	//clean everything temporaly needed
	free(myseq);
	free(agarder);
	for (i = 0; i < nseq; i++)
	{free(mesSeq[i].seq); free(mesSeq[i].name);}
	free(mesSeq);
//	print_distmat(my_mat);fflush(stdout);
	return my_mat;
}



/*------------------------------------------------------*/
/*
	Takes a distance file as an input (phylip format)
	Return a struc with a distance matrix
*/
struct DistanceMatrix read_distmat(  FILE *f_in , float ts_tv, char *ledir, FILE *fres)
{

	int a, b, c, d, matinf = 0, nbr = 0;
	int letter;                     /* used to read the file */
	long gpos;


	struct DistanceMatrix my_mat;   /* store distance matrix, names and matrix size */



	my_mat.n = 0;
	my_mat.names = NULL;
	my_mat.dist = NULL;
	my_mat.ratio_ts_tv = ts_tv;



	fscanf( f_in, "%ld", &my_mat.n);          /* phylip give the sequence number on the first line */

	while ( (letter = fgetc(f_in)) != '\n');  /* get rid of remaining blanks --if any-- and \n */



	/*
		Get memory
	*/

	my_mat.names = (char **)malloc( (size_t) sizeof(char *)*my_mat.n );
	if ( ! my_mat.names )fprintf(fres, "read_distmat:cannot allocate my_mat.names, bye\n"), fclose (fres), exit_properly(ledir);

	/*	for(a=0;a<my_mat.n; a++){
			my_mat.names[a] = (char *)malloc( (size_t) sizeof(char)*SIZE_NAME_DIST );
			if( ! my_mat.names[a] )
				fprintf(stderr, "read_distmat: cannot allocate my_mat.names[%d], bye",a), exit(4);
		}*/

	my_mat.dist = (double **)malloc( (size_t) sizeof(double *)*my_mat.n );
	if ( ! my_mat.dist )fprintf(fres, "read_distmat:cannot allocate my_mat.dist, bye"), fclose (fres), exit_properly(ledir);
	for (a = 0; a < my_mat.n; a++) {
		my_mat.dist[a] = (double *)malloc( (size_t) sizeof(double) * my_mat.n );
		if ( ! my_mat.dist[a] )
			fprintf(fres, "read_distmat: cannot allocate my_mat.dist[%d], bye", a), fclose (fres), exit_properly(ledir);
	}

//	fprintf(stderr,"%d to read\n",my_mat.n);
	for (a = 0; a < my_mat.n; a++) {


		/*
			read name
		*/
		c = d = 0;

		gpos = ftell(f_in);
		while (( (letter = fgetc(f_in)) != ' ') && letter !='\t' && !feof(f_in))
			{d++;}
		if (feof(f_in)) fprintf(fres, "read_distmat: PB in dist file format\n" ), fclose (fres), exit_properly(ledir);
		fseek(f_in, gpos, SEEK_SET);
	//	fprintf(stderr,"-->%d\n",d);
		my_mat.names[a] = (char *)malloc( (size_t) sizeof(char) * (d + 1) );
		while ( (letter = fgetc(f_in)) != ' ' && letter !='\t') {

			if (c < d) {

				my_mat.names[a][c] = (char)letter;
				c++;
			}

		}
		my_mat.names[a][c] = 0;
//		fprintf(stderr,"%d %s\n",a,my_mat.names[a]);

		if (a == 0) //for 1st line see if the matrix is full or lower will crash if upper
		{
			c = 0;
			gpos = ftell(f_in);
			while ( ( (letter = fgetc(f_in)) != '\n') && !feof(f_in))
				if (letter == ' ' ||  letter == '\t') c++;
//			printf("%d trouves sur ligne\n",c);	
			if (c < my_mat.n - 1)
				{matinf = 1;
				fprintf(stderr,"mat INF found\n");}
			fseek(f_in, gpos, SEEK_SET);
		}
		/*
			read all dist
		*/

		if (matinf)
		{
//			printf("read matrix inf<BR>");
			for (b = 0; b < nbr; b++)
			{
				fscanf( f_in, "%lf", ( my_mat.dist[a] + b) );
				my_mat.dist[b][a] = my_mat.dist[a][b];
			}
			nbr++;
		}
		else
			for (b = 0; b < my_mat.n; b++)
				fscanf( f_in, "%lf", ( my_mat.dist[a] + b) );


		while ( (letter = fgetc(f_in)) != '\n' && (letter != EOF));
	}


fprintf(stderr,"DONE\n");
	return my_mat;

}









void free_distmat(  struct DistanceMatrix mat ) {


	int a;
	for (a = 0; a < mat.n; a++) {
		free(mat.dist[a]);
		free(mat.names[a]);
	}

	free(mat.dist);
	free(mat.names);
}

