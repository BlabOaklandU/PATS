/*
	  This file is part of proteinortho.
	  (C) 2009 Marcus Lechner
 
	  proteinortho is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation; either version 2, or (at your
	  option) any later version.
 
	  proteinortho is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
 
	  You should have received a copy of the GNU General Public License
	  along with proteinortho; see the file COPYING.  If not, write to the
	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	  Boston, MA 02111-1307, USA.	
*/

/*
 *  TreebuilderC - core part of proteinorthos treebulider
 *  generates a newicktree from proteinortho's affinity matrix
 *
 *  branchlabels show the number of common genes in the subtrees
 *  branchlengths represent the number of new common genes since the last node
 *
 *  @author Marcus Lechner, Lydia Steiner
 *  @email marcus@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig
 *  @version 1.10 (for PO5)
 *  @date 2016-03-16
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basic-types.h"

BOOL** fetch_matrix(Uint* ccs, Uint* viecher, char* filename, Uint size, char*** pnamen);
void builder(BOOL** matrix, Uint viecher, Uint ccs, char*** pnamen);
void merge (Uint* max, Uint css, Uint viecher, BOOL* away, BOOL** matrix, int** aff, char*** pnamen);
BOOL getmax(BOOL** matrix, BOOL* away, Uint viecher, Uint ccs, int** aff, char*** pnamen);
int sim(BOOL** matrix, Uint a, Uint b, Uint ccs);
char* myitoa(int zahl, int base);
int get_number(char* word);

int main(int argc, char** argv) {
	/* print error if no filename given */
	if (argc != 2) {
		fprintf(stderr,"Error: no filename given.\nUseage: tbc FILENAME\n\n");		
		return 1;
	} 

	/* init variables */
	/* No. of ccs */
	Uint ccs;									
	/* No. of species */
	Uint viecher;	
	/* "string"-array of names, needs to be a pointer for it gets reallocaded by fetch_matrix */	
	char*** pnamen = (char***)malloc(sizeof(char**));							
	/* Just a counter */							
	Uint i;									
	
	fprintf(stderr," >Loading data...");
	/* read data from file */	
	BOOL** matrix = fetch_matrix(&ccs, &viecher, argv[1], 10000, pnamen);
	fprintf(stderr,"done\n");

	fprintf(stderr," >Generating tree...");
	/* Main part: generate the tree from matrix */	
	builder(matrix,viecher,ccs,pnamen);
	fprintf(stderr,"done\n");


	/* substitute every "_" to a whitespace in output tree */	
	i = 0;
	while ((*pnamen)[0][i] != '\0') {
		if ((*pnamen)[0][i] == '_') {(*pnamen)[0][i] = ' ';}
		i++;
	}
	/* and print it to stdout */	
	printf("%s;\n",(*pnamen)[0]);

	/* free memory */
	for(i = 0; i < viecher; i++){free(matrix[i]);}
	free(matrix);
//	for(i = 0; i < viecher; i++){free((*pnamen)[i]);}
	/* all other entries where already freed */
	free((*pnamen)[0]);
	free(*pnamen);free(pnamen);

	return 0;
}

/*  
 *  generates the tree
 *  matrix 	- ccs as adjacent matrix (0/1) 
 *  viecher	- number of species
 *  ccs	- number of ccs
 *  pnamen	- species names ordered the same way as matrix
 *  returns: nothing, but pnamen[0] contains the tree afterwards
 */
void builder(BOOL** matrix, Uint viecher, Uint ccs, char*** pnamen) {
	/* affinity matrix */
	int** aff = (int**)malloc(sizeof(int*)*viecher);
	unsigned int i;
	/* allocate further memory for affinity matrix */
	for (i = 0; i<viecher; i++) {
		aff[i] = (int*)malloc(sizeof(int)*viecher);
		/* initiate every entry with -1 */
		memset(aff[i], -1, sizeof(int)*viecher);
	}
	
	/* array to remember merged species */
	BOOL away[viecher];
	/* initiate every entry with 0 */
	memset(away, 0, sizeof(BOOL)*viecher);

	/* build the tree as long as there are at least 2 species */
	while (getmax(matrix, away, viecher, ccs, aff, pnamen));

	/* pnamen[0] contains the tree now */

	/* free affinity matrix */
	for (i = 0; i<viecher; i++) {
		free(aff[i]);
	}
	free(aff);
}

/*  
 *  merges two species in adjacent matrix
 *  only the common features (= 1 enties) stay, rest becomes 0
 *  loss of one species is marked in away-array
 *  need for recalculation of affinity is marked in aff with -1
 *  besides pnamen for the tree is bulided as merged
 *  max	- contains to entries which should be merged (the one, with the maximum similarity)
 *  ccs	- number of ccs
 *  viecher	- number of species
 *  away	- vector which tells which entry of matrix has been removed (= marked as removed)
 *  matrix 	- ccs as adjacent matrix (0/1)
 *  aff	- affinity matrix
 *  pnamen	- species names ordered the same way as matrix (including merged steps)
 *  returns:  nothing, but all involved data-structures represent the merge afterwards correctly
 */
void merge (Uint* max, Uint ccs, Uint viecher, BOOL* away, BOOL** matrix, int** aff, char*** pnamen) {
	unsigned int i, counter = 0;	

	/* only keep components that are present in both species */
	/* foreach component */	
	for (i = 0; i < ccs; i++) {
		/* if species 1 has it */
		if (matrix[max[0]][i] == '1') {
			/* give it the value of species 2, */
			matrix[max[0]][i] = matrix[max[1]][i];
		}
		/* now there is only a one, if both had it */
	}

	/* counter = affinity of the merged species */
	counter = aff[max[0]][max[1]];

	/* mark the species as away, the entry will not be used again */
	away[max[1]] = 1;
	/* mark every entry within the affinity matrix as 'recaculate me' */
	for (i = 0; i < viecher; i++) {
		/* only the upper triangle matrix is used */
		aff[max[0]][i] = -1;
		/* vice versa not needed */
	}

	/* fetch the number of proteins the specie's subtrees have in common */
	int gemein_0 = get_number((*pnamen)[max[0]]);
	int gemein_1 = get_number((*pnamen)[max[1]]);

	/* transform the number to text 
	   the difference with their affinity is the new branch length (= distance) */
	char* length_0 = myitoa(gemein_0-counter,10);
	char* length_1 = myitoa(gemein_1-counter,10);

	/* transform the number of common genes (= affinity) to text */
	char* anz = myitoa(counter,10);
	/* calculate how many place is needed for the new string 
         sum of both species strings (= subtree) + new numbers + parenteses + : + \0 */
	int newsize = strlen((*pnamen)[max[0]])+strlen((*pnamen)[max[1]])+strlen(anz)+strlen(length_0)+strlen(length_1)+8;
	/* and malloc */
	char* baum = (char*)malloc(sizeof(char)*newsize);
	/* make it an empty string, so strcat can handle it right */
	baum[0] = '\0';

	/* bulid the new string (= merged subtrees) */
	baum = strcat(baum,"(");
	baum = strcat(baum,(*pnamen)[max[0]]);
	baum = strcat(baum,":");
	baum = strcat(baum,length_0);
	baum = strcat(baum,",");
	baum = strcat(baum,(*pnamen)[max[1]]);
	baum = strcat(baum,":");
	baum = strcat(baum,length_1);
	baum = strcat(baum,")");
//	baum = strcat(baum,")[");
//	baum = strcat(baum,anz);
//	baum = strcat(baum,"]");

	/* free the space both merged subtrees */	
	free((*pnamen)[max[0]]);
	free((*pnamen)[max[1]]);
	/* store the merged tree in place of the first subtree */	
	(*pnamen)[max[0]] = baum;

	/* free */
	free(length_0);
	free(length_1);
	/* pointer gotten from myitoa */
	free(anz);

}


/*  
 *  finds the two species that have the highest similarity and merges them
 *  matrix 	- ccs as adjacent matrix (0/1)
 *  away	- vector which tells which entry of matrix has been removed (= marked as removed)
 *  viecher	- number of species
 *  ccs	- number of ccs
 *  aff	- affinity matrix
 *  pnamen	- species names ordered the same way as matrix (including merged steps)
 *  returns:  1 if there was a pair, 0 if there is no species/cluster left to comare (= finished)
 *            pnamen[0] contains the complete tree afterwards
 */
BOOL getmax(BOOL** matrix, BOOL* away, Uint viecher, Uint ccs, int** aff, char*** pnamen) {
	Uint i, j;
	int k, max = -1;
	
	/* Storage for the best pair */
	Uint bestpair[2];
	bestpair[0] = 0;
	bestpair[1] = 0;

	/* foreach species 1 (line) */
	for (i = 0; i < viecher-1; i++) {
		/* skip the linies with away entries */
		if (away[i]) {continue;}

		/* foreach species 2 (row) */
		for (j = i+1; j < viecher; j++) {
			/* skip the rows with away entries */
			if (away[j]) {continue;}
	
			/* recalculate similarity if needed */
			if (aff[i][j] == -1) {aff[i][j] = sim(matrix,i,j,ccs);}
	
			/* if similarity is bigger than best until now */
			k = aff[i][j];
			if (max < k) {
				/* store it and the corresponding indices */
				max = k;
				bestpair[0] = i;
				bestpair[1] = j;
			}
		}
	} 
	
	/* if both entries are the same (as initiated), we are done */
	if (bestpair[0] == bestpair[1]) {
		/* get the last branch length of the tree */
		char* anz = myitoa(get_number((*pnamen)[0]),10);
		/* recalculate size of the tree */
		int newsize = strlen((*pnamen)[0])+strlen(anz)+2;
		/* realloc */
		(*pnamen)[0] = (char*)realloc((*pnamen)[0],sizeof(char)*newsize);
		/* add the information about the root's branchlength */
		(*pnamen)[0] = strcat((*pnamen)[0],":");
		(*pnamen)[0] = strcat((*pnamen)[0],anz);
		/* free myitoa's string */
		free(anz);
		/* return finished */
		return 0;
	}

	/* else, merge the pairs */
	merge(bestpair,ccs,viecher,away,matrix,aff, pnamen);

	/* return done */
	return 1;
}

/*
 *  calculates the similarity between two speices
 *  matrix 	- ccs as adjacent matrix (0/1)
 *  a		- species 1
 *  b		- species 2
 *  ccs	- number of ccs
 *  return	- their simialarity (= number of common proteins)
 */
int sim(BOOL** matrix, Uint a, Uint b, Uint ccs) {
	Uint i;
	int counter = 0;
	/* foreach connected component */
	for (i = 0; i < ccs; i++) {
		/* if both species have an entry, count */
		if (matrix[a][i] == '1' && matrix[b][i] == '1') counter++;
	}
	/* return the sum of common entries */
	return counter;
}

/*
 *  fetches the matrix from a file
 *  ccs	- number of ccs (will be defined here)
 *  viecher	- number of species (will be defined here)
 *  filename- file to fetch the data from
 *  size	- initialy allocated size for the the matrix
 *  pnamen	- species names ordered the same way as in the returned matrix  (will be defined here)
 *  return	- pointer to the matrix
 *
 *  file format:
 *  no. of species
 *  names of species, tab seperated
 *  linewise connected component entries with tab seperated 1 or 0
 */
BOOL** fetch_matrix(Uint* ccs, Uint* viecher, char* filename, Uint size, char*** pnamen) {
	Uint i;

	/* open 0/1 ccs matrix */
	FILE* datei = fopen (filename, "r");
	/* read first line = no. of species */
	fscanf (datei, "%d\n", viecher);
	/* allocate memory for species names */
	*pnamen = (char**) malloc(sizeof(char*)*(*viecher));
	/* read species names and genecount from 2nd line of file */
	for (i = 0; i < (*viecher); i++) {
		/* allocate memory for name */
		/* Attention: dont use strings longer than 255 characters! */
		(*pnamen)[i] = (char*) malloc(sizeof(char)*255);
		/* load name in allocated memory */		
		if (i < (*viecher)-1) {
			/* normal case: */		
			fscanf (datei, "%s ",(*pnamen)[i]);
		}
		else {	
			/* special case: last species name is followed by newline */
			fscanf (datei, "%s \n",(*pnamen)[i]);
		}		
		/* resize memory to what is needed */
		(*pnamen)[i] = (char*) realloc((*pnamen)[i],sizeof(char)*(strlen((*pnamen)[i])+1));
	}
	
	/* number of connected components is 0 at the beginnig */
	*ccs = 0;
	/* allocate memory for 0/1 ccs matrix */
	BOOL** matrix = (BOOL**) malloc(sizeof(BOOL*)*(*viecher));

	/* Rowwise alloc */
	for (i = 0; i < (*viecher); i++) {
		matrix[i] = (BOOL*) malloc(sizeof(BOOL)*size);
	}

	/* will contain the first char read from line */
	BOOL first;
	
	/* read ccs from file, line by line, start with the first char of each line */
	while(fscanf(datei,"%c ",&first) != EOF){
		/* if more ccs than allocated memory -> allocate more memory */
		if(*ccs >= size)	{
			size += 10000;
			for (i = 0; i < (*viecher); i++) {
				matrix[i] = (BOOL*) realloc(matrix[i],sizeof(BOOL)*size); // BOOL*
			}
		}
		/* set first number to what was found first */		
		matrix[0][*ccs] = first;
		/* and gather the rest one by one from file */		
		for (i = 1; i < (*viecher)-1; i++) {
			/* store char in row */					
			fscanf (datei, "%c ",&(matrix[i][*ccs]));
		}
		/* special case: last char is followed by newline */					
		fscanf (datei, "%c \n",&(matrix[(*viecher)-1][*ccs]));
		/* we have a new cc in memory */
		(*ccs)++;
	}
	/* file read throug, close it */
	fclose (datei);
	
	/* resize memory to what is needed for size is likely to be bigger than ccs */
	for (i = 0; i < (*viecher); i++) {
		matrix[i] = (BOOL*) realloc(matrix[i],sizeof(BOOL)*(*ccs));
	}
	
	/* return the matrix pointer */
	return matrix;
}

/*
 *  finds the last [number] entry within a string and converts it to an int
 *  word	- the string to scan
 *  return	- number within the last [number] pattern
 */
int get_number(char* word) {
	int i, start = -1, ende = -1;	
	/* find start and end of the number */
	for (i = strlen(word); i >= 0; i--) {
		if (word[i] == ']') ende = i;
		if (word[i] == '[') {
			start = i+1;
			break;
		}
	}
	/* allocate memory for the number alone */
	char zahl[ende-start+1];
	/* copy the numbe to the new string */
	for (i = 0; i < ende-start; i++) {
		zahl[i] = word[start+i];
	}
	/* add the char finish sign */
	zahl[ende-start] = '\0';
	/* return the number as int */	
	return atoi(zahl);
}

/* 
 *  some compilers do not support itoa, so here is something which works for us 
 *  converts positive numbers to the corresponding text-form as char-array
 *  Attention: up to 64 signs are supported, base between 1 and 16
 *  if 'x' is included in the output, sth. went wrong
 *  zahl	- number to convert (should be postive
 *  base	- base to witch number will be converted
 *  returns - char array representing number zahl to base base as text
 *            Attention: returnvalue needs to be freed separately
 */
char* myitoa(int zahl, int base) {
	int i = 0, k;
	char reverse[64];

	/* Zero is easy, also negative numbers are cought here */
	if(zahl == 0){
		char* zero = (char*) malloc(sizeof(char)*2);
		zero[0] = '0';
		zero[1] = '\0';
		return zero;
	}

	/* convert to a reverse-"string" by modulu-calculation */
	while (zahl > 0) {
		switch(zahl%base) {
	        case  0:	reverse[i] = '0';break;
	        case  1:	reverse[i] = '1';break;
	        case  2:	reverse[i] = '2';break;
	        case  3:	reverse[i] = '3';break;
	        case  4:	reverse[i] = '4';break;
	        case  5:	reverse[i] = '5';break;
	        case  6:	reverse[i] = '6';break;
	        case  7:	reverse[i] = '7';break;
	        case  8:	reverse[i] = '8';break;
	        case  9:	reverse[i] = '9';break;
		  case  10:	reverse[i] = 'A';break;
	        case  11:	reverse[i] = 'B';break;
	        case  12:	reverse[i] = 'C';break;
	        case  13:	reverse[i] = 'D';break;
	        case  14:	reverse[i] = 'E';break;
	        case  15:	reverse[i] = 'F';break;
	        case  16:	reverse[i] = 'G';break;
		  default:  reverse[i] = 'x';
	       }
		/* point one position further */		
		i++;
		zahl /= base;
	}
	/* it's actually one sign less */
	i--;
	/* allocate memory for output */
	char* string = (char*) malloc(sizeof(char)*(i+2));
	/* reverse to calculated char-array to have the right order of signs */
	for (k = i; k >= 0; k--) {
		string[i-k] = reverse[k];
	}
	/* it's a string, so finish it with \0 */
	string[i+1] = '\0';
	/* and return it */
	return string;
}
