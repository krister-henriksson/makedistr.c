
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#define False 0
#define True 1
#define EPS 1e-16
#define LINESIZE 20
#define WORDSIZE 30
#define DEBUG 0


int mystrlen(char *ins){
  int n=0;
  while (*ins++ != '\0') n++;
  return n;
  /* Note: *ins++ same as *(ins++) due to operator (precedence and) associativity.
     Read as:
     1. Denote value of expression inside parentheses with V.
     V equals current value of ins.
     ins is set to be incremented after complete expression has been evaluated,
     i.e. after (*ins++ != '\0').
     2. Get value residing at address *V.
  */
}

int getline(FILE *fp, char **buf, int *bufsize);

int get_substrings( char *line,
		    char *delims,
		    char ***p_arg,
		    int *p_narg_alloc,
		    int **p_sizes
		    );



int main(int argc, char **argv){
  int sel, Nfiles;
  int i, j, k, Nx;
  double xmin, xmax, norm, t, binw, *x, *ft, *err_ft;
  int Nbins;
  int col;
  char *flag, **filename;
  int normYES;
  double x1, x2;
  int x1_used, x2_used;
  int use_xlog=0;
  int use_ylog=0;
  int use_dylog=0;
  int lenmax, len;



  /* ###############################################################################
     Help on usage.
     ############################################################################### */
  if (argc < 2){
    printf("Purpose: Make distribution.\n");
    printf("Usage:\n");
    printf("     %s \n", argv[0]);
    printf("Flags:\n");
    printf("     -x1        Consider values >=x1. Default: not used, all values considered\n");
    printf("     -x2        Consider values <=x2. Default: not used, all values considered\n");
    printf("     -nn        Do not normalize.\n");
    printf("     -c col     Column containing data. First column is denoted 0, etc. Default: 0\n");
    printf("     -N Nbins   Number of bins. Default: 100\n");
    printf("     -d binw    Bin width.\n");
    printf("     -f files   Data files. Use rest of line.\n");
    printf("\n");
    printf("\n");


    return(0);
  }

  x1_used = x2_used = 0;
  x1 = x2 = 0.0;
  normYES = 1;
  col = 0;
  Nbins = 100;
  sel = 2;
  Nfiles = 0;


  /* ###############################################################################
     Read arguments.
     ############################################################################### */
  for (i=1; i<argc; i++){
    lenmax=len=0;
    len=mystrlen(argv[i]);
    if (len>=lenmax) lenmax=len;
  }
  flag = (char *)calloc(lenmax+1, sizeof(char));

  sel = 0;
  for (i=1; i<argc; i++){
    strcpy(flag, argv[i]);

    if (strcmp(flag, "-x1") == 0){
      sscanf(argv[i+1], "%lg", &x1); i++; x1_used=1;
    }
    if (strcmp(flag, "-x2") == 0){
      sscanf(argv[i+1], "%lg", &x2); i++; x2_used=1;
    }
    if (strcmp(flag, "-nn") == 0){
      normYES = 0;
    }
    if (strcmp(flag, "-c") == 0){
      sscanf(argv[i+1], "%d", &col); i++;
    }
    if (strcmp(flag, "-N") == 0){
      sscanf(argv[i+1], "%d", &Nbins); i++; sel = 2;
    }
    if (strcmp(flag, "-d") == 0){
      sscanf(argv[i+1], "%lg", &binw); i++; sel = 1;
    }
    
    if (strcmp(flag, "-f") == 0){
      Nfiles = argc-i-1;
      lenmax=0, len=0;
      for (j=0; j<Nfiles; j++){
	len=mystrlen(argv[i+1+j]);
	if (len>=lenmax) lenmax=len;
      }

      filename = (char **)calloc( Nfiles, sizeof(char*));
      for (j=0; j<Nfiles; j++)
	filename[j] = (char *)calloc( lenmax+1, sizeof(char));
      
      for (j=0; j<Nfiles; j++){
	sscanf(argv[i+1+j], "%s", filename[j]);
      }
      
    }
  }


  Nx = 1000;
  x = (double *)calloc( Nx, sizeof(double));


  /* ###############################################################################
     Read data.
     ############################################################################### */
  {
    char *line=NULL, **arg;
    int narg, narg_alloc=0, lineSize=0, *strSizes=NULL;
    FILE *fp=NULL;
    
    i = 0;
    for (j=0; j<Nfiles; j++) {
      fp = fopen(filename[j], "r");

      /* Read a line */
      while (getline(fp, &line, &lineSize) != EOF){

	if ((line[0]=='#') || (line[0]=='!')) continue;

	narg = get_substrings(line, " \t", &arg, &narg_alloc, &strSizes);
	if (narg==0) continue;

	i++;

	if (i > Nx){
	  Nx *= 1.5;
	  x = (double *)realloc( x, Nx * sizeof(double));
	}

	sscanf(arg[col], "%lg", &(x[i-1]));
      }
      fclose(fp);
    }

    Nx = i;
    x = (double *)realloc( x, Nx * sizeof(double));
    if (x == NULL){
      printf("Unable to reallocate memory to %d sizeof(double). Exiting.\n", Nx);
      fflush(stdout);
      exit(1);
    }

    free(line);
    for (i=0; i<narg_alloc; ++i) free(arg[i]); free(arg);
    free(strSizes);
  }


  /* ###############################################################################
     Get min and max data points.
     ############################################################################### */
  j = 0;
  k = 0;
  xmin = 0.0;
  xmax = xmin;
  for (i=0; i<Nx; i++) {
    if (j==0 || (j>0 && x[i]<xmin) ){
      xmin = x[i]; j++;
    }
    if (k==0 || (k>0 && x[i]>xmax) ){
      xmax = x[i]; k++;
    }
  }
  if (x1_used) xmin = x1;
  if (x2_used) xmax = x2;

  if (sel==1){
    Nbins = floor( (xmax - xmin )/ binw );
    while (Nbins * binw > xmax) Nbins--;
    if (Nbins==0) Nbins=1;
  }
  else {
    if (Nbins==0) Nbins=1;
    binw = (xmax-xmin) / Nbins;
  }

  ft = (double *)calloc(Nbins, sizeof(double));
  if (ft == NULL) {
    fprintf(stderr, "Memory allocation failed (ft). Exit.\n");
    exit(1);
  }

  err_ft = (double *)calloc(Nbins, sizeof(double));
  if (err_ft == NULL) {
    fprintf(stderr, "Memory allocation failed (err_ft). Exit.\n");
    exit(1);
  }

  for (i=0; i<Nbins; i++){
    ft[i] = 0.0;
  }

  /*
    Put into bins:
  */
  for (i=0; i<Nx; i++){
    if (x1_used && x[i]<xmin) continue;
    if (x2_used && x[i]>xmax) continue;

    if (x[i]<1e-5 || x[i]>1e5) use_xlog=1;

    j = floor( (x[i]-xmin)/binw );
    if (j < 0)      j = 0;
    if (j >= Nbins) j = Nbins-1;
    ft[j] += 1.0;
  }
  
  /* Estimate the uncertainty. */
  for (i=0; i<Nbins; i++){
    err_ft[i] = sqrt(ft[i]);
  }

  if (normYES==1){
    norm = 0.0;
    for (i=0; i<Nbins; i++){
      norm += ft[i];
    }
  }

  for (i=0; i<Nbins; i++){
    if (ft[i]<1e-5 || ft[i]>1e5){
      use_ylog=1; break;
    }
    if (err_ft[i]<1e-5 || err_ft[i]>1e5){
      use_dylog=1;
      break;
    }
  }

  for (i=0; i<Nbins; i++){
    t = xmin + (i + 0.5) * binw;
    if (normYES==1 && Nbins > 1){
      ft[i] /= norm;
      err_ft[i] /= norm;
    }    
    if (use_xlog) printf("%25.10e", t);
    else          printf("%25.10f", t);
    if (use_ylog) printf(" %25.10e", ft[i]);
    else          printf(" %25.10f", ft[i]);
    if (use_dylog) printf(" %25.10e\n", err_ft[i]);
    else           printf(" %25.10f\n", err_ft[i]);
  }

  free(flag);
  for (j=0; j<Nfiles; j++) free(filename[j]);
  free(filename);
  free(x);
  free(ft);
  free(err_ft);

  return(0);
}






/* #############################################################################

   Function: Obtain a line of text from a file.

   Arguments:
   fp        File pointer, assumed initialized.
   buf       Pointer to char array, contains text line when function returns.
   bufsize   Size (number of characters allocated for) of the array 'buf'.

   Return value: Last character read, as integer.

   ########################################################################## */

int getline(FILE *fp, char **buf, int *bufsize){
  int ch, n=0;

  if (*bufsize==0 || *buf==NULL){
    *bufsize = LINESIZE;
    *buf = (char *)calloc(*bufsize, sizeof(char));
  }

  while( 1 ){
    ch = fgetc(fp);
    n++;

    if (n+1 > *bufsize){
      /* Allocate more space. Plus 1 used because terminating nul '\0'
	 will be added at end of string. */
      *bufsize *= 1.5;
      *buf = (char *)realloc( *buf, *bufsize * sizeof(char));
    }

    if (ch=='\n' || ch==EOF) break;

    (*buf)[n-1] = ch;
  }
  (*buf)[n-1] = '\0';

#if DEBUG
  printf("Read line: %s\n", *buf); fflush(stdout);
#endif

  return ch;
}


/* #############################################################################

   Function: Split a string into substrings based on given delimiters.

   Arguments:
   line           Input string, assumed initialized.
   delims         String containing delimiters.
   *p_arg         Output array, whose elements are the substrings.
   *p_narg_alloc  Size of output array.
   *p_sizes       Output array containing the sizes of the substrings.

   Return value: Number of substrings found.

   ########################################################################## */

int get_substrings( char *line,
		    char *delims,
		    char ***p_arg,
		    int *p_narg_alloc,
		    int **p_sizes
		    ){

  char ch;
  int i, j, istart, linesize, delimssize;
  int iLineChar, isubstring, substringIndex;
  int isDelimiter, readingSubstring;


  /* If space is not allocated, then allocate some. */
  if (*p_narg_alloc == 0){
    /* The part of the array that is not allocated is refered to by
       'istart'. This is the starting index of the unallocated part. */
    istart = 0; /* Whole array is unallocated. */

    /* We put some guessed size for the number of substrings: */
    *p_narg_alloc = LINESIZE;
    *p_arg = (char **)calloc(*p_narg_alloc, sizeof(char *));

    /* Now we have an array of pointers to char. These pointers will be
       assigned the addresses of the substrings, which are allocated later.
       For book keeping we need to allocate the array '*p_sizes' which will
       keep track of the sizes of the substrings. */
    *p_sizes = (int *)calloc(*p_narg_alloc, sizeof(int));

    /* We put some guessed size for all substrings: */
    for(i=istart; i<*p_narg_alloc; i++)
      (*p_sizes)[i] = WORDSIZE;
    /* Space for the substrings can now be allocated: */
    for (i=istart; i<*p_narg_alloc; i++)
      (*p_arg)[i] = (char *)calloc((*p_sizes)[i], sizeof(char));

    istart += *p_narg_alloc;
  }
  
  /*
    

   */
  linesize         = mystrlen(line);
  delimssize       = mystrlen(delims);
  isubstring       = 0;
  substringIndex   = 0;
  readingSubstring = False;

  /* Purposefully read also the terminating nul '\0' character! */
  for (iLineChar=0; iLineChar <= linesize; ++iLineChar){

    ch = line[iLineChar];

    isDelimiter = False;
    for (j=0; j<delimssize; ++j){
      if (ch==delims[j]){
	isDelimiter=True;
	break;
      }
    }

#if DEBUG
    if (iLineChar<linesize)
      printf("Char: %c   isDelimiter: %d   readingSubstring: %d  isubstring: %d  substringIndex: %d\n",
	     ch, isDelimiter, readingSubstring, isubstring, substringIndex);
#endif

    if (isDelimiter==True || ch=='\0'){

      if (readingSubstring == True){
	/* If we are in the process of adding to a substring, we now need to
	   finish the current substring and prepare for next one. */
	(*p_arg)[isubstring][substringIndex-1] = '\0';
	isubstring++;
	substringIndex = 0;
	readingSubstring = False;
	
	/* Make sure we have space for additional substrings: */
	if (isubstring > *p_narg_alloc){
	  *p_narg_alloc *= 1.5;
	  *p_arg = (char **)realloc(*p_arg,
				    *p_narg_alloc * sizeof(char *));
	  
	  *p_sizes = (int *)realloc(*p_sizes, *p_narg_alloc * sizeof(int));
	  
	  /* We put some guessed size for all NEW substrings: */
	  for(i=istart; i<*p_narg_alloc; i++)
	    (*p_sizes)[i] = WORDSIZE;
	  /* Space for the NEW substrings can now be allocated: */
	  for (i=istart; i<*p_narg_alloc; i++)
	    (*p_arg)[i] = (char *)calloc((*p_sizes)[i], sizeof(char));
	  
	  istart += *p_narg_alloc;
	}
      }
    }
    else {
      /* The read character is not one of the delimiters. Add it to the
	 current substring. */
      substringIndex++;

      if (substringIndex+1 > (*p_sizes)[isubstring]){
	(*p_sizes)[isubstring] *= 1.5;
	(*p_arg)[isubstring] = (char *)realloc( (*p_arg)[isubstring],
						(*p_sizes)[isubstring] * sizeof(char));
      }

      (*p_arg)[isubstring][substringIndex-1] = ch;
      readingSubstring = True;
    }

  }

  return isubstring;
}

