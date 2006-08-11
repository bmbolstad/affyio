/*************************************************************
 ** 
 ** It should be noted that Laurent Gautier provided the 
 ** initial CEL file parsing code as part of the file
 ** read_cdf.c.  This code served us well for the 1.0-1.2
 ** releases of Bioconductor. He should be commended for his 
 ** fine code.
 **
 ** The following code became the default parsing code
 ** at the BioC 1.3 release.
 **
 ** Also note that this code was originally part of the affy package
 ** before being moved to affyio
 **
 ************************************************************/
/*************************************************************
 **
 ** file: read_abatch.c
 **
 ** aim: read in from 1st to nth chips of CEL data
 **
 ** Copyright (C) 2003-2006    B. M. Bolstad
 **
 ** Created on Jun 13, 2003
 **
 ** Notes:
 **
 ** The following assumptions are made about text CEL files.
 ** 
 ** 1. A CEL file has a series of sections in the order
 **  
 **  [CEL]
 **  [HEADER]
 **  [INTENSITY]
 **  [MASKS]
 **  [OUTLIERS]
 **
 ** 2. As part of opening the file we will check that
 **    the first characters of the file are "[CEL]"
 **
 ** 3. In the [HEADER] section we expect lines beginning
 **   3a. Cols=
 **   3b. Rows=
 **      3ab.1 Cols should appear before Rows
 **   3c. DatHeader=
 **      3c.1 On the DatHeader line there should appear a
 **          string with the final characters ".1sq". We
 **          will assume that this is the name of the 
 **          CDF file (trim off the ".1sq")
 **
 ** 4. In the [INTENSITY] section there should be 
 **   4a. A line beginning "CellHeader="
 **   4b. After this line there should be cols*rows probe 
 **       intensity lines. Each of these lines should have
 **      4.b.1 Five tokens.
 **      4.b.2 the first token is an integer to be treated as
 **            the x location
 **      4.b.3 the second token is an integer to be treated as
 **            the y location
 **      4.b.4 the third token is a floating point number (double)
 **            to be treated as the probe intensity value.
 **       
 ** 5. The [MASKS] and [OUTLIERS] sections will be treated similarly
 **   5a. We look for a line beginning NumberCells=
 **       this will be the number of Masked or Outlier CELS 
 **       that we will expect to see.
 **   5b. For each line of these sections we will expect to see
 **       the first two items are integers indicating the 
 **       X and Y locations of probes that should be set to NA
 **       if the user sets the right flags.
 **
 ** The implementation of parsing for binary files was based upon 
 ** file format information supplied from Affymetrix.
 **
 **
 ** History
 **
 ** Jun 13, 2003 - Initial version
 ** Jun 14, 2003 - Further implementation
 ** Jun 15, 2003 - testing.
 ** Jun 16, 2003 - Extra verbosity (user controlled).
 ** Jun 17, 2003 - MASKS, OUTLIERS
 **                Added mechanism for reading the header.
 **                this function called ReadHeader
 **                this can be used to check all files same
 **                as first file.
 ** Jun 19, 2003 - Naming of columns of intensity matrix done in
 **                C code 
 ** Jun 29, 2003 - remove some unnecessary cruft from tokenize
 **                routine
 ** Jun 30, 2003 - remove tokenize step from read_cel_intensities. 
 **                aim is to gain more speed.
 ** Jul 1, 2003  - deal with compressed files.
 **                To avoid ugly pre-processor constructs code
 **                for gz functions is seperate from that in
 **                text files.
 **                Made BUF_SIZE 1024
 ** Jul 3, 2003  - add missing pre-processor command
 ** Jul 4, 2003  - A function for reading PM or MM or both
 **                as individual matrices.
 **                Made all the internal functions static
 **                to prevent namespace pollution
 ** Oct 3, 2003  - fix an error in check_cel_file /check_gzcel_file
 **                cdffile was not being properly checked
 ** Oct 14, 2003 - fix a long standing memory leak.
 ** Oct 17, 2003 - integrate binary format reading, Make
 **                it possible to read in a mixture of
 **                binary, text and gzipped text CEL files
 ** Nov 16, 2003 - More work on binary file support (in particular
 **                work on endian issues for no ia32 platforms) 
 ** Nov 20, 2003 - Fix endian bug in readfloat32. Clean up the
 **                some of the documentation
 ** Dec 2, 2003 - Fix up fopen commands (hopefully) for a problem
 **               on w32 machines with text files.
 ** May 7, 2004 - make strncmp strncasecmp. No more problems
 **               with cdf name capitalization. 
 ** May 11, 2004 - fix error message in case when zlib not available
 **                and a gzipped file has been supplied.
 ** Aug 16, 2005 - fix bit flipping when more than one number read
 ** Nov 15, 2005 - ability to read in SD and npixels values
 ** Nov 30, 2005 - remove compress argument from functions where it
 **                appears. it is legacy and has not been used in 
 **                a great deal of time.
 ** Dec 14, 2005 - Added ReadHeaderDetailed: which attempts to 
 **                be more complete then ReadHeader in what it returns
 **                from the header information basically the [CEL] and [HEADER]
 **                sections in the text files and similar information contained in 
 **                the first section of the binary CEL file format
 ** Jan 27, 2006 - Complete ReadHeaderDetailed for supported formats.
 **                Add in a C level object for storing contents of a single
 **                CEL file
 ** May 31, 2006 - Fix some compiler warnings
 ** Jul 17, 2006 - Fix application of masks and outliers for binary cel files.
 ** Jul 21, 2006 - Binary parser checks for file truncation
 ** Aug 11, 2006 - Additional truncation checks for text files
 **
 *************************************************************/
 
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "stdlib.h"
#include "stdio.h"

#define HAVE_ZLIB 1

#if defined(HAVE_ZLIB)
#include <zlib.h>
#endif

#define BUF_SIZE 1024


/****************************************************************
 **
 ** A structure for holding full header information
 **
 **
 **
 ***************************************************************/

typedef struct{
  char *cdfName;
  int cols;
  int rows;
  int GridCornerULx,GridCornerULy;	/* XY coordinates of the upper left grid corner in pixel coordinates.*/
  int GridCornerURx,GridCornerURy;    	/* XY coordinates of the upper right grid corner in pixel coordinates.*/
  int GridCornerLRx,GridCornerLRy;	/* XY coordinates of the lower right grid corner in pixel coordinates.*/
  int GridCornerLLx,GridCornerLLy;      /* XY coordinates of the lower left grid corner in pixel coordinates.*/ 
  char *DatHeader;
  char *Algorithm;
  char *AlgorithmParameters;
} detailed_header_info;



/******************************************************************
 **
 ** A "C" level object designed to hold information for a
 ** single CEL file
 **
 ** These should be created using the function
 **
 ** read_cel_file()
 **
 **
 **
 *****************************************************************/

typedef struct{
  detailed_header_info header;
  
  /** these are for storing the intensities, the sds and the number of pixels **/
  double *intensities;
  double *stddev;
  double *npixels;

  /** these are for storing information in the masks and outliers section **/
  
  int nmasks;
  int noutliers;

  short *masks_x, *masks_y;
  short *outliers_x, *outliers_y;

} CEL;







/****************************************************************
 ****************************************************************
 **
 ** Code for spliting strings into tokens. 
 ** Not heavily used anymore
 **
 ***************************************************************
 ***************************************************************/

/***************************************************************
 **
 ** tokenset
 ** 
 ** char **tokens  - a array of token strings
 ** int n - number of tokens in this set.
 **
 ** a structure to hold a set of tokens. Typically a tokenset is
 ** created by breaking a character string based upon a set of 
 ** delimiters.
 **
 **
 **************************************************************/

typedef struct{
  char **tokens;
  int n;
} tokenset;



/******************************************************************
 **
 ** tokenset *tokenize(char *str, char *delimiters)
 **
 ** char *str - a string to break into tokens
 ** char *delimiters - delimiters to use in breaking up the line
 **
 **
 ** RETURNS a new tokenset
 **
 ** Given a string, split into tokens based on a set of delimitors
 **
 *****************************************************************/

static tokenset *tokenize(char *str, char *delimiters){

  int i=0;

  char *current_token;
  tokenset *my_tokenset = Calloc(1,tokenset);
  my_tokenset->n=0;
  
  my_tokenset->tokens = NULL;

  current_token = strtok(str,delimiters);
  while (current_token != NULL){
    my_tokenset->n++;
    my_tokenset->tokens = Realloc(my_tokenset->tokens,my_tokenset->n,char*);
    my_tokenset->tokens[i] = Calloc(strlen(current_token)+1,char);
    strcpy(my_tokenset->tokens[i],current_token);
    my_tokenset->tokens[i][(strlen(current_token))] = '\0';
    i++;
    current_token = strtok(NULL,delimiters);
  }

  return my_tokenset; 
}


/******************************************************************
 **
 ** int tokenset_size(tokenset *x)
 **
 ** tokenset *x - a tokenset
 ** 
 ** RETURNS the number of tokens in the tokenset 
 **
 ******************************************************************/

static int tokenset_size(tokenset *x){
  return x->n;
}


/******************************************************************
 **
 ** char *get_token(tokenset *x, int i)
 **
 ** tokenset *x - a tokenset
 ** int i - index of the token to return
 ** 
 ** RETURNS pointer to the i'th token
 **
 ******************************************************************/

static char *get_token(tokenset *x,int i){
  return x->tokens[i];
}

/******************************************************************
 **
 ** void delete_tokens(tokenset *x)
 **
 ** tokenset *x - a tokenset
 ** 
 ** Deallocates all the space allocated for a tokenset 
 **
 ******************************************************************/

static void delete_tokens(tokenset *x){
  
  int i;

  for (i=0; i < x->n; i++){
    Free(x->tokens[i]);
  }
  Free(x->tokens);
  Free(x);
}

/*******************************************************************
 **
 ** int token_ends_with(char *token, char *ends)
 ** 
 ** char *token  -  a string to check
 ** char *ends_in   - we are looking for this string at the end of token
 **
 **
 ** returns  0 if no match, otherwise it returns the index of the first character
 ** which matchs the start of *ends.
 **
 ** Note that there must be one additional character in "token" beyond 
 ** the characters in "ends". So
 **
 **  *token = "TestStr"
 **  *ends = "TestStr"   
 **  
 ** would return 0 but if 
 **
 ** ends = "estStr"
 **
 ** we would return 1.
 **
 ** and if 
 ** 
 ** ends= "stStr"
 ** we would return 2 .....etc
 **
 **
 ******************************************************************/

static int token_ends_with(char *token, char *ends_in){
  
  int tokenlength = strlen(token);
  int ends_length = strlen(ends_in);
  int start_pos;
  char *tmp_ptr;
  
  if (tokenlength <= ends_length){
    /* token string is too short so can't possibly end with ends */
    return 0;
  }
  
  start_pos = tokenlength - ends_length;
  
  tmp_ptr = &token[start_pos];

  if (strcmp(tmp_ptr,ends_in)==0){
    return start_pos;
  } else {
    return 0;
  }
}


/****************************************************************
 ****************************************************************
 **
 ** Code for dealing with text CEL files.
 **
 ***************************************************************
 ***************************************************************/

/****************************************************************
 **
 ** void ReadFileLine(char *buffer, int buffersize, FILE *currentFile)
 **
 ** char *buffer  - place to store contents of the line
 ** int buffersize - size of the buffer
 ** FILE *currentFile - FILE pointer to an opened CEL file.
 **
 ** Read a line from a file, into a buffer of specified size.
 ** otherwise die.
 **
 ***************************************************************/

static void ReadFileLine(char *buffer, int buffersize, FILE *currentFile){
  if (fgets(buffer, buffersize, currentFile) == NULL){
    error("End of file reached unexpectedly. Perhaps this file is truncated.\n");
  }  
}	  


/****************************************************************
 **
 ** FILE *open_cel_file(char *filename)
 **
 ** char *filename - name of file to open
 **
 **
 ** RETURNS a file pointer to the open file
 **
 ** this file will open the named file and check to see that the 
 ** first characters agree with "[CEL]" 
 **
 ***************************************************************/

static FILE *open_cel_file(char *filename){
  
  const char *mode = "r";
  FILE *currentFile = NULL; 
  char buffer[BUF_SIZE];

  currentFile = fopen(filename,mode);
  if (currentFile == NULL){
     error("Could not open file %s", filename);
  } else {
    /** check to see if first line is [CEL] so looks like a CEL file**/
    ReadFileLine(buffer, BUF_SIZE, currentFile);
    if (strncmp("[CEL]", buffer, 4) == 0) {
      rewind(currentFile);
    } else {
      error("The file %s does not look like a CEL file",filename);
    }
  }
  
  return currentFile;

}

/******************************************************************
 **
 ** void findStartsWith(FILE *my_file,char *starts, char *buffer)
 **
 ** FILE *my_file - an open file to read from
 ** char *starts - the string to search for at the start of each line
 ** char *buffer - where to place the line that has been read.
 **
 **
 ** Find a line that starts with the specified character string.
 ** At exit buffer should contain that line
 **
 *****************************************************************/


static void  findStartsWith(FILE *my_file,char *starts, char *buffer){

  int starts_len = strlen(starts);
  int match = 1;

  do {
    ReadFileLine(buffer, BUF_SIZE, my_file);
    match = strncmp(starts, buffer, starts_len);
  } while (match != 0);
}


/******************************************************************
 **
 ** void AdvanceToSection(FILE *my_file,char *sectiontitle, char *buffer)
 **
 ** FILE *my_file - an open file
 ** char *sectiontitle - string we are searching for
 ** char *buffer - return's with line starting with sectiontitle
 **
 **
 *****************************************************************/

static void AdvanceToSection(FILE *my_file,char *sectiontitle, char *buffer){
  findStartsWith(my_file,sectiontitle,buffer);
}


/******************************************************************
 ** 
 ** int check_cel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2)
 **
 ** char *filename - the file to read
 ** char *ref_cdfName - the reference CDF filename
 ** int ref_dim_1 - 1st dimension of reference cel file
 ** int ref_dim_2 - 2nd dimension of reference cel file
 **
 ** returns 0 if no problem, 1 otherwise
 **
 ** The aim of this function is to read the header of the CEL file
 ** in particular we will look for the rows beginning "Cols="  and "Rows="
 ** and then for the line DatHeader=  to scope out the appropriate cdf
 ** file. An error() will be flagged if the appropriate conditions
 ** are not met.
 **
 **
 ******************************************************************/

static int check_cel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2){

  int i;
  int dim1,dim2;

  FILE *currentFile; 
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  currentFile = open_cel_file(filename);
  

  AdvanceToSection(currentFile,"[HEADER]",buffer);
  findStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  dim1 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  dim2 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  if ((dim1 != ref_dim_1) || (dim2 != ref_dim_2)){
    error("Cel file %s does not seem to have the correct dimensions",filename);
  }
  
  
  findStartsWith(currentFile,"DatHeader",buffer);
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    if (strncasecmp(get_token(cur_tokenset,i),ref_cdfName,strlen(ref_cdfName)) == 0){
      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be of %s type",filename,ref_cdfName);
    }
  }
  delete_tokens(cur_tokenset);
  fclose(currentFile);

  return 0;
}

/************************************************************************
 **
 ** int read_cel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel intensities for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_cel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index;
  double cur_mean;
  FILE *currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_cel_file(filename);
  
  AdvanceToSection(currentFile,"[INTENSITY]",buffer);
  findStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    if (strlen(buffer) <=2){
      Rprintf("Warning: found an empty line where not expected in %s.\nThis means that there is a cel intensity missing from the cel file.\nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, i);
      break;
    }

    current_token = strtok(buffer," \t");
    if (current_token == NULL){
       Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_x = atoi(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_y = atoi(current_token);
    current_token = strtok(NULL," \t");  
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_mean = atof(current_token);

    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = cur_mean;
    /* delete_tokens(cur_tokenset); */
  }

  fclose(currentFile);

  if (i != rows){
    return 1;
  }


  return 0;
}


/************************************************************************
 **
 ** int read_cel_file_stddev(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel stddev for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_cel_file_stddev(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index;
  double cur_mean, cur_stddev;
  FILE *currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_cel_file(filename);
  
  AdvanceToSection(currentFile,"[INTENSITY]",buffer);
  findStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    if (strlen(buffer) <=2){
      Rprintf("Warning: found an empty line where not expected in %s.\n This means that there is a cel intensity missing from the cel file.\n Sucessfully read to cel intensity %d of %d expected\n", filename, i-1, i);
      break;
    }

    current_token = strtok(buffer," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_x = atoi(current_token);

    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_y = atoi(current_token);
    
    current_token = strtok(NULL," \t");
     if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_mean = atof(current_token);

    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_stddev = atof(current_token);


    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = cur_stddev;
    /* delete_tokens(cur_tokenset); */
  }

  fclose(currentFile);

  if (i != rows){
    return 1;
  }

  return 0;
}




/************************************************************************
 **
 ** int read_cel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel stddev for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_cel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index,cur_npixels;
  double cur_mean, cur_stddev;
  FILE *currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_cel_file(filename);
  
  AdvanceToSection(currentFile,"[INTENSITY]",buffer);
  findStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    if (strlen(buffer) <=2){
      Rprintf("Warning: found an empty line where not expected in %s.\n This means that there is a cel intensity missing from the cel file.\n Sucessfully read to cel intensity %d of %d expected\n", filename, i-1, i);
      break;
    }

    current_token = strtok(buffer," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_x = atoi(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    
    cur_y = atoi(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    
    cur_mean = atof(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_stddev = atof(current_token);
    
    current_token = strtok(NULL," \t");  
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    
    cur_npixels = atoi(current_token);
    
    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = (double)cur_npixels;
    /* delete_tokens(cur_tokenset); */
  }

  fclose(currentFile);
  
  if (i != rows){
    return 1;
  }

  return 0;
}














/****************************************************************
 **
 ** void apply_masks(char *filename, double *intensity, int chip_num, 
 **                   int rows, int cols,int chip_dim_rows, 
 **                   int rm_mask, int rm_outliers)
 **
 ** char *filename    - name of file to open
 ** double *intensity - matrix of probe intensities
 ** int chip_num - the index 0 ...n-1 of the chip we are dealing with
 ** int rows - dimension of the intensity matrix
 ** int cols - dimension of the intensity matrix
 ** int chip_dim_rows - a dimension of the chip
 ** int rm_mask - if true locations in the MASKS section are set NA
 ** int rm_outliers - if true locations in the OUTLIERS section are set NA
 **
 ** This function sets the MASK and OUTLIER probes to NA
 ** 
 **
 ****************************************************************/

static void apply_masks(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers){
  
  int i;
  int numcells, cur_x, cur_y, cur_index;
  FILE *currentFile;
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  if ((!rm_mask) && (!rm_outliers)){
    /* no masking or outliers */
    return;
  }
  
  currentFile = open_cel_file(filename);
  /* read masks section */
  if (rm_mask){

    AdvanceToSection(currentFile,"[MASKS]",buffer);
    findStartsWith(currentFile,"NumberCells=",buffer); 
    cur_tokenset = tokenize(buffer,"=");
    numcells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    findStartsWith(currentFile,"CellHeader=",buffer); 


     for (i =0; i < numcells; i++){
       ReadFileLine(buffer, BUF_SIZE, currentFile);
       
       
       cur_tokenset = tokenize(buffer," \t");
       cur_x = atoi(get_token(cur_tokenset,0));
       cur_y = atoi(get_token(cur_tokenset,1));
       
       cur_index = cur_x + chip_dim_rows*(cur_y);
       intensity[chip_num*rows + cur_index] = R_NaN;
       delete_tokens(cur_tokenset); 
     }
  }

  /* read outliers section */

  if (rm_outliers){
    
    AdvanceToSection(currentFile,"[OUTLIERS]",buffer);
    findStartsWith(currentFile,"NumberCells=",buffer);
    cur_tokenset = tokenize(buffer,"=");
    numcells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    findStartsWith(currentFile,"CellHeader=",buffer); 
    for (i = 0; i < numcells; i++){
      ReadFileLine(buffer, BUF_SIZE, currentFile);      
      cur_tokenset = tokenize(buffer," \t");
      cur_x = atoi(get_token(cur_tokenset,0));
      cur_y = atoi(get_token(cur_tokenset,1));
      
      cur_index = cur_x + chip_dim_rows*(cur_y);
      intensity[chip_num*rows + cur_index] = R_NaReal;
      delete_tokens(cur_tokenset); 
    }
  }
  
  fclose(currentFile);

}


/****************************************************************
 **
 ** static void get_masks_outliers(char *filename, 
 **                         int *nmasks, short **masks_x, short **masks_y, 
 **                         int *noutliers, short **outliers_x, short **outliers_y
 ** 
 ** This gets the x and y coordinates stored in the masks and outliers sections
 ** of the cel files.
 **
 ****************************************************************/

static void get_masks_outliers(char *filename, int *nmasks, short *masks_x, short *masks_y, int *noutliers, short *outliers_x, short *outliers_y){
  
  FILE *currentFile;
  char buffer[BUF_SIZE];
  int numcells, cur_x, cur_y; 
  tokenset *cur_tokenset;
  int i;


  currentFile = open_cel_file(filename);
  /* read masks section */
  

  AdvanceToSection(currentFile,"[MASKS]",buffer);
  findStartsWith(currentFile,"NumberCells=",buffer); 
  cur_tokenset = tokenize(buffer,"=");
  numcells = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  findStartsWith(currentFile,"CellHeader=",buffer); 
  
  *nmasks = numcells;

  masks_x = Calloc(numcells,short);
  masks_y = Calloc(numcells,short);


  for (i =0; i < numcells; i++){
    ReadFileLine(buffer, BUF_SIZE, currentFile);
    
    
    cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    masks_x[i] = (short)cur_x;
    masks_y[i] = (short)cur_y;
    
    delete_tokens(cur_tokenset); 
  }
  
  
  /* read outliers section */
    
  AdvanceToSection(currentFile,"[OUTLIERS]",buffer);
  findStartsWith(currentFile,"NumberCells=",buffer);
  cur_tokenset = tokenize(buffer,"=");
  numcells = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  findStartsWith(currentFile,"CellHeader=",buffer); 

  *noutliers = numcells;
  outliers_x = Calloc(numcells,short);
  outliers_y = Calloc(numcells,short);


  for (i = 0; i < numcells; i++){
    ReadFileLine(buffer, BUF_SIZE, currentFile);      
    cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    Rprintf("%d: %d %d   %d\n",i, cur_x,cur_y, numcells);
    outliers_x[i] = (short)cur_x;
    outliers_y[i] = (short)cur_y;
    
    delete_tokens(cur_tokenset); 
  }
  
  
  fclose(currentFile);




}


/*************************************************************************
 **
 ** char *get_header_info(char *filename, int *dim1, int *dim2)
 **
 ** char *filename - file to open
 ** int *dim1 - place to store Cols
 ** int *dim2 - place to store Rows
 **
 ** returns a character string containing the CDF name.
 **
 ** gets the header information (cols, rows and cdfname)
 **
 ************************************************************************/

static char *get_header_info(char *filename, int *dim1, int *dim2){
  
  int i,endpos;
  char *cdfName = NULL;
  FILE *currentFile; 
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  currentFile = open_cel_file(filename);

  AdvanceToSection(currentFile,"[HEADER]",buffer);
  findStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  *dim1 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  *dim2 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  
  findStartsWith(currentFile,"DatHeader",buffer);
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(cur_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */
      
      cdfName= Calloc(endpos+1,char);
      strncpy(cdfName,get_token(cur_tokenset,i),endpos);
      cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  delete_tokens(cur_tokenset);
  fclose(currentFile);
  return(cdfName);
}


/*************************************************************************
 **
 ** void get_detailed_header_info(char *filename, detailed_header_info *header_info)
 **
 ** char *filename - file to open
 ** detailed_header_info *header_info - place to store header information
 **
 ** reads the header information from a text cdf file (ignoring some fields
 ** that are unused).
 **
 ************************************************************************/

static void get_detailed_header_info(char *filename, detailed_header_info *header_info){

  int i,endpos;
  FILE *currentFile; 
  char buffer[BUF_SIZE];
  char *buffercopy;

  tokenset *cur_tokenset;

  currentFile = open_cel_file(filename);

  AdvanceToSection(currentFile,"[HEADER]",buffer);

  findStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  header_info->cols = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  header_info->rows = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"GridCornerUL",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerULx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerULy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"GridCornerUR",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerURx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerURy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);
  
  findStartsWith(currentFile,"GridCornerLR",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerLRx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerLRy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"GridCornerLL",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerLLx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerLLy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"DatHeader",buffer);
  /* first lets copy the entire string over */

  buffercopy =  Calloc(strlen(buffer)+1,char);
  strcpy(buffercopy,buffer);
  cur_tokenset = tokenize(buffercopy,"\r\n");
  header_info->DatHeader = Calloc(strlen(get_token(cur_tokenset,0))-8,char);
  strcpy(header_info->DatHeader,(get_token(cur_tokenset,0)+10));  /* the +10 is to avoid the starting "DatHeader=" */
  Free(buffercopy);
  delete_tokens(cur_tokenset);

  
  /* now pull out the actual cdfname */ 
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(cur_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */
      
      header_info->cdfName= Calloc(endpos+1,char);
      strncpy( header_info->cdfName,get_token(cur_tokenset,i),endpos);
       header_info->cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  delete_tokens(cur_tokenset);
  
  findStartsWith(currentFile,"Algorithm",buffer);
  cur_tokenset = tokenize(buffer,"=\r\n");
  header_info->Algorithm = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
  strcpy(header_info->Algorithm,get_token(cur_tokenset,1));
  

  delete_tokens(cur_tokenset);

  findStartsWith(currentFile,"AlgorithmParameters",buffer);
  cur_tokenset = tokenize(buffer,"=\r\n");
  header_info->AlgorithmParameters = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
  strcpy(header_info->AlgorithmParameters,get_token(cur_tokenset,1));
  
  fclose(currentFile);

}










/***************************************************************
 **
 ** int isTextCelFile(char *filename)
 **
 ** test whether the file is a valid text cel file
 ** 
 **
 **************************************************************/


static int isTextCelFile(char *filename){

  const char *mode = "r";

  FILE *currentFile= NULL; 
  char buffer[BUF_SIZE];

  currentFile = fopen(filename,mode);
  if (currentFile == NULL){
    error("Could not open file %s", filename);
  } else {
    /** check to see if first line is [CEL] so looks like a CEL file**/
    ReadFileLine(buffer, BUF_SIZE, currentFile);
    fclose(currentFile);
    if (strncmp("[CEL]", buffer, 4) == 0) {
      return 1;
    }
  }
  return 0;
}


/****************************************************************
 ****************************************************************
 **
 ** Code for GZ files starts here.
 **
 ***************************************************************
 ***************************************************************/



#if defined(HAVE_ZLIB)


/****************************************************************
 **
 ** void ReadgzFileLine(char *buffer, int buffersize, FILE *currentFile)
 **
 ** char *buffer  - place to store contents of the line
 ** int buffersize - size of the buffer
 ** FILE *currentFile - FILE pointer to an opened CEL file.
 **
 ** Read a line from a gzipped file, into a buffer of specified size.
 ** otherwise die.
 **
 ***************************************************************/

static void ReadgzFileLine(char *buffer, int buffersize, gzFile currentFile){
  if (gzgets( currentFile,buffer, buffersize) == NULL){
    error("End of gz file reached unexpectedly. Perhaps this file is truncated.\n");
  }  
}


/****************************************************************
 **
 ** FILE *open_gz_cel_file(char *filename)
 **
 ** char *filename - name of file to open
 **
 **
 ** RETURNS a file pointer to the open file
 **
 ** this file will open the named file and check to see that the 
 ** first characters agree with "[CEL]" 
 **
 ***************************************************************/

static gzFile open_gz_cel_file(char *filename){
  
  const char *mode = "rb";

  gzFile currentFile= NULL; 
  char buffer[BUF_SIZE];

  currentFile = gzopen(filename,mode);
  if (currentFile == NULL){
     error("Could not open file %s", filename);
  } else {
    /** check to see if first line is [CEL] so looks like a CEL file**/
    ReadgzFileLine(buffer, BUF_SIZE, currentFile);
    if (strncmp("[CEL]", buffer, 4) == 0) {
      gzrewind(currentFile);
    } else {
      error("The file %s does not look like a CEL file",filename);
    }
  }
  
  return currentFile;

}



/******************************************************************
 **
 ** void gzfindStartsWith(gzFile *my_file,char *starts, char *buffer)
 **
 ** FILE *my_file - an open file to read from
 ** char *starts - the string to search for at the start of each line
 ** char *buffer - where to place the line that has been read.
 **
 **
 ** Find a line that starts with the specified character string.
 ** At exit buffer should contain that line
 **
 *****************************************************************/


static void  gzfindStartsWith(gzFile my_file,char *starts, char *buffer){

  int starts_len = strlen(starts);
  int match = 1;

  do {
    ReadgzFileLine(buffer, BUF_SIZE, my_file);
    match = strncmp(starts, buffer, starts_len);
  } while (match != 0);
}


/******************************************************************
 **
 ** void gzAdvanceToSection(gzFile my_file,char *sectiontitle, char *buffer)
 **
 ** FILE *my_file - an open file
 ** char *sectiontitle - string we are searching for
 ** char *buffer - return's with line starting with sectiontitle
 **
 **
 *****************************************************************/

static void gzAdvanceToSection(gzFile my_file,char *sectiontitle, char *buffer){
  gzfindStartsWith(my_file,sectiontitle,buffer);
}



/******************************************************************
 ** 
 ** int check_gzcel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2)
 **
 ** char *filename - the file to read
 ** char *ref_cdfName - the reference CDF filename
 ** int ref_dim_1 - 1st dimension of reference cel file
 ** int ref_dim_2 - 2nd dimension of reference cel file
 **
 ** returns 0 if no problem, 1 otherwise
 **
 ** The aim of this function is to read the header of the CEL file
 ** in particular we will look for the rows beginning "Cols="  and "Rows="
 ** and then for the line DatHeader=  to scope out the appropriate cdf
 ** file. An error() will be flagged if the appropriate conditions
 ** are not met.
 **
 **
 ******************************************************************/

static int check_gzcel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2){

  int i;
  int dim1,dim2;

  gzFile currentFile; 
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  currentFile = open_gz_cel_file(filename);
  

  gzAdvanceToSection(currentFile,"[HEADER]",buffer);
  gzfindStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  dim1 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  dim2 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  if ((dim1 != ref_dim_1) || (dim2 != ref_dim_2)){
    error("Cel file %s does not seem to have the correct dimensions",filename);
  }
  
  
  gzfindStartsWith(currentFile,"DatHeader",buffer);
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    if (strncasecmp(get_token(cur_tokenset,i),ref_cdfName,strlen(ref_cdfName)) == 0){
      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be of %s type",filename,ref_cdfName);
    }
  }
  delete_tokens(cur_tokenset);
  gzclose(currentFile);

  return 0;
}



/************************************************************************
 **
 ** int read_gzcel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel intensities for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_gzcel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index;
  double cur_mean;
  gzFile currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_gz_cel_file(filename);
  
  gzAdvanceToSection(currentFile,"[INTENSITY]",buffer);
  gzfindStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadgzFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    current_token = strtok(buffer," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_x = atoi(current_token);
    current_token = strtok(NULL," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_y = atoi(current_token);
    current_token = strtok(NULL," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_mean = atof(current_token);

    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = cur_mean;
    /* delete_tokens(cur_tokenset); */
  }

  gzclose(currentFile);
  
  if (i != rows){
    return 1;
  }

  return 0;
}

/************************************************************************
 **
 ** int read_gzcel_file_stddev(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel intensities for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_gzcel_file_stddev(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index;
  double cur_mean, cur_stddev;
  gzFile currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_gz_cel_file(filename);
  
  gzAdvanceToSection(currentFile,"[INTENSITY]",buffer);
  gzfindStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadgzFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    current_token = strtok(buffer," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_x = atoi(current_token);
    current_token = strtok(NULL," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_y = atoi(current_token);
    current_token = strtok(NULL," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_mean = atof(current_token);
  
    current_token = strtok(NULL," \t"); 
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }

    cur_stddev = atof(current_token);
    
    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = cur_stddev;
    /* delete_tokens(cur_tokenset); */
  }

  gzclose(currentFile);
  
  if (i != rows){
    return 1;
  }


  return 0;
}


/************************************************************************
 **
 ** int read_gzcel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols)
 **
 ** char *filename - the name of the cel file to read
 ** double *intensity  - the intensity matrix to fill
 ** int chip_num - the column of the intensity matrix that we will be filling
 ** int rows - dimension of intensity matrix
 ** int cols - dimension of intensity matrix
 **
 ** returns 0 if successful, non zero if unsuccessful
 **
 ** This function reads from the specified file the cel npixels for that
 ** array and fills a column of the intensity matrix.
 **
 ************************************************************************/

static int read_gzcel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){
  
  int i, cur_x,cur_y,cur_index,cur_npixels;
  double cur_mean, cur_stddev;
  gzFile currentFile; 
  char buffer[BUF_SIZE];
  /* tokenset *cur_tokenset;*/
  char *current_token;

  currentFile = open_gz_cel_file(filename);
  
  gzAdvanceToSection(currentFile,"[INTENSITY]",buffer);
  gzfindStartsWith(currentFile,"CellHeader=",buffer);  
  
  for (i=0; i < rows; i++){
    ReadgzFileLine(buffer, BUF_SIZE,  currentFile);
    /* cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    cur_mean = atof(get_token(cur_tokenset,2)); */
    
    current_token = strtok(buffer," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_x = atoi(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_y = atoi(current_token);
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_mean = atof(current_token);
  
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_stddev = atof(current_token);
    
    current_token = strtok(NULL," \t");
    if (current_token == NULL){
      Rprintf("Warning: found an incomplete line where not expected in %s.\nThe CEL file may be truncated. \nSucessfully read to cel intensity %d of %d expected\n", filename, i-1, rows);
      break;
    }
    cur_npixels = atoi(current_token);

    cur_index = cur_x + chip_dim_rows*(cur_y);
    intensity[chip_num*rows + cur_index] = (double)cur_npixels;
    /* delete_tokens(cur_tokenset); */
  }

  gzclose(currentFile);
  
  if (i != rows){
    return 1;
  }

  return 0;
}



/****************************************************************
 **
 ** void gz_apply_masks(char *filename, double *intensity, int chip_num, 
 **                   int rows, int cols,int chip_dim_rows, 
 **                   int rm_mask, int rm_outliers)
 **
 ** char *filename    - name of file to open
 ** double *intensity - matrix of probe intensities
 ** int chip_num - the index 0 ...n-1 of the chip we are dealing with
 ** int rows - dimension of the intensity matrix
 ** int cols - dimension of the intensity matrix
 ** int chip_dim_rows - a dimension of the chip
 ** int rm_mask - if true locations in the MASKS section are set NA
 ** int rm_outliers - if true locations in the OUTLIERS section are set NA
 **
 ** This function sets the MASK and OUTLIER probes to NA
 ** 
 **
 ****************************************************************/

static void gz_apply_masks(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers){
  
  int i;
  int numcells, cur_x, cur_y, cur_index;
  gzFile currentFile;
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  if ((!rm_mask) && (!rm_outliers)){
    /* no masking or outliers */
    return;
  }
  
  currentFile = open_gz_cel_file(filename);
  /* read masks section */
  if (rm_mask){

    gzAdvanceToSection(currentFile,"[MASKS]",buffer);
    gzfindStartsWith(currentFile,"NumberCells=",buffer); 
    cur_tokenset = tokenize(buffer,"=");
    numcells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    gzfindStartsWith(currentFile,"CellHeader=",buffer); 


     for (i =0; i < numcells; i++){
       ReadgzFileLine(buffer, BUF_SIZE, currentFile);
       
       
       cur_tokenset = tokenize(buffer," \t");
       cur_x = atoi(get_token(cur_tokenset,0));
       cur_y = atoi(get_token(cur_tokenset,1));
       
       cur_index = cur_x + chip_dim_rows*(cur_y);
       intensity[chip_num*rows + cur_index] = R_NaN;
       delete_tokens(cur_tokenset); 
     }
  }

  /* read outliers section */

  if (rm_outliers){
    
    gzAdvanceToSection(currentFile,"[OUTLIERS]",buffer);
    gzfindStartsWith(currentFile,"NumberCells=",buffer);
    cur_tokenset = tokenize(buffer,"=");
    numcells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    gzfindStartsWith(currentFile,"CellHeader=",buffer); 
    for (i = 0; i < numcells; i++){
      ReadgzFileLine(buffer, BUF_SIZE, currentFile);      
      cur_tokenset = tokenize(buffer," \t");
      cur_x = atoi(get_token(cur_tokenset,0));
      cur_y = atoi(get_token(cur_tokenset,1));
      
      cur_index = cur_x + chip_dim_rows*(cur_y);
      intensity[chip_num*rows + cur_index] = R_NaReal;
      delete_tokens(cur_tokenset); 
    }
  }
  
  gzclose(currentFile);

}


/*************************************************************************
 **
 ** char *gz_get_header_info(char *filename, int *dim1, int *dim2)
 **
 ** char *filename - file to open
 ** int *dim1 - place to store Cols
 ** int *dim2 - place to store Rows
 **
 ** returns a character string containing the CDF name.
 **
 ** gets the header information (cols, rows and cdfname)
 **
 ************************************************************************/

static char *gz_get_header_info(char *filename, int *dim1, int *dim2){
  
  int i,endpos;
  char *cdfName = NULL;
  gzFile currentFile; 
  char buffer[BUF_SIZE];
  tokenset *cur_tokenset;

  currentFile = open_gz_cel_file(filename);

  gzAdvanceToSection(currentFile,"[HEADER]",buffer);
  gzfindStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  *dim1 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  *dim2 = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  
  gzfindStartsWith(currentFile,"DatHeader",buffer);
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(cur_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */
      
      cdfName= Calloc(endpos+1,char);
      strncpy(cdfName,get_token(cur_tokenset,i),endpos);
      cdfName[endpos] = '\0';

      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  delete_tokens(cur_tokenset);
  gzclose(currentFile);
  return(cdfName);
}




/*************************************************************************
 **
 ** char *gz_get_detailed_header_info(char *filename, detailed_header_info *header_info)
 **
 ** char *filename - file to open
 ** detailed_header_info *header_info - place to store header information
 **
 ** reads the header information from a gzipped text cdf file (ignoring some fields
 ** that are unused).
 **
 ************************************************************************/

static void gz_get_detailed_header_info(char *filename, detailed_header_info *header_info){

  int i,endpos;
  gzFile *currentFile; 
  char buffer[BUF_SIZE];
  char *buffercopy;

  tokenset *cur_tokenset;

  currentFile = open_gz_cel_file(filename);

  gzAdvanceToSection(currentFile,"[HEADER]",buffer);

  gzfindStartsWith(currentFile,"Cols",buffer);  
  cur_tokenset = tokenize(buffer,"=");
  header_info->cols = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"Rows",buffer);
  cur_tokenset = tokenize(buffer,"=");
  header_info->rows = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"GridCornerUL",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerULx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerULy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"GridCornerUR",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerURx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerURy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);
  
  gzfindStartsWith(currentFile,"GridCornerLR",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerLRx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerLRy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"GridCornerLL",buffer);
  cur_tokenset = tokenize(buffer,"= ");
  header_info->GridCornerLLx  = atoi(get_token(cur_tokenset,1));
  header_info->GridCornerLLy  = atoi(get_token(cur_tokenset,2));
  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"DatHeader",buffer);
  /* first lets copy the entire string over */

  buffercopy =  Calloc(strlen(buffer)+1,char);
  strcpy(buffercopy,buffer);
  cur_tokenset = tokenize(buffercopy,"\r\n");
  header_info->DatHeader = Calloc(strlen(get_token(cur_tokenset,0))-8,char);
  strcpy(header_info->DatHeader,(get_token(cur_tokenset,0)+10));  /* the +10 is to avoid the starting "DatHeader=" */
  Free(buffercopy);
  delete_tokens(cur_tokenset);

  
  /* now pull out the actual cdfname */ 
  cur_tokenset = tokenize(buffer," ");
  for (i =0; i < tokenset_size(cur_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(cur_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */
      
      header_info->cdfName= Calloc(endpos+1,char);
      strncpy( header_info->cdfName,get_token(cur_tokenset,i),endpos);
       header_info->cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(cur_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  delete_tokens(cur_tokenset);
  
  gzfindStartsWith(currentFile,"Algorithm",buffer);
  cur_tokenset = tokenize(buffer,"=\r\n");
  header_info->Algorithm = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
  strcpy(header_info->Algorithm,get_token(cur_tokenset,1));
  

  delete_tokens(cur_tokenset);

  gzfindStartsWith(currentFile,"AlgorithmParameters",buffer);
  cur_tokenset = tokenize(buffer,"=\r\n");
  header_info->AlgorithmParameters = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
  strcpy(header_info->AlgorithmParameters,get_token(cur_tokenset,1));
  
  gzclose(currentFile);

}




/****************************************************************
 **
 ** static void gz_get_masks_outliers(char *filename, 
 **                         int *nmasks, short **masks_x, short **masks_y, 
 **                         int *noutliers, short **outliers_x, short **outliers_y
 ** 
 ** This gets the x and y coordinates stored in the masks and outliers sections
 ** of the cel files. (for gzipped text CEL files)
 **
 ****************************************************************/

static void gz_get_masks_outliers(char *filename, int *nmasks, short *masks_x, short *masks_y, int *noutliers, short *outliers_x, short *outliers_y){
  
  gzFile *currentFile;
  char buffer[BUF_SIZE];
  int numcells, cur_x, cur_y;
  tokenset *cur_tokenset;
  int i;


  currentFile = open_gz_cel_file(filename);
  /* read masks section */
  

  gzAdvanceToSection(currentFile,"[MASKS]",buffer);
  gzfindStartsWith(currentFile,"NumberCells=",buffer); 
  cur_tokenset = tokenize(buffer,"=");
  numcells = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  gzfindStartsWith(currentFile,"CellHeader=",buffer); 
  
  *nmasks = numcells;

  masks_x = Calloc(numcells,short);
  masks_y = Calloc(numcells,short);


  for (i =0; i < numcells; i++){
    ReadgzFileLine(buffer, BUF_SIZE, currentFile);
    
    
    cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    masks_x[i] = (short)cur_x;
    masks_y[i] = (short)cur_y;
    
    delete_tokens(cur_tokenset); 
  }
  
  
  /* read outliers section */
    
  gzAdvanceToSection(currentFile,"[OUTLIERS]",buffer);
  gzfindStartsWith(currentFile,"NumberCells=",buffer);
  cur_tokenset = tokenize(buffer,"=");
  numcells = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  gzfindStartsWith(currentFile,"CellHeader=",buffer); 

  *noutliers = numcells;
  outliers_x = Calloc(numcells,short);
  outliers_y = Calloc(numcells,short);


  for (i = 0; i < numcells; i++){
    ReadgzFileLine(buffer, BUF_SIZE, currentFile);      
    cur_tokenset = tokenize(buffer," \t");
    cur_x = atoi(get_token(cur_tokenset,0));
    cur_y = atoi(get_token(cur_tokenset,1));
    /* Rprintf("%d: %d %d   %d\n",i, cur_x,cur_y, numcells); */
    outliers_x[i] = (short)cur_x;
    outliers_y[i] = (short)cur_y;
    
    delete_tokens(cur_tokenset); 
  }
  
  
  gzclose(currentFile);




}








#endif



/***************************************************************
 **
 ** int isgzTextCelFile(char *filename)
 **
 ** test whether the file is a valid gzipped text cel file
 ** 
 **
 **************************************************************/

static int isgzTextCelFile(char *filename){
  
#if defined HAVE_ZLIB
  const char *mode = "rb"; 
 gzFile currentFile = NULL; 
 char buffer[BUF_SIZE];
 currentFile = gzopen(filename,mode);
 if (currentFile == NULL){
   error("Could not open file %s", filename);
 } else {
   /** check to see if first line is [CEL] so looks like a CEL file**/
   ReadgzFileLine(buffer, BUF_SIZE, currentFile);
   gzclose(currentFile);    /* fixed by WH 28 Dec 2003 */
   if (strncmp("[CEL]", buffer, 4) == 0) {
     return 1;
   }
 }
#endif 
 return 0;
}


/***************************************************************
 ***************************************************************
 **
 ** Code for manipulating the cdfInfo
 **
 ***************************************************************
 ***************************************************************/

/*************************************************************************
 **
 ** static int CountCDFProbes(SEXP cdfInfo)
 **
 ** SEXP cdfInfo - a list of matrices, each containing matrix of PM/MM probe 
 **                indicies
 **
 ** returns the number of probes (PM)
 **
 **
 **
 **
 *************************************************************************/


static int CountCDFProbes(SEXP cdfInfo){

  int i;
  int n_probes = 0;
  int n_probesets =  GET_LENGTH(cdfInfo);

  for (i =0; i < n_probesets; i++){
    n_probes +=INTEGER(getAttrib(VECTOR_ELT(cdfInfo,i),R_DimSymbol))[0]; 
  }


  return n_probes;
}


/*************************************************************************
 **
 ** static void  storeIntensities(double *CurintensityMatrix,double *pmMatrix,
 **                               double *mmMatrix, int curcol ,int rows,int cols,
 **                               int chip_dim_rows,SEXP cdfInfo)
 **
 ** double *CurintensityMatrix 
 **
 **
 *************************************************************************/

static void storeIntensities(double *CurintensityMatrix, double *pmMatrix, double *mmMatrix, int curcol, int rows, int cols, int tot_n_probes, SEXP cdfInfo, int which){
  
  int i = 0,j=0, currow=0;
  int n_probes=0;
  int n_probesets = GET_LENGTH(cdfInfo);
  double *cur_index;

  SEXP curIndices;
  
  for (i=0; i < n_probesets; i++){
    curIndices = VECTOR_ELT(cdfInfo,i);
    n_probes = INTEGER(getAttrib(curIndices,R_DimSymbol))[0];
    cur_index = NUMERIC_POINTER(AS_NUMERIC(curIndices));
    
    for (j=0; j < n_probes; j++){
      if (which >= 0){
	pmMatrix[curcol*tot_n_probes + currow] =  CurintensityMatrix[(int)cur_index[j] - 1]; 
      }
      if (which <= 0){
	mmMatrix[curcol*tot_n_probes + currow] =  CurintensityMatrix[(int)cur_index[j+n_probes] - 1];
      }
      currow++;
    }
  }
}



/****************************************************************
 ****************************************************************
 **
 ** These is the code for reading binary CEL files. (Note
 ** not currently viable outside IA32)
 **
 ***************************************************************
 ***************************************************************/

typedef struct{
  int magic_number;
  int version_number;
  int cols;
  int rows;
  int n_cells;
  int header_len;
  char *header;
  int alg_len;
  char *algorithm;
  int alg_param_len;
  char *alg_param;
  int celmargin;
  unsigned int n_outliers;
  unsigned int n_masks;
  int n_subgrids;
  FILE *infile;
} binary_header;



typedef  struct{
    float cur_intens;
    float cur_sd;
    short npixels;
} celintens_record;


typedef struct{
  short x;
  short y;
} outliermask_loc;




/*************************************************************************
 **
 ** Code for reading from the binary files, doing bit flipping if
 ** necessary (on big-endian machines)
 **
 **
 ************************************************************************/


static size_t fread_int32(int *destination, int n, FILE *instream){

  size_t result;
  
  result = fread(destination,sizeof(int),n,instream);

#ifdef WORDS_BIGENDIAN 
  while( n-- > 0 ){
  /* bit flip since all Affymetrix binary files are little endian */
    
    *destination=(((*destination>>24)&0xff) | ((*destination&0xff)<<24) |
		  ((*destination>>8)&0xff00) | ((*destination&0xff00)<<8));  
    destination++;
  }
#endif
  return result;
}



static size_t fread_uint32(unsigned int *destination, int n, FILE *instream){


  size_t result;

  result = fread(destination,sizeof(unsigned int),n,instream);

  
#ifdef WORDS_BIGENDIAN
  while( n-- > 0 ){
    /* bit flip since all Affymetrix binary files are little endian */
    *destination=(((*destination>>24)&0xff) | ((*destination&0xff)<<24) |
		  ((*destination>>8)&0xff00) | ((*destination&0xff00)<<8));  
    destination++;
  } 
#endif
  return result;
}



static size_t fread_int16(short *destination, int n, FILE *instream){
   size_t result;

   result = fread(destination,sizeof(short),n,instream);

#ifdef WORDS_BIGENDIAN
   while( n-- > 0 ){
     /* bit flip since all Affymetrix binary files are little endian */
     *destination=(((*destination>>8)&0xff) | ((*destination&0xff)<<8));  
     destination++;
   }
#endif
   return result;

}

#ifdef WORDS_BIGENDIAN
static void swap_float_4(float *tnf4)              /* 4 byte floating point numbers */
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}
#endif


static size_t fread_float32(float *destination, int n, FILE *instream){

  size_t result;

  
  result = fread(destination,sizeof(float),n,instream);

#ifdef WORDS_BIGENDIAN 
  while( n-- > 0 ){
    swap_float_4(destination);
    destination++;
  }
#endif
  
  return result;
}

static size_t fread_char(char *destination, int n, FILE *instream){

  size_t result;
  
  result = fread(destination,sizeof(char),n,instream);
  
#ifdef WORDS_BIGENDIAN 
  /* Probably don't need to do anything for characters */

#endif

  return result;

}




/*************************************************************
 **
 ** int isBinaryCelFile(char *filename)
 **
 ** filename - Name of the prospective binary cel file
 **
 ** Returns 1 if we find the appropriate parts of the 
 ** header (a magic number of 64 followed by version number of 
 ** 4)
 **
 **
 **
 *************************************************************/

static int isBinaryCelFile(char *filename){

  FILE *infile;

  int magicnumber;
  int version_number;
  
  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }
  
  if (!fread_int32(&magicnumber,1,infile)){
    return 0;
  }
  
  if (!fread_int32(&version_number,1,infile)){
    return 0;
  }


  if (magicnumber != 64){
    return 0;
  }

  if (version_number != 4){
    return 0;
  }

  
  fclose(infile);
  return 1;
}


/*************************************************************
 **
 ** static void delete_binary_header(binary_header *my_header)
 **
 ** binary_header *my_header
 **
 ** frees memory allocated for binary_header structure
 **
 *************************************************************/

static void delete_binary_header(binary_header *my_header){

  Free(my_header->header);
  Free(my_header->algorithm);
  Free(my_header->alg_param);
  Free(my_header);
}


/*************************************************************
 **
 ** static binary_header *read_binary_header(char *filename, int return_stream, FILE *infile)
 **
 ** char *filename - name of binary cel file
 ** int return_stream - if 1 return the stream as part of the header, otherwise close the
 **              file at end of function.
 **
 *************************************************************/

static binary_header *read_binary_header(char *filename, int return_stream){  /* , FILE *infile){ */
  
  FILE *infile;

  binary_header *this_header = Calloc(1,binary_header);
  
  /* Pass through all the header information */
  
  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  
  if (!fread_int32(&(this_header->magic_number),1,infile)){
    error("The binary file %s does not have the appropriate magic number\n",filename);
    return 0;
  }
  
  if (this_header->magic_number != 64){
    error("The binary file %s does not have the appropriate magic number\n",filename);
    return 0;
  }
  
  if (!fread_int32(&(this_header->version_number),1,infile)){
    return 0;
  }

  if (this_header->version_number != 4){
    error("The binary file %s is not version 4. Cannot read\n",filename);
    return 0;
  }
  
  if (!fread_int32(&(this_header->cols),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
    return 0;
  }

  if (!fread_int32(&(this_header->rows),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  

  if (!fread_int32(&(this_header->n_cells),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  
  if (this_header->n_cells != (this_header->cols)*(this_header->rows)){
    error("The number of cells does not seem to be equal to cols*rows in %s.\n",filename);
  }

  
  if (!fread_int32(&(this_header->header_len),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }

  this_header->header = Calloc(this_header->header_len+1,char);
  
  if (!fread(this_header->header,sizeof(char),this_header->header_len,infile)){
    error("binary file corrupted? Could not read any further.\n");
  }
  
  if (!fread_int32(&(this_header->alg_len),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  
  this_header->algorithm = Calloc(this_header->alg_len+1,char);
  
  if (!fread_char(this_header->algorithm,this_header->alg_len,infile)){
    error("binary file corrupted? Could not read any further.\n");
  }
  
  if (!fread_int32(&(this_header->alg_param_len),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  
  this_header->alg_param = Calloc(this_header->alg_param_len+1,char);
  
  if (!fread_char(this_header->alg_param,this_header->alg_param_len,infile)){
    error("binary file corrupted? Could not read any further.\n");
  }
    
  if (!fread_int32(&(this_header->celmargin),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  
  if (!fread_uint32(&(this_header->n_outliers),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }
  
  if (!fread_uint32(&(this_header->n_masks),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  }

  if (!fread_int32(&(this_header->n_subgrids),1,infile)){
    error("Binary file corrupted? Could not read any further\n");
  } 


  if (!return_stream){
    fclose(infile);
  } else {
    this_header->infile = infile;
  }
  
  
  return this_header;




}

/*************************************************************
 **
 ** static char *binary_get_header_info(char *filename, int *dim1, int *dim2)
 **
 ** this function pulls out the rows, cols and cdfname
 ** from the header of a binary cel file
 **
 *************************************************************/

static char *binary_get_header_info(char *filename, int *dim1, int *dim2){
  

  char *cdfName =0;
  tokenset *my_tokenset;

  int i = 0,endpos;
  
  binary_header *my_header;


  my_header = read_binary_header(filename,0);

  *dim1 = my_header->cols;
  *dim2 = my_header->rows;

  my_tokenset = tokenize(my_header->header," ");
    
  for (i =0; i < tokenset_size(my_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(my_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */      
      cdfName= Calloc(endpos+1,char);
      strncpy(cdfName,get_token(my_tokenset,i),endpos);
      cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(my_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  
  delete_binary_header(my_header);
  delete_tokens(my_tokenset);
  return(cdfName);
  
}





/*************************************************************************
 **
 ** void binary_get_detailed_header_info(char *filename, detailed_header_info *header_info)
 **
 ** char *filename - file to open
 ** detailed_header_info *header_info - place to store header information
 **
 ** reads the header information from a binary cdf file (ignoring some fields
 ** that are unused).
 **
 ************************************************************************/





static void binary_get_detailed_header_info(char *filename, detailed_header_info *header_info){

  /* char *cdfName =0; */
  tokenset *my_tokenset;
  tokenset *temp_tokenset;

  char *header_copy;
  char *tmpbuffer;


  
  int i = 0,endpos;
  
  binary_header *my_header;
  

  my_header = read_binary_header(filename,0);


  header_info->cols = my_header->cols;
  header_info->rows = my_header->rows;


  header_info->Algorithm = Calloc(strlen(my_header->algorithm)+1,char);
  
  strcpy(header_info->Algorithm,my_header->algorithm);

  header_info->AlgorithmParameters = Calloc(strlen(my_header->alg_param)+1,char);
  strncpy(header_info->AlgorithmParameters,my_header->alg_param,strlen(my_header->alg_param)-1);
  

  /* Rprintf("%s\n\n\n",my_header->header); */


  header_copy = Calloc(strlen(my_header->header) +1,char);
  strcpy(header_copy,my_header->header);
  my_tokenset = tokenize(header_copy,"\n");

  /** Looking for GridCornerUL, GridCornerUR, GridCornerLR, GridCornerLL and DatHeader */


  for (i =0; i < tokenset_size(my_tokenset);i++){
    /* Rprintf("%d: %s\n",i,get_token(my_tokenset,i)); */
    if (strncmp("GridCornerUL",get_token(my_tokenset,i),12) == 0){
      tmpbuffer = Calloc(strlen(get_token(my_tokenset,i))+1,char);
      strcpy(tmpbuffer,get_token(my_tokenset,i));

      temp_tokenset = tokenize(tmpbuffer,"= ");
      header_info->GridCornerULx  = atoi(get_token(temp_tokenset,1));
      header_info->GridCornerULy  = atoi(get_token(temp_tokenset,2));
      delete_tokens(temp_tokenset);
      Free(tmpbuffer);
    }
    if (strncmp("GridCornerUR",get_token(my_tokenset,i),12) == 0){
      tmpbuffer = Calloc(strlen(get_token(my_tokenset,i))+1,char);
      strcpy(tmpbuffer,get_token(my_tokenset,i));

      temp_tokenset = tokenize(tmpbuffer,"= ");
      header_info->GridCornerURx  = atoi(get_token(temp_tokenset,1));
      header_info->GridCornerURy  = atoi(get_token(temp_tokenset,2));
      delete_tokens(temp_tokenset);
      Free(tmpbuffer);
    }
    if (strncmp("GridCornerLR",get_token(my_tokenset,i),12) == 0){
      tmpbuffer = Calloc(strlen(get_token(my_tokenset,i))+1,char);
      strcpy(tmpbuffer,get_token(my_tokenset,i));

      temp_tokenset = tokenize(tmpbuffer,"= ");
      header_info->GridCornerLRx  = atoi(get_token(temp_tokenset,1));
      header_info->GridCornerLRy  = atoi(get_token(temp_tokenset,2));
      delete_tokens(temp_tokenset);
      Free(tmpbuffer);
    }
    if (strncmp("GridCornerLL",get_token(my_tokenset,i),12) == 0){
      tmpbuffer = Calloc(strlen(get_token(my_tokenset,i))+1,char);
      strcpy(tmpbuffer,get_token(my_tokenset,i));

      temp_tokenset = tokenize(tmpbuffer,"= ");
      header_info->GridCornerLLx  = atoi(get_token(temp_tokenset,1));
      header_info->GridCornerLLy  = atoi(get_token(temp_tokenset,2));
      delete_tokens(temp_tokenset);
      Free(tmpbuffer);
    }
    if (strncmp("DatHeader",get_token(my_tokenset,i),9) == 0){
      header_info->DatHeader = Calloc(strlen(get_token(my_tokenset,i))+1, char);
      strcpy(header_info->DatHeader,(get_token(my_tokenset,i)+10));
    }
  }
    

  delete_tokens(my_tokenset);


  Free(header_copy);

  header_copy = Calloc(my_header->header_len +1,char);
  strcpy(header_copy,my_header->header);
  my_tokenset = tokenize(header_copy," ");
    
  for (i =0; i < tokenset_size(my_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(my_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */      
      header_info->cdfName= Calloc(endpos+1,char);
      strncpy(header_info->cdfName,get_token(my_tokenset,i),endpos);
      header_info->cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(my_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }
  

  delete_tokens(my_tokenset);
  delete_binary_header(my_header);
  Free(header_copy);


}














/***************************************************************
 **
 ** static int check_cel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2)
 ** 
 ** This function checks a binary cel file to see if it has the 
 ** expected rows, cols and cdfname
 **
 **************************************************************/

static int check_binary_cel_file(char *filename, char *ref_cdfName, int ref_dim_1, int ref_dim_2){



  char *cdfName =0;
  tokenset *my_tokenset;

  int i = 0,endpos;
  
  binary_header *my_header;


  my_header = read_binary_header(filename,0);  

  if ((my_header->cols != ref_dim_1) || (my_header->rows != ref_dim_2)){
    error("Cel file %s does not seem to have the correct dimensions",filename);
  }
  
  my_tokenset = tokenize(my_header->header," ");



    
  for (i =0; i < tokenset_size(my_tokenset);i++){
    /* look for a token ending in ".1sq" */
    endpos=token_ends_with(get_token(my_tokenset,i),".1sq");
    if(endpos > 0){
      /* Found the likely CDF name, now chop of .1sq and store it */      
      cdfName= Calloc(endpos+1,char);
      strncpy(cdfName,get_token(my_tokenset,i),endpos);
      cdfName[endpos] = '\0';
      
      break;
    }
    if (i == (tokenset_size(my_tokenset) - 1)){
      error("Cel file %s does not seem to be have cdf information",filename);
    }
  }

  if (strncasecmp(cdfName,ref_cdfName,strlen(ref_cdfName)) != 0){
    error("Cel file %s does not seem to be of %s type",filename,ref_cdfName);
  }

  
  delete_binary_header(my_header);
  delete_tokens(my_tokenset);
  Free(cdfName);



  return 0;
}



/***************************************************************
 **
 ** static int read_binarycel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows)
 **
 ** 
 ** This function reads binary cel file intensities into the data matrix
 **
 **************************************************************/

static int read_binarycel_file_intensities(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){

  int i=0, j=0;
  int cur_index;
  
  int fread_err=0;


  celintens_record *cur_intensity = Calloc(1,celintens_record);
  binary_header *my_header;

  my_header = read_binary_header(filename,1);
  
  for (i = 0; i < my_header->rows; i++){
    for (j =0; j < my_header->cols; j++){
      cur_index = j + my_header->rows*i; /* i + my_header->rows*j; */
      fread_err = fread_float32(&(cur_intensity->cur_intens),1,my_header->infile);
      fread_err+= fread_float32(&(cur_intensity->cur_sd),1,my_header->infile);
      fread_err+=fread_int16(&(cur_intensity->npixels),1,my_header->infile);
      if (fread_err < 3){
	fclose(my_header->infile);
	delete_binary_header(my_header);
	Free(cur_intensity);
	return 1;
      }
      fread_err=0;
      intensity[chip_num*my_header->n_cells + cur_index] = (double )cur_intensity->cur_intens;
    }
  }
  
  fclose(my_header->infile);
  delete_binary_header(my_header);
  Free(cur_intensity);
  return(0);
}




/***************************************************************
 **
 ** static int read_binarycel_file_stdev(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows)
 **
 ** 
 ** This function reads binary cel file stddev values into the data matrix
 **
 **************************************************************/

static int read_binarycel_file_stddev(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){

  int i=0, j=0;
  int cur_index;

  int fread_err=0;
  
  celintens_record *cur_intensity = Calloc(1,celintens_record);
  binary_header *my_header;

  my_header = read_binary_header(filename,1);
  
  for (i = 0; i < my_header->rows; i++){
    for (j =0; j < my_header->cols; j++){
      cur_index = j + my_header->rows*i; /* i + my_header->rows*j; */
      fread_err = fread_float32(&(cur_intensity->cur_intens),1,my_header->infile);
      fread_err+= fread_float32(&(cur_intensity->cur_sd),1,my_header->infile);
      fread_err+= fread_int16(&(cur_intensity->npixels),1,my_header->infile);
      if (fread_err < 3){
	fclose(my_header->infile);
	delete_binary_header(my_header);
	Free(cur_intensity);
	return 1;
      }
      fread_err=0;
      intensity[chip_num*my_header->n_cells + cur_index] = (double )cur_intensity->cur_sd;
    }
  }
  
  fclose(my_header->infile);
  delete_binary_header(my_header);
  Free(cur_intensity);
  return(0);
}




/***************************************************************
 **
 ** static int read_binarycel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows)
 **
 ** 
 ** This function reads binary cel file npixels values into the data matrix
 **
 **************************************************************/

static int read_binarycel_file_npixels(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows){

  int i=0, j=0;
  int cur_index;

  int fread_err=0;
 
  celintens_record *cur_intensity = Calloc(1,celintens_record);
  binary_header *my_header;

  my_header = read_binary_header(filename,1);
  
  for (i = 0; i < my_header->rows; i++){
    for (j =0; j < my_header->cols; j++){
      cur_index = j + my_header->rows*i; /* i + my_header->rows*j; */
      fread_err = fread_float32(&(cur_intensity->cur_intens),1,my_header->infile);
      fread_err+= fread_float32(&(cur_intensity->cur_sd),1,my_header->infile);
      fread_err+= fread_int16(&(cur_intensity->npixels),1,my_header->infile);  
      if (fread_err < 3){
	fclose(my_header->infile);
	delete_binary_header(my_header);
	Free(cur_intensity);
	return 1;
      }
      fread_err=0;
      intensity[chip_num*my_header->n_cells + cur_index] = (double )cur_intensity->npixels;
    }
  }
  
  fclose(my_header->infile);
  delete_binary_header(my_header);
  Free(cur_intensity);
  return(0);
}
















/***************************************************************
 **
 ** static void binary_apply_masks(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers)
 **
 ** 
 **
 **************************************************************/

static void binary_apply_masks(char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers){
  
  int i=0;

  int cur_index;

  int sizeofrecords;

  outliermask_loc *cur_loc= Calloc(1,outliermask_loc);
  binary_header *my_header;

  my_header = read_binary_header(filename,1);
  
  sizeofrecords = 2*sizeof(float) + sizeof(short); /* sizeof(celintens_record) */
  
  fseek(my_header->infile,my_header->n_cells*sizeofrecords,SEEK_CUR);

  if (rm_mask){
    for (i =0; i < my_header->n_masks; i++){
      fread_int16(&(cur_loc->x),1,my_header->infile);
      fread_int16(&(cur_loc->y),1,my_header->infile);
      cur_index = (int)cur_loc->x + my_header->rows*(int)cur_loc->y; 
      /* cur_index = (int)cur_loc->y + my_header->rows*(int)cur_loc->x; */
      /*   intensity[chip_num*my_header->rows + cur_index] = R_NaN; */
      intensity[chip_num*rows + cur_index] =  R_NaN;
      

    }
  } else {
    fseek(my_header->infile,my_header->n_masks*sizeof(cur_loc),SEEK_CUR);

  }

  if (rm_outliers){
    for (i =0; i < my_header->n_outliers; i++){
      fread_int16(&(cur_loc->x),1,my_header->infile);
      fread_int16(&(cur_loc->y),1,my_header->infile);
      cur_index = (int)cur_loc->x + my_header->rows*(int)cur_loc->y; 
      /* intensity[chip_num*my_header->n_cells + cur_index] = R_NaN; */
      intensity[chip_num*rows + cur_index] =  R_NaN;
    }
  } else {
    fseek(my_header->infile,my_header->n_outliers*sizeof(cur_loc),SEEK_CUR);
  }
  
  fclose(my_header->infile);
  delete_binary_header(my_header);
 
  Free(cur_loc);

}

/****************************************************************
 **
 ** static void gz_get_masks_outliers(char *filename, 
 **                         int *nmasks, short **masks_x, short **masks_y, 
 **                         int *noutliers, short **outliers_x, short **outliers_y
 ** 
 ** This gets the x and y coordinates stored in the masks and outliers sections
 ** of the cel files. (for gzipped binary CEL files)
 **
 ****************************************************************/

static void binary_get_masks_outliers(char *filename, int *nmasks, short *masks_x, short *masks_y, int *noutliers, short *outliers_x, short *outliers_y){

  
  int i=0;



  outliermask_loc *cur_loc= Calloc(1,outliermask_loc);
  binary_header *my_header;

  my_header = read_binary_header(filename,1);

  fseek(my_header->infile,my_header->n_cells*sizeof(celintens_record),SEEK_CUR);
 

  *nmasks = my_header->n_masks;
  masks_x = Calloc(my_header->n_masks,short);
  masks_y = Calloc(my_header->n_masks,short);

  for (i =0; i < my_header->n_masks; i++){
    fread_int16(&(cur_loc->x),1,my_header->infile);
    fread_int16(&(cur_loc->y),1,my_header->infile);
    masks_x[i] = (cur_loc->x);
    masks_y[i] = (cur_loc->y);
  }


  *noutliers = my_header->n_outliers;
  outliers_x = Calloc(my_header->n_outliers,short);
  outliers_y = Calloc(my_header->n_outliers,short);
  


  for (i =0; i < my_header->n_outliers; i++){
    fread_int16(&(cur_loc->x),1,my_header->infile);
    fread_int16(&(cur_loc->y),1,my_header->infile);
    outliers_x[i] = (cur_loc->x);
    outliers_y[i] = (cur_loc->y);


  }
      

  delete_binary_header(my_header);
 
  Free(cur_loc);
  

}

/****************************************************************
 ****************************************************************
 **
 ** This is the code that interfaces with R
 **
 ***************************************************************
 ***************************************************************/

/************************************************************************
 **
 **  SEXP read_abatch(SEXP filenames, SEXP compress,  
 **                   SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, 
 **                   SEXP ref_cdfName)
 **
 ** SEXP filenames - an R list of filenames to read
 ** SEXP compress  - logical flag TRUE means files are *.gz
 ** SEXP rm_mask   - if true set MASKS  to NA
 ** SEXP rm_outliers - if true set OUTLIERS to NA
 ** SEXP rm_extra    - if true  overrides rm_mask and rm_outliers settings
 ** SEXP ref_cdfName - the reference CDF name to check each CEL file against 
 ** SEXP ref_dim     - cols/rows of reference chip
 ** SEXP verbose     - if verbose print out more information to the screen
 **
 ** RETURNS an intensity matrix with cel file intensities from
 ** each chip in columns
 **
 ** this function will read in all the cel files in a affybatch.
 ** this function will stop on possible errors with an error() call.
 **
 ** The intensity matrix will be allocated here. It will be given
 ** column names here. the column names that it will be given here are the 
 ** filenames.
 **
 *************************************************************************/

SEXP read_abatch(SEXP filenames, SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, SEXP ref_cdfName, SEXP ref_dim, SEXP verbose){
  
  int i; 
  
  int n_files;
  int ref_dim_1, ref_dim_2;

  char *cur_file_name;
  char *cdfName;
  double *intensityMatrix;

  SEXP intensity,names,dimnames;

  ref_dim_1 = INTEGER(ref_dim)[0];
  ref_dim_2 = INTEGER(ref_dim)[1];
  
  n_files = GET_LENGTH(filenames);
  
  PROTECT(intensity = allocMatrix(REALSXP, ref_dim_1*ref_dim_2, n_files));
  
  cdfName = CHAR(STRING_ELT(ref_cdfName,0));
  intensityMatrix = NUMERIC_POINTER(AS_NUMERIC(intensity));



  /* before we do any real reading check that all the files are of the same cdf type */

  for (i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    if (isTextCelFile(cur_file_name)){
      if (check_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
    } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      if (check_gzcel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
      
#else
       error("Compress option not supported on your platform\n");
#endif
     } else if (isBinaryCelFile(cur_file_name)){
      
       if (check_binary_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	 error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
       }
     } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
     }
  }

  /* 
     Now read in each of the cel files, one by one, filling out the columns of the intensity matrix.
  */
  
  for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (asInteger(verbose)){
	Rprintf("Reading in : %s\n",cur_file_name);
      }
      if (isTextCelFile(cur_file_name)){
	read_cel_file_intensities(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      read_gzcel_file_intensities(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
#else
      error("Compress option not supported on your platform\n");
#endif
    } else if (isBinaryCelFile(cur_file_name)){
	if (read_binarycel_file_intensities(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1)){
	  error("It appears that the file %s is corrupted.\n",cur_file_name);
	}
    } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
    }

  }
  

  /* Now lets go through all the files  filling in masks etc */


  if (asInteger(rm_mask) || asInteger(rm_outliers) || asInteger(rm_extra)){
    for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (isTextCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB	
	if (asInteger(rm_extra)){
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
#else
	error("Compress option not supported on your platform\n");
#endif
      } else if (isBinaryCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
      }
      



    }
  }
  
  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,n_files));
  for ( i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    SET_VECTOR_ELT(names,i,mkChar(cur_file_name));
  }
  SET_VECTOR_ELT(dimnames,1,names);
  setAttrib(intensity, R_DimNamesSymbol, dimnames);
  

  UNPROTECT(3);
  
  return intensity;  
}

/*************************************************************************
 **
 ** SEXP ReadHeader(SEXP filename)
 **
 ** SEXP filename - name of the file to Read.
 **
 ** RETURNS a List containing CDFName, Rows and Cols dimensions.
 ** 
 ** This function reads the HEADER of the CEL file, determines the
 ** CDF name and ROW,COL dimensions
 **
 *************************************************************************/

SEXP ReadHeader(SEXP filename){

  int ref_dim_1=0, ref_dim_2=0;

  char *cur_file_name;
  char *cdfName=0;

  SEXP headInfo;
  SEXP name;
  SEXP cel_dimensions;
  /* SEXP cel_dimensions_names; */
    
  PROTECT(cel_dimensions= allocVector(INTSXP,2));
  PROTECT(headInfo = allocVector(VECSXP,2));

  cur_file_name = CHAR(VECTOR_ELT(filename,0));
  


  /* check for type text, gzipped text or binary then ReadHeader */

  if (isTextCelFile(cur_file_name)){
    cdfName = get_header_info(cur_file_name, &ref_dim_1,&ref_dim_2);
  } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
    cdfName = gz_get_header_info(cur_file_name, &ref_dim_1,&ref_dim_2);
#else
    error("Compress option not supported on your platform\n");
#endif
  } else if (isBinaryCelFile(cur_file_name)){
    cdfName = binary_get_header_info(cur_file_name, &ref_dim_1,&ref_dim_2);
  } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
  }
  
  PROTECT(name = allocVector(STRSXP,1));
  SET_VECTOR_ELT(name,0,mkChar(cdfName));

  INTEGER(cel_dimensions)[0] = ref_dim_1;   /* This is cols */
  INTEGER(cel_dimensions)[1] = ref_dim_2;   /* this is rows */

  SET_VECTOR_ELT(headInfo,0,name);
  SET_VECTOR_ELT(headInfo,1,cel_dimensions);
  
  UNPROTECT(3);

  return headInfo;

}



/*************************************************************************
 **
 ** SEXP ReadHeaderDetailed(SEXP filename)
 **
 ** SEXP filename - name of the file to Read.
 **
 ** RETURNS a List containing CDFName, Rows and Cols dimensions plus more detailed header information
 **
 ** 
 ** This function reads the HEADER of the CEL file
 **
 *************************************************************************/


SEXP ReadHeaderDetailed(SEXP filename){

  SEXP HEADER;
  SEXP tmp_sexp;


  char *cur_file_name;
  detailed_header_info header_info;

  PROTECT(HEADER = allocVector(VECSXP,9)); /* return as a list */


  cur_file_name = CHAR(VECTOR_ELT(filename,0));
 

  if (isTextCelFile(cur_file_name)){
    get_detailed_header_info(cur_file_name,&header_info);
  } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
    gz_get_detailed_header_info(cur_file_name,&header_info);
#else
    error("Compress option not supported on your platform\n");
#endif
  } else if (isBinaryCelFile(cur_file_name)){
    binary_get_detailed_header_info(cur_file_name,&header_info);
  } else {
#if defined HAVE_ZLIB
    error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
    error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
  }

  /* Rprintf("%s\n",header_info.cdfName); */

  /* Copy everything across into the R data structure */
  
  PROTECT(tmp_sexp = allocVector(STRSXP,1));
  SET_VECTOR_ELT(tmp_sexp,0,mkChar(header_info.cdfName));
  SET_VECTOR_ELT(HEADER,0,tmp_sexp);
  UNPROTECT(1);
  PROTECT(tmp_sexp= allocVector(INTSXP,2));
  INTEGER(tmp_sexp)[0] = header_info.cols;   /* This is cols */
  INTEGER(tmp_sexp)[1] = header_info.rows;   /* this is rows */
  SET_VECTOR_ELT(HEADER,1,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp= allocVector(INTSXP,2));
  INTEGER(tmp_sexp)[0] = header_info.GridCornerULx;   
  INTEGER(tmp_sexp)[1] = header_info.GridCornerULy;   
  SET_VECTOR_ELT(HEADER,2,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp= allocVector(INTSXP,2));
  INTEGER(tmp_sexp)[0] = header_info.GridCornerURx;   
  INTEGER(tmp_sexp)[1] = header_info.GridCornerURy;   
  SET_VECTOR_ELT(HEADER,3,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp= allocVector(INTSXP,2));
  INTEGER(tmp_sexp)[0] = header_info.GridCornerLRx;   
  INTEGER(tmp_sexp)[1] = header_info.GridCornerLRy;   
  SET_VECTOR_ELT(HEADER,4,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp= allocVector(INTSXP,2));
  INTEGER(tmp_sexp)[0] = header_info.GridCornerLLx;   
  INTEGER(tmp_sexp)[1] = header_info.GridCornerLLy;   
  SET_VECTOR_ELT(HEADER,5,tmp_sexp);
  UNPROTECT(1);
   
  PROTECT(tmp_sexp = allocVector(STRSXP,1));
  SET_VECTOR_ELT(tmp_sexp,0,mkChar(header_info.DatHeader));
  SET_VECTOR_ELT(HEADER,6,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp = allocVector(STRSXP,1));
  SET_VECTOR_ELT(tmp_sexp,0,mkChar(header_info.Algorithm));
  SET_VECTOR_ELT(HEADER,7,tmp_sexp);
  UNPROTECT(1);

  PROTECT(tmp_sexp = allocVector(STRSXP,1));
  SET_VECTOR_ELT(tmp_sexp,0,mkChar(header_info.AlgorithmParameters));
  SET_VECTOR_ELT(HEADER,8,tmp_sexp);
  UNPROTECT(1);

  Free(header_info.Algorithm);
  Free(header_info.AlgorithmParameters);
  Free(header_info.DatHeader);
  Free(header_info.cdfName);

  UNPROTECT(1);
  return HEADER;
}

/*************************************************************************
 **
 ** SEXP read_probeintensities(SEXP filenames, SEXP compress,  SEXP rm_mask, 
 **                            SEXP rm_outliers, SEXP rm_extra, SEXP ref_cdfName, 
 **                            SEXP ref_dim, SEXP verbose, SEXP cdfInfo)
 **
 ** 
 ** SEXP filenames - an R list of filenames to read
 ** SEXP compress  - logical flag TRUE means files are *.gz
 ** SEXP rm_mask   - if true set MASKS  to NA
 ** SEXP rm_outliers - if true set OUTLIERS to NA
 ** SEXP rm_extra    - if true  overrides rm_mask and rm_outliers settings
 ** SEXP ref_cdfName - the reference CDF name to check each CEL file against
 ** SEXP ref_dim     - cols/rows of reference chip
 ** SEXP verbose     - if verbose print out more information to the screen
 ** SEXP cdfInfo     - locations of probes and probesets
 ** SEXP which       - Indicate whether PM, MM or both are required
 **
 ** returns an R list either one or two elements long. each element is a matrix
 ** either PM or MM elements
 **
 **
 ** This function reads probe intensites into PM and MM matrices. No 
 ** affybatch style matrix is created, Instead the order of the probes is
 ** dependent on the the information given in cdfInfo. cdfInfo is a list
 ** of matrices, each matricies has two columns. The first is assumed to
 ** be PM indices, the second column is assumed to be MM indices.
 **
 **
 *************************************************************************/
 

SEXP read_probeintensities(SEXP filenames,  SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, SEXP ref_cdfName, SEXP ref_dim, SEXP verbose, SEXP cdfInfo,SEXP which){

    
  int i; 
  
  int n_files;
  int ref_dim_1, ref_dim_2;
  int which_flag;  /* 0 means both, 1 means PM only, -1 means MM only */

  int num_probes;

  char *cur_file_name;
  char *cdfName;
  double *CurintensityMatrix, *pmMatrix=0, *mmMatrix=0;

  SEXP PM_intensity= R_NilValue, MM_intensity= R_NilValue, Current_intensity, names, dimnames;
  SEXP output_list,pmmmnames;
  
  if (strcmp(CHAR(STRING_ELT(which,0)),"pm") == 0){
    which_flag= 1;
  } else if (strcmp(CHAR(STRING_ELT(which,0)),"mm") == 0){
    which_flag = -1;
  } else {
    which_flag = 0;
  }



  ref_dim_1 = INTEGER(ref_dim)[0];
  ref_dim_2 = INTEGER(ref_dim)[1];
  
  n_files = GET_LENGTH(filenames);
  
  /* We will read in chip at a time */
  PROTECT(Current_intensity = allocMatrix(REALSXP, ref_dim_1*ref_dim_2, 1));
 
  cdfName = CHAR(STRING_ELT(ref_cdfName,0));
  CurintensityMatrix = NUMERIC_POINTER(AS_NUMERIC(Current_intensity));
  
  /* First check headers of cel files */

  /* before we do any real reading check that all the files are of the same cdf type */

  for (i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    if (isTextCelFile(cur_file_name)){
      if (check_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
    } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      if (check_gzcel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
      
#else
       error("Compress option not supported on your platform\n");
#endif
     } else if (isBinaryCelFile(cur_file_name)){
      
       if (check_binary_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	 error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
       }
     } else {  
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
     }
  }
  

  /* Lets count how many probes we have */
  
  num_probes = CountCDFProbes(cdfInfo);

  if (which_flag >= 0){
    PROTECT(PM_intensity = allocMatrix(REALSXP,num_probes,n_files));
    pmMatrix = NUMERIC_POINTER(AS_NUMERIC(PM_intensity));
  }

  if (which_flag <= 0){
    PROTECT(MM_intensity = allocMatrix(REALSXP,num_probes,n_files));
    mmMatrix = NUMERIC_POINTER(AS_NUMERIC(MM_intensity));
  }

  if (which_flag < 0){
    pmMatrix = NULL;
  }
  
  if (which_flag > 0){
    mmMatrix = NULL;
  }
  
  /* now lets read them in and store them in the PM and MM matrices */

  
  
  for (i=0; i < n_files; i++){ 
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    if (asInteger(verbose)){
      Rprintf("Reading in : %s\n",cur_file_name);
    }
    if (isTextCelFile(cur_file_name)){
      if(read_cel_file_intensities(cur_file_name,CurintensityMatrix, 0, ref_dim_1*ref_dim_2, n_files,ref_dim_1) !=0){
	error("The CEL file %s was corrupted. Data not read.\n",cur_file_name);
      }
      storeIntensities(CurintensityMatrix,pmMatrix,mmMatrix,i,ref_dim_1*ref_dim_2, n_files,num_probes,cdfInfo,which_flag);
    } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      if(read_gzcel_file_intensities(cur_file_name,CurintensityMatrix, 0, ref_dim_1*ref_dim_2, n_files,ref_dim_1)!=0){
	error("The CEL file %s was corrupted. Data not read.\n",cur_file_name);
      }
      storeIntensities(CurintensityMatrix,pmMatrix,mmMatrix,i,ref_dim_1*ref_dim_2, n_files,num_probes,cdfInfo,which_flag);
#else
      error("Compress option not supported on your platform\n");
#endif
    } else if (isBinaryCelFile(cur_file_name)){
      read_binarycel_file_intensities(cur_file_name,CurintensityMatrix, 0, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
      storeIntensities(CurintensityMatrix,pmMatrix,mmMatrix,i,ref_dim_1*ref_dim_2, n_files,num_probes,cdfInfo,which_flag);
    } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
    }
    
  }
  

  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,n_files));
  for ( i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    SET_VECTOR_ELT(names,i,mkChar(cur_file_name));
  }
  SET_VECTOR_ELT(dimnames,1,names);
  if (which_flag >=0){
    setAttrib(PM_intensity, R_DimNamesSymbol, dimnames);
  } 
  
  if (which_flag <=0){
    setAttrib(MM_intensity, R_DimNamesSymbol, dimnames);
  }
  
  if (which_flag == 0){
    PROTECT(output_list = allocVector(VECSXP,2));
    SET_VECTOR_ELT(output_list,0,PM_intensity);
    SET_VECTOR_ELT(output_list,1,MM_intensity);
    PROTECT(pmmmnames = allocVector(STRSXP,2));
    SET_VECTOR_ELT(pmmmnames,0,mkChar("pm"));
    SET_VECTOR_ELT(pmmmnames,1,mkChar("mm"));
  } else if (which_flag > 0){
    PROTECT(output_list = allocVector(VECSXP,1));
    SET_VECTOR_ELT(output_list,0,PM_intensity);
    PROTECT(pmmmnames = allocVector(STRSXP,1));
    SET_VECTOR_ELT(pmmmnames,0,mkChar("pm"));
  } else {
    PROTECT(output_list = allocVector(VECSXP,1));
    SET_VECTOR_ELT(output_list,0,MM_intensity);
    PROTECT(pmmmnames = allocVector(STRSXP,1));
    SET_VECTOR_ELT(pmmmnames,0,mkChar("mm"));
  }

  setAttrib(output_list,R_NamesSymbol,pmmmnames);
  
  if (which_flag != 0){
    UNPROTECT(6);
  } else {
    UNPROTECT(7);
  }
  return(output_list);

}









/************************************************************************
 **
 **  SEXP read_abatch_stddev(SEXP filenames, SEXP compress,  
 **                   SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, 
 **                   SEXP ref_cdfName)
 **
 ** SEXP filenames - an R list of filenames to read
 ** 
 ** SEXP rm_mask   - if true set MASKS  to NA
 ** SEXP rm_outliers - if true set OUTLIERS to NA
 ** SEXP rm_extra    - if true  overrides rm_mask and rm_outliers settings
 ** SEXP ref_cdfName - the reference CDF name to check each CEL file against 
 ** SEXP ref_dim     - cols/rows of reference chip
 ** SEXP verbose     - if verbose print out more information to the screen
 **
 ** RETURNS an intensity matrix with cel file stddev from
 ** each chip in columns
 **
 ** this function will read in all the cel files in a affybatch.
 ** this function will stop on possible errors with an error() call.
 **
 ** The intensity matrix will be allocated here. It will be given
 ** column names here. the column names that it will be given here are the 
 ** filenames.
 **
 *************************************************************************/

SEXP read_abatch_stddev(SEXP filenames,  SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, SEXP ref_cdfName, SEXP ref_dim, SEXP verbose){
  
  int i; 
  
  int n_files;
  int ref_dim_1, ref_dim_2;

  char *cur_file_name;
  char *cdfName;
  double *intensityMatrix;

  SEXP intensity,names,dimnames;

  ref_dim_1 = INTEGER(ref_dim)[0];
  ref_dim_2 = INTEGER(ref_dim)[1];
  
  n_files = GET_LENGTH(filenames);
  
  PROTECT(intensity = allocMatrix(REALSXP, ref_dim_1*ref_dim_2, n_files));
  
  cdfName = CHAR(STRING_ELT(ref_cdfName,0));
  intensityMatrix = NUMERIC_POINTER(AS_NUMERIC(intensity));



  /* before we do any real reading check that all the files are of the same cdf type */

  for (i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    if (isTextCelFile(cur_file_name)){
      if (check_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
    } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      if (check_gzcel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
      
#else
       error("Compress option not supported on your platform\n");
#endif
     } else if (isBinaryCelFile(cur_file_name)){
      
       if (check_binary_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	 error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
       }
     } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
     }
  }

  /* 
     Now read in each of the cel files, one by one, filling out the columns of the intensity matrix.
  */
  
  for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (asInteger(verbose)){
	Rprintf("Reading in : %s\n",cur_file_name);
      }
      if (isTextCelFile(cur_file_name)){
	read_cel_file_stddev(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      read_gzcel_file_stddev(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
#else
      error("Compress option not supported on your platform\n");
#endif
    } else if (isBinaryCelFile(cur_file_name)){
      	if (read_binarycel_file_stddev(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1)){
	  error("It appears that the file %s is corrupted.\n",cur_file_name);
	}
    } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
    }

  }
  

  /* Now lets go through all the files  filling in masks etc */


  if (asInteger(rm_mask) || asInteger(rm_outliers) || asInteger(rm_extra)){
    for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (isTextCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB	
	if (asInteger(rm_extra)){
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
#else
	error("Compress option not supported on your platform\n");
#endif
      } else if (isBinaryCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
      }
      



    }
  }
  
  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,n_files));
  for ( i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    SET_VECTOR_ELT(names,i,mkChar(cur_file_name));
  }
  SET_VECTOR_ELT(dimnames,1,names);
  setAttrib(intensity, R_DimNamesSymbol, dimnames);
  

  UNPROTECT(3);
  
  return intensity;  
}










/************************************************************************
 **
 **  SEXP read_abatch_npixels(SEXP filenames, SEXP compress,  
 **                   SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, 
 **                   SEXP ref_cdfName)
 **
 ** SEXP filenames - an R list of filenames to read
 ** 
 ** SEXP rm_mask   - if true set MASKS  to NA
 ** SEXP rm_outliers - if true set OUTLIERS to NA
 ** SEXP rm_extra    - if true  overrides rm_mask and rm_outliers settings
 ** SEXP ref_cdfName - the reference CDF name to check each CEL file against 
 ** SEXP ref_dim     - cols/rows of reference chip
 ** SEXP verbose     - if verbose print out more information to the screen
 **
 ** RETURNS an intensity matrix with cel file stddev from
 ** each chip in columns
 **
 ** this function will read in all the cel files in a affybatch.
 ** this function will stop on possible errors with an error() call.
 **
 ** The intensity matrix will be allocated here. It will be given
 ** column names here. the column names that it will be given here are the 
 ** filenames.
 **
 *************************************************************************/

SEXP read_abatch_npixels(SEXP filenames,  SEXP rm_mask, SEXP rm_outliers, SEXP rm_extra, SEXP ref_cdfName, SEXP ref_dim, SEXP verbose){
  
  int i; 
  
  int n_files;
  int ref_dim_1, ref_dim_2;

  char *cur_file_name;
  char *cdfName;
  double *intensityMatrix;

  SEXP intensity,names,dimnames;

  ref_dim_1 = INTEGER(ref_dim)[0];
  ref_dim_2 = INTEGER(ref_dim)[1];
  
  n_files = GET_LENGTH(filenames);
  
  PROTECT(intensity = allocMatrix(REALSXP, ref_dim_1*ref_dim_2, n_files));
  
  cdfName = CHAR(STRING_ELT(ref_cdfName,0));
  intensityMatrix = NUMERIC_POINTER(AS_NUMERIC(intensity));



  /* before we do any real reading check that all the files are of the same cdf type */

  for (i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    if (isTextCelFile(cur_file_name)){
      if (check_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
    } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      if (check_gzcel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
      }
      
#else
       error("Compress option not supported on your platform\n");
#endif
     } else if (isBinaryCelFile(cur_file_name)){
      
       if (check_binary_cel_file(cur_file_name,cdfName, ref_dim_1, ref_dim_2)){
	 error("File %s does not seem to have correct dimension or is not of %s chip type.", cur_file_name, cdfName);
       }
     } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
     }
  }

  /* 
     Now read in each of the cel files, one by one, filling out the columns of the intensity matrix.
  */
  
  for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (asInteger(verbose)){
	Rprintf("Reading in : %s\n",cur_file_name);
      }
      if (isTextCelFile(cur_file_name)){
	read_cel_file_npixels(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB
      read_gzcel_file_npixels(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1);
#else
      error("Compress option not supported on your platform\n");
#endif
    } else if (isBinaryCelFile(cur_file_name)){
	if (read_binarycel_file_npixels(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1)){
	  error("It appears that the file %s is corrupted.\n",cur_file_name);
	}
    } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
    }

  }
  

  /* Now lets go through all the files  filling in masks etc */


  if (asInteger(rm_mask) || asInteger(rm_outliers) || asInteger(rm_extra)){
    for (i=0; i < n_files; i++){ 
      cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
      if (isTextCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else if (isgzTextCelFile(cur_file_name)){
#if defined HAVE_ZLIB	
	if (asInteger(rm_extra)){
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  gz_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
#else
	error("Compress option not supported on your platform\n");
#endif
      } else if (isBinaryCelFile(cur_file_name)){
	if (asInteger(rm_extra)){
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,1,1);
	} else {
	  binary_apply_masks(cur_file_name,intensityMatrix, i, ref_dim_1*ref_dim_2, n_files,ref_dim_1,asInteger(rm_mask),asInteger(rm_outliers));
	}
      } else {
#if defined HAVE_ZLIB
       error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",cur_file_name);
#else
       error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",cur_file_name);
#endif
      }
      



    }
  }
  
  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,n_files));
  for ( i =0; i < n_files; i++){
    cur_file_name = CHAR(VECTOR_ELT(VECTOR_ELT(filenames,i),0));
    SET_VECTOR_ELT(names,i,mkChar(cur_file_name));
  }
  SET_VECTOR_ELT(dimnames,1,names);
  setAttrib(intensity, R_DimNamesSymbol, dimnames);
  

  UNPROTECT(3);
  
  return intensity;  
}




/****************************************************************
 ****************************************************************
 **
 ** The following is for "C" code interfacing with CEL files
 ** in a manner other than above.
 **
 ***************************************************************
 ***************************************************************/


/************************************************************************
 **
 ** CEL *read_cel_file(char *filename)
 **  
 ** Reads the contents of the CEL file into a "CEL" structure.
 ** Currently slightly inefficient (should be reimplemented more
 ** cleanly later)
 **
 **
 ************************************************************************/

CEL *read_cel_file(char *filename){
  
  CEL *my_CEL;


  my_CEL = Calloc(1, CEL);

  
  /** First get the header information **/


  if (isTextCelFile(filename)){
    get_detailed_header_info(filename,&my_CEL->header);
  } else if (isgzTextCelFile(filename)){
#if defined HAVE_ZLIB
    gz_get_detailed_header_info(filename,&my_CEL->header);
#else
    error("Compress option not supported on your platform\n");
#endif
  } else if (isBinaryCelFile(filename)){
    binary_get_detailed_header_info(filename,&my_CEL->header);
  } else {
#if defined HAVE_ZLIB
    error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",filename);
#else
    error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",filename);
#endif
  }


  /*** Now lets allocate the space for intensities, stdev, npixels ****/

  my_CEL->intensities = Calloc((my_CEL->header.cols)*(my_CEL->header.rows),double);
  my_CEL->stddev = Calloc((my_CEL->header.cols)*(my_CEL->header.rows),double);
  my_CEL->npixels = Calloc((my_CEL->header.cols)*(my_CEL->header.rows),double);

  if (isTextCelFile(filename)){
    read_cel_file_intensities(filename,my_CEL->intensities, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);	
    read_cel_file_stddev(filename,my_CEL->stddev, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
    read_cel_file_npixels(filename,my_CEL->npixels, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
  }  else if (isgzTextCelFile(filename)){
#if defined HAVE_ZLIB
    read_gzcel_file_intensities(filename,my_CEL->intensities, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);	
    read_gzcel_file_stddev(filename,my_CEL->stddev, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
    read_gzcel_file_npixels(filename,my_CEL->npixels, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
#else
    error("Compress option not supported on your platform\n");
#endif
  } else if (isBinaryCelFile(filename)){
    read_binarycel_file_intensities(filename,my_CEL->intensities, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);	
    read_binarycel_file_stddev(filename,my_CEL->stddev, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
    read_binarycel_file_npixels(filename,my_CEL->npixels, 0, (my_CEL->header.cols)*(my_CEL->header.rows), 1,my_CEL->header.cols);
  } else {
#if defined HAVE_ZLIB
    error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",filename);
#else
    error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",filename);
#endif
  }

  /*** Now add masks and outliers ***/

  if (isTextCelFile(filename)){
    get_masks_outliers(filename, &(my_CEL->nmasks), (my_CEL->masks_x), (my_CEL->masks_y), &(my_CEL->noutliers), (my_CEL->outliers_x), (my_CEL->outliers_y));
  } else if (isgzTextCelFile(filename)){
#if defined HAVE_ZLIB
    gz_get_masks_outliers(filename, &(my_CEL->nmasks), (my_CEL->masks_x), (my_CEL->masks_y), &(my_CEL->noutliers), (my_CEL->outliers_x), (my_CEL->outliers_y));
#else
    error("Compress option not supported on your platform\n");
#endif 
  } else if (isBinaryCelFile(filename)){



  } else {
#if defined HAVE_ZLIB
    error("Is %s really a CEL file? tried reading as text, gzipped text and binary\n",filename);
#else
    error("Is %s really a CEL file? tried reading as text and binary. The gzipped text format is not supported on your platform.\n",filename);
#endif
  }


  return my_CEL;


}


/**************************************************************************
 **
 **
 ** Read a single CEL file into an R list structure
 ** Mostly just for testing the above
 **
 **************************************************************************/




SEXP R_read_cel_file(SEXP filename){

  char *cur_file_name = CHAR(VECTOR_ELT(filename,0));

  CEL *myCEL =read_cel_file(cur_file_name);


  return filename;

}
