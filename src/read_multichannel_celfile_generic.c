/*************************************************************
 **
 ** file: read_multichannel_celfile_generic.c
 **
 ** Written by B. M. Bolstad <bmb@bmbolstad.com>
 **
 ** Aim is to read in Affymetrix CEL files in the
 ** "Command Console Generic Data" File Format
 ** This format is sometimes known as the Calvin format
 **
 ** As with other file format functionality in affyio
 ** gzipped files are accepted.
 **
 ** The implementation here is based upon openly available 
 ** file format information. The code here is not dependent or based
 ** in anyway on that in the Fusion SDK.
 **
 **
 ** History
 ** Sept 3, 2007 -Initial version
 ** Sept 9, 2007 - fix compiler warnings
 ** Oct 11, 2007 - fix missing DatHeader problem
 ** Feb 11, 2008 - add #include for inttypes.h in situations that stdint.h might not exist
 ** Feb 13, 2008 - fix problems with generic_get_detailed_header_info(), gzgeneric_get_detailed_header_info()
 ** May 18, 2009 - Add Ability to extract scan date from CEL file header
 ** May 25, 2010 - Multichannel CELfile support adapted from single channel parser
 ** Sep 4, 2017 - change gzFile* to gzFile
 **
 *************************************************************/
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#include <wchar.h>

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "read_generic.h"
#include "read_celfile_generic.h"
#include "read_multichannel_celfile_generic.h"
#include "read_abatch.h"

int isGenericMultiChannelCelFile(const char *filename){

  FILE *infile;
  generic_file_header file_header;
  generic_data_header data_header;
  
  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }

  if (!read_generic_file_header(&file_header,infile)){
    fclose(infile);
    return 0;
  }

  if (!read_generic_data_header(&data_header,infile)){
    Free_generic_data_header(&data_header);
    fclose(infile);
    return 0;
  }
  
  if (strcmp(data_header.data_type_id.value, "affymetrix-calvin-multi-intensity") !=0){
    Free_generic_data_header(&data_header);
 


   fclose(infile);
    return 0;
  }
  Free_generic_data_header(&data_header);
  
  fclose(infile);
  return 1;
}

static int compare_AWSTRING_Intensity(AWSTRING string){

  int rv = 0;
  if (string.len > 0){
    char *temp = Calloc(string.len+1,char);
    wcstombs(temp, string.value, string.len);
    
    rv = strcmp(temp,"Intensity");
    
    Free(temp);
  }
  return rv;
}



int multichannel_determine_number_channels(const char *filename){
  
  int j=0;
  int returnvalue = 0;
  
  
  FILE *infile;
  
  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;
  
  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);

  do {
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup;
    for (j=0; j < my_data_group.n_data_sets; j++){
      read_generic_data_set(&my_data_set,infile);
      if (!compare_AWSTRING_Intensity(my_data_set.data_set_name)){
	returnvalue++;
        break;
      }
      read_generic_data_set_rows(&my_data_set,infile); 
      Free_generic_data_set(&my_data_set);
    }
    Free_generic_data_group(&my_data_group);
    fseek(infile,next_group,SEEK_SET);
  } while (next_group > 0);	
  
  fclose(infile);
  Free_generic_data_header(&my_data_header);
  
  return(returnvalue);

}






char *multichannel_determine_channel_name(const char *filename, int channelindex){
  
  int k=0;
  char *returnvalue = 0;
  
  FILE *infile;
  
  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;
  
    uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);

  while (k < channelindex){
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    fseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  read_generic_data_group(&my_data_group,infile);
  if (my_data_group.data_group_name.len > 0){
    returnvalue = Calloc(my_data_group.data_group_name.len+1,char);
    wcstombs(returnvalue, my_data_group.data_group_name.value, my_data_group.data_group_name.len);
  } 
  Free_generic_data_group(&my_data_group);
  fclose(infile);
  Free_generic_data_header(&my_data_header);
  
  return(returnvalue);

}



/***************************************************************
 **
 ** static int read_binarycel_file_intensities(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows)
 **
 ** 
 ** This function reads binary cel file intensities into the data matrix
 **
 **************************************************************/

int read_genericcel_file_intensities_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex){

  int i=0, k=0;
  
  FILE *infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);

  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    fseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  /* Now the actual channel of intensities */
  read_generic_data_group(&my_data_group,infile);
  read_generic_data_set(&my_data_set,infile); 
  read_generic_data_set_rows(&my_data_set,infile);

  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((float *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_group(&my_data_group);
  fclose(infile);
  Free_generic_data_header(&my_data_header);
 
  return(0);
}





int read_genericcel_file_stddev_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex){

  int i=0, k=0;
    
  FILE *infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);

  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    fseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;	
  }

  read_generic_data_group(&my_data_group,infile);
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  read_generic_data_set(&my_data_set,infile); 
  read_generic_data_set_rows(&my_data_set,infile); 
  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((float *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  fclose(infile);


  return(0);





}



int read_genericcel_file_npixels_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows,  int channelindex){

  int i=0, k=0;
  
  FILE *infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);


  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){   
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    fseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  read_generic_data_group(&my_data_group,infile);

  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  

  read_generic_data_set(&my_data_set,infile); 
  read_generic_data_set_rows(&my_data_set,infile); 
  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((short *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  fclose(infile);


  return(0);




}




void generic_get_masks_outliers_multichannel(const char *filename, int *nmasks, short **masks_x, short **masks_y, int *noutliers, short **outliers_x, short **outliers_y, int channelindex){

  int i=0, k=0;
    
  FILE *infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      
    }
  

  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);


  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){  
    read_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    fseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  read_generic_data_group(&my_data_group,infile);

  /* passing the intensities */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* passing by the stddev */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  
  /* passing by the npixels */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* Now lets go for the "Outlier" */
  read_generic_data_set(&my_data_set,infile); 
   
  
  *noutliers = my_data_set.nrows;

  *outliers_x = Calloc(my_data_set.nrows,short); 
  *outliers_y = Calloc(my_data_set.nrows,short);
  
  read_generic_data_set_rows(&my_data_set,infile); 
  
  for (i=0; i < my_data_set.nrows; i++){
    (*outliers_x)[i] = ((short *)my_data_set.Data[0])[i];
    (*outliers_y)[i] = ((short *)my_data_set.Data[1])[i];
  }
  
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);

  /* Now lets go for the "Mask" */
  read_generic_data_set(&my_data_set,infile); 
   
  *nmasks = my_data_set.nrows;

  *masks_x = Calloc(my_data_set.nrows,short); 
  *masks_y = Calloc(my_data_set.nrows,short);
  
  
  read_generic_data_set_rows(&my_data_set,infile); 
  for (i=0; i < my_data_set.nrows; i++){
    (*outliers_x)[i] = ((short *)my_data_set.Data[0])[i];
    (*outliers_y)[i] = ((short *)my_data_set.Data[1])[i];
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  fclose(infile);
  
}





void generic_apply_masks_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers,  int channelindex){


  int i=0;
  int cur_index;
  
  short cur_x, cur_y;


  int nrows;
  int size;

  FILE *infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;
  nvt_triplet *triplet;
  AffyMIMEtypes cur_mime_type;

  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      
    }
 

  
  read_generic_file_header(&my_header, infile);
  read_generic_data_header(&my_data_header, infile);
  read_generic_data_group(&my_data_group,infile);

    
  triplet =  find_nvt(&my_data_header,"affymetrix-cel-rows");
  cur_mime_type = determine_MIMETYPE(*triplet);
  decode_MIME_value(*triplet,cur_mime_type, &nrows, &size);
  


  /* passing the intensities */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* passing by the stddev */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  
  /* passing by the npixels */
  read_generic_data_set(&my_data_set,infile); 
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* Now lets go for the "Outlier" */
  read_generic_data_set(&my_data_set,infile); 
  
  if (rm_outliers){
    read_generic_data_set_rows(&my_data_set,infile); 
    for (i=0; i < my_data_set.nrows; i++){
      cur_x = ((short *)my_data_set.Data[0])[i];
      cur_y = ((short *)my_data_set.Data[1])[i];
      cur_index = (int)cur_x + nrows*(int)cur_y; 
      intensity[chip_num*rows + cur_index] =  R_NaN;
    }
  }
  
  fseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);

  /* Now lets go for the "Mask" */
  read_generic_data_set(&my_data_set,infile); 
  if (rm_mask){
    read_generic_data_set_rows(&my_data_set,infile); 
    for (i=0; i < my_data_set.nrows; i++){
      cur_x = ((short *)my_data_set.Data[0])[i];
      cur_y = ((short *)my_data_set.Data[1])[i];
      cur_index = (int)cur_x + nrows*(int)cur_y; 
      intensity[chip_num*rows + cur_index] =  R_NaN;
    }
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  fclose(infile);
  
}

/*******************************************************************************************************
 *******************************************************************************************************
 **
 ** Code below supports gzipped command console format  MultiChannel CEL files
 **
 *******************************************************************************************************
 *******************************************************************************************************/


int isgzGenericMultiChannelCelFile(const char *filename){

  gzFile infile;
  generic_file_header file_header;
  generic_data_header data_header;
  
  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }

  if (!gzread_generic_file_header(&file_header,infile)){
    gzclose(infile);
    return 0;
  }

  if (!gzread_generic_data_header(&data_header,infile)){
    Free_generic_data_header(&data_header);
    gzclose(infile);
    return 0;
  }
  
  if (strcmp(data_header.data_type_id.value, "affymetrix-calvin-multi-intensity") !=0){
    Free_generic_data_header(&data_header);
 


   gzclose(infile);
    return 0;
  }
  Free_generic_data_header(&data_header);
  
  gzclose(infile);
  return 1;
}




/* basic idea is to count how many datagroups have a dataset called "Intensity" */



int gzmultichannel_determine_number_channels(const char *filename){
  
  int j=0;
 
  int returnvalue = 0;
   
  gzFile infile;
  
  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;
  
  generic_data_set my_data_set;

  uint32_t next_group =1;  


  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);
 
	
 
  do {
      gzread_generic_data_group(&my_data_group,infile);    
      next_group = my_data_group.file_position_nextgroup;	
      for (j=0; j < my_data_group.n_data_sets; j++){
         gzread_generic_data_set(&my_data_set,infile);
         if (!compare_AWSTRING_Intensity(my_data_set.data_set_name)){
	   returnvalue++;
 	   break;
       }
       gzread_generic_data_set_rows(&my_data_set,infile); 
       Free_generic_data_set(&my_data_set);
    }
    Free_generic_data_group(&my_data_group);
    gzseek(infile,next_group,SEEK_SET);
  }  while (next_group > 0);	
  
  gzclose(infile);
  Free_generic_data_header(&my_data_header);
  return(returnvalue);

}


char *gzmultichannel_determine_channel_name(const char *filename, int channelindex){
  
  int k=0;

  char *returnvalue = 0;
  
  gzFile infile;
  
  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;
  
  uint32_t next_group =1;  

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);

  while (k < channelindex){
    gzread_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    gzseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  gzread_generic_data_group(&my_data_group,infile);
  if (my_data_group.data_group_name.len > 0){
    returnvalue = Calloc(my_data_group.data_group_name.len+1,char);
    wcstombs(returnvalue, my_data_group.data_group_name.value, my_data_group.data_group_name.len);
  } 
  Free_generic_data_group(&my_data_group);
  gzclose(infile);
  Free_generic_data_header(&my_data_header);
  
  return(returnvalue);

}




int gzread_genericcel_file_intensities_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex){

  int i=0, k=0;
    
  gzFile infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);

  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){
    gzread_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    gzseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  /* Now the actual channel of intensities */
  gzread_generic_data_group(&my_data_group,infile);
  gzread_generic_data_set(&my_data_set,infile); 
  gzread_generic_data_set_rows(&my_data_set,infile);

  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((float *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_group(&my_data_group);
  gzclose(infile);
  Free_generic_data_header(&my_data_header);
 
  return(0);
}





int gzread_genericcel_file_stddev_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows,  int channelindex){

  int i=0, k=0;
    
  gzFile infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);

  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){    
    gzread_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    gzseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }

  gzread_generic_data_group(&my_data_group,infile);
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  gzread_generic_data_set(&my_data_set,infile); 
  gzread_generic_data_set_rows(&my_data_set,infile); 
  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((float *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  gzclose(infile);


  return(0);



}




int gzread_genericcel_file_npixels_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows,  int channelindex){

  int i=0, k=0;
  
  gzFile infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      return 0;
    }
  

  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);


  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){
    gzread_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    gzseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  gzread_generic_data_group(&my_data_group,infile);

  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  

  gzread_generic_data_set(&my_data_set,infile); 
  gzread_generic_data_set_rows(&my_data_set,infile); 
  for (i =0; i < my_data_set.nrows; i++){
    intensity[chip_num*my_data_set.nrows + i] = (double)(((short *)my_data_set.Data[0])[i]);
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  gzclose(infile);

  return(0);

}





void gzgeneric_get_masks_outliers_multichannel(const char *filename, int *nmasks, short **masks_x, short **masks_y, int *noutliers, short **outliers_x, short **outliers_y,  int channelindex){

  int i=0, k=0;
    
  gzFile infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;

  uint32_t next_group =1;  

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      
    }
  

  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);


  /* skip merrily through the file (optimise this with file pointers later) */
  while (k < channelindex){
    gzread_generic_data_group(&my_data_group,infile); 
    next_group = my_data_group.file_position_nextgroup; 
    gzseek(infile,next_group,SEEK_SET);
    Free_generic_data_group(&my_data_group);
    k++;
  }
  gzread_generic_data_group(&my_data_group,infile);

  /* passing the intensities */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* passing by the stddev */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  
  /* passing by the npixels */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* Now lets go for the "Outlier" */
  gzread_generic_data_set(&my_data_set,infile); 
   
  
  *noutliers = my_data_set.nrows;

  *outliers_x = Calloc(my_data_set.nrows,short); 
  *outliers_y = Calloc(my_data_set.nrows,short);
  
  gzread_generic_data_set_rows(&my_data_set,infile); 
  
  for (i=0; i < my_data_set.nrows; i++){
    (*outliers_x)[i] = ((short *)my_data_set.Data[0])[i];
    (*outliers_y)[i] = ((short *)my_data_set.Data[1])[i];
  }
  
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);

  /* Now lets go for the "Mask" */
  gzread_generic_data_set(&my_data_set,infile); 
   
  *nmasks = my_data_set.nrows;

  *masks_x = Calloc(my_data_set.nrows,short); 
  *masks_y = Calloc(my_data_set.nrows,short);
  
  
  gzread_generic_data_set_rows(&my_data_set,infile); 
  for (i=0; i < my_data_set.nrows; i++){
    (*outliers_x)[i] = ((short *)my_data_set.Data[0])[i];
    (*outliers_y)[i] = ((short *)my_data_set.Data[1])[i];
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);

  gzclose(infile);
}





void gzgeneric_apply_masks_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers,  int channelindex){


  int i=0;
  int cur_index;
  
  short cur_x, cur_y;

  int nrows;
  int size;

  gzFile infile;

  generic_file_header my_header;
  generic_data_header my_data_header;
  generic_data_group my_data_group;

  generic_data_set my_data_set;
  nvt_triplet *triplet;
  AffyMIMEtypes cur_mime_type;

  if ((infile = gzopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s\n",filename);
      
    }
 

  
  gzread_generic_file_header(&my_header, infile);
  gzread_generic_data_header(&my_data_header, infile);
  gzread_generic_data_group(&my_data_group,infile);

    
  triplet =  find_nvt(&my_data_header,"affymetrix-cel-rows");
  cur_mime_type = determine_MIMETYPE(*triplet);
  decode_MIME_value(*triplet,cur_mime_type, &nrows, &size);
  


  /* passing the intensities */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* passing by the stddev */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
  
  /* passing by the npixels */
  gzread_generic_data_set(&my_data_set,infile); 
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);
 
  /* Now lets go for the "Outlier" */
  gzread_generic_data_set(&my_data_set,infile); 
  
  if (rm_outliers){
    gzread_generic_data_set_rows(&my_data_set,infile); 
    for (i=0; i < my_data_set.nrows; i++){
      cur_x = ((short *)my_data_set.Data[0])[i];
      cur_y = ((short *)my_data_set.Data[1])[i];
      cur_index = (int)cur_x + nrows*(int)cur_y; 
      intensity[chip_num*rows + cur_index] =  R_NaN;
    }
  }
  
  gzseek(infile, my_data_set.file_pos_last, SEEK_SET); 
  Free_generic_data_set(&my_data_set);

  /* Now lets go for the "Mask" */
  gzread_generic_data_set(&my_data_set,infile); 
  if (rm_mask){
    gzread_generic_data_set_rows(&my_data_set,infile); 
    for (i=0; i < my_data_set.nrows; i++){
      cur_x = ((short *)my_data_set.Data[0])[i];
      cur_y = ((short *)my_data_set.Data[1])[i];
      cur_index = (int)cur_x + nrows*(int)cur_y; 
      intensity[chip_num*rows + cur_index] =  R_NaN;
    }
  }
  Free_generic_data_set(&my_data_set);
  Free_generic_data_header(&my_data_header);
  Free_generic_data_group(&my_data_group);
  
  gzclose(infile);
  
}



