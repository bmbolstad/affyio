#ifndef READ_MULTICHANNEL_CELFILE_GENERIC_H
#define READ_MULTICHANNEL_CELFILE_GENERIC_H

#include "read_abatch.h"

int isGenericMultiChannelCelFile(const char *filename);
int read_genericcel_file_intensities_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
int read_genericcel_file_stddev_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
int read_genericcel_file_npixels_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
void generic_get_masks_outliers_multichannel(const char *filename, int *nmasks, short **masks_x, short **masks_y, int *noutliers, short **outliers_x, short **outliers_y, int channelindex);
void generic_apply_masks_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers, int channelindex);
int multichannel_determine_number_channels(const char *filename);
char *multichannel_determine_channel_name(const char *filename, int channelindex);

int isgzGenericMultiChannelCelFile(const char *filename);
int gzread_genericcel_file_intensities_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
int gzread_genericcel_file_stddev_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
int gzread_genericcel_file_npixels_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int channelindex);
void gzgeneric_get_masks_outliers_multichannel(const char *filename, int *nmasks, short **masks_x, short **masks_y, int *noutliers, short **outliers_x, short **outliers_y, int channelindex);
void gzgeneric_apply_masks_multichannel(const char *filename, double *intensity, int chip_num, int rows, int cols,int chip_dim_rows, int rm_mask, int rm_outliers, int channelindex);
int gzmultichannel_determine_number_channels(const char *filename);
char *gzmultichannel_determine_channel_name(const char *filename, int channelindex);




#endif
