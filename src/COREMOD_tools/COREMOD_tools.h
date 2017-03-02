#ifndef _TOOLS_H
#define _TOOLS_H



int init_COREMOD_tools();


int COREMOD_TOOLS_mvProcCPUset(const char *csetname);

int create_counter_file(const char *fname, long NBpts);

int bubble_sort(double *array, long count);

void qs_float(float *array, long left, long right);
void qs_long(long *array, long left, long right);
void qs_double(double *array, long left, long right);
void qs_ushort(unsigned short *array, long left, long right);

void quick_sort_float(float *array, long count);
void quick_sort_long(long *array, long count);
void quick_sort_double(double *array, long count);
void quick_sort_ushort(unsigned short *array, long count);

void qs3(double *array, double *array1, double *array2, long left, long right);

void qs3_double(double *array, double *array1, double *array2, long left, long right);

void quick_sort3(double *array, double *array1, double *array2, long count);
void quick_sort3_float(float *array, float *array1, float *array2, long count);
void quick_sort3_double(double *array, double *array1, double *array2, long count);

void qs2l(double *array, long *array1, long left, long right);

void quick_sort2l(double *array, long *array1, long count);

void quick_sort2l_double(double *array, long *array1, long count);

void quick_sort3ll_double(double *array, long *array1, long *array2, long count);

int lin_regress(double *a, double *b, double *Xi2, double *x, double *y, double *sig, int nb_points);

int replace_char(char *content, char cin, char cout);

int read_config_parameter_exists(const char *config_file, const char *keyword);

int read_config_parameter(const char *config_file, const char *keyword, char *content);

float read_config_parameter_float(const char *config_file, const char *keyword);

long read_config_parameter_long(const char *config_file, const char *keyword);

int read_config_parameter_int(const char *config_file, const char *keyword);

long file_number_lines(const char *file_name);

FILE* open_file_w(const char *filename);

FILE* open_file_r(const char *filename);

int write_1D_array(double *array, long nbpoints, const char *filename);

int read_1D_array(double *array, long nbpoints, const char *filename);

int tp(const char *word);

int read_int_file(const char *fname);

int write_int_file(const char *fname, int value);

int write_float_file(const char *fname, float value);


int COREMOD_TOOLS_imgdisplay3D(const char *IDname, long step);

long COREMOD_TOOLS_statusStat(const char *IDstat_name, long indexmax);


#endif








































