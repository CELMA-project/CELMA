#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_



#define  JUST_WRITE_INI_FILE  1 
#define  READ_INI_FILE  0 


#ifndef USAGE
#define USAGE {\
fprintf(stderr,"usage : %s {-F|-R|-B|-I} name number [-N name number] [-P ini-filename] [-H]  [-<name> <value>] \n",argv[0]);\
fprintf(stderr,"\t\t -I ini-filename (for start from .ini file)...(usually 3rd step)\n");\
fprintf(stderr,"\t\t -F normal start from a start file with questions asked (for compatibility with real old startups)\n");\
fprintf(stderr,"\t\t -R restart of a run from an existing data file\n");\
fprintf(stderr,"\t\t -B (batch) start without questions asked.\n");\
fprintf(stderr,"\t\t -H write a sample ini-file %s_sample.ini (usually 1st step\n",para->codename );\
fprintf(stderr,"\t\t -P to translate [P]hysical parameters in ini-file to numerical (ie. normalized) parameter file (usually 2nd step)\n\n");\
fprintf(stderr,"Look at sample inifile >%s_sample.ini< in this directory.\n\n",para->codename);\
fprintf(stderr,"\t\t-D Option runs directly with build in parameters, for simple Debugging\n\n");\
exit(EXIT_FAILURE);\
}

#endif

#ifdef __cplusplus
extern "C"
{
#endif

void FUtils_IniStructure(HDF_DS *d,PARA *p);


/* User functions */

int FUtils_ReadArguments(int ,char **,HDF_DS *,PARA *, void (*CalcParameters)( PARA *,HDF_DS *));

void FUtils_SetupParaStructure(HDF_DS *data, PARA *para);
void FUtils_IniStructure(HDF_DS *d,PARA *p);
int FUtils_AllocateReadMemory(HDF_DS *data);

void FUtils_Read3DFieldbyName(double ***f_0,const char *name,HDF_DS *data,PARA *para);
void FUtils_Read3DFieldbyNumber(double ***f_0,int number,HDF_DS *data,PARA *para);
void FUtils_Read2DFieldbyName(double **f_0,const char *name,HDF_DS *data,PARA *para);
void FUtils_Read2DFieldbyNumber(double **f_0,int number,HDF_DS *data,PARA *para);


void FUtils_Write2dSpace(double **feld,const char *a,const char *name,long num,HDF_DS *ds,PARA *p, int create);
void FUtils_Write3dSpace(double ***a,const char *fname,const char *name,long num,HDF_DS *ds,PARA *p,int create);

/*  Local Functions */ 
/* Should only be called from functions in libcommon  */

int   FUtilsInt_MakeHDFFilename(char *filename, const char *name,int32 number);
void  FUtilsInt_SetupDataRange(HDF_DS *data,PARA *para);

void FUtilsInt_ReadHDF4Attributes(const char *name,int32 number,HDF_DS *data, PARA *para);
int FUtilsInt_WriteHDF4Attributes(const char *name,int32 number,HDF_DS *data,PARA *para);


int   FUtilsInt_ReadFileInfo(const char *name,int32 number,HDF_DS *data,char *Dateiname);
void  FUtilsInt_WriteHDF4(const char *name,int32 number,HDF_DS *data,PARA *para,const char *feldname);
int   FUtilsInt_ReadHDF4ByName( char *name,int32 number,const char *sdsname, HDF_DS *data,PARA *para);
int   FUtilsInt_ReadHDF4ByNumber( char *name,int32 number,int numsds, HDF_DS *data,PARA *para);
int   FUtilsInt_ReadHDF4Coordinates(const char *name,int32 number, HDF_DS *data,PARA *para);
void  FUtilsInt_AddToRunDatabase(HDF_DS *d, PARA *p, char **argv);
#ifdef __cplusplus
}
#endif


#endif /* _FILE_UTILS_H_ */


