/* HDF -Helper Functions                   
 * These are some small functions to allow easier access to HDF SD data-sets 
 * 29.08.1997 VN
 */

#include <assert.h>
#include <hdf.h> 
#include <df.h>
#include <mfhdf.h>
#include <unistd.h> 
#include <string.h>
#include <hdf_helper.h>

#undef DEBUG 



#ifdef DEBUG
#define COMM(a) a
#define BUGREPORT fprintf(stderr,"%s compiled %s %s, Line %i\n",__FILE__,__DATE__,__TIME__,__LINE__);
#else
#define COMM(a) 
#define BUGREPORT
#endif

/* Reads lcount double values with name name into field result
 * Nonexisting values default to default and are explicitly set
 * 
 * Input Parameters: SDS id id
 * name of varaible to read
 * pointer to hold result on output
 * lcount on input as number of variables in result
 * default is the user supplied default value
 */



void SDread_variable(int32 id,const char *name,double *result,int32 lcount,double def)
{
int32 indexl,num_type,count,j;
char tmpstr[8048];
char namel[256];
char buffer_c[8048];
void *buffer;


sprintf(namel,"%s",name);
buffer = &buffer_c[0];
for(j=0;j<lcount;j++) result[j] = def;


 COMM(fprintf(stderr,"%s %d:Try variable >%s< with %d realizations\n",__FILE__,__LINE__,name,(int)lcount);)
indexl = SDfindattr(id,name);
if(indexl < 0) 
  {
    COMM(fprintf(stderr,"%s %d: Variable %s not found. Using default >%f<.\n",__FILE__,__LINE__,name,def);)
    return;
  }

assert(SDattrinfo(id,indexl,tmpstr,&num_type,&count)+1);

if(count*sizeof(num_type) > 8000*sizeof(char))
 {
   fprintf(stderr,"%s %d: Variable %s has to many (%d ) realisations. Using default >%f<.\n",__FILE__,__LINE__,name,(int)count,def);
    return;
  }

if(count > lcount)
 {
   fprintf(stderr,"Warning: %s %d: Variable %s has more (%d) realisations than asked for>%d<.\n",__FILE__,__LINE__,name,(int)count,(int)lcount);
  }


assert(SDreadattr(id,indexl,(VOID*)buffer)+1);

BUGREPORT;

/* Copy lcount of elements over */
for(j=0;j<MIN(lcount,count);j++)
  {
    switch (num_type) 
      {
          case DFNT_FLOAT32:
              result[j] = (double) (((float*)buffer)[j]);
              break;
          case DFNT_FLOAT64:
              result[j] = (double) (((double*)buffer)[j]);
              break;
          case DFNT_INT8:
              result[j] = (double) (((int8*)buffer)[j]);
              break;
          case DFNT_INT16:
              result[j] = (double) (((int16*)buffer)[j]);
              break;
          case DFNT_INT32:
              result[j] = (double) (((int32*)buffer)[j]);
              break;    
      }
  }

 COMM(for(j=0;j<MIN(lcount,count);j++) fprintf(stderr,"%s %d: Variable %s[%d] = %f\n",__FILE__,__LINE__,name,(int)j,result[j]););
}
 
/*********************************************************************/

/* Reads a string with name name into field result
 * If string name is not found the string default is returned 
 *
 * Input Parameters: SDS id id
 * name of string to read
 * pointer to hold result on output
 * lcount holds length of array result 
 * contains default string to be copied into result
 *
 */


void SDread_string(int32 id,const char *name,char *result,int32 lcount,const char *def)
{
#define STRLENMAX 8000
int32 indexl,num_type,count;
char tmpstr[STRLENMAX+1];
char buffer_c[STRLENMAX+1];
void *buffer;


buffer = (void*) &buffer_c[0];

strncpy(tmpstr,name,STRLENMAX);
strncpy(result,def,lcount);

/* Check if name exists in HDF file */
indexl = SDfindattr(id,tmpstr);
if(indexl < 0) 
  {
      COMM(fprintf(stderr,"%s %d: %s not found, use default: >%s<\n",__FILE__,__LINE__,name,def););     
      return;
  }


/* Get info about variable */
assert(SDattrinfo(id,indexl,tmpstr,&num_type,&count)+1);

 if(count<STRLENMAX)
   {
   assert(SDreadattr(id,indexl,(VOID*)buffer)+1);
   
   ((char*)buffer)[count]='\0';
   if(num_type == DFNT_CHAR8)   
   {
       strncpy(result, ((char*)buffer),lcount);
       result[MIN(count,lcount)]='\0'; 
   }
   
   else    COMM(fprintf(stderr,"%s %d: %s is not DFNT_CHAR8 but %d, use default\n",__FILE__,__LINE__,name,(int)num_type););
   }
 else  COMM(fprintf(stderr,"%s %d: %s is to long to read: %d, use default %s\n",__FILE__,__LINE__,name,(int)count,def););

 COMM(fprintf(stderr, "%s %d:  Read %s: >%s< with lcount = %d\n",__FILE__,__LINE__,name,result,(int)lcount);)

}
