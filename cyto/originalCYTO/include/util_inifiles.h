/*****************************************/
/* Reading from files                    */
    
    
double Util_Str2Double(char *str);
double Util_ReadValueByName(const char *name, char *file, double def,double hadval,const char *desc,int isfloat,int rw);
int Util_ReadIniFile(HDF_DS *data,PARA *para,char *ininame,int rw);
char* Util_GetDesc(const char *full,const char *search,const char *def,char *result);
  
