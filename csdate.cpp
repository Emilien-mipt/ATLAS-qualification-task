#include "TDatime.h"

// function which returns starting day on the plot
inline TDatime day0() 
{
    return TDatime(2009,6,1,21,00,00);  // starting date on the plot
}

// function which returns number of seconds since starting day on the plot
inline Double_t cstime(int day, int month, int year, int hour=21, int min=0, int sec=0)
{
    TDatime d1(2000+year,month,day,hour,min,sec); 
    return d1.Convert()-day0().Convert();
}
