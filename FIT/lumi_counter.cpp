#include "TDatime.h"
#include <fstream>
#include <iostream>

inline unsigned int lumi_time(int day, int month, int year, int hour=21, int min=0, int sec=0)
{
    TDatime d1(2000+year,month,day,hour,min,sec);
    return d1.Convert();
}


float lumi_counter(int input_day, int input_month, int input_year){

  double t = lumi_time(input_day, input_month, input_year);
  t *= 1000;

  double t0 = 0, t1 = 0;
  double l0 = 0, l1 = 0;

  std::ifstream input("ATLAS_lumi_compressed.txt");

  if(input){
    do{
      t0 = t1;
      l0 = l1;
      input >> t1 >> l1;
    }
    while(t > t1);
  }
  else{
    std::cout << "Can't open the file" << std::endl;
  }
#ifdef DEBUG
  std::cout << std::endl;
  std::cout << t0 << " " << l0 << std::endl;
  std::cout << t << std::endl;
  std::cout << t1 << " " << l1 << std::endl;
#endif
  return l0;

}
