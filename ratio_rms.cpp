#include <cmath>
#define sq(a,b) sqrt(a*a + b*b)
//extern double RMS_errors[48][13];
float ratio_rms(int part, int ip, int row1, int row2)  {
       float RMS_errors[4][48][13] = {
          #include "LBA.rms"
          #include "LBC.rms"
          #include "EBA.rms"
          #include "EBC.rms"
       };
/*       cout << "TEST: Output for errors before filling for PMT " << ip << ": " << "and for partition " << partition << endl;
       for(int row=0; row<13; row++){
          cout << RMS_errors[ip][row] << ", ";
       }  
     cout << endl;*/
       cout << "Result error for PMT " << ip+1 << ", for row " << row1+1 << " and for row "<< row2+1 << " is "  << sq(RMS_errors[part][ip][row2], RMS_errors[part][ip][row2]) << endl;
       return sq(RMS_errors[part][ip][row2], RMS_errors[part][ip][row2]);
       //return result;
}
