// g++ -std=c++17 -Wall -I$EIC_SHELL_PREFIX/include -L$EIC_SHELL_PREFIX/lib SimSiPM.cpp -o SimSiPM -lsipm
#include <sipm/SiPM.h>

using namespace std;
using namespace sipm;

int main()
{
  SiPMProperties myProperties;
  myProperties.setSize(6.);
  myProperties.setPitch(15.);
  SiPMSensor mySensor(myProperties);

  vector<double> times = {12.460000, 12.570000, 12.500000, 12.790000, 12.850000, 12.600000, 12.680000, 12.550000, 12.620000, 12.620000, 12.750000, 12.870000, 12.750000};
  mySensor.resetState();
  mySensor.addPhotons(times);
  mySensor.runEvent();

  SiPMAnalogSignal mySignal = mySensor.signal();
  for(int i=0; i<mySignal.size(); i++)
    cout << mySignal[i] << ", ";
  double integral = mySignal.integral(10, 250, 0.5);
  cout << "Integral = " << integral << endl;

  return 0;
}
