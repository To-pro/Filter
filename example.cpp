#include <Arduino.h>
#define interval_t 1000
// #define PI 3.141592653589793


// #include "my_filter.h"
#include <Filter.h>

Filter LPF;
Filter HPF;
Filter LPF2;
Filter APD;
Command_Pre_Filter CPF;
Moving_Average_Filter MAF;

IntervalTimer Controller_Timer;

//variable
long int t_pre = 0;
long int t_p = 0;

double t = 0.0;
int t_idx=0;
double frequency1=3;
double frequency2=30;
double frequency3=150;
double frequency4=300;
double omega1=2*PI*frequency1;
double omega2=2*PI*frequency2;
double omega3=2*PI*frequency3;
double omega4=2*PI*frequency4;
double x0=sin(omega1*t);

double x=5*sin(omega1*t)+1*sin(omega2*t)+2*sin(omega3*t)+1*sin(omega4*t);
double y=0.0;
double u=0;




void ISR_Routing(){
  t_pre = micros();
  t_idx++;
  t=(double)t_idx*interval_t/1000000.0;
  double x0=sin(omega1*t);
  double dx0 = cos(omega1*t)*omega1;
  double x=5*sin(omega1*t)+1*sin(omega2*t)+2*sin(omega3*t)+1*sin(omega4*t);
  
  // double x0=0.001*t;
  // double dx0 = 0.001;

  // LPF2.Update_Filter(x);
  // HPF.Update_Filter(x);
  // y = LPF2.Get_Filtered_Data();
  // y = HPF.Get_Filtered_Data();
  MAF.Update_Filter(x);
  y = MAF.Get_Filtered_Data();

  // APD.Update_Filter(x0);
  // y = APD.Get_Filtered_Data();

  Serial.print(x,8);
  Serial.print("  ");
  Serial.print(y,8);
  // Serial.print("  ");
  // Serial.print(dx0,8);
  // Serial.print("  ");
  // Serial.print(t_p);
  // Serial.print("  ");
  Serial.print("\n");

  // t_p = micros()-t_pre;

  // if(t<=5){
  //   u=1;
  // }  
  // else if(t<=10){
  //   u=0;
  //   // Controller_Timer.end();
  // }
  // else{
  //   t_idx=0;
  // } 
  // CPF.Update_Filter(u);
  // double *y = CPF.Get_command();
  // Serial.print(u,8);
  // Serial.print("  ");
  // Serial.print(y[0],8);
  // Serial.print("  ");
  // Serial.print(y[1],8);
  // Serial.print("  ");
  // Serial.print(y[2],8);
  // Serial.print("  ");
  // Serial.print(y[3],8);
  // Serial.print("  ");
  // Serial.print(y[4],8);
  // Serial.print("\n");
}
void setup() {
  Serial.begin(115200);
  delay(1000);
  Serial.print("Setup\n");

  LPF.Init(ORDER_1_LOW_PASS,BACKWARD_METHOD,30);
  HPF.Init(ORDER_1_HIGH_PASS,BACKWARD_METHOD,300);

  LPF2.Init(ORDER_2_LOW_PASS,BACKWARD_METHOD,30,50);
  APD.Init(APPROXIMATED_DIFFERENTIATOR,BACKWARD_METHOD,300);

  MAF.Init(10);

  int n = 5;
  CPF.Init(0.8,n);
  // double *A_coefficient=CPF.Get_A_coefficient();
  // for (int i = 0;i<(n+1);i++)
  // {
  //   Serial.print(A_coefficient[i],4);
  //   Serial.print("  ");
  // }
  // Serial.print("\n");
  // while(1){};

  Controller_Timer.begin(ISR_Routing,interval_t);
  Controller_Timer.priority(255); 
}

void loop() {
  
}