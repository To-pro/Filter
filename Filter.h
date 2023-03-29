#ifndef FILTER_H
#define FILTER_H

#ifndef interval_t
#define interval_t 1000 //間隔時間(us)
#endif

#ifndef pi
#define pi 3.141592653589793
#endif
enum FILTER_TYPE
{
    ORDER_1_LOW_PASS,
    ORDER_1_HIGH_PASS,
    APPROXIMATED_DIFFERENTIATOR,
    ORDER_2_LOW_PASS,
    ODER_2_HIGH_PASS,
    ORER_1_COMPLEMENTARY,
    ODER_2_BUTTERWORTH_HPF,
    ODER_2_NOTCH_FILTER,
    ODER_2_BANDPASS_FILTER,
};
enum DISCRETE_METHOD
{
    BACKWARD_METHOD,
    BILINEAR_METHOD
};
class Filter
{
public:

    void Init(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc);
    void Init(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc1, double set_fc2);

    void Init(double fc_set);
    void reset();
    void Update_Filter(double present_sensor_data);
    double Get_Filtered_Data(void){
        return filtered_data[0];
    };

    void Get_Z_Domain_coefficient(double *coefficient);

protected:
    double alpha = 0.0;
    double tau = 0.0;
    double fc1 = 1.0;
    double fc2 = 1.0;

    double filtered_data[3] = {0.0};
    double sensor_data[3] = {0.0};

    int error_flag = 0; // 1:表示有問題(輸入不正確)

    // filter's z transfer function coefficient
    double alpha1, alpha2;
    double beta0, beta1, beta2;

    // filter's s transfer function coefficient
    double a0, a1, a2;
    double b0, b1, b2;
    void Compute_Z_Transform_Coefficient(DISCRETE_METHOD Discrete_method);

    // some filter's transfer function use
    double omega_n = 0.0;
    double damping_ratio = 0.0;

};

class Command_Pre_Filter
{
public:
    void Init(double fc_set,int n_set);
    // void Init(double fc_set,int n_set,double x0);
    void Update_Filter(double present_sensor_data);
    void reset();
    double* Get_command( ){
        return ptr_x;
    };
    double* Get_A_coefficient(){
        return ptr_A_coefficient;
    };
    
protected:
    double T = (double)interval_t / 1000000.0;
    int n = 1;
    double fc = 0.0;
    double tau = 0.0;
    double *ptr_A_coefficient;
    double *ptr_x;
    double *ptr_x_pre;
};

class Moving_Average_Filter
{
public:
    void Init(int n_set);
    void Update_Filter(double present_sensor_data);
    double Get_Filtered_Data(void){
        return filtered_data[0];
    };
protected:
    int n=0;
    double filtered_data[2]={0.0};
    double *ptr_x;
};

#endif
