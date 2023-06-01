#include "Filter.h"

void Filter::Init(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc)
    {
        switch (Filter_type)
        {
        case FILTER_TYPE::ORDER_1_LOW_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = 0.0;
            b0 = 1;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ORDER_1_HIGH_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = tau;
            b0 = 0.0;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::APPROXIMATED_DIFFERENTIATOR:
            fc1 = set_fc;
            tau = 1.0 / (2.0 * pi * fc1);
            // alpha = 1 / (tau + (interval_t / 1000000.0));
            b2 = 0.0;
            b1 = 1.0;
            b0 = 0.0;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ORDER_2_LOW_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = 0.0;
            b0 = 1;
            a2 = tau * tau;
            a1 = 2 * tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ODER_2_HIGH_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = tau * tau;
            a1 = 4 * pi * fc1;
            a0 = 4.0 * pi * pi * fc1 * fc1;
            break;
        // case FILTER_TYPE::ORER_1_COMPLEMENTARY:

        //     break;
        case FILTER_TYPE::ODER_2_BUTTERWORTH_HPF:
            fc1 = set_fc;
            omega_n = 2 * pi * fc1;
            damping_ratio = 0.707;
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = 1;
            a1 = 2 * omega_n;
            a0 = omega_n * omega_n;
            break;
        default:
            error_flag = 1;
            break;
        }
        Compute_Z_Transform_Coefficient(Discrete_method);
}

void Filter::Init(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc, double set_dt)
    {
        switch (Filter_type)
        {
        case FILTER_TYPE::ORDER_1_LOW_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = 0.0;
            b0 = 1;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ORDER_1_HIGH_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = tau;
            b0 = 0.0;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::APPROXIMATED_DIFFERENTIATOR:
            fc1 = set_fc;
            tau = 1.0 / (2.0 * pi * fc1);
            // alpha = 1 / (tau + (interval_t / 1000000.0));
            b2 = 0.0;
            b1 = 1.0;
            b0 = 0.0;
            a2 = 0.0;
            a1 = tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ORDER_2_LOW_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 0.0;
            b1 = 0.0;
            b0 = 1;
            a2 = tau * tau;
            a1 = 2 * tau;
            a0 = 1.0;
            break;
        case FILTER_TYPE::ODER_2_HIGH_PASS:
            fc1 = set_fc;
            tau = 1 / (2 * pi * fc1);
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = tau * tau;
            a1 = 4 * pi * fc1;
            a0 = 4.0 * pi * pi * fc1 * fc1;
            break;
        // case FILTER_TYPE::ORER_1_COMPLEMENTARY:

        //     break;
        case FILTER_TYPE::ODER_2_BUTTERWORTH_HPF:
            fc1 = set_fc;
            omega_n = 2 * pi * fc1;
            damping_ratio = 0.707;
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = 1;
            a1 = 2 * omega_n;
            a0 = omega_n * omega_n;
            break;
        default:
            error_flag = 1;
            break;
        }
        Compute_Z_Transform_Coefficient(Discrete_method,set_dt);
}

void Filter::Init2fc(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc1, double set_fc2)
    {
        switch (Filter_type)
        {
        case FILTER_TYPE::ORDER_2_LOW_PASS:
            fc1 = set_fc1;
            fc2 = set_fc2;
            b2 = 0.0;
            b1 = 0.0;
            b0 = 4.0 * pi * pi * fc1 * fc2;
            a2 = 1.0;
            a1 = 2.0 * pi * (fc1 + fc2);
            a0 = 4.0 * pi * pi * fc1 * fc2;
            break;
        case FILTER_TYPE::ODER_2_HIGH_PASS:
            fc1 = set_fc1;
            fc2 = set_fc2;
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = 1.0;
            a1 = 2.0 * pi * (fc1 + fc2);
            a0 = 4.0 * pi * pi * fc1 * fc2;
            break;
        case FILTER_TYPE::ODER_2_NOTCH_FILTER:
            fc1 = set_fc1; // center of the filtering frequency
            fc2 = set_fc2; // Bandwidth of the filtering frequency
            b2 = 1.0;
            b1 = 0.0;
            b0 = 4 * pi * pi * fc1 * fc1;
            a2 = 1.0;
            a1 = 2.0 * pi * fc2;
            a0 = 4.0 * pi * pi * fc1 * fc1;
            b2 = 0.0;
            b1 = 2 * damping_ratio * omega_n;
            b0 = 0.0;
            a2 = 1.0;
            a1 = 2 * damping_ratio * omega_n;
            a0 = omega_n * omega_n;
            break;
        case FILTER_TYPE::ODER_2_BANDPASS_FILTER:
            fc1 = set_fc1;
            fc2 = set_fc2;
            omega_n = 2 * pi * fc1;
            damping_ratio = fc2 / (2 * fc1);
            break;
        default:
            error_flag = 1;
            break;
        }
        Compute_Z_Transform_Coefficient(Discrete_method);
};

void Filter::Init2fc(FILTER_TYPE Filter_type, DISCRETE_METHOD Discrete_method, double set_fc1, double set_fc2, double set_dt)
    {
        switch (Filter_type)
        {
        case FILTER_TYPE::ORDER_2_LOW_PASS:
            fc1 = set_fc1;
            fc2 = set_fc2;
            b2 = 0.0;
            b1 = 0.0;
            b0 = 4.0 * pi * pi * fc1 * fc2;
            a2 = 1.0;
            a1 = 2.0 * pi * (fc1 + fc2);
            a0 = 4.0 * pi * pi * fc1 * fc2;
            break;
        case FILTER_TYPE::ODER_2_HIGH_PASS:
            fc1 = set_fc1;
            fc2 = set_fc2;
            b2 = 1.0;
            b1 = 0.0;
            b0 = 0.0;
            a2 = 1.0;
            a1 = 2.0 * pi * (fc1 + fc2);
            a0 = 4.0 * pi * pi * fc1 * fc2;
            break;
        case FILTER_TYPE::ODER_2_NOTCH_FILTER:
            fc1 = set_fc1; // center of the filtering frequency
            fc2 = set_fc2; // Bandwidth of the filtering frequency
            b2 = 1.0;
            b1 = 0.0;
            b0 = 4 * pi * pi * fc1 * fc1;
            a2 = 1.0;
            a1 = 2.0 * pi * fc2;
            a0 = 4.0 * pi * pi * fc1 * fc1;
            b2 = 0.0;
            b1 = 2 * damping_ratio * omega_n;
            b0 = 0.0;
            a2 = 1.0;
            a1 = 2 * damping_ratio * omega_n;
            a0 = omega_n * omega_n;
            break;
        case FILTER_TYPE::ODER_2_BANDPASS_FILTER:
            fc1 = set_fc1;
            fc2 = set_fc2;
            omega_n = 2 * pi * fc1;
            damping_ratio = fc2 / (2 * fc1);
            break;
        default:
            error_flag = 1;
            break;
        }
        Compute_Z_Transform_Coefficient(Discrete_method,set_dt);
};

void Filter::Reset()
{
    sensor_data[0]=0.0;
    sensor_data[1]=0.0;
    sensor_data[2]=0.0;
}


void Filter::Update_Filter(double present_sensor_data)
{
    sensor_data[2] = sensor_data[1];
    sensor_data[1] = sensor_data[0];
    sensor_data[0] = present_sensor_data;
    filtered_data[2] = filtered_data[1];
    filtered_data[1] = filtered_data[0];
    filtered_data[0] = -alpha2 * filtered_data[2] - alpha1 * filtered_data[1] + beta2 * sensor_data[2] + beta1 * sensor_data[1] + beta0 * sensor_data[0];
    // filtered_data[0] = alpha * (sensor_data[0]-sensor_data[1] + tau * filtered_data[1]);
}

void Filter::Compute_Z_Transform_Coefficient(DISCRETE_METHOD Discrete_method)
{
    alpha1 = 0.0;
    alpha2 = 0.0;
    beta0 = 0.0;
    beta1 = 0.0;
    beta2 = 0.0;
    double T = interval_t / 1000000.0;
    switch (Discrete_method)
    {
    case DISCRETE_METHOD::BACKWARD_METHOD:
        alpha1 = -(2 * a2 + T * a1) / (a0 * T * T + a1 * T + a2);
        alpha2 = a2 / (a0 * T * T + a1 * T + a2);
        beta0 = (b0 * T * T + b1 * T + b2) / (a0 * T * T + a1 * T + a2);
        beta1 = -(2 * b2 + T * b1) / (a0 * T * T + a1 * T + a2);
        beta2 = b2 / (a0 * T * T + a1 * T + a2);
        break;
    case DISCRETE_METHOD::BILINEAR_METHOD:
        alpha1 = -(-2 * a0 * T * T + 8 * a2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        alpha2 = (a0 * T * T - 2 * a1 * T + 4 * a2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta0 = (b0 * T * T + 2 * b1 * T + 4 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta1 = -(-2 * b0 * T * T + 8 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta2 = (b0 * T * T - 2 * b1 * T + 4 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        break;
    default:
        error_flag = 1;
        break;
    }
}

void Filter::Compute_Z_Transform_Coefficient(DISCRETE_METHOD Discrete_method,double set_dt)
{
    alpha1 = 0.0;
    alpha2 = 0.0;
    beta0 = 0.0;
    beta1 = 0.0;
    beta2 = 0.0;
    double T = set_dt;
    switch (Discrete_method)
    {
    case DISCRETE_METHOD::BACKWARD_METHOD:
        alpha1 = -(2 * a2 + T * a1) / (a0 * T * T + a1 * T + a2);
        alpha2 = a2 / (a0 * T * T + a1 * T + a2);
        beta0 = (b0 * T * T + b1 * T + b2) / (a0 * T * T + a1 * T + a2);
        beta1 = -(2 * b2 + T * b1) / (a0 * T * T + a1 * T + a2);
        beta2 = b2 / (a0 * T * T + a1 * T + a2);
        break;
    case DISCRETE_METHOD::BILINEAR_METHOD:
        alpha1 = -(-2 * a0 * T * T + 8 * a2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        alpha2 = (a0 * T * T - 2 * a1 * T + 4 * a2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta0 = (b0 * T * T + 2 * b1 * T + 4 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta1 = -(-2 * b0 * T * T + 8 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        beta2 = (b0 * T * T - 2 * b1 * T + 4 * b2) / (a0 * T * T + 2 * a1 * T + 4 * a2);
        break;
    default:
        error_flag = 1;
        break;
    }
}

void Command_Pre_Filter::Init(double fc_set,int n_set){
    fc = fc_set;
    n=n_set;
    tau = 1 / (2 * pi * fc);
    ptr_A_coefficient = new double[n+1];
    ptr_x = new double[n];
    ptr_x_pre = new double[n]; 
    double tau_pow_temp =1.0;

    for (int i = 0; i < n; i++)
    {
        if(i==0){
            ptr_A_coefficient[0]=1.0;
        }else{
            ptr_A_coefficient[i] = ptr_A_coefficient[i-1]*((double)((double)(n+1-i)/(double)i));
        }
        ptr_x[i]=0.0;
        ptr_x_pre[i]=0.0;
    }
    for (int i = 1; i <= n; i++)
    {
        tau_pow_temp=tau_pow_temp/tau;
        ptr_A_coefficient[n-i] = ptr_A_coefficient[n-i]*tau_pow_temp;
    }
    ptr_A_coefficient[n]=tau_pow_temp;
}

void Command_Pre_Filter::Reset()
{
    for (int i = 0; i < (n-1); i++)
    {
        ptr_x[i] = 0.0;
        ptr_x_pre[i] = 0.0;
    }
}

void Command_Pre_Filter::Reset(double x_0)
{
    ptr_x[0] = x_0;
    ptr_x_pre[0] = x_0;
    for (int i = 1; i < (n-1); i++)
    {
        ptr_x[i] = 0.0;
        ptr_x_pre[i] = 0.0;
    }
}

void Command_Pre_Filter::Update_Filter(double present_sensor_data){
    for (int i = 0; i < (n-1); i++)
    {
        ptr_x[i] = ptr_x_pre[i+1]*T+ptr_x_pre[i];
    }
    
    ptr_x[n-1] =0.0;
    for (int i = 0; i < n; i++)
    {
        ptr_x[n-1] = ptr_x[n-1] - ptr_A_coefficient[i]*ptr_x_pre[i];
    }
    ptr_x[n-1] = ptr_x_pre[n-1] + (ptr_x[n-1]+ ptr_A_coefficient[n]*present_sensor_data)*T;

    for (int i = 0; i < n; i++)
    {
        ptr_x_pre[i] = ptr_x[i];
    }
}

void Moving_Average_Filter::Init(int n_set){
    n = n_set;
    ptr_x = new double[n];
    for (int i = 0; i < n; i++)
    {
        ptr_x[i]=0.0;
    }
    filtered_data[0] = 0.0;
    filtered_data[1] = 0.0;
}

void Moving_Average_Filter::Update_Filter(double present_sensor_data){
    filtered_data[0] = filtered_data[1] + (present_sensor_data-ptr_x[n-1])/n;
    for (int i = n-1; i > 0; i--)
    {
        ptr_x[i]=ptr_x[i-1];
    }
    filtered_data[1]=filtered_data[0];
    ptr_x[0]=present_sensor_data;
}

void Moving_Average_Filter::Reset()
{
    for (int i = 0; i < n; i++)
    {
        ptr_x[i]=0.0;
    }
    filtered_data[0] = 0.0;
    filtered_data[1] = 0.0;
}
