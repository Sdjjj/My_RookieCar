#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <cmath>
using namespace std;


class QuinticPolynomial
{
public:
    double t;

    //初始位置参数
    double xp0, xv0, xa0;
    double yp0, yv0, ya0;
    
    //目标点位置参数
    double xp1, xv1, xa1;
    double yp1, yv1, ya1;

    //x方向上多项式系数
    double m0, m1, m2, m3, m4, m5;

    //y方向上多项式系数
    double n0, n1, n2, n3, n4, n5;



    //计算横向距离x关于时间t的五次多项式的系数
    QuinticPolynomial cal_coefficient_x(double xp0_input, double xv0_input, double xa0_input,
    double xp1_input, double xv1_input, double xa1_input, double t_input)
    {
        xp0 = xp0_input;
        xv0 = xv0_input;
        xa0 = xa0_input;

        xp1 = xp1_input;
        xv1 = xv1_input;
        xa1 = xa1_input;

        t = t_input;

        //定义矩阵
        Eigen::Matrix<double, 6, 6> A;
        Eigen::VectorXd X(6);
        Eigen::VectorXd B(6);

        A << 1, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 2, 0 , 0, 0 ,
              1, t, pow(t, 2), pow(t, 3), pow(t, 4), pow(t, 5),
              0, 1, 2*t, 3*pow(t, 2), 4*pow(t, 3), 5*pow(t, 4),
              0, 0, 2, 6*t, 12*pow(t, 2), 20*pow(t, 3);
        
        B << xp0, xv0, xa0, xp1, xv1, xa1;

        X = A.colPivHouseholderQr().solve(B);
        m0 = xp0;
        m1 = xv0;
        m2 = xa0/2;
        m3 = X[3];
        m4 = X[4];
        m5 = X[5];
    }

    double direction_x_distance(double t)
    {
        return m0+m1*t+m2*pow(t, 2)+m3*pow(t, 3)+m4*pow(t, 4)+m5*pow(t, 5);
    }

    double direction_x_velocity(double t)
    {
        return m1+2*m2*t+3*m3*pow(t, 2)+4*m4*pow(t, 3)+5*m5*pow(t, 4);
    }

    double direction_x_acceleration(double t)
    {
        return 2*m2+6*m3*t+12*m4*pow(t, 2)+20*m5*pow(t, 3);
    }

    //计算纵向距离y关于时间t的五次多项式的系数
    QuinticPolynomial cal_coefficient_y(double yp0_input, double yv0_input, double ya0_input,
    double yp1_input, double yv1_input, double ya1_input, double t_input)
    {
        yp0 = yp0_input;
        yv0 = yv0_input;
        ya0 = ya0_input;

        yp1 = yp1_input;
        yv1 = yv1_input;
        ya1 = ya1_input;

        t = t_input;
        
        Eigen::Matrix<double, 6, 6> C;
        Eigen::VectorXd Y(6);
        Eigen::VectorXd D(6);

        C << 1, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 0, 0,
             0, 0, 2, 0, 0, 0 ,
             1, t, pow(t, 2), pow(t, 3), pow(t, 4), pow(t, 5),
             0, 1, 2*t, 3*pow(t, 2), 4*pow(t, 3), 5*pow(t, 4),
             0, 0, 2, 6*t, 12*pow(t, 2), 20*pow(t, 3);

        D << yp0, yv0, ya0, yp1, yv1, ya1;   

        Y = C.colPivHouseholderQr().solve(D);

        n0 = yp0;
        n1 = yv0;
        n2 = ya0/2;
        n3 = Y[3];
        n4 = Y[4];
        n5 = Y[5];
    }

    double direction_y_distance(double t)
    {
        return n0+n1*t+n2*pow(t, 2)+n3*pow(t, 3)+n4*pow(t, 4)+n5*pow(t, 5);
    }

    double direction_y_velocity(double t)
    {
        return n1+2*n2*t+3*n3*pow(t, 2)+4*n4*pow(t, 3)+5*n5*pow(t, 4);
    }

    double direction_y_acceleration(double t)
    {
        return 2*n2+6*n3*t+12*n4*pow(t, 2)+20*n5*pow(t, 3);
    }

    double cal_thetar(double direction_x_velocity, double direction_y_velocity)
    {
        return atan2(direction_y_velocity, direction_x_velocity);
    }

    double cal_kappar(double direction_x_velocity, double direction_x_acceleration, 
    double direction_y_velocity, double direction_y_acceleration)
    {
        double input1 = direction_y_velocity*direction_x_acceleration - direction_x_velocity*direction_y_acceleration;
        double molecule = fabs(input1);

        double input2 = pow(pow(direction_x_velocity, 2) + pow(direction_y_velocity, 2), 3);
        double denominator = sqrt(input2);

        return molecule/denominator;
    }
};
