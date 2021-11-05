#include "quintcPolynomial.h"
#include "quinticType.h"
#include <iostream>
#include <cmath>
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include <vector>
#include <tf/tf.h>
using namespace std;

/**************               定义整车参数        ***************/
//轮胎侧偏刚度
double cf = -65494.663, cr = -115494.663;

//前后悬架载荷
double mass_fl = 500, mass_fr = 500, mass_rl = 520, mass_rr = 520;
double mass_front = mass_fl + mass_fr;
double mass_rear = mass_rl + mass_rr;
double m = mass_front + mass_rear;

//轴距
double L = 3.8; 
//前轴中心到质心的距离
double a = L * (1.0 - mass_front / m);
//后轴中心到质心的距离
double b = L * (1.0 - mass_rear / m);
//车辆绕z轴转动的转动惯量
double Iz = std::pow(a, 2) * mass_front + std::pow(b, 2) * mass_rear;
//轮胎最大转角(rad)
double wheel_max_degree = 0.6;


//采样时间
#define dt 0.1

// t-t0经历的时间
double T = 50;

double xend = 50.0;
double yend = 10.0;



//预瞄时间
#define ts 0.5

// 从ros中暂时读取的信息
//纵向车速vx（质心侧偏角很小，故纵向车速vx约等于质心速度v.vx应从传感器处获得，此处为防止分母为0 ，故设初值）
double vx = 0.01;
//横向车速（质心侧偏角很小，故纵向车速vy近似等于0）
double vy = 0;

//定义publisher的全局变量，便于在回调函数中更新参数
ros::Publisher LQR; //发布转角及线速度
ros::Publisher car_trajectory_pub; //小车的运动轨迹

nav_msgs::Path car_trajectory_msgs; //用于发布小车运动轨迹存储信息对象





/***********   五次多项式：计算x和y方向上的距离，速度，加速度      *******************/
quinticType pathpoint;  //由于只计算一次路径规划，故设为全局变量
void cal_path_number()
{   
    //quinticType path_temp; //需要多次计算路径规划时，设为局部变量

    QuinticPolynomial QPx;
    QPx.cal_coefficient_x(0.0, 0.0, 0.0, xend,
                           0.0, 0.0, T);
    //给出x方向上的初始和目标点横向状态参数：x0, vx0, ax0, x1,vx1, ax1
    //QPx.cal_coefficient_x(0, 0, 0, 10, 0, 0, time);

    QuinticPolynomial QPy;
    //给出y方向上的初始和目标点横向状态参数：y0, vy0, ay0, y1,vy1, ay1
    QPy.cal_coefficient_y(0.0, 0.0, 0.0, yend,
                           0.0, 0.0, T, xend);

    for(double t = 0; t < T; t += dt)
    {
        double x = QPx.calc_point_x(t);
        double xd = QPx.calc_point_xd(t);
        double xdd = QPx.calc_point_xdd(t);

        pathpoint.x.push_back(x);
        pathpoint.vx.push_back(xd);
        pathpoint.ax.push_back(xdd);

        double y_x_t = QPy.calc_point_y_x(x);
        double y_x_d = QPy.calc_point_y_x_d(x);
        double y_x_t_d = QPy.calc_point_y_t_d(y_x_d, xd);

        double y_x_dd = QPy.calc_point_y_x_dd(x);
        double y_x_t_dd = QPy.calc_point_y_t_dd(y_x_dd, xd, y_x_d, xdd);

        pathpoint.y.push_back(y_x_t);
        pathpoint.vy.push_back(y_x_t_d);
        pathpoint.ay.push_back(y_x_t_dd);

        // 压入曲率
        pathpoint.kappar.push_back(QPy.calc_point_k(y_x_dd, y_x_d));

    }

    int num = pathpoint.x.size();  
    for (int i = 0; i < num; i++) 
    {
    double dy = pathpoint.y[i + 1] - pathpoint.y[i];
    double dx = pathpoint.x[i + 1] - pathpoint.x[i];
    pathpoint.thetar.push_back(QPy.calc_point_thetar(dy, dx));
    }

}

/************        LQR 计算模块      ******************/
Eigen::Vector4d cal_LQR(double vx)
{
    Eigen::Vector4d K;
    Eigen::Matrix<double, 4, 4> A;
    Eigen::Vector4d B;

    A << 0, 1, 0, 0,
         0, (cf + cr) / (m * vx), -(cf + cr) / m, (a * cf - b * cr)/(m * vx), 
         0, 0, 0, 1, 
         0, (a * cf - b * cr) / (Iz * vx), -(a * cf - b * cr) / Iz, (a * a * cf + b * b * cr) / (Iz * vx);

    B << 0, -cf / m, 0, -a * cf / Iz;

    //将矩阵A 和 B 离散化       
    double td = 0.001; 
    Eigen::Matrix4d E;
    E.setIdentity(4, 4);

    Eigen::Matrix4d Ad; 
    Ad = (E - td * 0.5 * A).inverse() * (E + td * 0.5 * A);

    Eigen::Vector4d Bd;
    Bd = B * td;

    // Q为半正定的状态加权矩阵，为对角阵，Q元素变大意味着希望跟踪偏差能够快速趋于零
    Eigen::Vector4d diagonal_Q(20, 1, 1, 1);
    Eigen::Matrix4d Q = diagonal_Q.asDiagonal();
    
    //R为正定的控制加权矩阵，也为对角阵，R元素变大意味着希望控制输入能够尽可能小
    Eigen::Matrix<double, 1, 1> R;
    R << 15.0;

    //设置最大循环次数
    int MaxLoop = 500;

    //设置P_old 与P_new 之差的阈值
    double threshold = 10e-10;

    //考虑当vx=0时候的奇异性
    if (vx = 0)
    {
        K << 0, 0, 0, 0;
    }
    
    else
    {
        //迭代求解黎卡提方程
        Eigen::Matrix4d P_old = Q;

        for (size_t i = 0; i < MaxLoop; i++)
        {
            Eigen::Matrix4d P_new = Ad.transpose() * P_old * Ad - Ad.transpose() * P_old * Bd * 
                                    (R + Bd.transpose() * P_old * Bd).inverse() *  Bd.transpose() * P_old * Ad + Q;

            // apollo 计算范数的方式
            if (fabs((P_new - P_old).maxCoeff()) < threshold)
            {
                P_old = P_new;
                break;
            }
            
            P_old = P_new;
        }

        //求出反馈控制增益矩阵K
        K = (R + Bd.transpose() * P_old * Bd).inverse() * Bd.transpose() * P_old * Ad;

    }

    return K;
} 

/************          计算err,投影点曲率kr         ******************/
Eigen::VectorXd cal_err_kr(double x_now, double y_now, Eigen::Vector3d eular)
{   
    //计算预瞄后的x, y，预瞄时间为ts
    double x_pre = x_now + vx * ts * cos(eular[2]) - vy * ts * sin(eular[2]);
    double y_pre = y_now + vy * ts * cos(eular[2]) + vx * ts * sin(eular[2]);

    int length = pathpoint.x.size();

    cout << "length= " << length << endl;
    double temp1 =x_pre - pathpoint.x.at(0);
    //double temp1 =x_now - pathpoint.x.at(0);
    double temp2 =y_pre - pathpoint.y.at(0);
    //double temp2 =y_now - pathpoint.y.at(0);
    double d_min = pow(temp1, 2) + pow(temp2, 2);
    int index = 0;
    for (size_t i = 0; i < length; i++)
    {
        temp1 = x_pre - pathpoint.x.at(i);
        //double temp3 = x_now - pathpoint.x.at(i);
        temp2= y_pre - pathpoint.y.at(i);
        //double temp4= y_now - pathpoint.y.at(i);
        double d_temp = pow(temp1, 2) + pow(temp2, 2);

        if (d_temp < d_min)
        {
            d_min = d_temp;
            index = i;
        }
        
    }

    cout << "index = " << index << endl;

    //计算出最近规划点处的切线和法线方向向量tor和nor,误差的距离向量d_err
    Eigen::Vector2d tor;
    Eigen::Vector2d nor;
    Eigen::Vector2d d_err;

    tor << cos(pathpoint.thetar[index]) , sin(pathpoint.thetar[index]);
    nor << -sin(pathpoint.thetar[index]), cos(pathpoint.thetar[index]);
    d_err << x_pre - pathpoint.x[index], y_pre - pathpoint.y[index];
    //d_err << x_now - pathpoint.x[index], y_now - pathpoint.y[index];
     
    //计算 phi_dot   
    double phi_dot = vx * pathpoint.kappar[index];

    //计算预瞄的 phi,预瞄的时间为ts
    double phi = eular[2] + phi_dot * ts;
    //double phi = eular[2];
    
    //计算ed
    double ed = nor.transpose() * d_err;

    //计算es
    double es = tor.transpose() * d_err;

    //计算投影点处的thetar
    double projection_point_thetar = pathpoint.thetar[index] + pathpoint.kappar[index] * es;

    //计算ed的导数 ed_dot
     double ed_dot = vy * cos(phi - projection_point_thetar) + vx * sin(phi - projection_point_thetar);

    //计算e_phi，由于e_phi较小，sin近似等于e_phi,但sin在[0,2pai]可以表示正负，即方向
    double e_phi = sin(phi - projection_point_thetar);
    //double e_phi = phi - projection_point_thetar;

    //计算s_dot
    double s_dot_fenzi = vx * cos(phi - projection_point_thetar) - vy * sin(phi - projection_point_thetar);
    double s_dot = s_dot_fenzi / (1-pathpoint.kappar[index] * ed);

    //计算 e_phi_dot
    double e_phi_dot = phi_dot - pathpoint.kappar[index] * s_dot;

    //计算 kr
    double kr = pathpoint.kappar[index];

    Eigen::VectorXd err_kr(5, 1);
    err_kr[0] = ed;
    err_kr[1] = ed_dot;
    err_kr[2] = e_phi;
    err_kr[3] = e_phi_dot;
    err_kr[4] = kr;

    return err_kr;

}

/************        前馈控制计算模块      ******************/
double cal_delta_err(Eigen::Vector4d LQR_K, Eigen::VectorXd err_kr, double vx)
{
    double K3 = LQR_K[2];
    double kr = err_kr[4];

    double delta_err = kr * (a + b - b * K3 - m * pow(vx, 2)/(a + b) * 
                       (b / cf + a * K3 / cr - a / cr));

    
    return delta_err;
}

/************        最终控制计算模块      ******************/
double cal_delta(Eigen::VectorXd err_kr, Eigen::Vector4d LQR_K, double delta_err)
{
    Eigen::Vector4d err;
    err << err_kr[0], err_kr[1], err_kr[2], err_kr[3];

    double delta = -LQR_K.transpose() * err + delta_err;

    cout << "angel_unlimit:" << delta << endl;

    if (delta > wheel_max_degree)
    {
        delta = wheel_max_degree;
    }
    else if (delta < -wheel_max_degree)
    {
        delta = -wheel_max_degree;
    }

    return delta;
}

/************        更新从 Gazebo 小车仿真处拿到底盘小车速度vx      ******************/
void velocityCallBack(const geometry_msgs::TwistStamped &carWaypoint)
{
    vx = carWaypoint.twist.linear.x;

}

/************        更新 Gazebo 小车的位置信息并发布小车行驶的轨迹点     ******************/
void poseCallBack(const geometry_msgs::PoseStamped &nowPose)
{
    //将小车的四元数转换为PRY欧拉角
    tf::Quaternion quat;
    tf::quaternionMsgToTF(nowPose.pose.orientation, quat);
    Eigen::Vector3d RPY;
    tf::Matrix3x3(quat).getRPY(RPY[0], RPY[1], RPY[2]);


    //当前小车的xy坐标
    double x_now = nowPose.pose.position.x;
    double y_now = nowPose.pose.position.y;

    cout << "x_now:" << x_now << "    y_now:" << y_now << endl; 

    //计算前轮转角
    Eigen::VectorXd err_kr =cal_err_kr(x_now, y_now, RPY);
    Eigen::Vector4d LQR_K = cal_LQR(vx);
    double delta_err = cal_delta_err(LQR_K, err_kr, vx);
    double angle = cal_delta(err_kr, LQR_K, delta_err);

    cout << "angle:" << angle << endl;

    //初始化 geometry_msgs::Twist 类型的消息，并取别称为 vel_msg
    geometry_msgs::Twist vel_msg;

    //给 vel_msg 中的线速度 linear.x 和角度 anger.z 分别赋值
    vel_msg.linear.x = 5;
    vel_msg.angular.z = angle;

    //将计算好的转角及线速度值发布
    LQR.publish(vel_msg);
    

    //构建一个存储小车行驶路径点的信息对象car_now_pose
    geometry_msgs::PoseStamped car_now_pose_msg;
    car_now_pose_msg.pose.position.x = x_now;
    car_now_pose_msg.pose.position.y = y_now;

    car_now_pose_msg.pose.orientation.x = nowPose.pose.orientation.x;
    car_now_pose_msg.pose.orientation.y = nowPose.pose.orientation.y;
    car_now_pose_msg.pose.orientation.z = nowPose.pose.orientation.z;
    car_now_pose_msg.pose.orientation.w = nowPose.pose.orientation.w;

    car_now_pose_msg.header.frame_id = "world";
    car_now_pose_msg.header.stamp = ros::Time::now();

    //将含有时间戳等信息的小车当前轨迹点信息存入car_trajectory_msgs <navgs_msgs/path>中
    car_trajectory_msgs.poses.push_back(car_now_pose_msg);
    //发布小车轨迹点话题，后续可在rviz里订阅
    car_trajectory_pub.publish(car_trajectory_msgs);



}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "LQR");
    ros::NodeHandle n;

    //计算五次多项式轨迹
    cal_path_number();

    LQR = n.advertise<geometry_msgs::Twist>("/smart/cmd_vel", 20);

    

/**************************************************************/
    //发布计算出的规划轨迹
    ros::Publisher planning_traj_pub;
    planning_traj_pub = n.advertise<nav_msgs::Path>("planningTrajctoryPath", 20, true);

    nav_msgs::Path planning_traj_msg; //用于存储规划轨迹信息对象
    planning_traj_msg.header.frame_id = "world";
    planning_traj_msg.header.stamp = ros::Time::now();

    geometry_msgs::PoseStamped pose_planning; //用于存储planning_traj_msg对象中geometry_msgs::PoseStamped消息格式的信息的对象
    pose_planning.header.frame_id = "world";
    pose_planning.header.stamp = ros::Time::now();

    int len = pathpoint.x.size();
    for (size_t i = 0; i < len; i++)
    {
        pose_planning.pose.position.x = pathpoint.x[i];
        pose_planning.pose.position.y = pathpoint.y[i];
        pose_planning.pose.position.z = 0.0;

        pose_planning.pose.orientation.x = 0.0;
        pose_planning.pose.orientation.y = 0.0;
        pose_planning.pose.orientation.z = 0.0;
        pose_planning.pose.orientation.w = 0.0;

        planning_traj_msg.poses.push_back(pose_planning);
    }

    planning_traj_pub.publish(planning_traj_msg);
/**************************************************************/

/**************************************************************/
    //发布小车的运动轨迹
    car_trajectory_pub = n.advertise<nav_msgs::Path>("car_trajectory", 20, true);
    car_trajectory_msgs.header.frame_id = "world";
    car_trajectory_msgs.header.stamp = ros::Time::now();

/**************************************************************/

    ros::Subscriber car_velocity = n.subscribe("/smart/velocity", 20, velocityCallBack);
    ros::Subscriber car_pose = n.subscribe("/smart/rear_pose", 20, poseCallBack); 

    ros::spin();
    
    return 0;

}


