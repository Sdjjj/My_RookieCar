# My_RookieCar
采用Gazebo里搭建的小车模型，对ROS下C++写的基于车辆动力学模型的纯横向LQR控制算法在五次多项式路径下进行验证

仿真步骤：

1、使用前请安装Gazebo关节控制相关依赖
```
sudo apt-get install -y ros-kinetic-gazebo-ros-control
sudo apt-get install -y ros-kinetic-ros-control ros-kinetic-ros-controllers
sudo apt-get install -y ros-kinetic-gazebo-ros-control
```

2、编译ROS工作空间：

在/src所在目录下，终端输入：
   
        catkin_make

3、启动Gazebo下的空白世界并加载小车模型

新建终端，输入以下命令： 

        source /devel/setup.sh                 
        roslaunch car_model spawn_car.launch
        
4、启动rviz

rviz显示的配置文件为 myLQR.rviz ，请在Rviz里添加该配置文件

新建终端，输入以下命令：

        source /devel/setup.sh
        rivz      
        
5、启动LQR算法

新建终端，输入以下命令：


       source /devel/setup.sh
       rosrun road horizontalLQR

        
调整Q和R以及预瞄时间值，可改变不同速度下的跟踪效果以及车辆摆振效果
红色轨迹为生成的五次多项式轨迹，绿色轨迹为小车行驶轨迹
跟踪效果如图所示

![0}DW%)A%((8B2P4MKY~WHM9](https://user-images.githubusercontent.com/75204388/140480619-02f021b3-5c20-4e29-8003-5fae0d44fe73.png)


