#include <vector>
#include <eigen3/Eigen/Dense>
using namespace std;

class quinticType
{
public:
    vector<float> x;
    vector<float> vx;
    vector<float> ax;

    vector<float> y;
    vector<float> vy;
    vector<float> ay;

    vector<float> thetar;
    vector<float> kappar;

    Eigen::Vector2d tor;
    Eigen::Vector2d nor;

};