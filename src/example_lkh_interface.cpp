#include <ros/ros.h>
#include <lkh_tsp_solver/lkh_interface.h>

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "example_lkh_interface");
    ros::NodeHandle nh;

    int result = solveTSPLKH("/home/ahmad/personal_ws/src/lkh_tsp_solver/resource/pr2392.par");
    ROS_INFO("LKH TSP Result: %d", result);

    return 0;
}
