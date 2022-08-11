#include <ros/ros.h>
#include <lkh_tsp_solver/lkh_interface.h>
#include <Eigen/Eigen>
#include <random_numbers/random_numbers.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "example_lkh_interface");
    ros::NodeHandle nh;
    ros::Publisher pub;
    pub = nh.advertise<nav_msgs::Path>("/tour", 1, true);

    const int n = 7;
    Eigen::Vector3d pts[n];
    Eigen::MatrixXd cost_mat(n, n);
    pts[1] = Eigen::Vector3d(2, 1, 0);
    pts[0] = Eigen::Vector3d(10, 1, 0);
    pts[3] = Eigen::Vector3d(18, 1, 0);
    pts[2] = Eigen::Vector3d(26, 7, 0);
    pts[4] = Eigen::Vector3d(17, 9, 0);
    pts[6] = Eigen::Vector3d(8, 9, 0);
    pts[5] = Eigen::Vector3d(2, 9, 0);

    random_numbers::RandomNumberGenerator rng(1);

    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            /*
                The cost consist of two components (distance and angle).
                The distance is the Euclidean distance between the two points.
                The Euler angles (Roll, Pitch, and Yaw is calculated and the
                absolute value of the angles are added together).
            */
            Eigen::Vector3d v = pts[i] - pts[j];
            double d_cost = v.norm();
            double a_cost = fabs(atan2(v.y(), v.x())) * 180 / M_PI * 0; // TODO: 
            ROS_INFO("(%f,%f,%f) -> (%f,%f,%f): %f(dCost) + %f(aCost) = %f", pts[i].x(), pts[i].y(), pts[i].z(), pts[j].x(), pts[j].y(), pts[j].z(), j, d_cost, a_cost, d_cost + a_cost);
            // v = v.eulerAngles(0, 1, 2);
            // cost += (fabs(v.x()) + fabs(v.y()) + fabs(v.z()));
            cost_mat(i, j) = cost_mat(j, i) = (d_cost + a_cost);
        }
    }
    // std::cout << cost_mat << std::endl;
    LKH_TSP_Solver tsp("/home/ahmad/catkin_ws/src/lkh_tsp_solver/resource/", "test");
    tsp.setCostMatrix(cost_mat);
    ROS_INFO("LKH TSP Result: %d", tsp.solve());
    std::vector<int> result = tsp.getTour();

    nav_msgs::Path msg;
    msg.header.frame_id = "map";
    msg.header.stamp = ros::Time::now();
    msg.poses.resize(result.size());
    for (auto idx : result)
    {
        geometry_msgs::PoseStamped pt;
        pt.header = msg.header;
        pt.pose.orientation.w = 1;
        pt.pose.position.x = pts[idx].x();
        pt.pose.position.y = pts[idx].y();
        pt.pose.position.z = pts[idx].z();
        msg.poses.push_back(pt);
        std::cout << idx+1 << " -> " << pts[idx].x() << " " << pts[idx].y() << " " << pts[idx].z() << " " << std::endl;
    }
    pub.publish(msg);
    ros::Rate loop_rate(10);
    while (ros::ok())
    {
        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}
