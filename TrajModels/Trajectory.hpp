#pragma once
#include <vector>
#include <memory>
#include <initializer_list>
#include <stdexcept>

/*
    Abstract model of a Trajectory
	Template class Point needs to support the following functions: 
		- std::string  Point::to_string()
		- static float Point::dist(Point& a, Point& b)
*/ 
template<class Point>
class Trajectory
{
private:

	std::vector<std::shared_ptr<Point>> plist;

public: 

	Trajectory(std::initializer_list<std::shared_ptr<Point>> list) : plist(list) {};
    explicit Trajectory(std::vector<std::shared_ptr<Point>> list) : plist(list) {};
	
	Point& at(unsigned int i);
	unsigned int length();

    std::string to_string();

};


template<class Point>
Point& Trajectory<Point>::at(unsigned int i) {

	if (i >= plist.size()) {
		std::throw_with_nested(std::invalid_argument("Trajectory Index out of Bounds!"));
	}

	return *(plist.at(i));
}

template<class Point>
unsigned int Trajectory<Point>::length() {
	return plist.size();
}

template<class Point>
std::string Trajectory<Point>::to_string(){
    std::string traj_string;

    for(auto p_ptr : plist){
        traj_string += (*p_ptr).to_string() + " ";
    }

    return traj_string;
}