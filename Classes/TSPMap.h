#pragma once

#include <vector>
#include "City.h"
using namespace std;

class TSPMap {
private:
	vector<Vec2> mv_cities;
public:
	TSPMap() {}

	void setCitys(vector<City*>& cities) {
		mv_cities.clear();
		for (int i = 0, n = cities.size(); i < n; ++i) {
			mv_cities.push_back(cities.at(i)->getPosition());
		}
	}

	int numberOfCities() { return mv_cities.size(); }
	Vec2 PosOfCity(int i) { return mv_cities[i];}

	//calculate the square length of the route
	double lengthOfRouteSquare(vector<int>& route) {
		double sum = 0;
		auto &c = mv_cities;
		auto &r = route;
		int n = route.size();
		for (int i = 0; i < n-1; ++i) {
			sum += distSquare(c[ r[i] ], c[ r[i + 1] ]);
		}
		sum += distSquare(c[ r[n - 1] ], c[ r[0] ]);
		return sum;
	}

	//square length between point a and b
	inline double distSquare(Vec2& a, Vec2& b) { return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y); }

};