/*
This program randomly forms polygons with edges from 3--10
Note: These polygons are not just a mixture of random coordnates, they are coordnates that can form a shape when connected in order
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <algorithm>

using namespace std;

struct Point {
    float x, y;

    // Overloading the < operator to make sorting easier.
    bool operator < (const Point &p) const {
        return x < p.x || (x == p.x && y < p.y);
    }
};

int orientation(Point p, Point q, Point r) {
    float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;  // collinear
    return (val > 0) ? 1 : 2;  // clockwise or counterclockwise
}

// Function to compute the convex hull of a set of points.
vector<Point> convexHull(vector<Point> points) {
    int n = points.size();

    // There must be at least 3 points.
    if (n < 3) return {};

    vector<Point> hull;

    // Find the leftmost point.
    int l = 0;
    for (int i = 1; i < n; i++)
        if (points[i].x < points[l].x)
            l = i;

    int p = l, q;
    do {
        hull.push_back(points[p]);
        q = (p + 1) % n;

        for (int i = 0; i < n; i++) {
            if (orientation(points[p], points[i], points[q]) == 2)
                q = i;
        }
        p = q;
    } while (p != l);

    return hull;
}

int main() {
	freopen("in.txt","w",stdout);
	cout<<700<<' '<<700<<endl;
	cout<<fixed<<setprecision(2);
    srand(time(0));

    int numPoints = 10; 
    float maxWidth = 100.0;
    float maxHeight = 100.0;

    int numPolygons = 57;  // change this to generate a different number of polygons
	cout<<numPolygons<<endl;
    for (int j = 0; j < numPolygons; j++) {
        vector<Point> points;

        for (int i = 0; i < numPoints; i++) {
            Point p = { maxWidth * (float)rand() / RAND_MAX, maxHeight * (float)rand() / RAND_MAX };
            points.push_back(p);
        }

        vector<Point> hull = convexHull(points);

        cout << hull.size() << "\n";
        for (const auto& point : hull) {
            cout <<point.x << " " << point.y << "\n";
        }
        cout << "\n";
    }

    return 0;
}

