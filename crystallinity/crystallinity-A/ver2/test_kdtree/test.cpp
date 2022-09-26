#include <iostream>
#include <array>
#include <vector>
#include <kd-tree/kdtree.h>
#include <boost/progress.hpp>
#include <time.h>


class MyPoint : public std::array<double, 3>
{
public:

	// dimension of space (or "k" of k-d tree)
	// KDTree class accesses this member
	static const int DIM = 3;

	// the constructors
	MyPoint() {}
	MyPoint(double x, double y, double z)
	{ 
		(*this)[0] = x;
		(*this)[1] = y;
		(*this)[2] = z;
	}
};


std::vector<int> straight(
        const MyPoint &query,
        std::vector<MyPoint> &points,
        const double &radius
){
    double distance, radrad;
    radrad = radius * radius;
    std::vector<int> idxes;

    for (int i = 0; i < points.size(); i++){
        distance = 0.;
        for (int j = 0; j < 3; j++){
            distance += (query[j] - points[i][j])*(query[j] - points[i][j]);
        }
        if (distance <= radrad) idxes.push_back(i);
    }
    return idxes;
}


int main(int argc, char **argv)
{
    std::cout << "=====================================" << std::endl;
    std::cout << "     TEST FOR NEIGHBORHOOD QUERY" << std::endl;
    std::cout << "=====================================" << std::endl;
	const int seed = argc > 1 ? std::stoi(argv[1]) : 0;
	srand(seed);

	const double width = 100;
	const double height = 100;
	const double depth = 100;

	// generate points
	const int npoints = 30000;
	const double radius = 5;
	std::vector<MyPoint> points(npoints);
	for (int i = 0; i < npoints; i++)
	{
		const double x = (double)rand() / RAND_MAX * width;
		const double y = (double)rand() / RAND_MAX * height;
		const double z = (double)rand() / RAND_MAX * height;
		points[i] = MyPoint(x, y, z);
	}

	MyPoint query;
    std::cout << "dim: 3\nnum particles: " << npoints << std::endl;

    clock_t start, end, time;
    std::vector<std::vector<int> > num_neighbors(npoints, std::vector<int>(2));

    std::cout << "--- KBTree ---" << std::endl;
    start = clock();

    boost::progress_display show_progress(npoints);
    kdt::KDTree<MyPoint> kdtree(points);
    for (int i = 0; i < npoints; i++){
        query = points[i];
        const std::vector<int> radIndices = kdtree.radiusSearch(query, radius);
        num_neighbors[i][0] = radIndices.size();
        //for (int i = 0; i < radIndices.size(); i++){
        //    int idx = radIndices[i];
        //    std::cout << idx << ": ";
        //    std::cout << points[idx][0] << ", "
        //        << points[idx][1] << ", "
        //        << points[idx][2] << std::endl;
        //}
        ++show_progress;
    }
    end = clock();
    time = (double)(end - start) / (double)CLOCKS_PER_SEC * 1000.;
    std::cout << "--------------\n";
    std::cout << "time: " << time << " [ms]\n\n" << std::endl;

    std::cout << "--- straight ---" << std::endl;
    boost::progress_display show_progress2(npoints);
    start = clock();
    for (int i = 0; i < npoints; i++){
        query = points[i];
        const std::vector<int> idxes = straight(query, points, radius);
        num_neighbors[i][1] = idxes.size();
        //for (int i = 0; i < idxes.size(); i++){
        //    int idx = idxes[i];
        //    std::cout << idx << ": ";
        //    std::cout << points[idx][0] << ", "
        //        << points[idx][1] << ", "
        //        << points[idx][2] << std::endl;
        //}
        ++show_progress2;
    }
    end = clock();
    time = (double)(end - start) / (double)CLOCKS_PER_SEC * 1000.;
    std::cout << "--------------\n";
    std::cout << "time: " << time << " [ms]\n\n" << std::endl;
}
