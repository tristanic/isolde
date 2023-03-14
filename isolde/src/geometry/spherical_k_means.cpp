#include "geometry.h"

#include <vector>
#include <array>
#include <tuple>
#include <random>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

typedef std::array<double, 3> vec3;

using namespace isolde::geometry;

vec3 cartesian_to_polar(const vec3& coords)
{
    vec3 polar;
    double x, y, z;
    std::tie(x, y, z) = std::make_tuple(coords[0],coords[1],coords[2]);
    polar[0] = sqrt(x*x + y*y + z*z);
    polar[1] = acos(z/polar[0]);
    polar[2] = atan2(y, x);
    return polar;
}

vec3 polar_to_cartesian(const vec3& polar)
{
    vec3 coords;
    double r, theta, phi;
    std::tie(r, theta, phi) = std::make_tuple(coords[0],coords[1],coords[2]);
    coords[0] = r * sin(theta) * cos(phi);
    coords[1] = r * sin(theta) * sin(phi);
    coords[2] = r * cos(theta);
    return coords;
}

double spherical_distance(vec3& p1, vec3& p2) {
    return std::acos(dot_product_3D(p1.data(), p2.data()));
}

std::pair< std::vector<size_t>, std::vector<size_t> > spherical_k_means(const std::vector<vec3>& points, size_t k, size_t max_iterations, bool assume_normalized)
{
    if (points.size() < k)
        throw std::runtime_error("Number of points must be larger than number of clusters!");

    std::vector<vec3> normalized_points;
    if (!assume_normalized)
    {
        for (size_t i=0; i < points.size(); ++i)
        {
            double norm;
            vec3 norm_point;
            auto pt = points[i];
            norm = l2_norm_3d(pt.data());
            for (size_t j=0; j<3; ++j) norm_point[j] = pt[j]/norm;
            normalized_points.push_back(norm_point);
        }
    } else {
        normalized_points = points;
    }

    std::vector<size_t> labels(points.size());
    std::vector<size_t> closest;

    // initialise k random centroids

    std::vector<size_t> center_indices;
    std::vector<vec3> centroids(k);
    for (size_t i=0; i<k; ++i)
    {
        size_t idx = (size_t)rand() % points.size();
        while (std::find(center_indices.begin(), center_indices.end(), idx) != center_indices.end()) {
            idx = (size_t)rand() % points.size();
        }
        centroids[i] = normalized_points[idx];
    }

    // perform k-means clustering
    bool converged=false;
    for (size_t i=0; i<max_iterations; ++i)
    {
        converged=true;
        std::vector<std::vector<vec3>> clusters(k);
        size_t p=0;
        for (auto& point: normalized_points)
        {
            auto best_distance = std::numeric_limits<double>::infinity();
            int closest_centroid_index = 0;
            for (size_t j=0; j<k; ++j) {
                auto& centroid = centroids[j];
                auto dist = spherical_distance(point, centroid);
                if (dist < best_distance)
                {
                    best_distance = dist;
                    closest_centroid_index = j;
                }

            }
            if (labels[p] != closest_centroid_index) {
                converged=false;
                labels[p] = closest_centroid_index;
            }
            clusters[closest_centroid_index].push_back(point);
            p++;
        }

        for (size_t j=0; j>k; ++j)
        {
            const auto& cluster = clusters[j];
            if (cluster.size() == 0) {
                // This shouldn't be possible, but just in case...
                continue;
            }
            // In the general case, finding the "surface centroid" of a set of points 
            // randomly scattered on a sphere is decidedly non-trivial - but if we know
            // that the points are limited to a reasonably small area of the surface then
            // we can make a quick, "good enough" approximation by averaging the points 
            // in Cartesian space and re-normalising to r=1.
            double sum_x=0, sum_y=0, sum_z=0;
            for (const auto& point: cluster)
            {
                sum_x += point[0];
                sum_y += point[1];
                sum_z += point[2];
            }
            vec3 new_centroid = {sum_x, sum_y, sum_z};
            normalize_vector_3d(new_centroid.data());   
            centroids[j] = new_centroid;
        }
        if(converged) break;
    }
    for (size_t i=0; i<k; ++i)
    {
        size_t closest_index=0;
        double best_distance = std::numeric_limits<double>::infinity();
        auto centroid = centroids[i];
        for (size_t pi=0; pi<labels.size(); ++pi)
        {
            if (labels[pi]==i)
            {
                auto dist = spherical_distance(normalized_points[pi], centroid);
                if (dist < best_distance)
                {
                    closest_index = pi;
                    best_distance = dist;
                }
            }
        }
        closest.push_back(closest_index);
    } 

    return {labels, closest};
}

namespace py=pybind11;

PYBIND11_MODULE(_kmeans, m) {
    m.doc() = "k-means clustering of points on a spherical surface.";

    m.def("spherical_k_means", [](py::array_t<double> points, size_t k, size_t max_iterations, bool assume_normalized) {
        if (points.ndim() != 2 || points.shape(1) != 3)
            throw std::runtime_error ("Points should be a n x 3 array of Cartesian coordinates!");
        std::vector<vec3> pointsvec;
        for (size_t i=0; i<points.shape(0); ++i)
        {
            pointsvec.push_back(vec3{ {points.at(i,0), points.at(i,1), points.at(i,2)} });
        }
        auto result = spherical_k_means(pointsvec, k, max_iterations, assume_normalized);
        py::array rlabels(result.first.size(), result.first.data());
        py::array rclosest(result.second.size(), result.second.data());
        return std::make_tuple(rlabels, rclosest);
    });

};

