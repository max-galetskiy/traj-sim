#include "EuclideanPoint3D.hpp"
#include <cmath>

/*
    Converts GPS Coordinates to Cartesian coordinates (with 1 unit being 1km or 1m depending on the meters boolean)
 */
EuclideanPoint3D::EuclideanPoint3D(std::pair<float, float> longLat, bool meters) {

    // Convert to radians
    float longitude = (longLat.first * M_PI)/180;
    float latitude = (longLat.second * M_PI)/180;

    float earth_radius = meters ? 6371008 : 6371;

    x = earth_radius * cos(latitude) * cos(longitude);
    y = earth_radius * cos(latitude) * sin(longitude);
    z = earth_radius * sin(latitude);
}

std::string EuclideanPoint3D::to_string() {
    return "(" + std::to_string(x) + " , " + std::to_string(y) + " , " + std::to_string(z) + ")";
}

float EuclideanPoint3D::dist(EuclideanPoint3D &a, EuclideanPoint3D &b) {
    return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

fpair EuclideanPoint3D::circleLineIntersection(EuclideanPoint3D &center, EuclideanPoint3D &l1, EuclideanPoint3D &l2,
                                                                                  float r, int i) {

    // This approach follows the algorithm described by Wikipedia: https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    // Edge case l1 === l2
    if(l1.x == l2.x && l1.y == l2.y){
        return {i, i + 1};
    }

    /*
            Combining the equations dist(x,c)^2 = r ^2 and x = l1 + t * l2
            Leads to: d^2 * (u * u) + 2 d * (u * (o - c)) + (o - c)^2 - r^2 = 0
            Which follows the form of a quadratic formula in d
     */

    auto u = EuclideanPoint3D(l2.x - l1.x, l2.y - l1.y, l2.z - l1.z);
    auto c_to_l1 = EuclideanPoint3D(l1.x - center.x, l1.y - center.y, l1.z - center.z);

    float c_l1_norm = c_to_l1.x * c_to_l1.x + c_to_l1.y * c_to_l1.y + c_to_l1.z * c_to_l1.z;

    float a = u.x * u.x + u.y * u.y + u.z * u.z;
    float b = (u.x * c_to_l1.x + u.y * c_to_l1.y + u.z * c_to_l1.z); // normally there would be a factor of two, but in the simplification this falls away
    float c = c_l1_norm - r * r;

    float discriminant = b * b - a * c;

    // For Floating Point comparison to 0
    float threshold = 1e-8;

    if(discriminant < 0){
        return FreeSpaceCell::null_pair;
    }
    else if(discriminant < threshold){

        float t = -b/a;

        if(t < 0 || t > 1){
            return FreeSpaceCell::null_pair;
        }

        return {t + i, t + i};
    }

    float t1 = (-b + std::sqrt(discriminant))/a;
    float t2 = (-b - std::sqrt(discriminant))/a;

    // Handle over/underflow
    if((t1 < 0 && t2 < 0) || (t1 > 1 && t2 > 1)){   // Both intersection points are outside of line on same sides
        return FreeSpaceCell::null_pair;
    }

    t1 = (t1 > 1) ? 1 : ((t1 < 0) ? 0 : t1);
    t2 = (t2 > 1) ? 1 : ((t2 < 0) ? 0 : t2);

    return {std::min(t1,t2) + i, std::max(t1,t2) + i};

}

/*
    Again, creates Free Space cell [p1,p2] x [q1,q2]
 */
FreeSpaceCell EuclideanPoint3D::createFreeSpace(EuclideanPoint3D &p1, EuclideanPoint3D &q1, EuclideanPoint3D &p2, EuclideanPoint3D &q2, float epsilon, int i, int j, std::shared_ptr<DistCache<EuclideanPoint3D>> cache) {

    // Left   (i,j) -> (i + 1, j)
    auto left = EuclideanPoint3D::circleLineIntersection(q1, p1, p2, epsilon, i);

    // Bottom (i,j) -> (i, j + 1)
    auto bottom = EuclideanPoint3D::circleLineIntersection(p1, q1, q2, epsilon, j);

    // Right  (i,j + 1) -> (i + 1, j + 1)
    auto right = EuclideanPoint3D::circleLineIntersection(q2, p1, p2, epsilon, i);

    // Top    (i + 1,j) -> (i + 1, j + 1)
    auto top = EuclideanPoint3D::circleLineIntersection(p2, q1, q2, epsilon, j);

    return {i, j, fs_intervals_per_side, left, bottom, top, right};
}

/*
        Adjusts the right and top boundaries of a given Free Space cell (to refer to the reachable space and not the free space)
 */
void EuclideanPoint3D::propagateFreeSpace(EuclideanPoint3D &p1, EuclideanPoint3D &q1, EuclideanPoint3D &p2, EuclideanPoint3D &q2, FreeSpaceCell &cell,
                                          float epsilon, std::shared_ptr<DistCache<EuclideanPoint3D>> cache) {

    // If cell has no left and no bottom boundaries
    if(cell.getLeft().first == FreeSpaceCell::null_pair && cell.getBottom().first == FreeSpaceCell::null_pair){
        cell.setRight(FreeSpaceCell::null_pair);
        cell.setTop(FreeSpaceCell::null_pair);
        return;
    }

    // Create right and top boundaries (if they do not exist already)
    if(cell.getRight().first == FreeSpaceCell::null_pair){
        auto right = EuclideanPoint3D::circleLineIntersection(q2, p1, p2, epsilon, cell.i);
        cell.setRight(right);
    }
    if(cell.getTop().first == FreeSpaceCell::null_pair){
        auto top = EuclideanPoint3D::circleLineIntersection(p2, q1, q2, epsilon, cell.j);
        cell.setTop(top);
    }

    // If cell has both left and bottom boundaries, then there are no restrictions on the right and top boundaries
    if(cell.getLeft().first != FreeSpaceCell::null_pair && cell.getRight().first != FreeSpaceCell::null_pair){
        return;
    }

    // Bounds
    float i_min = -INFINITY;
    float j_min = -INFINITY;

    // Propagate from left
    if(cell.getLeft().first != FreeSpaceCell::null_pair){
        i_min = cell.getLeft().first.first;
    }

    // Propagate from bottom
    if(cell.getBottom().first != FreeSpaceCell::null_pair){
        j_min = cell.getBottom().first.first;
    }

    // Adjust bounds
    if(cell.getRight().first != FreeSpaceCell::null_pair){
        float r_beg = cell.getRight().first.first;
        float r_end = cell.getRight().first.second;

        r_beg = std::max(r_beg, i_min);

        cell.setRight({r_beg, r_end});
    }
    if(cell.getTop().first != FreeSpaceCell::null_pair){
        float t_beg = cell.getTop().first.first;
        float t_end = cell.getTop().first.second;

        t_beg = std::max(t_beg, j_min);

        cell.setTop({t_beg, t_end});
    }

}

float EuclideanPoint3D::dist_to_edge(EuclideanPoint3D &p, EuclideanPoint3D &e1, EuclideanPoint3D &e2) {

    auto edge_vec = EuclideanPoint3D(e2.x - e1.x, e2.y - e1.y, e2.z - e1.z);
    auto p_e1_vec = EuclideanPoint3D(p.x - e1.x, p.y - e1.y, p.z - e1.z);

    float projection = edge_vec.x * p_e1_vec.x + edge_vec.y * p_e1_vec.y + edge_vec.z * p_e1_vec.z;
    float edge_length = edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y + edge_vec.z * edge_vec.z;

    // Normalize projection
    projection = projection / (edge_length * edge_length);

    if(projection <= 0){ // closest point is e1
        return EuclideanPoint3D::dist(e1, p);
    }
    else if(projection >= 1) { // closest point is e2
        return EuclideanPoint3D::dist(e2, p);
    }else{ // closest point is in between e1 and e2
        auto closest_point = EuclideanPoint3D(e1.x + edge_vec.x * projection, e1.y + edge_vec.y * projection, e1.z + edge_vec.z * projection);
        return EuclideanPoint3D::dist(closest_point, p);
    }

}

float EuclideanPoint3D::dist_of_bisector(EuclideanPoint3D &p1, EuclideanPoint3D &p2, EuclideanPoint3D &e1, EuclideanPoint3D &e2) {

    float numerator = (p1.x - p2.x) * (p1.x + p2.x - 2 * e1.x * (e2.x - e1.x)) +
                      (p1.y - p2.y) * (p1.y + p2.y - 2 * e1.y * (e2.y - e1.y)) +
                      (p1.z - p2.z) * (p1.z + p2.z - 2 * e1.z * (e2.z - e1.z));

    float denominator = 2 * ((e2.x - e1.x) * (p1.x - p2.x) + (e2.y - e1.y) * (p1.y - p2.y) + (e2.z - e1.z) * (p1.z - p2.z));

    if(denominator == 0){
        return -1;
    }

    float t = numerator/denominator;

    if(t < 0 || t > 1){
        return -1;
    }

    auto intersection = EuclideanPoint3D(e1.x + t * (e2.x - e1.x), e1.y + t * (e2.y - e1.y), e1.z + t * (e2.z - e1.z));
    return EuclideanPoint3D::dist(p1, intersection);
}

void EuclideanPoint3D::computeFreeSpaceSide(EuclideanPoint3D &p1, EuclideanPoint3D &q1, EuclideanPoint3D &p2, EuclideanPoint3D &q2, FreeSpaceCell &cell, float epsilon, side s, std::shared_ptr<DistCache<EuclideanPoint3D>> cache) {
    switch (s) {
        case BOTTOM:
            cell.setBottom(EuclideanPoint3D::circleLineIntersection(p1, q1, q2, epsilon, cell.j));
            return;
        case LEFT:
            cell.setLeft(EuclideanPoint3D::circleLineIntersection(q1, p1, p2, epsilon, cell.i));
            return;
        case RIGHT:
            cell.setRight(EuclideanPoint3D::circleLineIntersection(q2, p1, p2, epsilon, cell.i));
            return;
        case TOP:
            cell.setTop(EuclideanPoint3D::circleLineIntersection(p2, q1, q2, epsilon, cell.j));
            return;
    }
}


