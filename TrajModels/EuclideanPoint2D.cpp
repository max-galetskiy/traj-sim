#include "EuclideanPoint2D.hpp"
#include <cmath>
#include <iostream>

// Projects LongLat Points into Euclidean space ("meters" decides whether m or km is the wanted unit)
EuclideanPoint2D::EuclideanPoint2D(std::pair<float, float> longLat, bool meters = true) {
    // Classical Mercator Projection (works well enough for local areas)
    float longitude = (longLat.first * M_PI)/180;
    float latitude = (longLat.second * M_PI)/180;

    float earth_radius = meters ? 6371008 : 6371;

    x = longitude * earth_radius ;
    y = 0.5f * logf((1 + sinf(latitude))/(1 - sinf(latitude))) * earth_radius;
}

std::string EuclideanPoint2D::to_string() {
	return "(" + std::to_string(x) + " , " + std::to_string(y) + ")";
}


float EuclideanPoint2D::dist(EuclideanPoint2D& a, EuclideanPoint2D& b)
{
	return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

/*
    Calculates the intersection of the circle and line, returns not the intersection points but the timestamps of the line, where
            - i      is the time at l1
            - i + 1  is the time at l2
    Returns a null pair if intersection does not exist, otherwise returns timestamps of the segment contained within the circle
*/
fpair EuclideanPoint2D::circleLineIntersection(EuclideanPoint2D &center, EuclideanPoint2D &l1, EuclideanPoint2D &l2, float r, int i) {

    // This approach follows the algorithm described by Wolfram MathWorld: https://mathworld.wolfram.com/Circle-LineIntersection.html

    // Edge case l1 === l2
    if(l1.x == l2.x && l1.y == l2.y){
        return {i,i+1};
    }

    // Translate lines so that the circle's center becomes the origin
    auto tl1 = EuclideanPoint2D(l1.x - center.x, l1.y - center.y);
    auto tl2 = EuclideanPoint2D(l2.x - center.x, l2.y - center.y);

    // For Floating Point comparison to 0
    float threshold = 1e-8;

    float x_diff = tl2.x - tl1.x;
    float y_diff = tl2.y - tl1.y;
    float diff = std::sqrt(x_diff * x_diff + y_diff * y_diff);
    float determinant = tl1.x * tl2.y - tl2.x * tl1.y;

    float discriminant = r*r * diff*diff - determinant*determinant;

    // No intersection
    if(discriminant < 0){
        return FreeSpaceCell::null_pair;
    }
    else{

        float signum = (y_diff < 0) ? -1 : 1;
        float root = std::sqrt(discriminant);

        // Since we only care about timestamps, we will only compute
        // the x points of intersections (and use y points in case the x difference equals zero)

        // One intersection point
        if(discriminant < threshold){

            float t;

            if(tl2.x != tl1.x){
                float x = (determinant * y_diff)/(diff*diff);
                t = (x - tl1.x) / (tl2.x - tl1.x);
            }else{
                float y = (-determinant * x_diff)/(diff*diff);
                t = (y - tl1.y) / (tl2.y - tl1.y);
            }

            // Handle over/underflow
            if(t < 0 || t > 1){
                return FreeSpaceCell::null_pair;
            }

            return {t + i, t + i};
        }
        // Two intersection points
        else{

            float t1;
            float t2;

            if(tl1.x != tl2.x){
                float root_term = signum * x_diff * root;
                float x1 = (determinant * y_diff + root_term)/(diff*diff);
                float x2 = (determinant * y_diff - root_term)/(diff*diff);
                t1 = (x1 - tl1.x) / (tl2.x - tl1.x);
                t2 = (x2 - tl1.x) / (tl2.x - tl1.x);

            }
            else{
                float root_term = std::abs(y_diff) * root;
                float y1 = (-determinant * x_diff + root_term)/(diff * diff);
                float y2 = (-determinant * x_diff - root_term)/(diff * diff);
                t1 = (y1 - tl1.y) / (tl2.y - tl1.y);
                t2 = (y2 - tl1.y) / (tl2.y - tl1.y);

            }


            // Handle over/underflow
            if((t1 < 0 && t2 < 0) || (t1 > 1 && t2 > 1)){   // Both intersection points are outside of line on same sides
                return FreeSpaceCell::null_pair;
            }

            t1 = (t1 > 1) ? 1 : ((t1 < 0) ? 0 : t1);
            t2 = (t2 > 1) ? 1 : ((t2 < 0) ? 0 : t2);

            return {std::min(t1,t2) + i, std::max(t1,t2) + i};
        }
    }

}

/*
     Creates the Free Space cell [p1,p2] x [q1,q2] (if this does not mean anything to you, refer to https://doi.org/10.1142/S0218195995000064)
 */
FreeSpaceCell EuclideanPoint2D::createFreeSpace(EuclideanPoint2D &p1, EuclideanPoint2D &q1, EuclideanPoint2D &p2, EuclideanPoint2D &q2, float epsilon, int i, int j, std::shared_ptr<DistCache<EuclideanPoint2D>> cache) {

    // Left   (i,j) -> (i + 1, j)
    auto left = EuclideanPoint2D::circleLineIntersection(q1, p1, p2, epsilon, i);

    // Bottom (i,j) -> (i, j + 1)
    auto bottom = EuclideanPoint2D::circleLineIntersection(p1, q1, q2, epsilon, j);

    // Right  (i,j + 1) -> (i + 1, j + 1)
    auto right = EuclideanPoint2D::circleLineIntersection(q2, p1, p2, epsilon, i);

    // Top    (i + 1,j) -> (i + 1, j + 1)
    auto top = EuclideanPoint2D::circleLineIntersection(p2, q1, q2, epsilon, j);

    return {i, j, fs_intervals_per_side, left, bottom, top, right};
}

/*
        Adjusts the right and top boundaries of a given Free Space cell to refer to the reachable space (and not just the free space)
 */
void EuclideanPoint2D::propagateFreeSpace(EuclideanPoint2D &p1, EuclideanPoint2D &q1, EuclideanPoint2D &p2, EuclideanPoint2D &q2, FreeSpaceCell &cell,float epsilon, std::shared_ptr<DistCache<EuclideanPoint2D>> cache) {

        // If cell has no left and no bottom boundaries
        if(cell.getLeft().first == FreeSpaceCell::null_pair && cell.getBottom().first == FreeSpaceCell::null_pair){
            cell.setRight(FreeSpaceCell::null_pair);
            cell.setTop(FreeSpaceCell::null_pair);
            return;
        }

        // Create right and top boundaries (if they do not exist already)
        if(cell.getRight().first == FreeSpaceCell::null_pair){
            auto right = EuclideanPoint2D::circleLineIntersection(q2, p1, p2, epsilon, cell.i);
            cell.setRight(right);
        }
        if(cell.getTop().first == FreeSpaceCell::null_pair){
            auto top = EuclideanPoint2D::circleLineIntersection(p2, q1, q2, epsilon, cell.j);
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

/*
        Computes given Free Space cell side
 */
void EuclideanPoint2D::computeFreeSpaceSide(EuclideanPoint2D& p1, EuclideanPoint2D& q1, EuclideanPoint2D& p2, EuclideanPoint2D& q2, FreeSpaceCell& cell, float epsilon, side s, std::shared_ptr<DistCache<EuclideanPoint2D>> cache){

    switch (s) {
        case BOTTOM:
            cell.setBottom(EuclideanPoint2D::circleLineIntersection(p1, q1, q2, epsilon, cell.j));
            return;
        case LEFT:
            cell.setLeft(EuclideanPoint2D::circleLineIntersection(q1, p1, p2, epsilon, cell.i));
            return;
        case RIGHT:
            cell.setRight(EuclideanPoint2D::circleLineIntersection(q2, p1, p2, epsilon, cell.i));
            return;
        case TOP:
            cell.setTop(EuclideanPoint2D::circleLineIntersection(p2, q1, q2, epsilon, cell.j));
            return;
    }
}

/*
    Calculates the distance between p and the edge/segment defined by [e1,e2]
 */
float EuclideanPoint2D::dist_to_edge(EuclideanPoint2D &p, EuclideanPoint2D &e1, EuclideanPoint2D &e2) {

    auto edge_vec = EuclideanPoint2D(e2.x - e1.x, e2.y - e1.y);
    auto p_e1_vec = EuclideanPoint2D(p.x - e1.x, p.y - e1.y);

    float projection = edge_vec.x * p_e1_vec.x + edge_vec.y * p_e1_vec.y;
    float edge_length = edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y;

    if(edge_length == 0){
        return -1;
    }

    // Normalize projection
    projection = projection / (edge_length * edge_length);

    if(projection <= 0){ // closest point is e1
        return EuclideanPoint2D::dist(e1, p);
    }
    else if(projection >= 1) { // closest point is e2
        return EuclideanPoint2D::dist(e2, p);
    }else{ // closest point is in between e1 and e2
        auto closest_point = EuclideanPoint2D(e1.x + edge_vec.x * projection, e1.y + edge_vec.y * projection);
        return EuclideanPoint2D::dist(closest_point, p);
    }

}

/*
    Computes the distance of p1/p2 to the intersection of p1 and p2s bisector with the edge e1e2
 */
float EuclideanPoint2D::dist_of_bisector(EuclideanPoint2D &p1, EuclideanPoint2D &p2, EuclideanPoint2D &e1, EuclideanPoint2D &e2) {

    /*

            Origin of this formula:
            We aim to solve the equation dist(p1,X) = dist(p2,X)
            with X being the intersection of the bisector between p1 and p2 and the line spanning from e1 to e2
            Thus X has the form  e1 + t * (e2 - e1)

            Thus, one needs to solve the equation:
            dist(p1,X) = dist(p2,X)
            <=>
            (p1.x - x.x)^2 + (p1.y - x.y)^2 = (p2.x - x.x)^2 + (p2.y - x.y)^2

            After a lot of arithmetic fiddling, you will then end up with the solution
            t = numerator/denominator
            With the equations for numerator and denominator displayed below

     */

    float numerator = (p1.x - p2.x) * (p1.x + p2.x - 2 * e1.x * (e2.x - e1.x)) + (p1.y - p2.y) * (p1.y + p2.y - 2 * e1.y * (e2.y - e1.y));
    float denominator = 2 * ((e2.x - e1.x) * (p1.x - p2.x) + (e2.y - e1.y) * (p1.y - p2.y));

    if(denominator == 0){
        return -1;
    }

    float t = numerator/denominator;

    if(t < 0 || t > 1){
        return -1;
    }

    auto intersection = EuclideanPoint2D(e1.x + t * (e2.x - e1.x), e1.y + t * (e2.y - e1.y));
    return EuclideanPoint2D::dist(p1, intersection);

}