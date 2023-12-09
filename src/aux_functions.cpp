
#include <Rcpp.h>
using namespace Rcpp;              
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
//


//' @name compute_partial_arc_lengths
//' @title Compute partial arc lengths
//' @description Computes the partial arc lengths from a matrix of coordinates of a SpatialLine
//' @param coords [nx2 matrix] Matrix of the points of the lines
// [[Rcpp::export]]
Eigen::MatrixXd compute_partial_arc_lengths(Eigen::MatrixXd coords) {
    
    Eigen::MatrixXd out_mat(coords.rows(),3);
    Eigen::VectorXd arc_lengths(coords.rows());
    
    int i;

    // arc_lengths(0) = 0;
    // for(i=0; i<coords.rows()-1; i++){
    //     Eigen::VectorXd v = coords.row(i+1) - coords.row(i);
    //     arc_lengths(i+1) = arc_lengths(i) + v.norm();
    // } 
    
    // previous version
    // arc_lengths(coords.rows()-1) = 0; 
    // for(i = coords.rows()-1; i > 0; i--){
    //     Eigen::VectorXd v = coords.row(i) - coords.row(i-1);
    //     arc_lengths(i-1) = arc_lengths(i) + v.norm();
    // }  
    
    // Trapezoidal rule
    arc_lengths(coords.rows()-1) = 0; 
    Eigen::VectorXd v = coords.row(coords.rows()-1) - coords.row(coords.rows()-2);
    arc_lengths(coords.rows()-2) = arc_lengths(coords.rows()-1) + v.norm();  
    for(i = coords.rows()-2; i > 0; i--){
        Eigen::VectorXd v1 = coords.row(i) - coords.row(i-1);
        Eigen::VectorXd v2 = coords.row(i+1) - coords.row(i);
        v2(1) = v2(1)/v2(0)*V1(0);
        v2(0) = v1(0);
        arc_lengths(i-1) = arc_lengths(i) + (v1.norm() + v2.norm())/2;
    } 
    
    
    out_mat.col(0) = coords.col(0);
    out_mat.col(1) = coords.col(1);
    out_mat.col(2) = arc_lengths;

    return(out_mat);

}
