#ifndef IGL_TUFTED_COVER_LAPLACIAN_H
#define IGL_TUFTED_COVER_LAPLACIAN_H
#include <Eigen/Core>
#include <Eigen/Sparse>
// Given a possibly non-manifold mesh, construct the tufted cover Laplacian
// and mass matrix
//
// Inputs:
//   V     #V by 3 list of vertex positions
//   F     #F by 3 list of triangle indices into the rows of V
//   delta intrinsic mollification tolerance (default: 0.000001)
// Outputs:
//   L  #V by #V sparse Laplacian matrix
//   B  #V by #V sparse mass matrix
void tufted_cover_laplacian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  double delta,
  Eigen::SparseMatrix<double> & L,
  Eigen::SparseMatrix<double> & B);

void tufted_cover_laplacian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L,
  Eigen::SparseMatrix<double> & B);
#endif
