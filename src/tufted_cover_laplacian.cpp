#include "tufted_cover_laplacian.h"
#include "tufted_cover.h"
#include "flip_to_delaunay.h"
#include <igl/edge_lengths.h>
#include <igl/cotmatrix_intrinsic.h>
#include <igl/massmatrix_intrinsic.h>

void tufted_cover_laplacian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L,
  Eigen::SparseMatrix<double> & B)
{
  tufted_cover_laplacian(V,F,0.000001,L,B);
}


void tufted_cover_laplacian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  double delta,
  Eigen::SparseMatrix<double> & L,
  Eigen::SparseMatrix<double> & B)
{
  Eigen::MatrixXi Ftilde, Gtilde;
  tufted_cover(V,F,Ftilde,Gtilde);

  Eigen::MatrixXd l, lout;
  igl::edge_lengths(V,Ftilde,l);

  // intrinsic mollification
  if (delta > 0)
  {
    double epsilon = 0;
    for (int i=0; i<l.rows(); i++)
    for (int j=0; j<3; j++)
    {
      double value = delta - l(i,(j+1)%3) - l(i,(j+2)%3) + l(i,j);
      if (value > epsilon)
        epsilon = value;
    }
    l = l.array() + epsilon;
  }

  Eigen::MatrixXi Fout, Gout;
  flip_to_delaunay(Ftilde,Gtilde,l, Fout,Gout,lout);

  igl::cotmatrix_intrinsic(lout,Fout,L);
  L *= 0.5;

  igl::massmatrix_intrinsic(lout,Fout,igl::MASSMATRIX_TYPE_DEFAULT,B);
  B *= 0.5;
}