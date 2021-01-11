#include "tufted_cover.h"
#include "flip_to_delaunay.h"
#include "tufted_cover_laplacian.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/edge_lengths.h>
#include <iostream>

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  igl::read_triangle_mesh(argc>1?argv[1]:"../data/nonmanifold.off",V,F);

	Eigen::MatrixXi Ftilde, Gtilde, Fout, Gout;
	tufted_cover(V,F,Ftilde,Gtilde);

  Eigen::MatrixXd lin, lout;
  igl::edge_lengths(V,Ftilde,lin);
  flip_to_delaunay(Ftilde,Gtilde,lin,Fout,Gout,lout);

  Eigen::SparseMatrix<double> L,B;
  tufted_cover_laplacian(V,F,L,B);

  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
  V,v      view input mesh
  C,c      view tufted cover
  X,x      view intrinsic Delaunay triangulation
)";

  const auto show_input_mesh = [&]()
  {
    viewer.data().clear();
	  viewer.data().set_mesh(V,F);
  };

  const auto show_tufted_cover = [&]()
  {
    viewer.data().clear();
	  viewer.data().set_mesh(V,Ftilde);
  };

  const auto show_delaunay = [&]()
  {
    viewer.data().clear();
	  viewer.data().set_mesh(V,Fout);
  };

  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key,int)
  {
    switch(key)
    {
      case 'V':
      case 'v':
        show_input_mesh();
        return true;
      case 'C':
      case 'c':
      	show_tufted_cover();
        return true;
      case 'X':
      case 'x':
      	show_delaunay();
        return true;
    }
    return false;
  };

  show_input_mesh();
  viewer.launch();
}