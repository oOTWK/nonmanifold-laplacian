#include "tufted_cover.h"
#include <igl/per_face_normals.h>
#include <vector>

void tufted_cover(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & Ftilde,
  Eigen::MatrixXi & Gtilde)
{
  struct face {int index, side, orientation;};
  struct edge_faces {
    struct f_elmt {face f[2]; double sort_key;};
    std::vector<f_elmt> faces;
    Eigen::Vector3d edge_vector;
  };
  std::vector<edge_faces> edges;

  const int n = V.rows();
  const int m = F.rows();

  Ftilde.resize(2*m,3);
  Gtilde.resize(2*m,3);

  Ftilde.topLeftCorner(m,3) = F;
  Ftilde.col(0).tail(m) = F.col(1);
  Ftilde.col(1).tail(m) = F.col(0);
  Ftilde.col(2).tail(m) = F.col(2);

  Eigen::MatrixXd N;
  igl::per_face_normals(V,F,N);

  Eigen::MatrixXi eid = Eigen::MatrixXi::Zero(n,n);

  for (int c=0, x=0; x<m; x++)
  for (int y=0; y<3; y++)
  {
    int i = F(x,y);
    int j = F(x,(y+1)%3);

    if (i > j)
      std::swap(i,j);

    if (eid(i,j) == 0)
    {
      eid(i,j) = ++c;
      edges.push_back({{}, V.row(j)-V.row(i)});
    }
  }

  // collect incident faces for each edge
  for (int x=0; x<m; x++)
  {
    int f_ind[2] = {x,x+m};
    for (int y=0; y<3; y++)
    {
      int i = F(x,y);
      int j = F(x,(y+1)%3);
      int side[2] {y,y};
      if (y == 1)
        side[1] = 2;
      else if (y == 2)
        side[1] = 1;

      int ori;
      if (i > j)
      {
        std::swap(i,j);
        ori = 1;
      }
      else ori = 0;

      struct edge_faces & ef = edges[eid(i,j)-1];

      // compute sort_key for cyclice ordering of faces
      double sort_key;
      if (ef.faces.empty())
        sort_key = 0;
      else
      {
        Eigen::Vector3d first_n = N.row((ef.faces[0].f[0].index)%3);
        Eigen::Vector3d n = N.row(x);

        sort_key = first_n.dot(n);

        if ((first_n.cross(n)).dot(ef.edge_vector) > 0)
          sort_key += 2.0;
        else
          sort_key *= -1;
      }

      ef.faces.push_back({{{f_ind[ori], side[ori], 0},
                          {f_ind[(ori+1)%2], side[(ori+1)%2], 1}},
                          sort_key});
    }
  }


  auto sort_func = [](const auto & f1, const auto & f2) { return (f1.sort_key > f2.sort_key); };

  // glue faces in circular order
  for (auto e : edges)
  {
    std::sort(e.faces.begin()+1, e.faces.end(), sort_func); 

    face * f = & (e.faces[0].f[0]);
    for (int i=1; i<e.faces.size(); i++)
    {
      face * g1 = & e.faces[i].f[0];
      face * g2 = & e.faces[i].f[1];

      if (f->orientation == g1->orientation)
        std::swap(g1,g2);

      Gtilde(f->index,f->side) = g1->index;
      Gtilde(g1->index,g1->side) = f->index;
      f = g2;
    }

    face * g = & e.faces[0].f[1];
    Gtilde(f->index,f->side) = g->index;
    Gtilde(g->index,g->side) = f->index;
  }
}

