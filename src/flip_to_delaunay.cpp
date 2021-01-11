#include "flip_to_delaunay.h"
#include <vector>
#include <iostream>

void flip_to_delaunay(
  const Eigen::MatrixXi & Fin,
  const Eigen::MatrixXi & Gin,
  const Eigen::MatrixXd & lin,
  Eigen::MatrixXi & Fout,
  Eigen::MatrixXi & Gout,
  Eigen::MatrixXd & lout)
{
  // store an edge with its two incident faces; Fout(f1,e1) == Fout(f2,e2)
  typedef struct { int f1, f2, e1, e2; } edge;
  std::vector<edge> edges_to_check; // stack or LIFO

  const int m = Fin.rows();

  Fout = Fin;
  Gout = Gin;
  lout = lin.array().square();


  auto get_edge_id = [](int f, int u, int v, const auto & Fout) {
    if (u == Fout(f,0) && v == Fout(f,1))
      return 0;
    if (u == Fout(f,1) && v == Fout(f,2))
      return 1;
    if (u == Fout(f,2) && v == Fout(f,0))
      return 2;
    std::cout << "Error flip_to_delaunay(): can't find edge!" << std::endl;
    return -1;
  };

  auto get_cot = [](double a, double b, double c) {
    double s = b + c - a;
    return s / sqrt(4*b*c - s*s);
  };


  // put all edges to stack to check
  for (int f1=0; f1<m; f1++)
  for (int e1=0; e1<3; e1++)
  {
    int f2 = Gout(f1,e1);
    int e2 = get_edge_id(f2, Fout(f1,(e1+1)%3), Fout(f1,e1), Fout);
    edges_to_check.push_back({ f1, f2, e1, e2 });
  }

  // main loop: check intrinsic Delaunay, and flip edge
  while (!edges_to_check.empty())
  {
    edge e = edges_to_check.back();
    edges_to_check.pop_back();

    // if edge is already checked (updated), skip
    if ( e.f1 == e.f2 
      || Gout(e.f1,e.e1) != e.f2
      || Gout(e.f2,e.e2) != e.f1
      || Fout(e.f1,e.e1) != Fout(e.f2,(e.e2+1)%3)
      || Fout(e.f1,(e.e1+1)%3) != Fout(e.f2,e.e2) )
      continue;

    // check intrinsic Delaunay
    double a1 = lout(e.f1,(e.e1+2)%3);
    double b1 = lout(e.f1,(e.e1+1)%3);
    double c1 = lout(e.f1,e.e1);
    double a2 = lout(e.f2,(e.e2+2)%3);
    double b2 = lout(e.f2,(e.e2+1)%3);
    double c2 = lout(e.f2,e.e2);

    if (get_cot(a1,b1,c1) + get_cot(a2,b2,c2) >= 0)
      continue;


    /***   flip edge   ***/
    // update Fout
    int i = Fout(e.f1,e.e1);
    int j = Fout(e.f1,(e.e1+1)%3);
    int k = Fout(e.f1,(e.e1+2)%3);
    int m = Fout(e.f2,(e.e2+2)%3);

    Fout(e.f1,0) = m;
    Fout(e.f1,1) = k;
    Fout(e.f1,2) = i;
    Fout(e.f2,0) = k;
    Fout(e.f2,1) = m;
    Fout(e.f2,2) = j;

    // update Gout
    int f_c = Gout(e.f1,(e.e1+1)%3);
    int f_d = Gout(e.f1,(e.e1+2)%3);
    int f_e = Gout(e.f2,(e.e2+1)%3);
    int f_f = Gout(e.f2,(e.e2+2)%3);

    Gout(e.f1,0) = e.f2;
    Gout(e.f1,1) = f_d;
    Gout(e.f1,2) = f_e;
    Gout(e.f2,0) = e.f1;
    Gout(e.f2,1) = f_f;
    Gout(e.f2,2) = f_c;
    Gout(f_c,get_edge_id(f_c,k,j,Fout)) = e.f2;
    Gout(f_e,get_edge_id(f_e,m,i,Fout)) = e.f1;

    // update l
    double l_ij = lout(e.f1,(e.e1+2)%3);
    double l_ik = lout(e.f1,(e.e1+1)%3);
    double l_kj = lout(e.f1,e.e1);
    double l_im = lout(e.f2,e.e2);
    double l_jm = lout(e.f2,(e.e2+1)%3);

    // compute new edge length
    double s, t;
    s = l_ij + l_ik - l_kj;
    t = l_ij * l_ik;
    double cos_a = 0.5*s/sqrt(t);
    double sin_a = 0.5*sqrt((4*t - s*s)/t);
    s = l_ij + l_im - l_jm;
    t = l_ij * l_im;
    double cos_b = 0.5*s/sqrt(t);
    double sin_b = 0.5*sqrt((4*t - s*s)/t);
    double cos = cos_a*cos_b - sin_a*sin_b;
    double l_km = l_im + l_ik - 2*sqrt(l_im*l_ik)*cos;

    lout(e.f1,0) = l_ik;
    lout(e.f1,1) = l_im;
    lout(e.f1,2) = l_km;
    lout(e.f2,0) = l_jm;
    lout(e.f2,1) = l_kj;
    lout(e.f2,2) = l_km;

    // put edges to stack
    if (e.f1 != (f_d=Gout(e.f1,1)))
      edges_to_check.push_back({ e.f1, f_d, 1, get_edge_id(f_d,i,k,Fout) });
    if (e.f1 != (f_e=Gout(e.f1,2)))
      edges_to_check.push_back({ e.f1, f_e, 2, get_edge_id(f_e,m,i,Fout) });
    if (e.f2 != (f_f=Gout(e.f2,1)))
      edges_to_check.push_back({ e.f2, f_f, 1, get_edge_id(f_f,j,m,Fout) });
    if (e.f2 != (f_c=Gout(e.f2,2)))
      edges_to_check.push_back({ e.f2, f_c, 2, get_edge_id(f_c,k,j,Fout) });
  }

  lout = lout.cwiseSqrt();
}