#include "generate_mesh.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace akantu;

void generateMesh(Mesh& mesh, Real length, Real x0, Real xf, UInt N) {
  Real h_uniform = length / (static_cast<Real>(N) * 0.5);

  // write the gmsh file
  std::ofstream geo_file("bar.geo");
  if (!geo_file.is_open()) {
    throw std::runtime_error("Failed to open file for writing: bar.geo");
  }
  geo_file << "Point(1) = { " << x0 << ", 0, 0, " << h_uniform << " };\n"

           << "Point(2) = { " << xf << ", 0, 0, " << h_uniform << " };\n"

           << "Point(3) = { " << xf << ", " << h_uniform << ", 0, " << h_uniform
           << " };\n"

           << "Point(4) = { " << x0 << ", " << h_uniform << ", 0, " << h_uniform
           << " };\n\n"

           << "Line(1) = {1,2};\n"
           << "Line(2) = {2,3};\n"
           << "Line(3) = {3,4};\n"
           << "Line(4) = {4,1};\n\n"

           << "Line Loop(1) = {1,2,3,4};\n"
           << "Plane Surface(1) = {1};\n"

           << "Physical Surface(1) = {1};\n"
           << "Physical Line(\"Yblocked\") = {1,3};\n"
           << "Physical Line(\"right\") = {2};\n"
           << "Physical Line(\"left\") = {4};\n\n"

           << "Transfinite Surface {1} Left;\n";

  geo_file.close();

  // Generate the mesh using Gmsh
  int return_code = std::system("gmsh -2 -order 1 -o bar.msh bar.geo");
  if (return_code != 0) {
    throw std::runtime_error("Gmsh failed to generate the mesh");
  }
  std::cout << "Info    : Mesh generated" << std::endl;

  mesh.read("bar.msh");
}