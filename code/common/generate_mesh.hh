#pragma once

#include <aka_common.hh>
#include <mesh.hh>

/// @brief  Generate a 1D mesh with N elements between x0 and xf, with a total
/// length of length.
/// @param length length of the mesh
/// @param x0 right boundary of the mesh
/// @param xf left boundary of the mesh
/// @param N number of triangular elements in the mesh
/// @return mesh object
void generateMesh(akantu::Mesh& mesh,
                  akantu::Real length,
                  akantu::Real x0,
                  akantu::Real xf,
                  akantu::UInt N);