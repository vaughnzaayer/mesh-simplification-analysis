#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/quadric_error_simplification.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"

#include "geometrycentral/surface/direction_fields.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <Eigen/Dense>
#include <fmt/base.h>
#include <fmt/format.h>


using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<RichSurfaceMeshData> rich_data;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Structs
struct QueueEntry {
  double cost;
  std::pair<Vertex, Vertex> pair;
  Vector3 target;
  bool operator>(const QueueEntry& other) const { return cost > other.cost; }
};

struct MeshIteration {
  int iter;
  ManifoldSurfaceMesh mesh;
  VertexPositionGeometry geometry;
  bool operator>(const MeshIteration& other) const { return iter > other.iter; }
};

// Some algorithm parameters
float param1 = 1.0;
double tolerance = 1.0;
size_t target;

// Meta information
int qed_iter_count = 0;
bool recordOutput = true;
std::string outPath = "~/MeshAnalysis/QEDIterations/";

// Create containers for mesh data
std::unique_ptr<VertexData<Eigen::Matrix4d>> vertexQuadrics;



// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
                     std::to_string(param1));

  geometry->requireVertexGaussianCurvatures();
  psMesh->addVertexScalarQuantity("curvature",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);
}

void countVertices() {
  for (Vertex v : mesh->vertices()) {
    std::cout << "Vertex " << v << " has degree " << v.degree() << std::endl;
  }
}

//TODO: Write a separate function for computing the optimal target and cost of each vertex pair
// Vector3 computeTargetVertex(std::pair<Vertex, Vertex> p) {
//   Eigen::Vector4d tmp;
//   Vector3 target;
//   Eigen::Matrix4d Q = (*vertexQuadrics)[p.first] + (*vertexQuadrics)[p.second];
//   if (std::abs(Q.determinant()) > 1e-9) {
//     tmp = Q.inverse() * Eigen::Vector4d(0,0,0,1);
//     target = Vector3{tmp[0], tmp[1], tmp[2]};
//   } else {
//     target = (geometry->vertexPositions[p.first] + geometry->vertexPositions[p.second]) * 0.5;
//   }
//   return target;
// }

// double computeCost(std::pair<Vertex, Vertex> p, Vector3 v) {
//   Eigen::Matrix4d Q = (*vertexQuadrics)[p.first] + (*vertexQuadrics)[p.second];
//   Eigen::Vector4d tmp = Eigen::Vector4d(v[0], v[1], v[2], 1);
//   return (tmp.transpose() * Q * tmp)(0,0);
// }

// std::vector<std::vector<size_t>> getNonManifoldFaces(Vertex v1, Vertex v2) {
//   std::vector<std::vector<size_t>> newFaceIndices;
//   for (Face f : mesh->faces()) {
//     std::vector<size_t> faceVerts;
//     bool isDegenerate = false;
    
//     for (Vertex fv : f.adjacentVertices()) {
//         size_t idx = (fv == v2) ? v1.getIndex() : fv.getIndex();
//         faceVerts.push_back(idx);
//     }
    
//     if (faceVerts.size() == 3) {
//         if (faceVerts[0] == faceVerts[1] || 
//             faceVerts[1] == faceVerts[2] || 
//             faceVerts[2] == faceVerts[0]) {
//             isDegenerate = true; 
//         }
//     }
    
//     if (!isDegenerate) {
//         newFaceIndices.push_back(faceVerts);
//     }
//   }
//   return newFaceIndices;
// }

// void collapsePair(QueueEntry p) {
//   Vertex v1 = p.pair.first;
//   Vertex v2 = p.pair.second;

//   // Compute new connectivity
//   std::vector<std::vector<size_t>> newFaces = getNonManifoldFaces(v1, v2);

//   // Move v1 to target position
//   geometry->vertexPositions[v1] = p.target;

//   // Make a new mesh with potentially non-manifold topology
//   std::unique_ptr<SurfaceMesh> newMesh;
//   std::unique_ptr<VertexPositionGeometry> newGeo;
  
//   std::vector<Vector3> positions;
//   positions.resize(mesh->nVerticesCapacity());

//   for (Vertex v : mesh->vertices()) {
//     positions[v.getIndex()] = geometry->inputVertexPositions[v];
//   }

//   std::tie(newMesh, newGeo) = makeSurfaceMeshAndGeometry(newFaces, positions);
// }



// void QEMsimplification() {
//   geometry->requireFaceNormals();

//   // Compute Q for each vertex 
//   for (Vertex v : mesh->vertices()) {
//     Eigen::Matrix4d quad = Eigen::Matrix4d::Zero();
//     for (Face f : v.adjacentFaces()) {
//       Vector3 normal = geometry->faceNormals[f];
//       double a = normal[0];
//       double b = normal[1];
//       double c = normal[2];
//       double d = dot(normal, geometry->vertexPositions[v]);
//       Eigen::Vector4d plane(a, b, c, d);
//       quad += (plane * plane.transpose());
//     }
//     (*vertexQuadrics)[v] = quad;
//   }

//   // Find all valid pairs
//   std::vector<std::pair<Vertex, Vertex>> valid;
//   for (Vertex v : mesh->vertices()) {
//     for (Vertex u : mesh->vertices()) {
//       if (u != v) {
//         double dist = norm(geometry->vertexPositions[v] - geometry->vertexPositions[u]);
//         for (Vertex nbr : v.adjacentVertices()) {
//           if (nbr == u) {valid.push_back(std::pair<Vertex, Vertex>(u,v));}
//         }
//         if (dist <= tolerance) {valid.push_back(std::pair<Vertex, Vertex>(u,v));}
//       }
//     }
//   }

//   // Compute optimal contraction target and cost, add to min. heap
//   std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> heap;
//   for (std::pair<Vertex, Vertex> p : valid) {
//     Vector3 v = computeTargetVertex(p);
//     double cost = computeCost(p, v);
//     heap.push(QueueEntry{cost, p, v});
    
//   }

//   // Keep track of the different resolutions

//   // Simplificiation loop
//   while (mesh->nVertices() < target && !heap.empty()) {
//     // Pop the top element from the queue
//     QueueEntry top = heap.top();
//     heap.pop();

//     // Check 
    




//   }



// }

// Save current mesh data
void saveMesh() {
  // Save obj data
  auto name = psMesh->getName();
  std::filesystem::create_directories(outPath);
  writeSurfaceMesh(*mesh, *geometry, fmt::format("{0}{1}_qed_iteration_{2}.obj", outPath, name, qed_iter_count));
}



// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  if (ImGui::Button("do work")) {
    doWork();
  }

  ImGui::SliderFloat("QED tolerance", &param1, 0., 5.);

  if (ImGui::Button("Run QED Simplification")) {
    if (recordOutput) { saveMesh(); }
    qed_iter_count += 1;
    quadricErrorSimplify(*mesh, *geometry, static_cast<double>(param1));
    auto tmp_ps_mesh = polyscope::registerSurfaceMesh(fmt::format("QED Mesh Iteration {}", qed_iter_count), 
      geometry->vertexPositions, 
      mesh->getFaceVertexList());
  }
  
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
  args::Positional<bool> recordOutput(parser, "record output", "Record output to \'~/MeshAnalysis/QEDIterations/\'.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();  

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));
  vertexQuadrics = std::make_unique<VertexData<Eigen::Matrix4d>>(*mesh);

  // Save first iteration of the mesh

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  VertexData<Vector3> vBasisY(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
    vBasisY[v] = geometry->vertexTangentBasis[v][1];
  }

  auto vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
