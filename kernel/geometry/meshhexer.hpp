// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_MESHHEXER) || defined(DOXYGEN)

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/string.hpp>

namespace FEAT::Geometry
{
  struct Gap
  {
    double diameter;
    double confidence;

    Gap(double d, double c) : diameter(d), confidence(c)
    {
    }
  };

  class MeshHexerSurfaceMesh
  {
    void* _surface_mesh;

    friend MeshHexerSurfaceMesh load_from_file(const String& filename, bool triangulate);

    explicit MeshHexerSurfaceMesh(void* sm) noexcept : _surface_mesh(sm)
    {
    }

  public:
    ~MeshHexerSurfaceMesh();

    MeshHexerSurfaceMesh(const MeshHexerSurfaceMesh&) = delete;
    MeshHexerSurfaceMesh& operator=(const MeshHexerSurfaceMesh&) = delete;

    MeshHexerSurfaceMesh(MeshHexerSurfaceMesh&&) = default;
    MeshHexerSurfaceMesh& operator=(MeshHexerSurfaceMesh&&) = default;

    /**
     * \brief Check if the surface mesh describes a closed surface
     *
     * \returns True, if there are no border edges in the surface mesh
     */
    bool is_closed() const;

    /**
     * \brief Check if the winding order of the faces of the surface mesh is consistent
     *
     * \return True, if either all faces of the surface mesh are wound clockwise or all faces of the surface mesh are
     * wound counter-clockwise
     */
    bool is_wound_consistently() const;

    /**
     * \brief Check if the surface normals of all faces of the surface mesh point towards the unbounded side
     *
     * \pre The surface mesh is closed, call SurfaceMesh::is_closed() to check
     * \pre The surface mesh is wound consistently, call SurfaceMesh::is_wound_consistently() to check
     *
     * \returns True, if the surface normals of all faces of the surface mesh point towards the unbounded side, i.e.
     * outside.
     */
    bool is_outward_oriented() const;

    /**
     * \brief Compute the smallest aspect ratio of any face of the surface mesh
     *
     * \returns The minimal aspect ratio of any face of the surface mesh
     */
    double minimal_aspect_ratio() const;

    /**
     * \brief Compute the largest aspect ratio of any face of the surface mesh
     *
     * \returns The maximum aspect ratio of any face of the surface mesh
     */
    double maximal_aspect_ratio() const;

    /**
     * \brief Compute Gap candidates of the surface mesh.
     *
     * Gaps try to capture the diameter of the volume described by the surface mesh.
     * The diameter is measured from the barycentric centers of all faces of the surface mesh.
     * For each face of the surface mesh the largest sphere that just touches the center of
     * the face and one other point of the surface mesh is computed.
     * We call this sphere the maximal inscribed sphere (MIS) and use the diameter of that sphere
     * as the diameter at that face.
     *
     * \verbatim
     * │                 │
     * │     ∙∙∙∙∙∙∙     │
     * │  ∙∙∙       ∙∙∙  │
     * │ ∙             ∙ │
     * │∙    Diameter   ∙│
     * │∙---------------∙│
     * │∙               ∙│
     * │ ∙             ∙ │
     * │  ∙∙∙       ∙∙∙  │
     * │     ∙∙∙∙∙∙∙     │
     * │                 │
     * \endverbatim
     *
     * Note that not all these MIS span gaps that correspond to the intuitive diameter of the
     * surface mesh at that point. Badly reconstructed meshes for example might contain many
     * small "valleys" on surfaces that are supposed to be "smooth"
     *
     * \verbatim
     *              ∙∙∙∙∙∙∙
     *            ∙∙       ∙∙
     *\         ∙           ∙         /
     *  \\      ∙             ∙      //
     *    \\    ∙             ∙    //
     *      \\  ∙             ∙  //
     *        \\ ∙           ∙ //
     *          \\∙∙       ∙∙//
     *            \\∙∙∙∙∙∙∙//
     *              \\   //
     *                \ /
     * \endverbatim
     *
     * Each gap is given a confidence score in the range [0, 1] that shows how likely that gap
     * corresponds to a real diameter. That score is determined from:
     * - the aspect ratios of involved faces
     * - self-intersections of the surface mesh
     * - relative sizes of involved faces
     * - similarity of normals at the touching points of the MIS
     * - topological distance (distance across the surface between the touching points of the MIS)
     * - dihedral angles near the touching points of the MIS
     *
     * Generally gaps with a score above 0.95 have a good change of corresponding to real diameters.
     *
     * \returns A vector containing one gap candiate for each face of the surface mesh
     */
    std::vector<Gap> gaps();

    /**
     * \brief Compute the smallest Gap of the surface mesh
     *
     * See SurfaceMesh::gaps() for details on how gaps area measured and scored.
     *
     * \returns The smallest gap, by diameter, of all gaps with a score above 0.95,
     * or the smallest gap among the top 10 percent of all gaps, if no such gaps exist.
     */
    Gap min_gap();

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] x_min Minimum x coordinate of bounding box
     * \param[in] y_min Minimum y coordinate of bounding box
     * \param[in] z_min Minimum z coordinate of bounding box
     * \param[in] x_max Maximum x coordinate of bounding box
     * \param[in] y_max Maximum y coordinate of bounding box
     * \param[in] z_max Maximum z coordinate of bounding box
     * \param[in] levels Number of levels in intended mesh hierarchy
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     */
    template<typename DataType_>
    ConformalMesh<Shape::Hexahedron, 3, DataType_> fbm_mesh(
      DataType_ x_min,
      DataType_ y_min,
      DataType_ z_min,
      DataType_ x_max,
      DataType_ y_max,
      DataType_ z_max,
      std::uint64_t levels);

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] levels Number of levels in intended mesh hierarchy
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     */
    template<typename DataType_>
    ConformalMesh<Shape::Hexahedron, 3, DataType_> fbm_mesh(std::uint64_t levels);

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] x_min Minimum x coordinate of bounding box
     * \param[in] y_min Minimum y coordinate of bounding box
     * \param[in] z_min Minimum z coordinate of bounding box
     * \param[in] x_max Maximum x coordinate of bounding box
     * \param[in] y_max Maximum y coordinate of bounding box
     * \param[in] z_max Maximum z coordinate of bounding box
     * \param[in] levels Number of levels in intended mesh hierarchy
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     *
     * \returns A conformal base mesh and suggested refinement levels
     */
    template<typename DataType_>
    std::pair<ConformalMesh<Shape::Hexahedron, 3, DataType_>, SubdivisionLevels> prepare_fbm_mesh(
      DataType_ x_min,
      DataType_ y_min,
      DataType_ z_min,
      DataType_ x_max,
      DataType_ y_max,
      DataType_ z_max,
      std::uint64_t levels);

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] levels Number of levels in intended mesh hierarchy
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     *
     * \returns A conformal base mesh and suggested refinement levels
     */
    template<typename DataType_>
    std::pair<ConformalMesh<Shape::Hexahedron, 3, DataType_>, SubdivisionLevels> prepare_fbm_mesh(std::uint64_t levels);

    /**
     * \brief Write the surface mesh to disk
     *
     * \param[in] filename Filename to write to. Must end in .ply or .vtu
     *
     * Writes the surface mesh to disk as a .ply or .vtu file. The written file
     * contains mesh properties that have been calculated as intermediate results,
     * such as maximal inscribed spheres, or topological distances.
     * These properties will be reused if the written mesh is read again.
     *
     * \returns A result indicating success or containing an error message
     */
    void write_to_file(const std::string& filename);
  }; // MeshHexerSurfaceMesh

  /**
   * \brief Load a surface mesh from a file
   *
   * \param[in] filename     Path to mesh file
   * \param[in] triangulate  If set, input mesh will automatically be triangulated
   *
   * Supported file formats are
   * - .off
   * - .obj
   * - .stl
   * - .ply
   * - .ts
   * - .vtp
   * - .vtu
   */
  MeshHexerSurfaceMesh load_from_file(const String& filename, bool triangulate = false);
} // namespace FEAT::Geometry

#endif
