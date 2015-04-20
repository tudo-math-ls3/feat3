#pragma once
#ifndef KERNEL_GEOMETRY_MESH_SMOOTHER_HPP
#define KERNEL_GEOMETRY_MESH_SMOOTHER_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/lafem/dense_vector.hpp>

namespace FEAST
{
  namespace Geometry
  {

    /**
     * \brief Baseclass for mesh optimisation algorithms
     *
     * This abstract class is the baseclass for all mesh optimisation algorithms, which can be direct, variational
     * or something entirely different.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \author Jordi Paul
     *
     */
    template<typename DataType_, typename MemType_, typename TrafoType_>
    class MeshSmoother
    {
      public:
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVector<MemType, DataType> VectorType;

        /// Coordinates of the vertices, as they get changed in the optimisation process
        VectorType _coords[MeshType::world_dim];

      public:
        /// The underlying transformation
        const TrafoType& _trafo;
        /// The mesh for the underlying transformation
        MeshType& _mesh;

      public:
        /// Constructor
        explicit MeshSmoother(const TrafoType& trafo_) :
          _trafo (trafo_),
          _mesh ( const_cast<MeshType&>(_trafo.get_mesh()) )
          {
            for(int d = 0; d < MeshType::world_dim; ++d)
              _coords[d]= std::move(VectorType(_mesh.get_num_entities(0)));

          }

        virtual ~MeshSmoother()
        {
        }

        /// \brief Initialises parts of the MeshSmoother not set in in the constructor
        virtual void init()
        {
          get_coords();
        }

        /// \brief Gets the coordinates from the underlying mesh and saves them in _coords.
        virtual void get_coords()
        {
          const typename MeshType::VertexSetType& vertex_set = _mesh.get_vertex_set();

          for(Index i(0); i < _mesh.get_num_entities(0); ++i)
          {
            for(int d(0); d < MeshType::world_dim; ++d)
              _coords[d](i,DataType(vertex_set[i][d]));
          }
        }

        /// \brief Sets the coordinates in the underlying mesh to _coords.
        virtual void set_coords()
        {
          typename MeshType::VertexSetType& vertex_set = _mesh.get_vertex_set();

          for(Index i(0); i < _mesh.get_num_entities(0); ++i)
          {
            for(int d(0); d < MeshType::world_dim; ++d)
              vertex_set[i][d] = _coords[d](i);
          }
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() = 0;

        /// \brief Prepares the mesh optimiser for application
        virtual void prepare() = 0;

    }; // class MeshSmoother


  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_MESH_SMOOTHER_HPP
