/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 * Generally, Peter wanted to take a deeper look at the base mesh implementation.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_BM_CELL_DATA_VALIDATION_HPP
#define KERNEL_BM_CELL_DATA_VALIDATION_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/cell_data.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/vertex.hpp>

/// FEAST namespace
namespace FEAST
{
  /// BaseMesh namespace comprising base mesh specific code
  namespace BaseMesh
  {
    /**
    * \brief class for validating vertex neighbours
    *
    * \tparam cell_space_dim_
    * cell and space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in
    * a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note
    * Q: Why are there three classes CellDataValidationVertNeigh, CellDataValidationEdgeNeigh and
    *    CellDataValidationFaceNeigh? Can they not be united (the code is quite similar)?
    * A: 1) That is a consequence of that fact, that we explicitely store vertices, edges and faces instead of
    *    putting all in one array of general "items".
    *    2) If we only had to compare subcells here, we maybe could use the most general class BaseMeshItem (being
    *    parent of vertices, edges and faces). However, we also have to access children/parents of those subcells,
    *    and that is not possible if we only handle them as BaseMeshItem. (Note that for vertices we don't need
    *    parent/child information.)
    *    3) The code differs for vertices on the one hand and edges/faces on the other hand: For the latter we
    *    eventually have to traverse the history (parent/child relations) for validating neighbourhood, which is not
    *    necessary for vertices.
    *    4) Constructs like 'if(cell_space_dim == 1) then ... else if(cell_space_dim == 2) ...' are not possible:
    *    For cell_space_dim = 3, e.g., one would use the function num_faces() which is simply not available for
    *    cell_space_dim < 3, i.e. the code would not compile.
    *    5) Of course, one could modify the whole class concept, such that the code duplication below can be avoided,
    *    but that make the whole thing even more complicated than it already is. So, my current approach is to accept
    *    the code duplication below rather than blowing up the class hierarchy/relations or using too abstract concepts
    *    for handling/storing subcells.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellDataValidationVertNeigh
    {

    private:

    public:

      /**
      * \brief validates whether vertex neighbours are set correctly
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c, std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidationVertNeigh::validate()");
        ASSERT(cell_space_dim_ >= 1, "Cell + space dimension " + stringify(cell_space_dim_) + " must be at least 1.");
        try
        {
          std::string s;
          for(unsigned char ivertex(0) ; ivertex < c->num_vertices() ; ++ivertex)
          {
            // get vertex neighbours at vertex ivertex
            std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*> neighs
              = c->neighbours_item(SDIM_VERTEX, ivertex);
            for(unsigned int ineigh = 0 ; ineigh < neighs.size() ; ineigh++)
            {
              // check whether ineigh-th neighbour is active
              if(!neighs[ineigh]->active())
              {
                // if not, throw an error
                s = "Cell " + neighs[ineigh]->print_index() + ", being the vertex-neighbour of cell " + c->print_index()
                    + " (at vertex " + stringify((int)ivertex) +  ", " + stringify(ineigh)
                    + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: Deactivated until update of neighbourhood has been implemented! Instead, only print error message.
//                throw new InternalError(s);
                std::cerr << "ERROR: " << s;
                stream << "ERROR: " << s;
              }
              else // neighs[ineigh]->active()
              {
                // Get the vertex neighbours of the current neighbour.
                // This is an array of vectors, where the ivertex_in_neigh-th entry of the array is the vector of
                // neighbours at the ivertex_in_neigh-th vertex.
                std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>*
                  neighs_of_neigh = neighs[ineigh]->neighbours_subdim(SDIM_VERTEX);

                // init flag indicating whether neighbour has been found
                bool neighbour_found = false;
                for(unsigned char ivertex_in_neigh = 0 ; ivertex_in_neigh < neighs[ineigh]->num_vertices() ;
                    ivertex_in_neigh++)
                {
                  for(unsigned int nn = 0 ; nn < neighs_of_neigh[ivertex_in_neigh].size() ; nn++)
                  {
                    if(neighs_of_neigh[ivertex_in_neigh][nn] == c)
                    {
                      // debug output
                      s = "Vertex-neighbourhood between cell " + c->print_index() + " (at vertex "
                          + stringify((int)ivertex) + ", " + stringify(ineigh) + "-th pos.) and cell "
                          + neighs[ineigh]->print_index() + " (vertex = " + stringify((int)ivertex_in_neigh)
                          + ", " + stringify(nn) + "-th pos.) found!\n";
                      stream << s;
                      neighbour_found = true;

                      // check whether the two neighbours really share the same vertex
                      if (c->vertex(ivertex) != neighs[ineigh]->vertex(ivertex_in_neigh))
                      {
                        s = "Vertex-neighbours cell " + c->print_index() + " (at vertex " + stringify((int)ivertex)
                            + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                            + " (vertex = " + stringify((int)ivertex_in_neigh) + ", " + stringify(nn)
                            + "-th pos.) do not share the same vertex! Aborting program!\n";
                        throw new InternalError(s);
                      }
                      break;
                    }
                  }  // for 0 <= nn < neighs_of_neigh[ivertex_in_neigh].size()
                  if(neighbour_found)
                  {
                    break;
                  }
                } // for 0 <= ivertex_in_neigh < neighs[ineigh]->num_vertices()
                if(!neighbour_found)
                {
                  s = "No vertex-neighbour between cell " + c->print_index() + " (at vertex " + stringify((int)ivertex)
                      + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                      + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= ivertex < c->num_vertices()
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      } // validate(...)
    }; // class CellDataValidationVertNeigh


    /**
    * \brief class for validating edge neighbours
    *
    * \tparam cell_space_dim_
    * cell and space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in
    * a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note See note in class CellDataValidationVertNeigh.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellDataValidationEdgeNeigh
    {

    private:

    public:

      /**
      * \brief validates whether edge neighbours are set correctly
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c, std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidationEdgeNeigh::validate()");
        ASSERT(cell_space_dim_ >= 2, "Cell + space dimension " + stringify(cell_space_dim_) + " must be at least 2.");
        try
        {
          std::string s;
          for(unsigned char iedge(0) ; iedge < c->num_edges() ; ++iedge)
          {
            // get edge neighbours at edge iedge
            std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*> neighs
              = c->neighbours_item(SDIM_EDGE, iedge);
            for(unsigned int ineigh = 0 ; ineigh < neighs.size() ; ineigh++)
            {
              // check whether ineigh-th neighbour is active
              if(!neighs[ineigh]->active())
              {
                // if not, throw an error
                s = "Cell " + neighs[ineigh]->print_index() + ", being the edge-neighbour of cell "
                    + c->print_index() + " (at edge " + stringify((int)iedge) +  ", " + stringify(ineigh)
                    + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: Deactivated until update of neighbourhood has been implemented! Instead, only print error message.
//                throw new InternalError(s);
                std::cerr << "ERROR: " << s;
                stream << "ERROR: " << s;
              }
              else // neighs[ineigh]->active()
              {
                // Get the edge neighbours of the current neighbour.
                // This is an array of vectors, where the iedge_in_neigh-th entry of the array is the vector of
                // neighbours at the iedge_in_neigh-th edge.
                std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>*
                  neighs_of_neigh = neighs[ineigh]->neighbours_subdim(SDIM_EDGE);

                // init flag indicating whether neighbour has been found
                bool neighbour_found = false;
                for(unsigned char iedge_in_neigh = 0 ; iedge_in_neigh < neighs[ineigh]->num_edges() ; iedge_in_neigh++)
                {
                  for(unsigned int nn = 0 ; nn < neighs_of_neigh[iedge_in_neigh].size() ; nn++)
                  {
                    if(neighs_of_neigh[iedge_in_neigh][nn] == c)
                    {
                      // debug output
                      s = "Edge-neighbourhood between cell " + c->print_index() + " (at edge " + stringify((int)iedge)
                          + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                          + " (edge = " + stringify((int)iedge_in_neigh) + ", " + stringify(nn) + "-th pos.) found!\n";
                      stream << s;
                      neighbour_found = true;

                      // We have to check whether the refinement levels of the two edges differ (which means, that the
                      // one with the smaller ref. level is the parent of the other one). If this is the case, then we
                      // have to "go up" in the hierarchy of the edge with the greater ref. level until the levels (and
                      // the edges) are equal.
                      Cell<1, cell_space_dim_, world_dim_>* edge_c = c->edge(iedge);
                      Cell<1, cell_space_dim_, world_dim_>* edge_neigh = neighs[ineigh]->edge(iedge_in_neigh);
                      char level_diff = edge_c->refinement_level() - edge_neigh->refinement_level();
                      while(level_diff < 0)
                      {
                        // Ref. level of c's edge is smaller than that of the neighbour's edge.
                        // Get parent of the latter.
                        edge_neigh = edge_neigh->parent();
                        level_diff++;
                      }
                      while(level_diff > 0)
                      {
                        // Ref. level of c's edge is greater than that of the neighbour's edge.
                        // Get parent of the former.
                        edge_c = edge_c->parent();
                        level_diff--;
                      }
                      ASSERT(level_diff == 0, "Level difference " + stringify(level_diff) + " must be 0.");
                      ASSERT(edge_c->refinement_level() == edge_neigh->refinement_level(), "Edge ref. level "
                             + stringify(edge_c->refinement_level()) + " must equal ref. level "
                             + stringify(edge_neigh->refinement_level()) + " of neighbour edge.");
                      // check whether the two neighbours really share the same edge
                      if (edge_c != edge_neigh)
                      {
                        s = "Edge-neighbours cell " + c->print_index() + " (at edge " + stringify((int)iedge)
                            + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                            + " (edge = " + stringify((int)iedge_in_neigh) + ", " + stringify(nn)
                            + "-th pos.) do not share the same edge! Aborting program!\n";
                        throw new InternalError(s);
                      }
                      break;
                    }
                  }  // for 0 <= nn < neighs_of_neigh[iedge_in_neigh].size()
                  if(neighbour_found)
                  {
                    break;
                  }
                } // for 0 <= iedge_in_neigh < neighs[ineigh]->num_edges()
                if(!neighbour_found)
                {
                  s = "No edge-neighbour between cell " + c->print_index() + " (at edge " + stringify((int)iedge)
                      + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                      + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= iedge < c->num_edges()
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      } // validate(...)
    }; // class CellDataValidationEdgeNeigh


    /**
    * \brief class for validating face neighbours
    *
    * \tparam cell_space_dim_
    * cell and space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in
    * a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \note See note in class CellDataValidationVertNeigh.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellDataValidationFaceNeigh
    {

    private:

    public:

      /**
      * \brief validates whether face neighbours are set correctly
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c, std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidationFaceNeigh::validate()");
        ASSERT(cell_space_dim_ >= 3, "Cell + space dimension " + stringify(cell_space_dim_) + " must be at least 3.");
        try
        {
          std::string s;
          for(unsigned char iface(0) ; iface < c->num_faces() ; ++iface)
          {
            // get face neighbours at face iface
            std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*> neighs
              = c->neighbours_item(SDIM_FACE, iface);
            for(unsigned int ineigh = 0 ; ineigh < neighs.size() ; ineigh++)
            {
              // check whether ineigh-th neighbour is active
              if(!neighs[ineigh]->active())
              {
                // if not, throw an error
                s = "Cell " + neighs[ineigh]->print_index() + ", being the face-neighbour of cell " + c->print_index()
                    + " (at face " + stringify((int)iface) +  ", " + stringify(ineigh)
                    + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: Deactivated until update of neighbourhood has been implemented! Instead, only print error message.
//                throw new InternalError(s);
                std::cerr << "ERROR: " << s;
                stream << "ERROR: " << s;
              }
              else // neighs[ineigh]->active()
              {
                // Get the face neighbours of the current neighbour.
                // This is an array of vectors, where the iface_in_neigh-th entry of the array is the vector of
                // neighbours at the iface_in_neigh-th face.
                std::vector<Cell<cell_space_dim_, cell_space_dim_, world_dim_>*>*
                  neighs_of_neigh = neighs[ineigh]->neighbours_subdim(SDIM_FACE);

                // init flag indicating whether neighbour has been found
                bool neighbour_found = false;
                for(unsigned char iface_in_neigh = 0 ; iface_in_neigh < neighs[ineigh]->num_faces() ; iface_in_neigh++)
                {
                  for(unsigned int nn = 0 ; nn < neighs_of_neigh[iface_in_neigh].size() ; nn++)
                  {
                    if(neighs_of_neigh[iface_in_neigh][nn] == c)
                    {
                      // debug output
                      s = "Face-neighbourhood between cell " + c->print_index() + " (at face " + stringify((int)iface)
                          + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                          + " (face = " + stringify((int)iface_in_neigh) + ", " + stringify(nn) + "-th pos.) found!\n";
                      stream << s;
                      neighbour_found = true;

                      // We have to check whether the refinement levels of the two faces differ (which means, that the
                      // one with the smaller ref. level is the parent of the other one). If this is the case, then we
                      // have to "go up" in the hierarchy of the face with the greater ref. level until the levels (and
                      // the faces) are equal.
                      Cell<2, cell_space_dim_, world_dim_>* face_c = c->face(iface);
                      Cell<2, cell_space_dim_, world_dim_>* face_neigh = neighs[ineigh]->face(iface_in_neigh);
                      char level_diff = face_c->refinement_level() - face_neigh->refinement_level();
                      while(level_diff < 0)
                      {
                        // Ref. level of c's face is smaller than that of the neighbour's face.
                        // Get parent of the latter.
                        face_neigh = face_neigh->parent();
                        level_diff++;
                      }
                      while(level_diff > 0)
                      {
                        // Ref. level of c's face is greater than that of the neighbour's face.
                        // Get parent of the former.
                        face_c = face_c->parent();
                        level_diff--;
                      }
                      ASSERT(level_diff == 0, "Level difference " + stringify(level_diff) + " must be 0.");
                      ASSERT(face_c->refinement_level() == face_neigh->refinement_level(), "Face ref. level "
                             + stringify(face_c->refinement_level()) + " must equal ref. level "
                             + stringify(face_neigh->refinement_level()) + " of neighbour face.");
                      // check whether the two neighbours really share the same face
                      if (face_c != face_neigh)
                      {
                        s = "Face-neighbours cell " + c->print_index() + " (at face " + stringify((int)iface) + ", "
                            + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (face = "
                            + stringify((int)iface_in_neigh) + ", " + stringify(nn)
                            + "-th pos.) do not share the same face! Aborting program!\n";
                        throw new InternalError(s);
                      }
                      break;
                    }
                  }  // for 0 <= nn < neighs_of_neigh[iface_in_neigh].size()
                  if(neighbour_found)
                  {
                    break;
                  }
                } // for 0 <= iface_in_neigh < neighs[ineigh]->num_faces()
                if(!neighbour_found)
                {
                  s = "No face-neighbour between cell " + c->print_index() + " (at face " + stringify((int)iface)
                      + ", " + stringify(ineigh) + "-th pos.) and cell " + neighs[ineigh]->print_index()
                      + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= iface < c->num_faces()
        }
        catch(Exception& e)
        {
          ErrorHandler::exception_occured(e);
        }
      } // validate(...)
    }; // class CellDataValidationFaceNeigh


    /**
    * \brief class for validating cell-specific data like neighbourhood information
    *
    * This general class is empty, only specialisations for cell_dim_ = space_dim_ are implemented.
    *
    * \tparam cell_dim_
    * cell dimension (must be <= space dimension)
    *
    * \tparam space_dim_
    * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellDataValidation
    {

    private:

    public:

      /**
      * \brief dummy function called by cells with dimension smaller than space dimension
      *
      * \param[in] c
      * dummy cell
      *
      * \param[in,out] stream
      * dummy stream
      */
      static void validate_neighbourhood(
        Cell<cell_dim_, space_dim_, world_dim_> const *,
        std::ostream&)
      {
        // do nothing here
      }
    };


    /**
    * \brief specialisation of class CellDataValidation for cell_dim_ = space_dim_ = 3
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataValidation<3, 3, world_dim_>
    {

    private:

    public:

      /**
      * \brief validates neighbourhood of 3D cells
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate_neighbourhood(
        Cell<3, 3, world_dim_> const * c,
        std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidation::validate_neighbourhood()");
        CellDataValidationVertNeigh<3, world_dim_>::validate(c, stream);
        CellDataValidationEdgeNeigh<3, world_dim_>::validate(c, stream);
        CellDataValidationFaceNeigh<3, world_dim_>::validate(c, stream);
      }
    };


    /**
    * \brief specialisation of class CellDataValidation for cell_dim_ = space_dim_ = 2
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataValidation<2, 2, world_dim_>
    {

    private:

    public:

      /**
      * \brief validates neighbourhood of 2D cells
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate_neighbourhood(
        Cell<2, 2, world_dim_> const * c,
        std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidation::validate_neighbourhood()");
        CellDataValidationVertNeigh<2, world_dim_>::validate(c, stream);
        CellDataValidationEdgeNeigh<2, world_dim_>::validate(c, stream);
      }
    };


    /**
    * \brief specialisation of class CellDataValidation for cell_dim_ = space_dim_ = 1
    *
    * \tparam world_dim_
    * world dimension (determines the number of coordinates)
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataValidation<1, 1, world_dim_>
    {

    private:

    public:

      /**
      * \brief validates neighbourhood of 1D cells
      *
      * \param[in] c
      * cell whose neighbourhood is to be validated
      *
      * \param[in,out] stream
      * stream validation info is written into
      */
      static void validate_neighbourhood(
        Cell<1, 1, world_dim_> const * c,
        std::ostream& stream)
      {
        CONTEXT("BaseMesh::CellDataValidation::validate_neighbourhood()");
        CellDataValidationVertNeigh<1, world_dim_>::validate(c, stream);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BM_CELL_DATA_VALIDATION_HPP
