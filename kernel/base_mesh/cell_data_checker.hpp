#pragma once
#ifndef KERNEL_BM_CELL_DATA_CHECKER_HPP
#define KERNEL_BM_CELL_DATA_CHECKER_HPP 1

// includes, system
#include <iostream> // for std::ostream
#include <cassert>  // for assert()
#include <vector>   // for std::vector

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/base_mesh/cell_data.hpp>
#include <kernel/base_mesh/cell.hpp>
#include <kernel/base_mesh/vertex.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    /**
    * \brief class for checking vertex neighbours
    *
    * \note
    * Q: Why are there three classes CellDataCheckerVertexNeighbours, CellDataCheckerEdgeNeighbours and
    *    CellDataCheckerFaceNeighbours? Can they not be united (the code is quite similar)?
    * A: 1) That is a consequence of that fact, that we explicitely store vertices, edges and faces instead of
    *    putting all in one array of general "items".
    *    2) If we only had to compare subcells here, we maybe could use the most general class BaseMeshItem (being
    *    parent of vertices, edges and faces). However, we also have to access children/parents of those subcells,
    *    and that is not possible if we only handle them as BaseMeshItem. (Note that for vertices we don't need
    *    parent/child information.)
    *    3) The code differs for vertices on the one hand and edges/faces on the other hand: For the latter we
    *    eventually have to traverse the history (parent/child relations) for checking neighbourhood, which is not
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
    class CellDataCheckerVertexNeighbours
    {
    private:
    public:
      /// checks whether vertex neighbours are set correctly
      static void check(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c)
      {
        assert(cell_space_dim_ >= 1);
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
                s = "Cell " + neighs[ineigh]->print_index() + ", being the vertex-neighbour of cell "
                    + c->print_index() + " (at vertex " + StringUtils::stringify((int)ivertex) +  ", "
                    + StringUtils::stringify(ineigh) + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: Deactivated until update of neighbourhood has been implemented! Instead, only print error message.
//                throw new InternalError(s);
                std::cerr << "ERROR: " << s;
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
    + StringUtils::stringify((int)ivertex) + ", " + StringUtils::stringify(ineigh)
    + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (vertex = "
    + StringUtils::stringify((int)ivertex_in_neigh) + ", " + StringUtils::stringify(nn)
    + "-th pos.) found!\n";
std::cout << s;
                      neighbour_found = true;

                      // check whether the two neighbours really share the same vertex
                      if (c->vertex(ivertex) != neighs[ineigh]->vertex(ivertex_in_neigh))
                      {
                        s = "Vertex-neighbours cell " + c->print_index() + " (at vertex "
                            + StringUtils::stringify((int)ivertex) + ", " + StringUtils::stringify(ineigh)
                            + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (vertex = "
                            + StringUtils::stringify((int)ivertex_in_neigh) + ", " + StringUtils::stringify(nn)
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
                  s = "No vertex-neighbour between cell " + c->print_index() + " (at vertex "
                      + StringUtils::stringify((int)ivertex) + ", " + StringUtils::stringify(ineigh)
                      + "-th pos.) and cell " + neighs[ineigh]->print_index() + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= ivertex < c->num_vertices()
        }
        catch(InternalError* e)
        {
          std::cerr << e->message() << std::endl;
          exit(1);
        }
      } // check(...)
    }; // class CellDataCheckerVertexNeighbours


    /**
    * \brief class for checking edge neighbours
    *
    * \note See note in class CellDataCheckerVertexNeighbours.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellDataCheckerEdgeNeighbours
    {
    private:
    public:
      /// checks whether edge neighbours are set correctly
      static void check(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c)
      {
        assert(cell_space_dim_ >= 2);
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
                    + c->print_index() + " (at edge " + StringUtils::stringify((int)iedge) +  ", "
                    + StringUtils::stringify(ineigh) + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: temporarily deactivated!
//                throw new InternalError(s);
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
s = "Edge-neighbourhood between cell " + c->print_index() + " (at edge "
    + StringUtils::stringify((int)iedge) + ", " + StringUtils::stringify(ineigh)
    + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (edge = "
    + StringUtils::stringify((int)iedge_in_neigh) + ", " + StringUtils::stringify(nn)
    + "-th pos.) found!\n";
std::cout << s;
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
                      assert(level_diff == 0);
                      assert(edge_c->refinement_level() == edge_neigh->refinement_level());
                      // check whether the two neighbours really share the same edge
                      if (edge_c != edge_neigh)
                      {
                        s = "Edge-neighbours cell " + c->print_index() + " (at edge "
                            + StringUtils::stringify((int)iedge) + ", " + StringUtils::stringify(ineigh)
                            + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (edge = "
                            + StringUtils::stringify((int)iedge_in_neigh) + ", " + StringUtils::stringify(nn)
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
                  s = "No edge-neighbour between cell " + c->print_index() + " (at edge "
                      + StringUtils::stringify((int)iedge) + ", " + StringUtils::stringify(ineigh)
                      + "-th pos.) and cell " + neighs[ineigh]->print_index() + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= iedge < c->num_edges()
        }
        catch(InternalError* e)
        {
          std::cerr << e->message() << std::endl;
          exit(1);
        }
      } // check(...)
    }; // class CellDataCheckerEdgeNeighbours


    /**
    * \brief class for checking face neighbours
    *
    * \note See note in class CellDataCheckerVertexNeighbours.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_space_dim_,
      unsigned char world_dim_>
    class CellDataCheckerFaceNeighbours
    {
    private:
    public:
      /// checks whether face neighbours are set correctly
      static void check(Cell<cell_space_dim_, cell_space_dim_, world_dim_> const * c)
      {
        assert(cell_space_dim_ >= 3);
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
                s = "Cell " + neighs[ineigh]->print_index() + ", being the face-neighbour of cell "
                    + c->print_index() + " (at face " + StringUtils::stringify((int)iface) +  ", "
                    + StringUtils::stringify(ineigh) + "-th pos.), has children which must not be the case!\n";
                s += "It seems the neighbours of cell " + c->print_index() + " have not been updated correctly!\n";
// COMMENT_HILMAR: temporarily deactivated!
//                throw new InternalError(s);
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
s = "Face-neighbourhood between cell " + c->print_index() + " (at face "
    + StringUtils::stringify((int)iface) + ", " + StringUtils::stringify(ineigh)
    + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (face = "
    + StringUtils::stringify((int)iface_in_neigh) + ", " + StringUtils::stringify(nn)
    + "-th pos.) found!\n";
std::cout << s;
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
                      assert(level_diff == 0);
                      assert(face_c->refinement_level() == face_neigh->refinement_level());
                      // check whether the two neighbours really share the same face
                      if (face_c != face_neigh)
                      {
                        s = "Face-neighbours cell " + c->print_index() + " (at face "
                            + StringUtils::stringify((int)iface) + ", " + StringUtils::stringify(ineigh)
                            + "-th pos.) and cell " + neighs[ineigh]->print_index() + " (face = "
                            + StringUtils::stringify((int)iface_in_neigh) + ", " + StringUtils::stringify(nn)
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
                  s = "No face-neighbour between cell " + c->print_index() + " (at face "
                      + StringUtils::stringify((int)iface) + ", " + StringUtils::stringify(ineigh)
                      + "-th pos.) and cell " + neighs[ineigh]->print_index() + " found! Aborting program!\n";
                  throw new InternalError(s);
                }
              }
            } // for 0 <= ineigh < neighs.size()
          } // for 0 <= iface < c->num_faces()
        }
        catch(InternalError* e)
        {
          std::cerr << e->message() << std::endl;
          exit(1);
        }
      } // check(...)
    }; // class CellDataCheckerFaceNeighbours


    /**
    * \brief class for checking cell-specific data like neighbourhood information
    *
    * This general class is empty, only specialisations for cell_dim_ = space_dim_ are implemented.
    *
    * \author Hilmar Wobker
    */
    template<
      unsigned char cell_dim_,
      unsigned char space_dim_,
      unsigned char world_dim_>
    class CellDataChecker
    {
    private:
    public:
      /// dummy function called by cells with dimension smaller than space dimension
      static void check_neighbourhood(Cell<cell_dim_, space_dim_, world_dim_> const * c)
      {
      }
    };


    /**
    * \brief specialisation of class CellDataChecker for cell_dim_ = space_dim_ = 3
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataChecker<3, 3, world_dim_>
    {
    private:
    public:
      static void check_neighbourhood(Cell<3, 3, world_dim_> const * c)
      {
        CellDataCheckerVertexNeighbours<3, world_dim_>::check(c);
        CellDataCheckerEdgeNeighbours<3, world_dim_>::check(c);
        CellDataCheckerFaceNeighbours<3, world_dim_>::check(c);
      }
    };


    /**
    * \brief specialisation of class CellDataChecker for cell_dim_ = space_dim_ = 2
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataChecker<2, 2, world_dim_>
    {
    private:
    public:
      static void check_neighbourhood(Cell<2, 2, world_dim_> const * c)
      {
        CellDataCheckerVertexNeighbours<2, world_dim_>::check(c);
        CellDataCheckerEdgeNeighbours<2, world_dim_>::check(c);
      }
    };


    /**
    * \brief specialisation of class CellDataChecker for cell_dim_ = space_dim_ = 1
    *
    * \author Hilmar Wobker
    */
    template<unsigned char world_dim_>
    class CellDataChecker<1, 1, world_dim_>
    {
    private:
    public:
      static void check_neighbourhood(Cell<1, 1, world_dim_> const * c)
      {
        CellDataCheckerVertexNeighbours<1, world_dim_>::check(c);
      }
    };
  } // namespace BaseMesh
} // namespace FEAST

#endif // #define KERNEL_BM_CELL_DATA_CHECKER_HPP
