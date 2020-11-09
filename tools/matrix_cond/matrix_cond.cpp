// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////
//
// Matrix Condition Tool
// =====================
//
// This tool computes the extremal singular values of a matrix and its
// spectral condition number.
//
// \author Peter Zajac
//
///////////////////////////////////////////////////////////////////////////////


#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/diagonal_precond.hpp>
#include <kernel/solver/pcg.hpp>

#include <deque>
#include <vector>
#include <map>

namespace MatrixCond
{
  using namespace FEAT;
  using namespace FEAT::LAFEM;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef DenseVector<MemType, DataType, IndexType> VectorType;
  typedef SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  typedef NoneFilter<MemType, DataType, IndexType> FilterType;

  class MatrixATA
  {
  public:
    typedef MatrixCond::MemType MemType;
    typedef MatrixCond::DataType DataType;
    typedef MatrixCond::IndexType IndexType;
    typedef MatrixCond::VectorType VectorTypeL;
    typedef MatrixCond::VectorType VectorTypeR;

  private:
    const MatrixType& _mat_a;
    const MatrixType& _mat_at;
    mutable VectorType _vec_tmp;

  public:
    explicit MatrixATA(const MatrixType& mat_a, const MatrixType& mat_at) :
      _mat_a(mat_a), _mat_at(mat_at), _vec_tmp(mat_a.rows())
    {
    }

    VectorTypeL create_vector_l() const
    {
      return _mat_at.create_vector_l();
    }

    VectorTypeR create_vector_r() const
    {
      return _mat_a.create_vector_r();
    }

    Index rows() const
    {
      return _mat_at.rows();
    }

    Index columns() const
    {
      return _mat_a.columns();
    }

    void apply(VectorType& r, const VectorType& x) const
    {
      _mat_a.apply(_vec_tmp, x);
      _mat_at.apply(r, _vec_tmp);
    }

    void apply(VectorType& r, const VectorType& x, const VectorType& y, const DataType alpha = DataType(1)) const
    {
      _mat_a.apply(_vec_tmp, x);
      _mat_at.apply(r, _vec_tmp, y, alpha);
    }
  };

  bool import_matrix(MatrixType& matrix, const String& filename, FileMode file_mode)
  {
    try
    {
      switch(file_mode)
      {
      case FileMode::fm_mtx:
      case FileMode::fm_csr:
        matrix.read_from(file_mode, filename);
        return true;

      case FileMode::fm_ell:
        {
          SparseMatrixELL<MemType, DataType, IndexType> mtmp;
          mtmp.read_from(file_mode, filename);
          matrix.convert(mtmp);
          return true;
        }

      case FileMode::fm_coo:
        {
          SparseMatrixCOO<MemType, DataType, IndexType> mtmp;
          mtmp.read_from(file_mode, filename);
          matrix.convert(mtmp);
          return true;
        }

      case FileMode::fm_bm:
        {
          SparseMatrixBanded<MemType, DataType, IndexType> mtmp;
          mtmp.read_from(file_mode, filename);
          matrix.convert(mtmp);
          return true;
        }

      default:
        return false;
      }
    }
    catch(...)
    {
      return false;
    }

    return false;
  }

  class SinValBase
  {
  public:
    DataType value;
    DataType sv_tol;
    Index max_iter;
    Index num_iter;
    bool plot;
    bool want;
    bool have;

    SinValBase() :
      value(0.0),
      sv_tol(1E-6),
      max_iter(10000),
      num_iter(0),
      plot(true),
      want(true),
      have(false)
    {
    }
  };

  class MaxSinVal : public SinValBase
  {
  public:
    MaxSinVal() : SinValBase() {}

    template<typename Matrix_>
    int calc_max_sv(const Matrix_& matrix, VectorType& vec)
    {
      const DataType eps_sqr = Math::sqr(Math::eps<DataType>());

      // normalize initial vector
      value = vec.norm2();
      vec.scale(vec, DataType(1) / value);
      VectorType vec_tmp(vec.clone());

      // iterate power method
      for(num_iter = 1; num_iter <= max_iter; ++num_iter)
      {
        // appy mat-vec product
        matrix.apply(vec_tmp, vec);
        // compute new norm
        DataType sv_old = value;
        value = vec_tmp.norm2();
        // compute normalized vector
        vec.scale(vec_tmp, DataType(1) / value);
        // print new norm
        if(plot)
        {
          std::cout << "SV-Max Iter " << stringify(num_iter).pad_front(6) << ": " << stringify_fp_sci(value);
          // avoid division by zero
          if(eps_sqr*Math::abs(value-sv_old) < value)
            std::cout << " [" << stringify_fp_sci(Math::abs(value-sv_old)/value) << "]";
          std::cout << std::endl;
        }
        // check tolerance
        if(Math::abs(value-sv_old) < sv_tol*value)
        {
          // okay
          have = true;
          return int(num_iter);
        }
      }

      // no convergence
      return -1;
    }

    int calc_max_sv_spd(const MatrixType& matrix, VectorType& vec)
    {
      std::cout << "Computing maximal singular value by powering A..." << std::endl;
      return calc_max_sv(matrix, vec);
    }

    int calc_max_sv_ata(const MatrixType& matrix, const MatrixType& transpo, VectorType& vec)
    {
      std::cout << "Computing maximal singular value by powering A^T*A..." << std::endl;
      // create virtual A^T*A matrix
      MatrixATA mat_ata(matrix, transpo);
      int rtn = calc_max_sv(mat_ata, vec);
      if(have)
        value = Math::sqrt(value);
      return rtn;
    }
  }; // class MaxSinVal


  class MinSinVal : public SinValBase
  {
  public:
    Index pcg_maxiter;
    Index pcg_numiter;

    MinSinVal() :
      SinValBase(),
      pcg_maxiter(10000),
      pcg_numiter(0)
    {
    }

    bool calc_ata_diag(VectorType& vdiag, const MatrixType& matrix, const MatrixType& transpo)
    {
      // compute inverse main diagonal of A^T * A
      vdiag = VectorType(matrix.columns(), DataType(0));
      DataType* vdi = vdiag.elements();

      // fetch the transposed matrix arrays
      const IndexType* row_ptr_t = transpo.row_ptr();
      const DataType* data_t = transpo.val();

      DataType dmin(1E+99), dmax(0.0);
      for(Index i(0); i < transpo.rows(); ++i)
      {
        for(Index j(row_ptr_t[i]); j < row_ptr_t[i+1]; ++j)
        {
          vdi[i] += Math::sqr(data_t[j]);
        }
        dmin = Math::min(dmin, vdi[i]);
        dmax = Math::max(dmax, vdi[i]);
      }

      const DataType deps = Math::eps<DataType>();
      if(dmin < dmax*deps)
      {
        // matrix singular
        return false;
      }

      // invert diagonal entries
      for(Index i(0); i < transpo.rows(); ++i)
      {
        vdi[i] = DataType(1) / vdi[i];
      }

      return true;
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    int calc_min_sv(const Matrix_& matrix, const Filter_& filter, Precond_& precon, VectorType& vec)
    {
      const DataType eps_sqr = Math::sqr(Math::eps<DataType>());

      // create CG solver
      Solver::PCG<Matrix_, Filter_> solver(matrix, filter, precon);
      solver.init();
      solver.set_max_iter(pcg_maxiter);

      // normalize initial vector
      value = vec.norm2();
      vec.scale(vec, DataType(1) / value);
      VectorType vec_tmp(vec.clone());

      // iterate power method
      for(num_iter = 1; num_iter <= max_iter; ++num_iter)
      {
        // appy mat-vec product
        Solver::Status status = solver.apply(vec_tmp, vec);
        switch(status)
        {
        case Solver::Status::diverged:
        case Solver::Status::aborted:
        case Solver::Status::undefined:
          return -2;
        default:
          break;
        }
        pcg_numiter += solver.get_num_iter();

        // compute new norm
        DataType sv_old = value;
        value = vec_tmp.norm2();
        // compute normalized vector
        vec.scale(vec_tmp, DataType(1) / value);
        // print new norm
        if(plot)
        {
          std::cout << "SV-Min Iter " << stringify(num_iter).pad_front(6) << ": " << stringify_fp_sci(value);
          // avoid division by zero
          if(eps_sqr*Math::abs(value-sv_old) < value)
            std::cout << " [" << stringify_fp_sci(Math::abs(value-sv_old)/value) << "]";
          // get PCG statistics
          {
            Index cgit = solver.get_num_iter();
            DataType def1 = solver.get_def_final();
            std::cout <<  " {" << stringify(cgit).pad_front(5) << " : "  << stringify_fp_sci(def1) << "}";
          }
          std::cout << std::endl;
        }
        // check tolerance
        if(Math::abs(value-sv_old) < sv_tol*value)
        {
          // okay
          have = true;
          return int(num_iter);
        }
      }

      // no convergence
      return -1;
    }

    int calc_min_sv_spd(const MatrixType& matrix, const FilterType& filter, VectorType& vec)
    {
      std::cout << "Computing minimal singular value by powering PCG(A)..." << std::endl;
      auto precon = Solver::new_jacobi_precond(matrix, filter);
      int rtn = calc_min_sv(matrix, filter, precon, vec);
      if(rtn > 0)
        value = DataType(1) / value;
      return rtn;
    }

    int calc_min_sv_ata(const MatrixType& matrix, const MatrixType& transpo, const FilterType& filter, VectorType& vec)
    {
      std::cout << "Computing minimal singular value by powering PCD(A^T*A)..." << std::endl;

      // create virtual A^T*A matrix
      MatrixATA mat_ata(matrix, transpo);

      // compute inverse main diagonal of A^T * A
      VectorType vdiag;
      if(!calc_ata_diag(vdiag, matrix, transpo))
      {
        // matrix singular
        return -2;
      }

      // create a diagonal preconditioner
      auto precon = Solver::new_diagonal_precond(vdiag, filter);

      // compute minsv
      int rtn = calc_min_sv(mat_ata, filter, precon, vec);
      if(rtn > 0)
        value = DataType(1) / Math::sqrt(value);
      return rtn;
    }
  };

  void print_help()
  {
    std::cout << std::endl;
    std::cout << "USAGE: matrix-cond <filename> [<format>] [options...]" << std::endl;
    std::cout << std::endl;
    //            123456789-123456789-123456789-123456789-123456789-123456789-
    std::cout << "This tool can compute the extremal singular values and the" << std::endl;
    std::cout << "spectral condition number of a given matrix stored in the" << std::endl;
    std::cout << "file named <filename>." << std::endl;
    std::cout << std::endl;
    std::cout << "WARNING: This tool is experimental; use at own risk!" << std::endl;
    std::cout << std::endl;
    std::cout << "The maximal singular value is computed iteratively by" << std::endl;
    std::cout << "applying the power method onto the normal system A^T * A." << std::endl;
    std::cout << "The minimal singular value is computed iteratively by" << std::endl;
    std::cout << "applying the power method onto a Jacobi-preconditioned" << std::endl;
    std::cout << "PCG solver applied onto the normal system A^T * A." << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "The following options are available for this tool:" << std::endl;
    std::cout << std::endl;
    std::cout << "--spd" << std::endl;
    std::cout << "Assume that the matrix is symmetric and positive definite." << std::endl;
    std::cout << "In this case, the power method is applied onto the matrix A" << std::endl;
    std::cout << "itself rather than onto its normal system A^T * A." << std::endl;
    std::cout << "This usually speeds up the solution process." << std::endl;
    std::cout << std::endl;
    std::cout << "--no-svmin / --no-svmax" << std::endl;
    std::cout << "Do not compute the minimal or maximal singular values, resp." << std::endl;
    std::cout << std::endl;
    std::cout << "--svmin-noplot / --svmax-noplot" << std::endl;
    std::cout << "Do not plot the power method iterations to cout." << std::endl;
    std::cout << std::endl;
    std::cout << "--svmin-maxiter / --svmax-maxiter" << std::endl;
    std::cout << "Specifies the maximum number of iterations for the power" << std::endl;
    std::cout << "method for the computation the extremal singular values." << std::endl;
    std::cout << "If not specified, 10000 iterations are performed at most." << std::endl;
    std::cout << std::endl;
    std::cout << "--svmin-tol / --svmax-tol" << std::endl;
    std::cout << "Specifies the tolerance for the stopping criterion of the" << std::endl;
    std::cout << "power method. If not specified, 1E-6 is used." << std::endl;
    std::cout << std::endl;
    std::cout << "--svmin-read-vector <filename> [<format>] /" << std::endl;
    std::cout << "--svmax-read-vector <filename> [<format>]" << std::endl;
    std::cout << "Specifies the filename of the vector that is to be used as" << std::endl;
    std::cout << "the initial vector for the power method. If not given, the" << std::endl;
    std::cout << "power method start with the vector whose entries are all 1." << std::endl;
    std::cout << std::endl;
    std::cout << "--svmin-write-vector <filename> [<format>] /" << std::endl;
    std::cout << "--svmax-write-vector <filename> [<format>]" << std::endl;
    std::cout << "Writes out the eigenvectors to the minimal and maximal" << std::endl;
    std::cout << "singular values if the power method converged." << std::endl;
    std::cout << "If not specified, the vectors are not written out." << std::endl;
    std::cout << std::endl;
    std::cout << "--pcg-maxiter" << std::endl;
    std::cout << "Specifies the maximum number of PCG iterations for the" << std::endl;
    std::cout << "power method for the computation of the minimal singular" << std::endl;
    std::cout << "value. If not given, 10000 iterations are performed at most." << std::endl;
  }

  int main(int argc, char* argv[])
  {
    // create a matrix filemode map
    std::map<String, FileMode> mat_fm_map;
    mat_fm_map.insert(std::make_pair("mtx", FileMode::fm_mtx));
    mat_fm_map.insert(std::make_pair("ell", FileMode::fm_ell));
    mat_fm_map.insert(std::make_pair("csr", FileMode::fm_csr));
    mat_fm_map.insert(std::make_pair("coo", FileMode::fm_coo));
    mat_fm_map.insert(std::make_pair("bm", FileMode::fm_bm));

    // create a vector filemode map
    std::map<String, FileMode> vec_fm_map;
    vec_fm_map.insert(std::make_pair("exp", FileMode::fm_exp));
    vec_fm_map.insert(std::make_pair("dv" , FileMode::fm_dv));

    // create an argument parser
    SimpleArgParser args(argc, argv);

    // add all supported options
    args.support("spd");
    args.support("no-svmin");
    args.support("no-svmax");
    args.support("svmin-noplot");
    args.support("svmax-noplot");
    args.support("svmin-maxiter");
    args.support("svmax-maxiter");
    args.support("svmin-tol");
    args.support("svmax-tol");
    args.support("svmin-read-vector");
    args.support("svmax-read-vector");
    args.support("svmin-write-vector");
    args.support("svmax-write-vector");
    args.support("pcg-maxiter");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
      return 1;
    }

    // create objects
    MaxSinVal max_sv;
    MinSinVal min_sv;

    // vector i/o vars
    String svmax_vec_read_name;
    String svmin_vec_read_name;
    String svmax_vec_write_name;
    String svmin_vec_write_name;
    FileMode svmax_vec_read_fm(FileMode::fm_dv);
    FileMode svmin_vec_read_fm(FileMode::fm_dv);
    FileMode svmax_vec_write_fm(FileMode::fm_dv);
    FileMode svmin_vec_write_fm(FileMode::fm_dv);

    // check what to do today
    bool assume_spd = (args.check("spd") >= 0);
    min_sv.want = (args.check("no-svmin") < 0);
    max_sv.want = (args.check("no-svmax") < 0);
    min_sv.plot = (args.check("svmin-noplot") < 0);
    max_sv.plot = (args.check("svmax-noplot") < 0);
    args.parse("svmax-maxiter", max_sv.max_iter);
    args.parse("svmin-maxiter", min_sv.max_iter);
    args.parse("svmax-tol", max_sv.sv_tol);
    args.parse("svmin-tol", min_sv.sv_tol);
    args.parse("pcg-maxiter", min_sv.pcg_maxiter);
    bool read_svmax_vec = (args.parse("svmax-read-vector", svmax_vec_read_name, string_mapped(svmax_vec_read_fm, vec_fm_map)) > 0);
    bool read_svmin_vec = (args.parse("svmin-read-vector", svmin_vec_read_name, string_mapped(svmin_vec_read_fm, vec_fm_map)) > 0);
    bool write_svmax_vec = (args.parse("svmax-write-vector", svmax_vec_write_name, string_mapped(svmax_vec_write_fm, vec_fm_map)) > 0);
    bool write_svmin_vec = (args.parse("svmin-write-vector", svmin_vec_write_name, string_mapped(svmin_vec_write_fm, vec_fm_map)) > 0);

    // query matrix name and matrix format
    String matrix_name;
    FileMode matrix_mode(FileMode::fm_exp);

    // fetch number of skipped arguments
    int nsa = args.num_skipped_args();
    if(nsa <= 1)
    {
      // no skipped arguments; print help and exit
      print_help();
      return 0;
    }

    // treat second argument as input matrix name
    matrix_name = args.get_arg(1);

    // now let's check the matrix format
    String mat_frm;
    if(nsa > 2)
      mat_frm = args.get_arg(2);
    else
    {
      // fetch extension of filename
      mat_frm = matrix_name.split_by_string(".").back();
    }

    // try to map format
    if(!string_mapped_lookup(matrix_mode, mat_fm_map, mat_frm))
    {
      std::cerr << "ERROR: Failed to determine matrix file format for '" << matrix_name << "'" << std::endl;
      return 1;
    }

    std::cout << "Reading matrix '" << matrix_name << "'..." << std::endl;

    // let's import the matrix
    MatrixType matrix;
    if(!import_matrix(matrix, matrix_name, matrix_mode))
    {
      std::cerr << "ERROR: Failed to import matrix '" << matrix_name << "'" << std::endl;
      return 1;
    }

    std::cout << "Computing matrix transpose..." << std::endl;

    // now let's transpose the matrix
    MatrixType transpo = matrix.transpose();

    // create virtual A^T*A matrix
    MatrixATA mat_ata(matrix, transpo);

    // create a dummy filter
    FilterType filter;

    // compute maximal singular value
    if(max_sv.want)
    {
      std::cout << std::endl;

      // create initial vector
      VectorType vec_svmax(matrix.rows(), DataType(1));

      // read initial vector
      if(read_svmax_vec)
      {
        std::cout << "Reading initial vector from '" << svmax_vec_read_name << "'..." << std::endl;
        vec_svmax.read_from(svmax_vec_read_fm, svmax_vec_read_name);
      }

      // compute singular value by A^T*A
      int ret = 0;
      if(assume_spd)
        ret = max_sv.calc_max_sv_spd(matrix, vec_svmax);
      else
        ret = max_sv.calc_max_sv_ata(matrix, transpo, vec_svmax);

      // everything okay?
      if(ret > 0)
      {
        std::cout << "Power method converged after " << max_sv.num_iter << " iterations" << std::endl;
        if(write_svmax_vec)
        {
          std::cout << "Writing out vector to '" << svmax_vec_write_name << "'..." << std::endl;
          vec_svmax.write_out(svmax_vec_write_fm, svmax_vec_write_name);
        }
      }
      else if(ret == -1)
      {
        std::cout << "WARNING: Power method did not converge; maximal singular value may be garbage!" << std::endl;
      }
    }

    // compute minimal singular value
    if(min_sv.want)
    {
      std::cout << std::endl;

      // create initial vector
      VectorType vec_svmin(matrix.rows(), DataType(1));

      // read initial vector
      if(read_svmin_vec)
      {
        std::cout << "Reading initial vector from '" << svmin_vec_read_name << "'..." << std::endl;
        vec_svmin.read_from(svmin_vec_read_fm, svmin_vec_read_name);
      }

      // compute singular value by A^T*A
      int ret = 0;
      if(assume_spd)
        ret = min_sv.calc_min_sv_spd(matrix, filter, vec_svmin);
      else
        ret = min_sv.calc_min_sv_ata(matrix, transpo, filter, vec_svmin);

      // everything okay?
      if(ret > 0)
      {
        std::cout << "Power method converged after " << min_sv.num_iter << " iterations" << std::endl;
        if(write_svmin_vec)
        {
          std::cout << "Writing out vector to '" << svmin_vec_write_name << "'..." << std::endl;
          vec_svmin.write_out(svmin_vec_write_fm, svmin_vec_write_name);
        }
      }
      else if(ret == -1)
      {
        std::cout << "WARNING: Power method did not converge; minimal singular value may be garbage!" << std::endl;
      }
      else if(ret == -2)
      {
        std::cout << "ERROR: PCG solver did not converge; minimal singular value may be garbage!" << std::endl;
      }
    }

    std::cout << std::endl << "Summary:" << std::endl;
    if(max_sv.want)
    {
      std::cout << "Max SV PM Iter....: " << max_sv.num_iter << std::endl;
    }
    if(min_sv.want)
    {
      std::cout << "Min SV PM Iter....: " << min_sv.num_iter << std::endl;
      std::cout << "Min SV PCG Iter...: " << min_sv.pcg_numiter << std::endl;
    }
    if(max_sv.have)
      std::cout << "Max Singular Value: " << stringify_fp_sci(max_sv.value) << std::endl;
    if(min_sv.have)
      std::cout << "Min Singular Value: " << stringify_fp_sci(min_sv.value) << std::endl;
    if(max_sv.have && min_sv.have)
      std::cout << "Spectral Condition: " << stringify_fp_sci(max_sv.value / min_sv.value) << std::endl;

    // okay
    return 0;
  }
} // namespace MatrixCond

int main(int argc, char* argv[])
{
  return MatrixCond::main(argc, argv);
}
