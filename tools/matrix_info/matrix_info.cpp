///////////////////////////////////////////////////////////////////////////////
//
// Matrix Information Tool
// =======================
//
// This is a tool which performs a simple analysis of a given matrix.
// This tool certainly cannot complete with the analytical powers of tools
// like Matlab or Maple, but it provides a fast and simple way to get some
// basic information about a matrix.
//
// See the documentation page "The matrix-info Tool" in the FEAST documentation
// for detailed information about this tool and its usage.
//
// \author Peter Zajac
//
///////////////////////////////////////////////////////////////////////////////


#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

#include <deque>
#include <vector>
#include <map>

namespace MatrixInfo
{
  using namespace FEAST;
  using namespace FEAST::LAFEM;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef DenseVector<MemType, DataType, IndexType> VectorType;
  typedef SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;

  // padding length
  static constexpr std::size_t pad_len = 27;

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

  String percentify(double d)
  {
    int di = int(1000000.0 * d);
    int i0 = di / 10000;
    int i1 = ((di % 10000) + 5) / 10;
    return (stringify(i0).pad_front(2) + ".") + (stringify(i1).pad_front(3, '0')) + "%";
  }

  // decodes a bit-wise combined matrix content code
  String dec_content_code(int code)
  {
    // 0x1: has zero entries
    // 0x2: has positive entries
    // 0x4: has negative entries
    switch(code)
    {
    case 0: return "empty";
    case 1: return "zero";
    case 2: return "positive";
    case 3: return "non-negative";
    case 4: return "negative";
    case 5: return "non-positive";
    case 6: return "non-zero";
    case 7: return "mixed";
    default: return "-invalid-";
    }
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

    // create an argument parser
    SimpleArgParser args(argc, argv);

    // add all supported options
    args.support("nztol");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;
      return 1;
    }

    // query options
    DataType nz_tol = Math::sqrt(Math::eps<DataType>());
    int ia = args.parse("nztol", nz_tol);
    if(ia < 0)
    {
      std::cerr << "ERROR: Failed to parse '" << args.get_arg(-ia) << "' as non-zero tolenance" << std::endl;
      return 1;
    }

    // query matrix name and matrix format
    String matrix_name;
    FileMode matrix_mode(FileMode::fm_exp);

    // fetch number of skipped arguments
    int nsa = args.num_skipped_args();
    if(nsa <= 1)
    {
      // no skipped arguments
      std::cout << "USAGE: matrix-info <filename> [<format>] [options...]" << std::endl;
      std::cout << std::endl;
      std::cout << "See the doxygen documentation for details about this tool's" << std::endl;
      std::cout << "options and its usage." << std::endl;
      std::cout << std::endl;
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
      std::deque<String> ds;
      matrix_name.split_by_charset(ds, ".");
      mat_frm = ds.back();
    }

    // try to map format
    if(!string_mapped_lookup(matrix_mode, mat_fm_map, mat_frm))
    {
      std::cerr << "ERROR: Failed to determine matrix file format for '" << matrix_name << "'" << std::endl;
      return 1;
    }

    // let's import the matrix
    MatrixType matrix;
    if(!import_matrix(matrix, matrix_name, matrix_mode))
    {
      std::cerr << "ERROR: Failed to import matrix '" << matrix_name << "'" << std::endl;
      return 1;
    }

    // print matrix name
    std::cout << std::endl;
    std::cout << String("Matrix Filename").pad_back(pad_len, '.') << ": " << matrix_name << std::endl;

    // now let's transpose the matrix
    MatrixType transpo = matrix.transpose();

    // get the matrix sizes
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const Index nmin = Math::min(nrows, ncols);
    const Index nnze = matrix.used_elements();

    // get the matrix arrays
    const IndexType* row_ptr = matrix.row_ptr();
    const IndexType* col_idx = matrix.col_ind();
    const DataType* mat_data = matrix.val();

    // get the transposed matrix arrays
    const IndexType* row_ptr_t = transpo.row_ptr();
    const IndexType* col_idx_t = transpo.col_ind();
    const DataType* mat_data_t = transpo.val();

    // allocate two vectors
    std::vector<DataType> vsum_row(nrows, DataType(0));
    std::vector<DataType> vsum_col(ncols, DataType(0));
    std::vector<DataType> vmax_row(nrows, DataType(0));
    std::vector<DataType> vmax_col(ncols, DataType(0));
    std::vector<DataType> vabs_dia(nmin, DataType(0));
    std::vector<Index> vnze_row(nrows, Index(0));
    std::vector<Index> vnze_col(ncols, Index(0));
    std::vector<Index> diag_ptr(nmin, ~Index(0));
    std::vector<Index> row_degree_distribution_counts(10, Index(0));

    // initialise statistical values
    Index row_degree = 0;
    Index row_degree_idx = 0;
    Index row_bandw = 0;
    Index row_bandw_idx = 0;
    Index col_degree = 0;
    Index col_degree_idx = 0;
    Index col_bandw = 0;
    Index col_bandw_idx = 0;
    String row_degree_distribution = "";
    Index average_row_degree = 0;

    // various norms
    DataType norm_col_sum = DataType(0);
    DataType norm_row_sum = DataType(0);
    DataType norm_frobenius = DataType(0);
    DataType norm_max = DataType(0);

    Index norm_col_sum_idx = 0;
    Index norm_row_sum_idx = 0;
    Index norm_max_row_idx = 0;
    Index norm_max_col_idx = 0;

    // row/col dominance factory
    DataType row_dom = DataType(0);
    DataType col_dom = DataType(0);
    Index row_dom_idx = nmin;
    Index col_dom_idx = nmin;

    bool have_struct_diag = (nrows == ncols);

    // first, let's loop over the transposed matrix rows, which correspond to our matrix columns
    for(Index j(0); j < ncols; ++j)
    {
      // compute col degree
      Index my_col_degree = row_ptr_t[j+1] - row_ptr_t[j];

      // store col degree
      vnze_col[j] = my_col_degree;

      // update degree
      if(my_col_degree > col_degree)
      {
        col_degree = my_col_degree;
        col_degree_idx = j;
      }

      // compute bandwidth
      if(my_col_degree > 0)
      {
        Index my_col_bandw = col_idx_t[row_ptr_t[j+1]-1] - col_idx_t[row_ptr_t[j]] + 1;
        if(my_col_bandw > col_bandw)
        {
          col_bandw = my_col_bandw;
          col_bandw_idx = j;
        }
      }
    }

    // now loop over the matrix rows
    for(Index i(0); i < nrows; ++i)
    {
      // compute degree
      Index my_row_degree = row_ptr[i+1] - row_ptr[i];

      average_row_degree += my_row_degree;

      // store row degree
      vnze_row[i] = my_row_degree;

      // update degree
      if(my_row_degree > row_degree)
      {
        row_degree = my_row_degree;
        row_degree_idx = i;
      }


      // compute bandwidth
      if(my_row_degree > 0)
      {
        Index my_row_bandw = col_idx[row_ptr[i+1]-1] - col_idx[row_ptr[i]] + 1;
        if(my_row_bandw > row_bandw)
        {
          row_bandw = my_row_bandw;
          row_bandw_idx = i;
        }
      }

      bool my_struct_diag = false;

      // an loop over all non-zeros in the row
      for(Index k(row_ptr[i]); k < row_ptr[i+1]; ++k)
      {
        // get the column index
        Index j = col_idx[k];

        // get our absolute entry
        const DataType abs_ij = Math::abs(mat_data[k]);

        // did we find the diagonal entry?
        if(j == i)
        {
          my_struct_diag = true;
          diag_ptr[i] = k;
          vabs_dia[i] = abs_ij;
        }

        // update row/col sums
        vsum_row[i] += abs_ij;
        vsum_col[j] += abs_ij;

        // update row/col max
        vmax_row[i] = Math::max(vmax_row[i], abs_ij);
        vmax_col[j] = Math::max(vmax_col[j], abs_ij);

        // update frobenius norm
        norm_frobenius += Math::sqr(abs_ij);

        // update max norm
        if(abs_ij > norm_max)
        {
          norm_max = abs_ij;
          norm_max_row_idx = i;
          norm_max_col_idx = i;
        }
      }

      // do we have a structural diagonal entry?
      have_struct_diag = have_struct_diag && my_struct_diag;
    }

    average_row_degree /= nrows;

    // now loop again over the matrix rows
    // we need to know the max row degree to calculate row distribution
    for(Index i(0); i < nrows; ++i)
    {
      float percentage(float(vnze_row[i]) / float(row_degree) * float(10));
      percentage = Math::floor(percentage);
      percentage = Math::min(percentage, float(9));
      row_degree_distribution_counts.at(Index(percentage)) += Index(1);
    }
    // format distribution string
    for(Index i(0) ; i < Index(row_degree_distribution_counts.size()) ; ++i)
    {
      float percentage(float(row_degree_distribution_counts[i]) / float(nrows) * float(100));
      row_degree_distribution += stringify(Index(percentage));
      if (i < row_degree_distribution_counts.size() - 1)
        row_degree_distribution += " | ";
    }

    // compute frobenius norm
    norm_frobenius = Math::sqrt(norm_frobenius);

    // loop over row info
    for(Index i(0); i < nrows; ++i)
    {
      // update row-sum norm
      if(vsum_row[i] > norm_row_sum)
      {
        norm_row_sum = vsum_row[i];
        norm_row_sum_idx = i;
      }

      // do we have a numerical diagonal entry?
      if((i < nmin) && (vabs_dia[i] > (vsum_row[i] * nz_tol)))
      {
        // yes, so let's check if it is dominant
        DataType my_row_dom = (vsum_row[i] / vabs_dia[i]) - DataType(1);
        if(my_row_dom > row_dom)
        {
          row_dom = my_row_dom;
          row_dom_idx = i;
        }
      }
    }

    // loop over column info
    for(Index j(0); j < ncols; ++j)
    {
      // update col-sum norm
      if(vsum_col[j] > norm_col_sum)
      {
        norm_col_sum = vsum_col[j];
        norm_col_sum_idx = j;
      }

      // do we have a numerical diagonal entry?
      if((j < nmin) && (vabs_dia[j] > (vsum_col[j] * nz_tol)))
      {
        // yes, so let's check if it is dominant
        DataType my_col_dom = (vsum_col[j] / vabs_dia[j]) - DataType(1);
        if(my_col_dom > col_dom)
        {
          col_dom = my_col_dom;
          col_dom_idx = j;
        }
      }
    }

    // matrix content codes; see dec_content_code for details
    int mat_content = 0;
    int dia_content = 0;
    int upp_content = 0;
    int low_content = 0;

    // now let's perform another sweep over the matrix
    Index numerical_non_zero = 0;
    for(Index i(0); i < nrows; ++i)
    {
      // loop over the row
      for(Index k(row_ptr[i]); k < row_ptr[i+1]; ++k)
      {
        // fetch the column index
        Index j = col_idx[k];

        // fetch the entry and its absolute
        DataType a_ij = mat_data[k];
        DataType abs_ij = Math::abs(a_ij);

        // determine whether the entry is really non-zero
        // we do this by checking whether its absolute value against the maximum norm
        bool treat_as_zero = (abs_ij < nz_tol * norm_max);

        // increment numerical non-zero
        if(!treat_as_zero)
          ++numerical_non_zero;

        // compute content code
        int content_code = (treat_as_zero ? 0x1 : (a_ij > DataType(0) ? 0x2 : 0x4));

        // update matrix content code
        mat_content |= content_code;

        // update part matrix
        if(i < j)
          upp_content |= content_code;
        else if(i > j)
          low_content |= content_code;
        else
          dia_content |= content_code;
      }
    }

    bool sym_symbolic = false;
    bool sym_numeric = false;

    // is the matrix square?
    if(nrows == ncols)
    {

      // let's see if the matrix is symmetric
      sym_symbolic = true;

      // compare the row pointer arrays
      for(Index i(0); sym_symbolic && (i <= nrows); ++i)
      {
        sym_symbolic = sym_symbolic && (row_ptr[i] == row_ptr_t[i]);
      }

      // compare the col-index arrays
      for(Index i(0); sym_symbolic && (i < nnze); ++i)
      {
        sym_symbolic = sym_symbolic && (col_idx[i] == col_idx_t[i]);
      }

      // compare the data arrays
      sym_numeric = sym_symbolic;
      for(Index i(0); sym_numeric && (i < nnze); ++i)
      {
        sym_numeric = sym_numeric && (Math::abs(mat_data[i] - mat_data_t[i]) <  nz_tol * norm_max);
      }
    }

    // print some basic information
    std::cout << String("Non-Zero Tolerance").pad_back(pad_len, '.') << ": " << scientify(nz_tol) << std::endl;
    std::cout << String("Number of Rows").pad_back(pad_len, '.') << ": " << nrows << std::endl;
    std::cout << String("Number of Columns").pad_back(pad_len, '.') << ": " << ncols << std::endl;
    std::cout << String("Max Row Degree").pad_back(pad_len, '.') << ": " << row_degree << " (max in row " << row_degree_idx << ")" << std::endl;
    std::cout << String("Max Column Degree").pad_back(pad_len, '.') << ": " << col_degree << " (max in col " << col_degree_idx << ")" << std::endl;
    std::cout << String("Row Degree Distribution").pad_back(pad_len, '.') << ": " << row_degree_distribution << std::endl;
    std::cout << String("Average Row Degree").pad_back(pad_len, '.') << ": " << average_row_degree << std::endl;
    std::cout << String("Max Row Bandwidth").pad_back(pad_len, '.') << ": " << row_bandw << " (max in row " << row_bandw_idx << ")" << std::endl;
    std::cout << String("Max Column Bandwidth").pad_back(pad_len, '.') << ": " << col_bandw << " (max in col " << col_bandw_idx << ")" << std::endl;
    std::cout << String("Symbolical Non-Zeros").pad_back(pad_len, '.') << ": " << nnze << std::endl;
    std::cout << String("Numerical Non-Zeros").pad_back(pad_len, '.') << ": " << numerical_non_zero;
    if(nnze > 0)
      std::cout << " (" << percentify(double(numerical_non_zero) / double(nnze)) << ")";
    std::cout << std::endl;
    std::cout << String("Has Main Diag Structure").pad_back(pad_len, '.') << ": " << (have_struct_diag ? "yes" : "no") << std::endl;
    std::cout << String("Has Symmetric Structure").pad_back(pad_len, '.') << ": ";
    if(nrows != ncols)
      std::cout << "no (matrix is not square)" << std::endl;
    else
      std::cout << (sym_symbolic ? "yes" : "no") << std::endl;;
    std::cout << String("Has Symmetric Data").pad_back(pad_len, '.') << ": ";
    if(nrows != ncols)
      std::cout << "no (matrix is not square)" << std::endl;
    else
      std::cout << (sym_numeric ? "yes" : "no") << std::endl;;

    // plot content information
    std::cout << String("Matrix Content").pad_back(pad_len, '.') << ": " << dec_content_code(mat_content) << std::endl;
    std::cout << String("Diagonal Content").pad_back(pad_len, '.') << ": " << dec_content_code(dia_content) << std::endl;
    std::cout << String("Lower Content").pad_back(pad_len, '.') << ": " << dec_content_code(low_content) << std::endl;
    std::cout << String("Upper Content").pad_back(pad_len, '.') << ": " << dec_content_code(upp_content) << std::endl;

    // plot diagonal dominance factory
    std::cout << String("Row-Diag Dominance Factor").pad_back(pad_len, '.') << ": ";
    if(row_dom_idx >= nmin)
      std::cout << "-NA-" << std::endl;
    else
      std::cout << scientify(row_dom) << " (in row " << row_dom_idx << ")" << std::endl;
    std::cout << String("Col-Diag Dominance Factor").pad_back(pad_len, '.') << ": ";
    if(col_dom_idx >= nmin)
      std::cout << "-NA-" << std::endl;
    else
      std::cout << scientify(col_dom) << " (in col " << col_dom_idx << ")" << std::endl;

    // plot norms
    std::cout << String("Frobenius Norm").pad_back(pad_len, '.') << ": " << scientify(norm_frobenius) << std::endl;
    std::cout << String("Col-Sum Norm").pad_back(pad_len, '.') << ": " << scientify(norm_col_sum)
      << " (max in col " << norm_col_sum_idx << ")" << std::endl;
    std::cout << String("Row-Sum Norm").pad_back(pad_len, '.') << ": " << scientify(norm_row_sum)
      << " (max in row " << norm_row_sum_idx << ")" << std::endl;
    std::cout << String("Maximum Norm").pad_back(pad_len, '.') << ": " << scientify(norm_max)
      << " (max in row " << norm_max_row_idx << ", col " << norm_max_col_idx << ")" << std::endl;

    // okay
    return 0;
  }
} // namespace MatrixInfo

int main(int argc, char* argv[])
{
  return MatrixInfo::main(argc, argv);
}
