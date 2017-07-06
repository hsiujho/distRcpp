
#include <Rcpp.h>
using namespace Rcpp;

// generic function for gunifrac_distance
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline double gunifrac_distance(InputIterator1 begin1, InputIterator1 end1,
                                InputIterator2 begin2,
                                InputIterator3 begin3, double powv) {

  // value to return
  double up=0;
  double down=0;

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  InputIterator3 it3 = begin3;

  // for each input item
  while (it1 != end1) {

    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    double brlen = *it3++;

    // accumulate if appropirate
    if (d1 > 0 || d2 > 0){
      double d3 = d1+d2;
      double bl = brlen * pow(d3, powv);
      up += bl*std::abs(d1-d2) / d3;
      down += bl;
    }
  }
  return (up/down);
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct GufDistance : public Worker {

  // input matrix to read from
  const RMatrix<double> mat;
  const NumericVector brlen;
  const double powv;

  // output matrix to write to
  RMatrix<double> rmat;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  GufDistance(const NumericMatrix mat, const NumericVector brlen, const double powv, NumericMatrix rmat)
    : mat(mat), brlen(brlen), powv(powv), rmat(rmat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {

        //        NumericVector col1 = mat.column(i);
        //       NumericVector col2 = mat.column(j);
        RMatrix<double>::Column col1 = mat.column(i);
        RMatrix<double>::Column col2 = mat.column(j);

        double dist = gunifrac_distance(col1.begin(), col1.end(), col2.begin(), brlen.begin(),powv);
        rmat(i,j) = dist;
        rmat(j,i) = rmat(i,j);

      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix GunifracRcpp(NumericMatrix mat, NumericVector brlen, double powv) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat.ncol(), mat.ncol());

  // create the worker
  GufDistance gufDistance(mat,brlen,powv, rmat);

  // call it with parallelFor
  parallelFor(0, mat.ncol(), gufDistance);

  return rmat;
}
