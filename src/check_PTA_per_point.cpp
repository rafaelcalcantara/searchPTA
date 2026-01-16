#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix checkPTA(NumericVector pta_nodes, List points) {

  int k = pta_nodes.size();
  int n = points.size();

  LogicalMatrix out(n, k);

  for (int i = 0; i < n; i++) {
    NumericVector point = points[i];

    // Convert "point" to std::unordered_set for faster look-up
    std::unordered_set<double> set_i(point.begin(), point.end());

    for (int j = 0; j < k; j++) {
      out(i, j) = set_i.count(pta_nodes[j]) > 0;
    }
  }

  CharacterVector colnames = as<CharacterVector>(wrap(pta_nodes));

  CharacterVector rownames = as<CharacterVector>(wrap(points.names()));

  out.attr("dimnames") = List::create(rownames, colnames);


  return out;
}
