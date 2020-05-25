#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "config.h"
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]

const int max_nodes( const int & max_tree_depth ) {
  const int &result = 2*(pow(2, max_tree_depth))-1;
  return result;
}


class Tree_node {
public:
  arma::rowvec coef;
  arma::uvec filter;
};



const Tree_node vanilla_comparison( const arma::mat& X ){

  Tree_node result;

  int selection_var = round((X.n_cols-1)*randu());
  arma::vec splitting = X.col(selection_var);
  double comparison = splitting.min() + (splitting.max() - splitting.min()) * randu();

  arma::uvec filtering = find(splitting < comparison);

  arma::rowvec coef(2);

  coef(0) = selection_var;
  coef(1) = comparison;

  result.filter = filtering;
  result.coef = coef;

  return result;
}


const Tree_node extended_comparison( const arma::mat& X, const int& extension ){

  Tree_node result;

  arma::rowvec mins = min(X, 0);
  arma::rowvec maxes = max(X, 0);
  arma::vec intercepts = randu(X.n_cols);
  intercepts = (mins + (maxes - mins)).t() % intercepts;

  arma::vec norm_vec = randn(X.n_cols);
  arma::vec cols_to_smpl = linspace(0,(X.n_cols-1), X.n_cols);
  int n_extend = (X.n_cols - 1 - extension);
  cols_to_smpl = shuffle(cols_to_smpl);
  cols_to_smpl.set_size(n_extend);
  arma::uvec index = conv_to<uvec>::from(cols_to_smpl);
  norm_vec.elem(index).fill(0);

  arma::mat X_ = X;
  X_.each_row() +=  intercepts.t();
  arma::vec filter_selection = X_ * norm_vec;

  arma::uvec filtering = find(filter_selection < 0);

  arma::rowvec coef(norm_vec.n_elem + intercepts.n_elem);

  for(int i = 0; i < intercepts.n_elem; i++){
    coef(i) = intercepts(i);
  }
  int j = 0;
  for(int i = intercepts.n_elem; i < norm_vec.n_elem + intercepts.n_elem; i++){
    coef(i) = norm_vec(j);
    j++ ;
  }

  result.filter = filtering;
  result.coef = coef;

  return result;
}


const Tree_node split_cpp( const arma::mat& X, const int& extension,
                           const std::string& type = "vanilla" ){
  Tree_node result;

  if(type == "vanilla"){

    result = vanilla_comparison( X );
  }
  else if(type == "extended") {
    result = extended_comparison( X, extension );
  }
    return result;
}


class iTree{
   public:
     arma::mat X;
     arma::mat Id_mat;
     arma::mat pred_mat;
};


NumericVector negative_subset(NumericVector vec, NumericVector filter) {

  LogicalVector subsets (vec.length(), 0);

  for(int i = 0; i < filter.length(); i++){
    for(int j = 0; j< filter.length(); j++){
      if( filter[j] == i){
        subsets[ filter[i] ] = 1;
      }
    }
  }
  return vec[ !subsets ];
}


arma::uvec uvec_diff( arma::uvec first, arma::uvec second){

  NumericVector it_first(first.begin(), first.end());
  NumericVector it_second(second.begin(), second.end());

  NumericVector res = negative_subset(it_first, it_second);

  arma::uvec result(as<arma::uvec>(res));

  return result;
}

iTree init_iTree( const arma::mat& X,
                  const int& max_tree_depth,
                  const int& extension_level,
                  const std::string& type ){
  iTree new_tree;
  new_tree.X = X;

  int max_nodes_per_tree_depth = max_nodes(max_tree_depth);

  new_tree.Id_mat = mat( max_nodes_per_tree_depth, 5, fill::zeros );
  // "TerminalID", "Type","Size","Left", "Right"
  if(type == "vanilla"){
    new_tree.pred_mat = mat( max_tree_depth, 2, fill::zeros );
  }
  else if(type == "extended"){
    new_tree.pred_mat = mat(max_tree_depth, 2*X.n_cols, fill::zeros );
  }

  return new_tree;
}



void recursion( const arma::uvec& index,
                const int& max_depth,
                const int& node_index,
                const int& current_depth,
                iTree& used_tree,
                const int& extension_level,
                const std::string& type ){
// dont sample from columns which only contain duplicates
  arma::mat X = used_tree.X.rows(index);

  bool found_duplicates = false;
  int comp = X.n_cols;

  for( int i = 0; i < comp; i++ ){
    arma::vec un_id = unique( X.col(i));
   if( un_id.n_elem == 1){
     found_duplicates = true;
     break;
   }
  }

  if(current_depth >= max_depth || found_duplicates || index.n_elem == 1){
    used_tree.Id_mat(node_index, 2 ) = -1;
    used_tree.Id_mat(node_index, 3 ) = index.n_elem;
    return;
  }

  Tree_node res = split_cpp( X, extension_level, type );
  // "TerminalID", "Type","Size","Left", "Right"

  int nodeLeft =   2 * node_index ;
  int nodeRight =  2 * node_index + 1;
  used_tree.Id_mat(node_index, 3 ) = nodeLeft;
  used_tree.Id_mat(node_index, 4) = nodeRight;
  used_tree.Id_mat(node_index, 1) = 1;


  used_tree.pred_mat.row(current_depth) = res.coef;

  arma::uvec left_filter = res.filter;
  arma::uvec right_filter = uvec_diff(index, res.filter);

  recursion( left_filter, max_depth, nodeLeft, current_depth + 1,
             used_tree, extension_level, type );
  recursion( right_filter, max_depth, nodeRight, current_depth + 1,
             used_tree, extension_level, type );
}

// [[Rcpp::export]]
Rcpp::List exporter( const arma::mat& X_in, const int& ext = 10,
                     const int& max_depth = 10,
                     const std::string& type = "vanilla"){


  arma::vec almost_index = linspace(0,(X_in.n_rows-1), X_in.n_rows);
  arma::uvec actually_index = conv_to<uvec>::from(almost_index);

  iTree res = init_iTree(X_in, max_depth, ext, type );

  recursion(actually_index, max_depth, 0, 0, res, ext, type);

  return Rcpp::List::create( res.Id_mat, res.pred_mat  );
}


// [[Rcpp::export]]
Rcpp::List extended_tester(const arma::mat& X, const int& ext){

  Tree_node ins_res = extended_comparison(X, ext);
  iTree res =  init_iTree(X, 8, 1, "extended");

  // res.pred_mat.row(1) = ins_res.coef;

  return Rcpp::List::create( res.pred_mat );
}





