// updateW_discovR_grouplasso.cpp //

// rcpp script for conditional update of W for discovR

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double soft (double x, double lambda){	
	double x2;
	double result;
	double sgn;

	x2 = std::abs(x) - lambda;

	if (x2 < 0){ 
		x2 = 0;
	} 

	if (x < 0){
		sgn = -1;
	} else {
		sgn = 1;
	}

	result = sgn * x2;

	return result;
}

// [[Rcpp::export]]
double losscal (arma::mat Y, 
	arma::mat X, 
	arma::mat W, 
	arma::mat Px, 
	arma::mat Py, 
	double alpha, 
	arma::vec lasso_w, 
	double ridge_w,
	arma::vec lasso_y, 
	double ridge_y){

	arma::mat lasso_w_mat;
	arma::mat lasso_y_mat;
  	int r;
  	int R = lasso_w.size();
  	double pca_loss;
  	double reg_loss;
  	double result;

	lasso_w_mat = arma::zeros( W.n_rows, W.n_cols );

	for (r = 0; r < R; r++){
		lasso_w_mat(arma::span(), arma::span(r)).fill(lasso_w(r));
	}

	lasso_y_mat = arma::zeros( Py.n_rows, Py.n_cols );	
	
	for (r = 0; r < R; r++){
		lasso_y_mat(arma::span(), arma::span(r)).fill(lasso_y(r));
	}

	pca_loss = accu(square(X - X * W * trans(Px))) / accu(square(X));
	reg_loss = accu(square(Y - X * W * trans(Py))) / accu(square(Y));

    result = (alpha) * reg_loss +
      (1 - alpha) * pca_loss +
      accu(lasso_w_mat % abs(W)) +
      ridge_w * accu(square(W)) +
      accu(lasso_y_mat % abs(Py)) +
      ridge_y * accu(square(Py));
    
    return (result);

}




// [[Rcpp::export]]
Rcpp::List updateW_cpp_ridge (
	arma::mat X,
  	arma::mat Y,
	arma::mat W,
	arma::mat Px,
  	arma::mat Py,
	int R,
 	double alpha,
	arma::vec lasso_w,
	double ridge_w,
	arma::vec lasso_y,
	double ridge_y){

  	arma::mat W_new;
  	arma::mat W_r;

  	int J;
	arma::mat r_h;
	arma::mat s_h;
	arma::mat soft_y_h;
	arma::mat soft_x_h;
	arma::mat h_soft;
	double h_softened;

	arma::mat denom;
	double W_rh_new;
	arma::mat W_r_new;

	double loss0;
	double loss1;
	arma::vec loss_hist;
	int iter = 0;
  
    int r;

    // defining initial objects
    double SSY = accu(square(Y));
    double SSX = accu(square(X));

    
  	// initial loss calculation #
  	loss0 = losscal(Y, X, W, Px, Py, 
  				   alpha, lasso_w, ridge_w, 
                   lasso_y, ridge_y);
  

    loss_hist = arma::zeros( (R+1) * W.n_rows );
    
    loss_hist(iter) = loss0;
    
    iter = iter + 1;
  
    W_new = W;

    J = X.n_cols;

   	for (r = 0; r < R; r++){
   		Rcpp::checkUserInterrupt();

      	W_r = arma::vectorise(W_new.col(r));

		W_r_new = W_r;

      // elementwise coordinate descent
      for (int h = 0; h < J; h++){

		r_h = Y - X * W_new * trans(Py) + (X.col(h) * W_r(h,0)) * trans(arma::vectorise(Py(arma::span(), arma::span(r))));

		s_h = X * arma::vectorise(Px(arma::span(), arma::span(r))) - (X * arma::vectorise(W_new(arma::span(), arma::span(r))) - X.col(h) * W_r(h,0));

		soft_y_h = (2 * alpha / SSY) * (r_h * arma::vectorise(Py(arma::span(), arma::span(r))));

		soft_x_h = (2 * (1 - alpha) / SSX) * s_h;

		h_soft = trans(soft_y_h + soft_x_h) * X.col(h);

		double to_soft_h = h_soft(0,0);

		h_softened = soft(to_soft_h, lasso_w(r));

		denom = (((2 * alpha / SSY) * (trans(arma::vectorise(Py(arma::span(), arma::span(r)))) * arma::vectorise(Py(arma::span(), arma::span(r)))) +
				(2 * (1 - alpha) / SSX)) * (trans(X.col(h)) * X.col(h))) + 
				(2 * ridge_w);

		double denom_double = denom(0,0);

		W_rh_new = h_softened / denom_double;


        W_r_new(h,0) = W_rh_new;

		// updating the W_new matrix 
        W_new.col(r) = W_r_new;

        
        // loss increase checking
	  	loss1 = losscal(Y, X, W_new, Px, Py, 
	  		alpha, lasso_w, ridge_w,
	  		lasso_y, ridge_y);

	    // if loss increased:
      	if ((loss1 - loss0) > 1e-10){
        	Rcpp::stop("loss increase via coordinate descent");
      	}

      	// otherwise, save the loss and continue
      	loss_hist(iter,0) = loss1;

      	iter = iter + 1;
      
      	loss0 = loss1;

        W_r = arma::vectorise(W_new.col(r));


      }

    }

  Rcpp::List result;
  
  result["W_new"] = W_new;
  result["W_r"] = W_r;
  // result["r"] = r;
  // result["k"] = k;
  // result["block_index"] = block_index;
  // result["r_k"] = r_k;
  // result["s_k"] = s_k;
  // result["soft_y"] = soft_y;
  // result["soft_x"] = soft_x;
  // result["k_soft"] = k_soft;
  // result["k_soft_colsums"] = k_soft_colsums;
  result["h_softened"] = h_softened;
  // result["r_kh"] = r_kh;
  // result["s_kh"] = s_kh;
  // result["soft_y_h"] = soft_y_h;
  // result["soft_x_h"] = soft_x_h;
  // result["kh_soft"] = kh_soft;
  // result["kh_softened"] = kh_softened;
  // result["denom"] = denom;
  result["W_rh_new"] = W_rh_new;
  result["W_r_new"] = W_r_new;
  // result["loss1"] = loss1;
  result["loss_hist"] = loss_hist;

  return result;
}

