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
double grouplasso_norm (
	arma::mat x,
	Rcpp::List blockindex,
	int R,
 	arma::vec grouplasso){

	double ssq;
	arma::mat result; 
	int nblock;
	double l2norm;

	nblock = blockindex.size();
	
	int r;
	int i;

	l2norm = 0;

	for (r = 0; r < R; r++){
		for (i = 0; i < nblock; i++){
			Rcpp::NumericVector block;
			int block_min;
			int block_max;
			arma::mat inblock;
			double blocksize;

			block = blockindex(i);
			block_min = min(block);
			block_max = max(block);

			block_min = block_min - 1;
			block_max = block_max - 1;

			inblock = x(arma::span(block_min, block_max), r);

			ssq = sqrt(as_scalar(trans(inblock) * inblock));

			blocksize = block.size();

			l2norm = l2norm + grouplasso(r) * ssq * sqrt(blocksize);
		}
	}

	return l2norm;
}


// [[Rcpp::export]]
double losscal (arma::mat Y, 
	arma::mat X, 
	arma::mat W, 
	arma::mat Px, 
	arma::mat Py, 
	double alpha, 
	arma::vec lasso_w, 
	arma::vec grouplasso_w, 
	arma::vec lasso_y, 
	double ridge_y, 
	Rcpp::List blockindex){

	arma::mat lasso_w_mat;
	arma::mat lasso_y_mat;
  	int r;
  	int R = lasso_w.size();
  	double pca_loss;
  	double reg_loss;
  	double l2norm;
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

    l2norm = grouplasso_norm(W, blockindex, W.n_cols, grouplasso_w);
    
    result = (alpha) * reg_loss +
      (1 - alpha) * pca_loss +
      accu(lasso_w_mat % abs(W)) +
      l2norm +
      accu(lasso_y_mat % abs(Py)) +
      ridge_y * accu(square(Py));
    
    return (result);

}


// [[Rcpp::export]]
Rcpp::List updateW_cpp (
	arma::mat X,
  	arma::mat Y,
	arma::mat W,
	arma::mat Px,
  	arma::mat Py,
	int R,
 	double alpha,
	arma::vec lasso_w,
	arma::vec grouplasso_w,
	arma::vec lasso_y,
	double ridge_y,
	Rcpp::List blockindex){

  	arma::mat W_new;
  	double nblock;

    arma::mat Xk;
  	arma::mat Wk;
    arma::mat r_k;
    arma::mat s_k;
    arma::mat soft_y;
    arma::mat soft_x;
    arma::mat k_soft;
    arma::mat k_soft_colsums;
    arma::mat k_softened;
  	int J_k;
    Rcpp::IntegerVector block_index;
	arma::mat r_kh;
	arma::mat s_kh;
	arma::mat soft_y_h;
	arma::mat soft_x_h;
	arma::mat kh_soft;
	double kh_softened;
	arma::mat denom;
	double Wkh_new;
	arma::mat Wk_new;
	double loss0;
	double loss1;
	arma::vec loss_hist;
	int iter = 0;
  
    int r;
    int k;

    nblock = blockindex.size();

    // defining initial objects
    double SSY = accu(square(Y));
    double SSX = accu(square(X));

    
  	// initial loss calculation #
  	loss0 = losscal(Y, X, W, Px, Py, 
  				   alpha, lasso_w, grouplasso_w, 
                   lasso_y, ridge_y, blockindex);
  
    loss_hist = arma::zeros( (R+1) * (nblock+1) * W.n_rows );
    
    loss_hist(iter) = loss0;
    
    iter = iter + 1;
  
    W_new = W;

   	for (r = 0; r < R; r++){

   		for (k = 0; k < nblock; k++){
   			
      	 // specifying Xk
	   	 Rcpp::NumericVector block;
		  int block_min;
		  int block_max;

		  block = blockindex(k);
		  block_min = min(block);
		  block_max = max(block);

		  block_min = block_min - 1;
		  block_max = block_max - 1;

        Xk = X(arma::span(), arma::span(block_min, block_max));

		J_k = block.size();
      
        // specifying Wk
		Wk = arma::vectorise(W_new(arma::span(block_min, block_max), r));

		// calculating what should go in the soft thresholding
        r_k = Y - X * W_new * trans(Py) + Xk * Wk * trans(arma::vectorise(Py(arma::span(), arma::span(r))));

        s_k = X * arma::vectorise(Px(arma::span(), arma::span(r))) - (X * arma::vectorise(W_new(arma::span(), arma::span(r))) - Xk * Wk);

        soft_y = (2 * alpha / SSY) * (r_k * arma::vectorise(Py(arma::span(), arma::span(r))));
        
        soft_x = (2 * (1 - alpha) / SSX) * s_k;

        k_soft = arma::diagmat(soft_y + soft_x) * Xk;

		k_soft_colsums = sum(k_soft, 0); // dim = 0 does column sums
		
		// k_softened is the result of the softening (column vector)
		k_softened = trans(k_soft_colsums);


		// soft-thesholding for group k, via for-loop
        for (int i = 0; i < J_k; i++){

          double to_soft = k_soft_colsums(0, i);

          double softened = soft(to_soft, lasso_w(r));

          k_softened(i, 0) = softened;
      	}

      	// calculating the l2 norm
        double groupnorm = sqrt(accu(square(k_softened)));

        // blockwise sparsity 
        if (groupnorm < (grouplasso_w(r) * sqrt(J_k))){
        	
        	Wk_new = arma::zeros( Wk.n_rows, Wk.n_cols );

		    // update the entire W_new
		    W_new(arma::span(block_min, block_max),r) = Wk_new;

		    // loss increase checking
		  	loss1 = losscal(Y, X, W_new, Px, Py, 
  				   alpha, lasso_w, grouplasso_w, 
                   lasso_y, ridge_y, blockindex);
  
		      	// if loss increased:
          	if (loss1 > loss0){
            	Rcpp::stop("loss increase after block sparsified");
          	}

          	// otherwise, save the loss and continue
          	loss_hist(iter) = loss1;

          	iter = iter + 1;
          
          	loss0 = loss1;

        } else {
          
          Xk = X(arma::span(), arma::span(block_min, block_max));
          
          Wk = arma::vectorise(W_new(arma::span(block_min, block_max), r));

          // elementwise coordinate descent
          for (int h = 0; h < J_k; h++){

			r_kh = Y - X * W_new * trans(Py) + (Xk.col(h) * Wk(h,0)) * trans(arma::vectorise(Py(arma::span(), arma::span(r))));

			s_kh = X * arma::vectorise(Px(arma::span(), arma::span(r))) - (X * arma::vectorise(W_new(arma::span(), arma::span(r))) - Xk.col(h) * Wk(h,0));

			soft_y_h = (2 * alpha / SSY) * (r_kh * arma::vectorise(Py(arma::span(), arma::span(r))));

			soft_x_h = (2 * (1 - alpha) / SSX) * s_kh;

			kh_soft = trans(soft_y_h + soft_x_h) * Xk.col(h);

			double to_soft_kh = kh_soft(0,0);

			kh_softened = soft(to_soft_kh, lasso_w(r));

			denom = (((2 * alpha / SSY) * (trans(arma::vectorise(Py(arma::span(), arma::span(r)))) * arma::vectorise(Py(arma::span(), arma::span(r)))) +
					(2 * (1 - alpha) / SSX)) * (trans(Xk.col(h)) * Xk.col(h))) + 
					(grouplasso_w(r) * sqrt(J_k) / sqrt(accu(square(Wk))));

			double denom_double = denom(0,0);

			Wkh_new = kh_softened / denom_double;

			Wk_new = Wk;

            Wk_new(h,0) = Wkh_new;

 			// updating the W_new matrix 
            W_new(arma::span(block_min, block_max),r) = Wk_new;

            
            // loss increase checking
		  	loss1 = losscal(Y, X, W_new, Px, Py, 
  				   alpha, lasso_w, grouplasso_w, 
                   lasso_y, ridge_y, blockindex);
  
		    // if loss increased:
          	if ((loss1 - loss0) > 1e-10){
            	Rcpp::stop("loss increase via coordinate descent");
          	}

          	// otherwise, save the loss and continue
          	loss_hist(iter) = loss1;

          	iter = iter + 1;
          
          	loss0 = loss1;

            Wk = arma::vectorise(W_new(arma::span(block_min, block_max), r));


          }

        }

	  }
	}

  Rcpp::List result;
  
  result["W_new"] = W_new;
  // result["Wk"] = Wk;
  // result["r"] = r;
  // result["k"] = k;
  // result["block_index"] = block_index;
  // result["r_k"] = r_k;
  // result["s_k"] = s_k;
  // result["soft_y"] = soft_y;
  // result["soft_x"] = soft_x;
  // result["k_soft"] = k_soft;
  // result["k_soft_colsums"] = k_soft_colsums;
  // result["k_softened"] = k_softened;
  // result["r_kh"] = r_kh;
  // result["s_kh"] = s_kh;
  // result["soft_y_h"] = soft_y_h;
  // result["soft_x_h"] = soft_x_h;
  // result["kh_soft"] = kh_soft;
  // result["kh_softened"] = kh_softened;
  // result["denom"] = denom;
  // result["Wkh_new"] = Wkh_new;
  // result["Wk_new"] = Wk_new;
  // result["loss1"] = loss1;
  result["loss_hist"] = loss_hist;

  return result;
}

