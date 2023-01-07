// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>
//using boost::math::policies::policy;
//using boost::math::tools::newton_raphson_iterate;
//using boost::math::tools::halley_iterate; //
//using boost::math::tools::eps_tolerance; // Binary functor for specified number of bits.
//using boost::math::tools::bracket_and_solve_root;
//using boost::math::tools::toms748_solve;

#include <boost/math/special_functions/next.hpp> // For float_distance.
#include <tuple> // for std::tuple and std::make_tuple.
#include <boost/math/special_functions/cbrt.hpp> // For boost::math::cbrt.
#include <boost/math/special_functions/pow.hpp> // boost::math::pow<5,double>

using namespace Rcpp;
using namespace std;
using namespace boost::math::tools;


//' @useDynLib geneplot
//' @importFrom Rcpp sourceCpp

// Note that these two functions produce different results if given a vector of decimal numbers
// [[Rcpp::export]]
double double_sum (NumericVector nu) {
    //int a = sum(nu);
    return(sum(nu));
}
// [[Rcpp::export]]
int int_sum (NumericVector nu) {
    //int a = sum(nu);
    return(sum(nu));
}

// Function for calculating the probabilities for the DCM distribution at a
// single locus for a single population
NumericVector rcpp_calc_probs(NumericVector nu, bool leave_one_out=false) {
    // When calculating the sum of the nu values, must use a DOUBLE not an int,
    // because if using an int it will not add up any decimal parts of the nu values
    double n = sum(nu);
    int k = nu.length();

    // Calculate the probabilities of all genotypes at this locus, using
    // Hardy-Weinberg formula
    NumericVector probs(k+(k*(k-1)/2));
    int idx=0;
    for (int ii=0; ii < k; ii++) {
        for (int jj=ii; jj < k; jj++) {
            if (leave_one_out) {
                // For LOO version of SPDF saddlepoint, we want to remove copies
                // of the alleles from the observed + prior before calculating
                // the genotype probabilities. We remove as many copies as we
                // can, so if there aren't any copies in the observed data then
                // we don't remove any for calculating the LOO genotype
                // probability.
                //
                // Note that nu[idx] == 1 would mean that there were no alleles
                // of that type observed, and nu[idx] == 2 would mean fewer than
                // 2 alleles of that type observed, because the prior will
                // always be bigger than 0 so prior + nobserved will always be
                // bigger than 1 if nobserved is at least 1. That's why we only
                // remove 2 copies if nu > 2 and only remove 1 if nu > 1
                if (jj == ii) {
                    if (nu[ii] > 2) probs[idx] = (nu[ii]-2)*(nu[ii]-1)/((n-2)*(n-1));
                    else if (nu[ii] > 1) probs[idx] = (nu[ii]-1)*nu[ii]/((n-1)*n);
                    else probs[idx] = (nu[ii])*(nu[ii]+1)/(n*(n+1));
                } else {
                    if (nu[ii] > 1 && nu[jj] > 1) probs[idx] = 2*(nu[ii]-1)*(nu[jj]-1)/((n-2)*(n-1));
                    else if (nu[ii] > 1 && nu[jj] <= 1) probs[idx] = 2*(nu[ii]-1)*nu[jj]/((n-1)*n);
                    else if (nu[ii] <= 1 && nu[jj] > 1) probs[idx] = 2*nu[ii]*(nu[jj]-1)/((n-1)*n);
                    else probs[idx] = 2*nu[ii]*nu[jj]/(n*(n+1));
                }
            } else {
                if (jj == ii) {
                    probs[idx] = nu[ii]*nu[jj]/n/(n+1) + nu[ii]/n/(n+1);
                } else {
                    probs[idx] = 2*nu[ii]*nu[jj]/n/(n+1);
                }
            }
            idx++;
        }
    }

    // Rcout << "probs : " << probs << "\n";
    return probs;
}

// Function for calculating the probabilities for the DCM distribution at a
// single locus for one or two populations
// [[Rcpp::export]]
List rcpp_calc_dist(NumericMatrix nu, bool leave_one_out = false, bool differenceLGPs = false) {
    NumericVector values;
    NumericVector probs;

    // If nu has 2 rows, then treat each row separately - the first is the base
    // population nu values - for the probabilities, the second is the
    // comparison population nu values
    if (nu.nrow() > 1)
    {
        if (differenceLGPs) {
            NumericVector nu1 = nu.row(0); // Base population
            NumericVector nu2 = nu.row(1); // Comparison population, we want to know how well A individuals fit into B

            // Calculate the probabilities for populations A and B, then remove
            // any that are zero before calculating logs. Do this BEFORE doing
            // any sorting by size, so that the corresponding entries in the two
            // populations are still aligned with each other

            NumericVector probsA = rcpp_calc_probs(nu1);
            NumericVector probsA_LOO = rcpp_calc_probs(nu1,leave_one_out=leave_one_out);
            NumericVector probsB = rcpp_calc_probs(nu2);

            if (is_true(any(probsA_LOO == 0)) || is_true(any(probsB == 0)))
            {
                for (int idx=0; idx < nu1.length();) {
                    if (probsA_LOO[idx] == 0 || probsB[idx] == 0) {
                        probsA.erase(idx);
                        probsA_LOO.erase(idx);
                        probsB.erase(idx);
                    } else {
                        idx++;
                    }
                }
            }

            values = log(probsA_LOO) - log(probsB);
            probs = probsA; // the probs are the NON-LOO probs for popA
        } else {
            NumericVector nu1 = nu.row(0); // Base population
            NumericVector nu2 = nu.row(1); // Comparison population, we want to know how well B fits into A

            // Rcout << "nu1 : " << nu1 << "\n";
            // Rcout << "nu2 : " << nu2 << "\n";

            probs = rcpp_calc_probs(nu2);

            // Rcout << "probs from nu2 : " << probs << "\n";

            // Calculate the probabilities for population 1, then remove any
            // that are zero before calculating logs. Also remove the
            // corresponding entries from the pop2 probs -- do this BEFORE doing
            // any sorting by size, so that the corresponding entries in the two
            // populations are still aligned with each other
            NumericVector probs1 = rcpp_calc_probs(nu1);

            // Rcout << "probs from nu1 : " << probs1 << "\n";

            if (is_true(any(probs1 == 0)))
            {
                for (int idx=0; idx < nu1.length();) {
                    if (probs1[idx] == 0) {
                        probs1.erase(idx);
                        probs.erase(idx);
                    } else {
                        idx++;
                    }
                }
            }
            values = log(probs1);
            // Rcout << "values from nu1 : " << values << "\n";
        }

                // Rcout << "unsorted values : " << values << "\n";

        // sort the values in increasing order and the probabilities in the same order
        // First create a vector of indices
        IntegerVector sort_order = seq_along(values) - 1;
        // Then sort that vector by the values of y
        std::sort(sort_order.begin(), sort_order.end(), [&](int i, int j){return values[i] < values[j];});

        // Rcout << "sort order : " << sort_order << "\n";
        values = values[sort_order];

                // Rcout << "unsorted probs : " << probs << "\n";
        probs = probs[sort_order];
                // Rcout << "sorted probs : " << probs << "\n";
                // Rcout << "sorted values : " << values << "\n";
        //         Rcout << "sorting : " << sort_order << "\n";
        //         Rcout << "sorted probs : " << probs << "\n";
    }
    else
    {
        probs = rcpp_calc_probs(nu.row(0));

        if (leave_one_out) {
            NumericVector probs_for_values = rcpp_calc_probs(nu.row(0),leave_one_out=true);

            // remove any zero probabilities
            if (is_true(any(probs_for_values == 0))) {
                for (int idx=0; idx < probs.length();) {
                    if (probs_for_values[idx] == 0) {
                        probs_for_values.erase(idx);
                        probs.erase(idx);
                    } else {
                        idx++;
                    }
                }
            }

            // Rcout << "probs : " << probs << "\n";
            // Rcout << "probs_for_values : " << probs_for_values << "\n";

            // sort the probabilities in increasing order so that when the table
            // of frequencies is created, it will be in the same order as unique(probs)
            // First create a vector of indices
            IntegerVector sort_order = seq_along(probs) - 1;
            // Then sort that vector by the values of y
            std::sort(sort_order.begin(), sort_order.end(), [&](int i, int j){return probs[i] < probs[j];});
            probs = probs[sort_order];
            probs_for_values = probs_for_values[sort_order];

            values = log(probs_for_values);
        } else {
            // remove any zero probabilities
            if (is_true(any(probs == 0))) {
                // Remove genotypes with zero probability
                // Note: do NOT include an iterator in the for loop for erasing zero
                // elements because when we apply erase that will automatically increment it
                // so don't want to use ++it in that circumstance
                NumericVector::iterator it;
                for(it = probs.begin(); it <= probs.end();) {
                    if (*it == 0){
                        it = probs.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            // sort the probabilities in increasing order
            probs.sort();
            values = log(probs);
        }
    }

    // Also calculate different powers of the values
    NumericVector r0power_values(values.length(), 1);
    NumericVector r1power_values = values;
    NumericVector r2power_values = pow(values,2);
    NumericVector r3power_values = pow(values,3);

    List L = List::create(Named("probs") = probs, Named("values") = values,
                          Named("r0power_values") = r0power_values,
                          Named("r1power_values") = r1power_values,
                          Named("r2power_values") = r2power_values,
                          Named("r3power_values") = r3power_values);

    return L;
}

// Function for calculating the minimum value of the DCM distribution at a single locus
// [[Rcpp::export]]
double rcpp_calc_min_beta(NumericMatrix nu_mat) {

    // If nu values are provided for multiple populations, assume that the first
    // one determines the minimum possible value i.e. reduce nu to just that
    // population
    NumericVector nu = nu_mat.row(0);

    // Use a DOUBLE to store the sum of nu values, because if using an int it will
    // lose the decimal parts of the nu values
    double n = sum(nu);
    double min_val = 0;
    if (nu.length() > 1) {
        // if the least common allele type has a frequency that is less than
        // twice the frequency of the next least common allele, then the minimum
        // log-prob will be the log-prob of the homozygote of the least common
        // allele otherwise it will be the log-prob of the heterozygote of the
        // two least common alleles
        nu.sort();
        double min_nu = nu[0];
        double next_min_nu = nu[1];
        if (min_nu + 1 < 2*next_min_nu) {
            min_val = log(min_nu*(min_nu+1)/(n*(n+1)));
        } else {
            min_val = log(2*min_nu*next_min_nu/(n*(n+1)));
        }
    }
    return min_val;
}

// Function for calculating the maximum value of the DCM distribution at a single locus
// [[Rcpp::export]]
double rcpp_calc_max_beta(NumericMatrix nu_mat) {
    // If nu values are provided for multiple populations, assume that the first
    // one determines the maximum possible value i.e. reduce nu to just that
    // population
    NumericVector nu = nu_mat.row(0);

    // Use a DOUBLE to store the sum of nu values, because if using an int it will
    // lose the decimal parts of the nu values
    double n = sum(nu);
    double max_val = 0;
    if (nu.length() > 1) {
        nu.sort(true);
        // if the most common allele type has a frequency that is more than
        // twice the frequency of the next most common allele, then the maximum
        // log-prob will be the log-prob of the homozygote of the most common
        // allele; otherwise it will be the log-prob of the heterozygote of the
        // two most common alleles
        double max_nu = nu[0];
        double next_max_nu = nu[1];
        if ((max_nu + 1) > 2*next_max_nu) {
            max_val = log(max_nu*(max_nu+1)/(n*(n+1)));
        } else {
            max_val = log(2*max_nu*next_max_nu/(n*(n+1)));
        }
    }
    return max_val;
}

// Function for combining the DCM distributions from multiple loci
// [[Rcpp::export]]
List rcpp_calc_multi_locus_dist(List nutab, bool leave_one_out = false, bool differenceLGPs = false) {

    List dist (nutab.length());
    double min_dist = 0;
    double max_dist = 0;

    for (int ll=0; ll < nutab.length(); ll++) {
        // Copying the value of row 0 to the vector v
        NumericMatrix nu = nutab[ll];

        // Rcout << "leave_one_out : " << leave_one_out << "\n";
        // Rcout << "differenceLGPs : " << differenceLGPs << "\n";
        List loc_dist = rcpp_calc_dist(nu, leave_one_out=leave_one_out, differenceLGPs=differenceLGPs);
        dist[ll] = loc_dist;

        // when the values of the distribution are the differences between the
        // LGPs with respect to the two populations, or leave_one_out is being
        // used, then need to take the min and max from the calculated values
        // instead of usual min/max calc
        if (differenceLGPs || leave_one_out) {
            NumericVector values = loc_dist["values"];
            min_dist += min(values);
            max_dist += max(values);
        } else {
            // can simply add the locus values because they are logs of probabilities
            min_dist += rcpp_calc_min_beta(nu);
            max_dist += rcpp_calc_max_beta(nu);
        }

    }

    List L = List::create(_["dist"] = dist, _["min"] = min_dist, _["max"] = max_dist);

    return L;
}

// [[Rcpp::export]]
double rcpp_calc_multi_locus_K(List dist, double s) {

    double K = 0;

    for (int ll=0; ll < dist.length(); ll++) {

        List loc_dist = dist[ll];
        NumericVector probs = loc_dist["probs"];
        NumericVector values = loc_dist["values"];
        K += log(sum(exp(s*values)*probs));
    }

    return K;
}

// [[Rcpp::export]]
double rcpp_calc_multi_locus_K1(List dist, double s) {

    double K1 = 0;

    for (int ll=0; ll < dist.length(); ll++) {

        List loc_dist = dist[ll];
        NumericVector probs = loc_dist["probs"];
        NumericVector values = loc_dist["values"];

        // Don't try to calculate this quicker because even for one pop, if
        // leave_one_out = true then values != log(probs)
        NumericVector exp_svals = exp(s*values);
        K1 += sum(values*exp_svals*probs)/sum(exp_svals*probs);
    }

    return K1;
}

// [[Rcpp::export]]
double rcpp_calc_multi_locus_K2(List dist, double s) {

    double K2 = 0;
    double K2_component;
    double M0;
    double M1;
    double M2;
    NumericVector exp_svals;

    for (int ll=0; ll < dist.length(); ll++) {

        List loc_dist = dist[ll];
        NumericVector probs = loc_dist["probs"];
        NumericVector values = loc_dist["values"];
        NumericVector r2_values = loc_dist["r2power_values"];

        NumericVector exp_svals = exp(s*values);

        // Don't try to calculate this quicker because even for one pop, if
        // leave_one_out = true then values != log(probs)
        M0 = sum(exp_svals*probs);
        M1 = sum(values*exp_svals*probs);
        M2 = sum(r2_values*exp_svals*probs);

        // Rcout << "M0 : " << M0 << "\n";
        // Rcout << "M1 : " << M1 << "\n";
        // Rcout << "M2 : " << M2 << "\n";

        K2_component = M2/M0 - pow(M1/M0,2);

        // Rcout << "K2 component : " << K2_component << "\n";

        if (!std::isnan(K2_component) && K2_component < 0) K2_component = 0;

        K2 += K2_component;
    }

    return K2;
}
// [[Rcpp::export]]
double rcpp_calc_multi_locus_K3(List dist, double s) {

    double K3 = 0;
    double K3_component;
    double M0;
    double M1;
    double M2;
    double M3;

    for (int ll=0; ll < dist.length(); ll++) {

        List loc_dist = dist[ll];
        NumericVector probs = loc_dist["probs"];
        NumericVector values = loc_dist["values"];
        NumericVector r2_values = loc_dist["r2power_values"];
        NumericVector r3_values = loc_dist["r3power_values"];

        NumericVector exp_svals = exp(s*values);

        M0 = sum(exp_svals*probs);
        M1 = sum(values*exp_svals*probs);
        M2 = sum(r2_values*exp_svals*probs);
        M3 = sum(r3_values*exp_svals*probs);

        // Rcout << "M0 : " << M0 << "\n";
        // Rcout << "M1 : " << M1 << "\n";
        // Rcout << "M2 : " << M2 << "\n";
        // Rcout << "M3 : " << M3 << "\n";

        K3_component = M3/M0 - 3*M1*M2/pow(M0,2) + 2*pow(M1/M0,3);

        // Rcout << "K3 component : " << K3_component << "\n";

        K3 += K3_component;
    }

    return K3;
}

// [[Rcpp::export]]
double rcpp_calc_mu(List dist) {
    double mu = rcpp_calc_multi_locus_K1(dist, 0);
    return mu;
}

// [[Rcpp::export]]
double rcpp_calc_sigma(List dist) {
    double sigma = sqrt(rcpp_calc_multi_locus_K2(dist, 0));
    return sigma;
}

// [[Rcpp::export]]
double rcpp_calc_what(double x, double sh, List dist, bool use_x = false) {
    double wh_inner;
    double wh;

    if (use_x) {
        wh_inner = sh*x - rcpp_calc_multi_locus_K(dist, sh);
    } else {
        wh_inner = sh*rcpp_calc_multi_locus_K1(dist, sh) - rcpp_calc_multi_locus_K(dist, sh);
    }

    if (wh_inner < 0)
    {
        wh_inner = 1e-16;
    }
    if (sh < 0) {
        wh = -sqrt(2*wh_inner);
    } else {
        wh = sqrt(2*wh_inner);
    }

    return wh;
}

// [[Rcpp::export]]
double rcpp_calc_uhat(double sh, List dist) {
    return sh*sqrt(rcpp_calc_multi_locus_K2(dist, sh));
}

struct shat_functor_deriv
{ // Functor also returning 1st derivative.
    shat_functor_deriv(double const& to_find_root_of, List dist) : _x(to_find_root_of), _dist(dist)
    { // Constructor stores value to find root of, and distribution to use
        // for example: calling shat_functor_deriv(x) to use to get K1 root of x.
    }
    std::pair<double, double> operator()(double const& s)
    {
        // Return both f(s) and f'(s)
        double fs = rcpp_calc_multi_locus_K1(_dist, s) - _x;    // Difference (estimate K1(s) - value).
        double ds = rcpp_calc_multi_locus_K2(_dist, s);    // 1st derivative (K2(s)).
        return std::make_pair(fs, ds);   // 'return' both fs and ds.
    }
private:
    double _x;       // Store value to find the K1 root of
    List _dist;      // Store distribution to use
};

// [[Rcpp::export]]
double rcpp_calc_shat(double x, List dist)
{
    double lower = -50;
    double upper = 150;

    double f_lower = rcpp_calc_multi_locus_K1(dist, lower) - x;
    double f_upper = rcpp_calc_multi_locus_K1(dist, upper) - x;

    if (isnan(f_lower))
    {
        while(lower < upper - 10 && isnan(f_lower))
        {
            lower += 10;
            f_lower = rcpp_calc_multi_locus_K1(dist, lower) - x;
        }
    }

    if (isnan(f_upper))
    {
        while(upper > lower + 10 && isnan(f_upper))
        {
            upper -= 10;
            f_upper = rcpp_calc_multi_locus_K1(dist, upper) - x;
        }
    }

    if (isnan(f_lower) || isnan(f_upper)) {
        double result = nan("0");
        return result;
    }

    // Rcout << "lower : " << lower << "\n";
    // Rcout << "upper : " << upper << "\n";
    // Rcout << "f_lower : " << f_lower << "\n";
    // Rcout << "f_upper : " << f_upper << "\n";

    // return shat using 1st derivative of K1 and Newton_Raphson.
    double guess = 0;                                        // Rough guess is zero.
    double min = lower;                                      // Minimum possible value is lower
    double max = upper;                                      // Maximum possible value is upper
    const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type double.
    int get_digits = static_cast<int>(digits * 0.6);         // Accuracy doubles with each step, so stop when we have
    // just over half the digits correct.
    const std::uintmax_t maxit = 1000;
    std::uintmax_t it = maxit;
    double result = newton_raphson_iterate(shat_functor_deriv(x, dist), guess, min, max, get_digits, it);
    return result;
}

// Note that this function does NOT have logten option because it should only
// be called from the Rcpp functions
// rcpp_calc_fhat calls it AFTER converting x if logten = true
// rcpp_calc_Qhat calls it ONLY with logten = false
// So either way, it is only used with logten = false
// Also, do NOT export this function, because it doesn't have the logten conversion in it
double rcpp_calc_fhat_pt(double x, List dist, double min_dist, double max_dist, double mean_dist) {
    double sh;
    double fh;
    double K2;
    if (x >= max_dist || x <= min_dist) {
        fh = 0;
    } else {
        if (x == mean_dist) {
            sh = 0;
        } else {
            sh = rcpp_calc_shat(x, dist);
        }

        K2 = rcpp_calc_multi_locus_K2(dist, sh);
        if (isnan(K2)) fh = nan("0");
        else if (K2 == 0) fh = 0;
        else fh = exp(rcpp_calc_multi_locus_K(dist, sh) - sh*x)/sqrt(2*M_PI*K2);

        // If still infinite, resort to the approximation K1(sh) \approx x
        if (!isfinite(fh)) fh = exp(rcpp_calc_multi_locus_K(dist, sh) - sh*rcpp_calc_multi_locus_K1(dist, sh))/sqrt(2*M_PI*K2);
    }
    return fh;
}

// Function for calculating the saddlepoint approximation to the PDF
// [[Rcpp::export]]
NumericVector rcpp_calc_fhat(NumericVector x_vec, List dist,
                             double min_dist, double max_dist, double mean_dist,
                             bool logten = false) {

    NumericVector fh_vec (x_vec.length());

    // Create a copy of xvals in case we need to change the copy from log10 to log;
    // if we change the original vector it will change it by reference so it will
    // be changed outside the function. Need to use clone() to make it a deep
    // copy rather than just more references to the original.
    NumericVector x_vec_calc = clone(x_vec);
    if (logten) x_vec_calc = x_vec_calc/log10(exp(1));

    // Rcout << "max x : " << max(x_vec) << "\n";
    // Rcout << "max x calc : " << max(x_vec_calc) << "\n";

    double x;
    for(int idx=0; idx < x_vec_calc.length(); idx++) {
        x = x_vec_calc[idx];

        fh_vec[idx] = rcpp_calc_fhat_pt(x, dist, min_dist, max_dist, mean_dist);
    }

    return fh_vec;
}

// Note that this function does NOT have logten option because it should only
// be called from the Rcpp functions
// rcpp_calc_Fhat calls it AFTER converting x if logten = true
// rcpp_calc_qsearch_params calls it ONLY with logten = false
// So either way, it is only used with logten = false
// Also, do NOT export this function, because it doesn't have the logten conversion in it
double rcpp_calc_Fhat_pt(double x, List dist, double min_dist, double max_dist,
                         double mean_dist, bool use_x = false) {

    int nloci = dist.length();

    double sh, wh, uh;
    double K2, K3;
    double Fh;
    if (x >= max_dist) {
        Fh = 1;
    } else if (x <= min_dist) {
        Fh = 0;
    } else {
        // Within the vicinity of mean_dist, Fhat uses a different formula.
        // If this second formula is used only for x==mean_dist then Fhat is
        // smooth, but in fact due to the difficulty of accurately
        // calculating s_hat, w_hat and u_hat accurately in the vicinity of
        // mean_dist, it is necessary to use the Fhat(mean_dist) formula
        // also for x near mean_dist, to avoid major discontinuities in
        // Fhat. Although Fhat(mean_dist) is not precisely accurate for
        // x !=  mu, it is far more accurate than the other formula for Fhat
        // when applied to x close to mean_dist. How close to mean_dist you
        // can get before the numerical discontinuity arises depends on the
        // number of loci if (x != mean_dist)
        if (abs(x - mean_dist) > nloci*1E-5) {
            sh = rcpp_calc_shat(x, dist);

            if (isnan(sh)) {
                Fh = nan("0");
            } else {
                wh = rcpp_calc_what(x, sh, dist, use_x=use_x);
                uh = rcpp_calc_uhat(sh, dist);
                // Use the scalar form of pnorm and dnorm
                Fh = R::pnorm(wh, 0.0, 1.0, true, false) + R::dnorm(wh, 0.0, 1.0, false)*(1/wh - 1/uh);
            }
        } else {
            K3 = rcpp_calc_multi_locus_K3(dist,0);
            K2 = rcpp_calc_multi_locus_K2(dist,0);
            Fh = 0.5 + K3/(6*sqrt(2*M_PI)*pow(K2,1.5));
        }

        // Correct Fh if it is smaller than 0 or bigger than 1
        if (!isnan(Fh) && Fh < 0) {
            Fh = 0;
        } else if (!isnan(Fh) && Fh > 1) {
            Fh = 1;
        }
    }
    return Fh;
}

// Function for calculating the saddlepoint approximation to the CDF
// [[Rcpp::export]]
NumericVector rcpp_calc_Fhat(NumericVector x_vec, List dist,
                             double min_dist, double max_dist, double mean_dist,
                             bool logten = false, bool use_x = false)
{
    NumericVector Fh_vec (x_vec.length());

    // Create a copy of xvals in case we need to change the copy from log10 to log;
    // if we change the original vector it will change it by reference so it will
    // be changed outside the function.
    NumericVector x_vec_calc = clone(x_vec);
    if (logten) x_vec_calc = x_vec_calc/log10(exp(1));

    double x;

    for(int idx=0; idx < x_vec_calc.length(); idx++) {
        x = x_vec_calc[idx];
        Fh_vec[idx] = rcpp_calc_Fhat_pt(x, dist, min_dist, max_dist, mean_dist, use_x);
    }
    return Fh_vec;
}

// Function to calculate the upper limit of x values to use when searching for
// quantiles, determined by where the calculation of shat starts becoming less
// accurate, and the lower limit of x values, determined when K2 becomes NA
// [[Rcpp::export]]
List rcpp_calc_qsearch_params(List dist, double min_dist, double max_dist,
                              double sh_tol = 1e-4, double max_quantile = 0.995,
                              bool logten = false)
{
    double Fh_test;
    double mean_dist;
    bool found_max_shat_limit;
    double diff_max, diff_min, x_max, x_min, Fh_max, Fh_min;
    double sh_test, x_test;

    mean_dist = rcpp_calc_mu(dist);

    // Find the upper limit for x by taking diff away from the distribution max
    diff_max = 0.1;
    found_max_shat_limit = false;

    Fh_test = rcpp_calc_Fhat_pt(max_dist - diff_max, dist, min_dist, max_dist, mean_dist);

    // If already hitting NaNs for Fhat, try for bigger distance away from xmax.
    while(isnan(Fh_test) && diff_max < (max_dist - min_dist - 0.1)) {
        diff_max = diff_max + 0.1;
        Fh_test = rcpp_calc_Fhat_pt(max_dist-diff_max, dist, min_dist, max_dist, mean_dist);
    }

    // If not yet hitting NaNs for Fhat, try for smaller distance away from xmax.
    while(!found_max_shat_limit && diff_max > 1e-12 && !isnan(Fh_test))
    {
        // Rcout << "Trying for smaller distance from xmax\n";
        // Rcout << "Fh_test is " << Fh_test << "\n";

        sh_test = rcpp_calc_shat(max_dist-diff_max, dist);
        x_test = rcpp_calc_multi_locus_K1(dist, sh_test);
        if (abs(x_test - max_dist + diff_max) > sh_tol || Fh_test - max_quantile <= 0) found_max_shat_limit = true;
        else diff_max = diff_max/10;

        Fh_test = rcpp_calc_Fhat_pt(max_dist-diff_max, dist, min_dist, max_dist, mean_dist);
    }
    // go back to the previous diff value, either because for this one Fhat is
    // NA, or because we've found the point where x and K1(s(x)) diverge, so we
    // want to go back a step
    diff_max = diff_max*10;

    // # Find the lower limit by taking diffs away from the distribution min
    // diff_min = 0
    // Changed on 2017-02-07 because having diff_min = 0 is pointless since
    // Fhat(min_dist+diff_min) = Fhat(min_dist) which is automatically set to 0,
    // so we never look for higher diff_min.
    diff_min = 0.1;

    Fh_test = rcpp_calc_Fhat_pt(min_dist+diff_min, dist, min_dist, max_dist, mean_dist);

    // If already hitting NaNs for Fhat, try for bigger distance away from xmax.
    while((isnan(Fh_test) || Fh_test == 1) && diff_min < (max_dist - diff_max - min_dist - 0.01))
    {
        diff_min = diff_min + 0.01;
        Fh_test = rcpp_calc_Fhat_pt(min_dist+diff_min, dist, min_dist, max_dist, mean_dist);
    }

    if (diff_min > 0.1) {
        Rcout << "diff_min is " << diff_min << "\n";
    }

    x_max=max_dist-diff_max;
    x_min=min_dist+diff_min;

    // Need to calculate the upper and lower limits of Fhat for logten=FALSE,
    // because the xvalues tested have been in the space of log, not log10
    Fh_max = rcpp_calc_Fhat_pt(x_max, dist, min_dist, max_dist, mean_dist);
    Fh_min = rcpp_calc_Fhat_pt(x_min, dist, min_dist, max_dist, mean_dist);

    if (isnan(Fh_max)) Rcout << "Fh_max is nan\n";

    if (x_max > max_dist) Rcout << "x_max > max_dist\n";

    List L = List::create(Named("x_max") = x_max, Named("x_min") = x_min,
                          Named("Fh_max") = Fh_max, Named("Fh_min") = Fh_min,
                          Named("diff_max") = diff_max, Named("diff_min") = diff_min,
                          Named("dist") = dist, Named("mean_dist") = mean_dist,
                          Named("min_dist") = min_dist, Named("max_dist") = max_dist,
                          Named("logten") = logten,
                          Named("max_quantile") = max_quantile);
    return L;
}

struct Qhat_functor_deriv
{ // Functor also returning 1st derivative.
    Qhat_functor_deriv(double const& to_find_inverse_of, List dist,
                       double min_dist, double max_dist, double mean_dist) :
                       _y(to_find_inverse_of), _dist(dist),
                       _min_dist(min_dist), _max_dist(max_dist), _mean_dist(mean_dist)
    { // Constructor stores value to find root of, and distribution to use
        // for example: calling Qhat_functor_deriv(x) to use to get Fhat inverse at y.
    }
    std::pair<double, double> operator()(double const& x)
    {
        // Return both f(x) and f'(x).
        double fx = rcpp_calc_Fhat_pt(x, _dist, _min_dist, _max_dist, _mean_dist) - _y;   // Difference (estimate Fhat(x) - value).
        double dx = rcpp_calc_fhat_pt(x, _dist, _min_dist, _max_dist, _mean_dist);        // 1st derivative (fhat(x)).
        return std::make_pair(fx, dx);   // 'return' both fx and dx.
    }
private:
    double _y;         // Store value to find the Fhat inverse of
    List _dist;        // Store distribution to use
    double _min_dist;  // Store value of min_dist
    double _max_dist;  // Store value of max_dist
    double _mean_dist; // Store value of mean_dist
};


// [[Rcpp::export]]
double rcpp_calc_Qhat(double y, List Fhat_qsearch_params, List nutab,
                      bool logten = false)
{
    List dist = Fhat_qsearch_params["dist"];
    double mean_dist = Fhat_qsearch_params["mean_dist"];
    double min_dist = Fhat_qsearch_params["min_dist"];
    double max_dist = Fhat_qsearch_params["max_dist"];
    double Fh_min = Fhat_qsearch_params["Fh_min"];
    double Fh_max = Fhat_qsearch_params["Fh_max"];

    if (isnan(y) || isnan(Fh_max)) stop("Invalid value to find quantile for.");

    double Qh;

    if (y >= Fh_max) {
        Qh = max_dist;
        if (logten) Qh = Qh/log(10);
    } else if (y <= Fh_min) {// (y <= 0) was the old version, changed on 2017-02-07
        Qh = min_dist;
        if (logten) Qh = Qh/log(10);
    } else {
        if (isnan(Fh_max) || isnan(Fh_min))
        {
            Qh = nan("0");
        }
        else if (((Fh_max-y <= 0) & (Fh_min-y <= 0)) || ((Fh_max-y >= 0) & (Fh_min-y >= 0)))
        {
            Qh = nan("0");
            stop("Can't find quantile, CDF values same at both ends of search range.");
        }
        else
        {
            // return shat using 1st derivative of Fhat and Newton_Raphson.
            double guess = 0;                                        // Rough guess is zero.
            double min = Fhat_qsearch_params["x_min"];               // Minimum possible value is lower
            double max = Fhat_qsearch_params["x_max"];               // Maximum possible value is upper
            const int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type double.
            int get_digits = static_cast<int>(digits * 0.6);         // Accuracy doubles with each step, so stop when we have
            // just over half the digits correct.
            const std::uintmax_t maxit = 1000;
            std::uintmax_t it = maxit;

            try {
                double Qh = newton_raphson_iterate(
                    Qhat_functor_deriv(y, dist, min_dist, max_dist, mean_dist),
                    guess, min, max, get_digits, it);

                if (isnan(Qh)) stop("NaN quantile value.");

                if (Qh > max_dist) stop("quantile value beyond the upper limit of the distribution."); // Need to check this BEFORE converting via logten!

                if (logten) Qh = Qh/log(10);
                return Qh;
            } catch(std::exception &ex) {
                forward_exception_to_r(ex);
                Qh = nan("0");
            } catch(...) {
                ::Rf_error("c++ exception (unknown reason)");
            }
        }
    }

    return Qh;
}




/*** R
### NOTE: MANY of these tests require access to the low-level Rcpp functions,
### which will not be available unless you add the Rcpp::export tag above them!
# source("saddlepoint_spdf_master_from_geneplot_research.R")
# rcpp_calc_dist(matrix(c(1,2,0),nrow=1))
# calc.probs.func(c(1,2,0))
#
# library(testthat)
# nu1 = matrix(2)
# probs1 = rcpp_calc_probs(nu1)
# expect_equivalent(probs1,1)
# dist1 = rcpp_calc_dist(nu1)
# expect_equivalent(dist1$probs,1)
#
# nu2 = matrix(c(1.5,1.5), nrow=1)
# probs2 = rcpp_calc_probs(nu2)
# expect_equivalent(sort(probs2,decreasing=FALSE),c(1.5*2.5/3/4,1.5*2.5/3/4,2*1.5*1.5/3/4))
# # expect_equivalent(log(1.5*2.5/3/4),rcpp_calc_min_beta(nu2))
# # expect_equivalent(log(2*1.5*1.5/3/4),rcpp_calc_max_beta(nu2))
#
# nu3 = matrix(c(9.5,1.5),nrow=1)
# probs3 = rcpp_calc_probs(nu3)
# expect_equivalent(sort(probs3,decreasing=FALSE),c(1.5*2.5/11/12,2*1.5*9.5/11/12,9.5*10.5/11/12))
# # expect_equivalent(log(1.5*2.5/11/12),rcpp_calc_min_beta(nu3))
# # expect_equivalent(log(9.5*10.5/11/12),rcpp_calc_max_beta(nu3))
#
# nu4 = matrix(c(10+1/3,1/3,1/3),nrow=1)
# probs4 = rcpp_calc_probs(nu4)
# expect_equivalent(sort(probs4,decreasing=FALSE),c(2*1/3*1/3/11/12, 1/3*(1/3+1)/11/12, 1/3*(1/3+1)/11/12, 2*1/3*(10+1/3)/11/12, 2*1/3*(10+1/3)/11/12, (10+1/3)*(11+1/3)/11/12))
# # expect_equivalent(log(2*1/3*1/3/11/12),rcpp_calc_min_beta(nu4))
# # expect_equivalent(log((10+1/3)*(11+1/3)/11/12),rcpp_calc_max_beta(nu4))
#
# nu23 = rbind(nu2,nu3)
# dist23 = rcpp_calc_dist(nu23)
# dist23r = calc.probs.func(nu23)
# expect_equivalent(dist23$probs, dist23r['probs',])
# expect_equivalent(dist23$values, dist23r['values',])
#
# dist23 = rcpp_calc_dist(nu23, TRUE)
# dist23r = calc.probs.func(nu23, TRUE)
# expect_equivalent(dist23$probs, dist23r['probs',])
# expect_equivalent(dist23$values, dist23r['values',])
#
# dist23 = rcpp_calc_dist(nu23, differenceLGPs=TRUE)
# dist23r = calc.probs.func(nu23, differenceLGPs=TRUE)
# expect_equivalent(dist23$probs, dist23r['probs',])
# expect_equivalent(dist23$values, dist23r['values',])
#
# dist23 = rcpp_calc_dist(nu23, TRUE,TRUE)
# dist23r = calc.probs.func(nu23, TRUE,TRUE)
# expect_equivalent(dist23$probs, dist23r['probs',])
# expect_equivalent(dist23$values, dist23r['values',])
#
# dist_info = rcpp_calc_multi_locus_dist(list(nu2,nu3,nu4))
# expect_equivalent(sum(c(0,log(1.5*2.5/3/4),log(1.5*2.5/11/12),log(2*1/3*1/3/11/12))), dist_info$min)
# expect_equivalent(sum(c(0,log(2*1.5*1.5/3/4),log(9.5*10.5/11/12),log((10+1/3)*(11+1/3)/11/12))), dist_info$max)
#
# expect_equivalent(rcpp_calc_multi_locus_K(dist_info$dist, s=0), multi.K(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K(dist_info$dist, s=1), multi.K(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K1(dist_info$dist, s=0), multi.K1(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K1(dist_info$dist, s=1), multi.K1(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K2(dist_info$dist, s=0), multi.K2(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K2(dist_info$dist, s=1), multi.K2(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K3(dist_info$dist, s=0), multi.K3(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K3(dist_info$dist, s=1), multi.K3(dist_info$dist, s=1))
#
# ## Leave-one-out
# dist_info = rcpp_calc_multi_locus_dist(list(nu2,nu3,nu4), leave_one_out = TRUE)
# # dist_info = calc.multi.locus.probs.func(list(as.vector(nu2),as.vector(nu3),as.vector(nu4)), leave_one_out = TRUE)
# expect_equivalent(rcpp_calc_multi_locus_K(dist_info$dist, s=0), multi.K(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K(dist_info$dist, s=1), multi.K(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K1(dist_info$dist, s=0), multi.K1(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K1(dist_info$dist, s=1), multi.K1(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K2(dist_info$dist, s=0), multi.K2(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K2(dist_info$dist, s=1), multi.K2(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_multi_locus_K3(dist_info$dist, s=0), multi.K3(dist_info$dist, s=0))
# expect_equivalent(rcpp_calc_multi_locus_K3(dist_info$dist, s=1), multi.K3(dist_info$dist, s=1))
#
# expect_equivalent(rcpp_calc_mu(dist_info$dist), mu(dist_info$dist))
# expect_equivalent(rcpp_calc_sigma(dist_info$dist), sigma(dist_info$dist))
#
# expect_equivalent(rcpp_calc_what(1, 1, dist_info$dist, TRUE), what(1, 1, dist_info$dist, TRUE))
# expect_equivalent(rcpp_calc_what(1, 1, dist_info$dist, FALSE), what(1, 1, dist_info$dist, FALSE))
# expect_equivalent(rcpp_calc_uhat(1, dist_info$dist), uhat(1, dist_info$dist))
#
# shat(-2, dist_info$dist)
# rcpp_calc_shat(-2, dist_info$dist)
#
# Convert the Rcpp output into the appropriate format for R version of fhat
# multi_dist = list(dist=lapply(dist_info$dist, function(loc) {
#     print(loc)
#     mat = rbind(loc$values, loc$probs)
#     rownames(mat) = c("values","probs")
#     mat
# }), min=dist_info$min, max=dist_info$max)
#
# mean_dist = rcpp_calc_mu(dist_info$dist)
#
# x_vec = -13:1
#
# fhat(x_vec, multi_dist, logten=TRUE)
# rcpp_calc_fhat(x_vec, dist_info$dist, dist_info$min, dist_info$max, mean_dist, logten=TRUE)
# fhat(x_vec, multi_dist, logten=FALSE)
# rcpp_calc_fhat(x_vec, dist_info$dist, dist_info$min, dist_info$max, mean_dist, logten=FALSE)
#
# Fh = Fhat(x_vec, multi_dist, mean_dist, logten=TRUE)
# Fh_rcpp = rcpp_calc_Fhat(x_vec, dist_info$dist, dist_info$min, dist_info$max, mean_dist, logten=TRUE)
# expect_equivalent(Fh_rcpp, Fh)
# Fh = Fhat(x_vec, multi_dist, mean_dist, logten=FALSE)
# Fh_rcpp = rcpp_calc_Fhat(x_vec, dist_info$dist, dist_info$min, dist_info$max, mean_dist, logten=FALSE)
# expect_equivalent(Fh_rcpp, Fh)
# Fh = Fhat(x_vec, multi_dist, mean_dist, logten=FALSE, use.x.in.wh = TRUE)
# Fh_rcpp = rcpp_calc_Fhat(x_vec, dist_info$dist, dist_info$min, dist_info$max, mean_dist, use_x=TRUE)
# expect_equivalent(Fh_rcpp, Fh)
#
# ## two_pops = FALSE
# qparams = get.qsearch.params(multi_dist, dist_info$dist, s.tol=1e-4, logten=FALSE, max.quantile=0.995, twoPops=FALSE)
# qparams_rcpp = rcpp_calc_qsearch_params(dist_info$dist, dist_info$min, dist_info$max, sh_tol=1e-4,
#                                         max_quantile=0.995, logten=FALSE)
# expect_equivalent(qparams_rcpp$x_min, qparams$x.min)
# expect_equivalent(qparams_rcpp$x_max, qparams$x.max)
# expect_equivalent(qparams_rcpp$Fh_min, qparams$Fh.min)
# expect_equivalent(qparams_rcpp$Fh_max, qparams$Fh.max)
# expect_equivalent(qparams_rcpp$diff_min, qparams$diff.min)
# expect_equivalent(qparams_rcpp$diff_max, qparams$diff.max)
#
# ## two_pops = TRUE
# qparams = get.qsearch.params(multi_dist, dist_info$dist, s.tol=1e-4, logten=FALSE, max.quantile=0.995, twoPops=TRUE)
# qparams_rcpp = rcpp_calc_qsearch_params(dist_info$dist, dist_info$min, dist_info$max, sh_tol=1e-4,
#                                         max_quantile=0.995, logten=FALSE)
# expect_equivalent(qparams_rcpp$x_min, qparams$x.min)
# expect_equivalent(qparams_rcpp$x_max, qparams$x.max)
# expect_equivalent(qparams_rcpp$Fh_min, qparams$Fh.min)
# expect_equivalent(qparams_rcpp$Fh_max, qparams$Fh.max)
# expect_equivalent(qparams_rcpp$diff_min, qparams$diff.min)
# expect_equivalent(qparams_rcpp$diff_max, qparams$diff.max)
*/

