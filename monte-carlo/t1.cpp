#include <iostream>
#include <random>


double monte_carlo_call_price(const int &num_sims, const double &S, const double &K, const double &r, const double &v, const double &T)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> nd{0, 1};

    double S_adjust = S * exp(T * (r - 0.5 * v * v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i = 0; i < num_sims; i++)
    {
        double gauss_bm = nd(gen);
        S_cur = S_adjust * exp(sqrt(v * v * T) * gauss_bm);
        payoff_sum += std::max(S_cur - K, 0.0);
    }

    return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);
}


double monte_carlo_put_price(const int &num_sims, const double &S, const double &K, const double &r, const double &v, const double &T)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> nd{0, 1};

    double S_adjust = S * exp(T * (r - 0.5 * v * v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i = 0; i < num_sims; i++)
    {
        double gauss_bm = nd(gen);
        S_cur = S_adjust * exp(sqrt(v * v * T) * gauss_bm);
        payoff_sum += std::max(K - S_cur, 0.0);
    }

    return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);
}


int main()
{
    int num_sims = 10000;    // Number of simulated asset paths
    double S = 100.0;        // Option price
    double K = 100.0;        // Strike price
    double r = 0.05;         // Risk-free rate (5%)
    double v = 0.2;          // Volatility of the underlying (20%)
    double T = 1.0;          // One year until expiry

    double call = monte_carlo_call_price(num_sims, S, K, r, v, T);
    double put = monte_carlo_put_price(num_sims, S, K, r, v, T);

    std::cout << call << std::endl;
    std::cout << put << std::endl;
    
    return 0;
}
