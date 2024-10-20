#include <gmpxx.h>
#include <random>
#include <tuple>

std::tuple<mpz_class, mpz_class, mpz_class> extended_euclid(mpz_class& a, mpz_class& b);

void eratosthenes(std::vector<int>& primes, int n);

bool is_prime(gmp_randclass& rr, mpz_class& number, int steps);

long long log(mpz_class& C, long long num);