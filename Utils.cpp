#include <Utils.hpp>

std::tuple<mpz_class, mpz_class, mpz_class> extended_euclid(mpz_class& a, mpz_class& b) {
    mpz_class x = 1;
    mpz_class y = 0;
    mpz_class g = mpz_class(a);
    mpz_class r = 0;
    mpz_class s = 1;
    mpz_class t = mpz_class(b);
    while (t > 0) {
        mpz_class q = g / t;
        mpz_class u = x - q * r;
        mpz_class v = y - q * s;
        mpz_class w = g - q * t;
        x = mpz_class(r);
        y = mpz_class(s);
        g = mpz_class(t);
        r = mpz_class(u);
        s = mpz_class(v);
        t = mpz_class(w);
    }
    return {g, x, y};
}

void eratosthenes(std::vector<int>& primes, int n) {
    std::vector<int> m(n + 1, 0);
    for (int i = 2; i < n; i++) {
        if (m[i] == 0) {
            m[i] = i;
            primes.push_back(i);
        }
        for (std::vector<int>::size_type j = 0; j < primes.size() && primes[j] <= m[i] && 
            i * primes[j] <= n; j++) {
            m[i * primes[j]] = primes[j];
        }
    }
}

mpz_class fast_power_mod(mpz_class& a, mpz_class& u, mpz_class& mod) {
    if (u == 0) {
        return 1;
    }
    mpz_class half = u / 2;
    mpz_class sqrt = fast_power_mod(a, half, mod);
    if (u % 2 == 0) {
        return (sqrt * sqrt) % mod;
    } 
    return (a * sqrt * sqrt) % mod;
}

bool is_prime(gmp_randclass& rr, mpz_class& n, int steps) {
    while (steps > 0) {
        mpz_class a = rr.get_z_range(n - 1);
        if (a <= 2) {
            continue;
        }
        steps -= 1;
        int t = 0;
        mpz_class u = n - 1;
        while (u % 2 == 0) {
            u /= 2;
            t += 1;
        }
        std::vector<mpz_class> x = {fast_power_mod(a, u, n), 0};
        for (int i = 1; i < t + 1; i++) {
            int cur = i % 2;
            int prev = (i + 1) % 2;
            x[cur] = (x[prev] * x[prev]) % n;
            if (x[cur] % n == 1 && x[prev] % n != 1 && x[prev] % n != n - 1) {
                return false;
            }
        }
        if (x[t % 2] != 1) {
            return false;
        }
    }
    return true;
}

long long log(mpz_class& C, long long num) {
    mpf_class cf = mpf_class(C);
    mpz_class one = mpz_class(1);
    mpz_class to_log = C + 2 * mpz_class(floor(sqrt(cf))) + 1;
    long long deg = 0;
    while (one < to_log) {
        deg += 1;
        one *= mpz_class((int) num);
    }
    return deg;
}