#include <Ecm.hpp>
#include <cmath>

Lenstra_ECM::Lenstra_ECM() {
    complete = false;
    primes = {2, 3, 5};
    n = 4;
    result = 2;
    B = 5;
}

void Lenstra_ECM::set_uncompletness() {
    complete = false;
}

mpz_class Lenstra_ECM::get_result() {
    return result;
}

bool Lenstra_ECM::is_not_satisfied(gmp_randclass& rr, mpz_class& m) {
    if (is_prime(rr, m, 100)) {
        result = -1;
        return true;
    }
    if (m % 2 == 0) {
        result = 2;
        return true;
    }
    if (m % 3 == 0) {
        result = 3;
        return true;
    }
    return false;
}

std::variant<bool, Point> Lenstra_ECM::multiply_point(
    mpz_class& n, 
    mpz_class& a, 
    long long num, 
    std::variant<bool, Point>& point
) {
    long long steps = static_cast<long long>(std::log2(static_cast<double>(num)));
    std::variant<bool, Point> result = Point(0, 0, true);
    std::variant<bool, Point> current_pow = point;
    for (long long i = 0; i < steps + 1; i++) {
        if (num % 2 == 1) {
            result = add_points(n, a, result, current_pow);
            if (std::holds_alternative<bool>(result)) {
                return false;
            }
        }
        num /= 2;
        current_pow = add_points(n, a, current_pow, current_pow);
        if (std::holds_alternative<bool>(current_pow)) {
                return false;
        }
    }
    return result;
}

void normalize(mpz_class& x, mpz_class& n) {
    if (x < 0) {
        x += n;
    }
}

std::variant<bool, Point> Lenstra_ECM::add_points(
    mpz_class& n, 
    mpz_class& a, 
    std::variant<bool, Point>& point1, 
    std::variant<bool, Point>& point2
) {
    Point p1 = std::get<Point>(point1);
    Point p2 = std::get<Point>(point2);
    if (p1.get_is_o()) {
        return p2;
    } else if (p2.get_is_o()) {
        return p1;
    }
    mpz_class alpha;
    if (p1.get_x() == p2.get_x()) {
        if (p1.get_y() != p2.get_y()) {
            return Point(0, 0, true);
        } else {
            if (p1.get_y() == 0) {
                return Point(0, 0, true);
            } else {
                mpz_class input = 2 * p1.get_y();
                std::tuple<mpz_class, mpz_class, mpz_class> gcds = extended_euclid(input, n);
                mpz_class gcd = std::get<0>(gcds);
                if (gcd != 1 && gcd != n) {
                    result = gcd;
                    return false;
                } 
                auto inverse = std::get<1>(gcds);
                alpha = ((3 * p1.get_x() * p1.get_x() + a) * inverse) % n;
            }
        }
    } else {
        mpz_class input = p1.get_x() - p2.get_x();
        normalize(input, n);
        std::tuple<mpz_class, mpz_class, mpz_class> gcds = extended_euclid(input, n);
        mpz_class gcd = std::get<0>(gcds);
        if (gcd != 1 && gcd != n) {
            result = gcd;
            return false;
        }
        mpz_class inverse = std::get<1>(gcds);
        alpha = ((p1.get_y() - p2.get_y()) * inverse) % n;
    }
    normalize(alpha, n);
    mpz_class beta = (p1.get_y() - alpha * p1.get_x()) % n;
    normalize(beta, n);
    mpz_class x = (alpha * alpha - p1.get_x() - p2.get_x()) % n;
    normalize(x, n);
    mpz_class y = (-(alpha * x + beta)) % n;
    normalize(y, n);
    return Point(x, y);
}

bool Lenstra_ECM::try_ecm(
    mpz_class& n, 
    mpz_class& a, 
    int B, 
    mpz_class& C, 
    Point& point
) {
    std::variant<bool, Point> Q = point;
    for (int prime : primes) {
        if (complete) {
            return true;
        }
        if (prime > B) {
            break;
        }
        long long steps = log(C, prime);
        for (int j = 0; j < steps; j++) {
            Q = multiply_point(n, a, prime, Q);
            if (std::holds_alternative<bool>(Q)) {
                return true;
            }
        }
    }
    return false;
}

void Lenstra_ECM::factor(mpz_class& n, int B, mpz_class& C) {
    gmp_randclass rr(gmp_randinit_mt);
    std::mt19937 mt{std::random_device{}()};
    rr.seed(mt());
    if (is_not_satisfied(rr, n)) {
        return;
    }
    mtx.lock();
    if (B > this->B) {
        this->B = B;
        std::vector<int> prims = {};
        eratosthenes(prims, B);
        primes = prims;
    }
    mtx.unlock();
    while (true) {
        if (complete) {
            break;
        }
        mpz_class a = rr.get_z_range(n);
        mpz_class u = rr.get_z_range(n);
        mpz_class v = rr.get_z_range(n);
        mpz_class b = (v * v - u * u * u - a * u) % n;
        normalize(b, n);
        mpz_class D = (4 * a * a * a + 27 * b * b) % n;
        std::tuple<mpz_class, mpz_class, mpz_class> gcds = extended_euclid(D, n);
        mpz_class gcd = std::get<0>(gcds);
        if (gcd > 1 && gcd < n) {
            result = gcd;
            complete = true;
            break;
        } else if (gcd == n) {
            continue;
        }
        Point point = Point(u, v);
        if (try_ecm(n, a, B, C, point)) {
            complete = true;
            break;
        }
    }
}