#include <Point.cpp>
#include <variant>
#include <atomic>
#include <mutex>

class Lenstra_ECM {
public:
    Lenstra_ECM();
    void factor(mpz_class& n, int B, mpz_class& C);
    mpz_class get_result();
    void set_uncompletness();
private:
    std::vector<int> primes;
    mpz_class n;
    mpz_class result;
    int B;
    std::mutex mtx;
    std::atomic<bool> complete;

    bool is_not_satisfied(gmp_randclass& rr, mpz_class& m);
    std::variant<bool, Point> add_points(
        mpz_class& n, 
        mpz_class& a, 
        std::variant<bool, Point>& point1, 
        std::variant<bool, Point>& point2
    );
    std::variant<bool, Point> multiply_point(
        mpz_class& n, 
        mpz_class& a, 
        long long num, 
        std::variant<bool, Point>& point
    );
    bool try_ecm(
        mpz_class& n, 
        mpz_class& a, 
        int B, 
        mpz_class& C, 
        Point& point
    );
};