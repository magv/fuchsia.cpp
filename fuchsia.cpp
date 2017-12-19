#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <assert.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <tuple>

using namespace GiNaC;
using namespace std;

/* LOGGING
 * ============================================================
 *
 * As with everything in C++, you can "optimize" this piece of
 * code to e.g. use compile-only format string parsing, minimize
 * number of created functions, and so on, but at the expense
 * of kilolines of code, and your own sanity lost in the fight
 * versus byzantine template rules. Please don't.
 */

bool log_verbose = false;

static auto _log_starttime = chrono::steady_clock::now();
static auto _log_lasttime = chrono::steady_clock::now();
static int _log_depth = 0;

const char*
log_adv(const char *fmt)
{
    for (int i = 0; ; i++) {
        if (fmt[i] == '{') {
            cout.write(fmt, i);
            return fmt + i + 2;
        }
        if (fmt[i] == 0) {
            cout.write(fmt, i);
            return fmt + i;
        }
    }
}

/* This function is used to print objects into the log. Override it
 * for the data types you care about to modify their appearance.
 */
template<typename T> static inline void
log_format(ostream &o, const T &value)
{
    o << value;
}

void
log_print_start(const char *lvl)
{
    auto t = chrono::steady_clock::now();
    auto dt = chrono::duration_cast<chrono::duration<double>>(t - _log_starttime).count();
    cout << "\033[32m[" << lvl << " " << std::fixed << std::setprecision(4) << dt << "s +";
    cout << chrono::duration_cast<chrono::duration<double>>(t - _log_lasttime).count() << "s";
    for (int i = 0; i < _log_depth; i++) {
        cout << " *";
    }
    cout << "]\033[0m ";
    _log_lasttime = t;
}

template<typename T> const char *
log_print_one(const char *fmt, const T &value)
{
    fmt = log_adv(fmt);
    log_format(cout, value);
    return fmt;
}

void
log_print_end(const char *fmt)
{
    cout << fmt << endl;
}

struct _sequencehack {
    template<typename ...Args>
    _sequencehack(Args &&...) {}
};

template<typename ...Args> static void
log_fmt(const char *lvl, const char *fmt, const Args &...args)
{
    log_print_start(lvl);
    (void) _sequencehack {
        (fmt = log_print_one(fmt, args), 0)
        ...
    };
    log_print_end(fmt);
}

/* Log an debug message. These can be suppressed by setting
 * log_verbose to false.
 */
template<typename... Args> static inline void
logd(const char *fmt, const Args &...args)
{
    if (log_verbose) {
        log_fmt("dbg", fmt, args...);
    }
}

/* Log an information message.
 */
template<typename... Args> static inline void
logi(const char *fmt, const Args &...args)
{
    log_fmt("inf", fmt, args...);
}

/* Log an error message.
 */
template<typename... Args> static inline void
loge(const char *fmt, const Args &...args)
{
    log_fmt("err", fmt, args...);
}

template<typename F>
struct _scopeexithack {
    _scopeexithack(F f) : f(f) {}
    ~_scopeexithack() { f(); }
    F f;
};

/* Place this macro as the start of a function, and you'll get
 * log entries every time this function is entered and exited.
 *
 * As a general rule, if you're adding logging statements to a
 * function, add LOGME to its start as well.
 */
#define LOGME \
    const auto &__log_func = __func__; \
    logd("> {}()", __log_func); \
    const auto __log_t0 = _log_lasttime; \
    _log_depth++; \
    auto __log_f = [&]{ \
        _log_depth--; \
        if (log_verbose) { \
            auto t = chrono::steady_clock::now(); \
            auto dt = chrono::duration_cast<chrono::duration<double>>(t - __log_t0).count(); \
            logd("< {}(+{}s)",__log_func, dt); \
        } \
    }; \
    auto __log_s = _scopeexithack<decltype(__log_f)>(__log_f);

/* The error base
 * ============================================================
 */

class fuchsia_error : public exception {
    const char *message;
    public:
    fuchsia_error(const char *message) : message(message) {};
    virtual const char* what() const throw() { return message; }
};

/* Miscellaneous general utilities
 * ============================================================
 */

matrix
ex_to_matrix(const ex &m)
{
    return ex_to<matrix>(m.evalm());
}

matrix
identity_matrix(unsigned n)
{
    matrix id(n, n);
    for (unsigned i = 0; i < n; i++) {
        id(i, i) = 1;
    }
    return id;
}

matrix
matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs)
{
    return m.transpose().solve(vars.transpose(), rhs.transpose()).transpose();
}

/* Infinity object.
 */
struct infinity_s { infinity_s() {} };
typedef structure<infinity_s> infinity_t;

const ex infinity = infinity_t();

template<> void
infinity_t::print(const print_context & c, unsigned level) const
{
    if (is_a<print_tree>(c))
        inherited::print(c, level);
    c.s << "Infinity";
}

/* Load matrix from a stream in Mathematica format.
 */
pair<matrix, symtab>
load_matrix(istream &i, const symtab &table)
{
    parser reader(table);
    ex x = reader(i);
    // TODO: check that the input is indeed a matrix of rational
    // expressions, with no imaginary or irrational numbers; signal
    // errors if this is not the case.
    matrix m = ex_to<matrix>(lst_to_matrix(ex_to<lst>(x)));
    return make_pair(m, reader.get_syms());
}

/* Load matrix from a file in Mathematica format.
 */
pair<matrix, symtab>
load_matrix(const char *filename, const symtab &table)
{
    ifstream i(filename);
    return load_matrix(i, table);
}

/* Write matrix to a stream object in Mathematica format.
 */
void
save_matrix(ostream &f, const matrix &m)
{
    f << "{";
    for (unsigned i = 0; i < m.rows(); i++) {
        if (i != 0) f << ",\n";
        f << "{";
        for (unsigned j = 0; j < m.cols(); j++) {
            if (j != 0) f << ",";
            f << m(i, j);
        }
        f << "}";
    }
    f << "}";
}

/* Save matrix to a file in Mathematica format.
 */
void
save_matrix(const char *filename, const matrix &m)
{
    ofstream f(filename);
    save_matrix(f, m);
}

/* Divide one polynomial in x by another, return the quotient
 * and the remainder.
 */
pair<ex, ex>
poly_divmod(const ex &a, const ex &b, const symbol &x)
{
    ex bb = b.expand();
    int deg_b = bb.degree(x);
    if (deg_b == 0)
        return make_pair(a/b, ex(0));
    ex aa = a.expand();
    int deg_a = aa.degree(x);
    if (deg_a < deg_b)
        return make_pair(0, aa);
    vector<ex> va(deg_a + 1);
    vector<ex> vb(deg_b + 1);
    for (int i = 0; i <= deg_a; i++) va[i] = aa.coeff(x, i);
    for (int i = 0; i <= deg_b; i++) vb[i] = bb.coeff(x, i);
    ex bn = normal(vb[deg_b]);
    ex q = 0;
    for (int i = deg_a; i >= deg_b; i--) {
        ex k = normal(va[i]/bn);
        q += k*pow(x, i - deg_b);
        // This zero assignment is only here to clean up the memory earlier.
        va[i] = 0;
        for (int j = 0; j < deg_b; j++)
            va[j+i-deg_b] -= vb[j]*k;
    }
    ex r = 0;
    for (int i = 0; i < deg_b; i++) {
        r += normal(va[i])*pow(x, i);
    }
    return make_pair(q, r);
}

/* Decompose polynomial over rationals e(x) into {c_i(x)},
 * such that e(x) = \sum_i c_i(x) p(x)^i.
 */
lst
poly_decompose(const ex &e, const ex &p, const symbol &x)
{
    lst c;
    ex q = e;
    while (!q.is_zero()) {
        auto qr = poly_divmod(q, p, x);
        c.append(qr.second);
        // XXX: is normal() call needed here?
        q = qr.first;
    }
    return c;
}

/* Canonical form for rational expressions.
 */
ex
ratcan(const ex &e)
{
    ex nd = e.normal().numer_denom();
    return nd.op(0).expand()/nd.op(1).expand();
}

/* Canonical form for matrices of rational expressions.
 */
matrix
ratcan(const matrix &m)
{
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.nops(); i++) {
        r.let_op(i) = ratcan(m.op(i));
    }
    return r;
}

/* Normal form for matrices of rational expressions.
 *
 * This is similar to 'ratcan', but doesn't expand numerator and
 * denominator needlessly.
 */
matrix
normal(const matrix &m)
{
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.nops(); i++) {
        r.let_op(i) = normal(m.op(i));
    }
    return r;
}

/* Iterate through factors of e, call yield(f, k) for each
 * factor of the form f^k.
 *
 * Note that this function doesn't factor e itself, only iterates
 * through the factors already explicitly present.
 */
template <typename F> void
factor_iter(const ex &e, F yield)
{
    if (is_a<mul>(e)) {
        for (const auto &f : e) {
            if (is_a<power>(f)) {
                yield(f.op(0), ex_to<numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (is_a<power>(e)) {
            yield(e.op(0), ex_to<numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

/* Factor an expression and iterate through its factors the same
 * way factor_iter does it.
 *
 * This would be normally the same as:
 *     factor_iter(factor(...), ...)
 *
 * ... but GiNaC's factor(...) is known to expand its arguments
 * (leading to slowness and sometimes hangs), while this routine
 * doesn't. Consider this to be a bug fix.
 */
template<typename F> void
factor_and_iter(const ex &e, F yield)
{
    factor_iter(e,
        [&](const ex &f1, int k1) {
            factor_iter(factor(f1),
                [&](const ex &f2, int k2) {
                    yield(f2, k1*k2);
                });
        });
}

ex
factor_fixed(const ex &e)
{
    ex res = 1;
    factor_and_iter(e,
        [&](const ex &f, int k) {
            res *= pow(f, k);
        });
    return res;
}

/* Iterate through terms of e, call yield(t) for each one.
 */
template <typename F> void
term_iter(const ex &e, F yield)
{
    if (is_a<add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

/* Iterate through terms of a partial fraction decomposition
 * of a rational expression in x. Call yield(p, q, k) for
 * each p*q^k term.
 */
template<typename F> void
partial_fraction_iter(const ex &e, const symbol &x, F yield)
{
    ex nd = e.normal().numer_denom();
    ex numer = nd.op(0);
    ex denom = nd.op(1);
    auto qr = poly_divmod(numer, denom, x);
    ex q = qr.first;
    ex r = qr.second;
    // At this point e == q + r/denom.
    struct Term { ex f; int n; };
    vector<Term> factors;
    ex denomf = 1;
    exmap sym2ex;
    factor_iter(factor_fixed(denom),
        [&](const ex &f, int i) {
            // We'll replace coefficient expressions with newly
            // generated symbols, and will subs() them back in
            // near the end. This makes it possible to work with
            // larger expressions.
            ex frenamed = 0;
            for (int i = 0; i <= f.degree(x); i++) {
                ex c = f.coeff(x, i);
                if (!is_a<numeric>(c) && !is_a<symbol>(c)) {
                    symbol t;
                    sym2ex[t] = c;
                    frenamed += t*pow(x, i);
                } else {
                    frenamed += c*pow(x, i);
                }
            }
            denomf *= pow(frenamed, i);
            factors.push_back(Term{frenamed, i});
        });
    lst clist;
    ex eq = -r;
    // For each factor...
    for (const auto &f : factors) {
        int deg = f.f.degree(x);
        // For each power of factor...
        for (int k = 1; k <= f.n; k++) {
            ex pfnum = 0;
            // For each power of the parfrac term numerator...
            for (int l = 0; l < deg; l++) {
                symbol c;
                clist.append(c);
                pfnum += c*pow(x, l);
            }
            eq += expand(pfnum*(denomf/pow(f.f, k)));
        }
    }
    lst eqlist;
    int deg_lo = eq.ldegree(x);
    int deg_hi = eq.degree(x);
    for (int k = deg_lo; k <= deg_hi; k++) {
        eqlist.append(eq.coeff(x, k) == 0);
    }
    ex sol = lsolve(eqlist, clist, solve_algo::gauss);
    q = q.expand().collect(x);
    for (int k = q.ldegree(x); k <= q.degree(x); k++) {
        ex c = q.coeff(x, k);
        if (!c.is_zero()) yield(c, x, k);
    }
    // For each factor...
    unsigned idx = 0;
    for (const auto &f : factors) {
        int deg = f.f.degree(x);
        // For each power of factor...
        for (int k = 1; k <= f.n; k++) {
            ex pfnum = 0;
            // For each power of the parfrac term numerator...
            for (int l = 0; l < deg; l++) {
                assert(idx < sol.nops());
                pfnum += (sol.op(idx).rhs().subs(sym2ex))*pow(x, l);
                idx++;
            }
            if (!pfnum.is_zero())
                yield(pfnum, f.f.subs(sym2ex), -k);
        }
    }
}

/* Transform a rational expression into partial fraction form in x.
 */
ex
partial_fraction(const ex &e, const symbol &x)
{
    ex pf = 0;
    partial_fraction_iter(e, x, [&](const auto &p, const auto &q, int n) {
        pf += p*pow(q, n);
    });
    return pf;
}

/* Matrix in a partial fraction form:
 *     M(x) = \sum_i R_i (x-p_i)^{k_i},
 *     where p_i = 0, if k_i >= 0.
 *
 * Note that this class expects its p_i arguments to be in some
 * sort of canonical form. Honor this expectation. Or else.
 */
struct pfmatrix {
    typedef pair<ex, int> key;
    struct key_is_less {
        bool operator ()(const key &k1, const key &k2) const;
    };

    map<key, matrix, key_is_less> residues;
    unsigned nrows;
    unsigned ncols;
    symbol x;

    pfmatrix(unsigned nrows, unsigned ncols, const symbol &x);
    pfmatrix(const matrix &m, const symbol &x);
    matrix &operator ()(const ex &p, int k);
    matrix to_matrix();
    void normalize();
    // M += C*(x-pi)^ki
    void add(const matrix &C, const ex &p1, int k1);
    // M += C*(x-pi)^ki/(x-p2)
    void add_div(const matrix &C, const ex &p1, int k1, const ex &p2);
    // M += C*(x-pi)^ki*(x-p2)
    void add_mul(const matrix &C, const ex &p1, int k1, const ex &p2);
    // M += C*(x-p1)^k1*(x-p2)^k2
    void add_pow(const matrix &C, const ex &p1, int k1, const ex &p2, int k2);
    // M' = T^-1 M T
    pfmatrix with_constant_t(const matrix &T) const;
    // M' = L M R
    pfmatrix with_constant_t(const matrix &L, const matrix &R) const;
    pfmatrix with_balance_t(const matrix &P, const ex &x1, const ex &x2) const;
};

bool
pfmatrix::key_is_less::operator ()(const key &k1, const key &k2) const
{
    return (k1.second < k2.second) || ex_is_less()(k1.first, k2.first);
}

pfmatrix::pfmatrix(unsigned nrows, unsigned ncols, const symbol &x)
    : nrows(nrows), ncols(ncols), x(x)
{ }

pfmatrix::pfmatrix(const matrix &m, const symbol &x)
    : nrows(m.rows()), ncols(m.cols()), x(x)
{
    for (unsigned i = 0; i < nrows; i++) {
        for (unsigned j = 0; j < ncols; j++) {
            // Caveat emptor! This algorithm may place
            // (non-canonical) zero cells into Ci values,
            // and possible whole zero Ci matrices (canonical
            // and/or not) into this->residues.
            term_iter(m(i, j),
                [&](const ex &t) {
                    partial_fraction_iter(t, x,
                        [&](const auto &p, const auto &q, int k) {
                            int deg = q.degree(x);
                            if (k >= 0) {
                                // q == x
                                (*this)(0, k)(i, j) += p;
                            }
                            else if (deg == 1) {
                                // q == c0 + c1*x
                                ex c0 = q.coeff(x, 0);
                                ex c1 = q.coeff(x, 1);
                                (*this)(ratcan(-c0/c1), k)(i, j) += p*pow(c1, k);
                            }
                            else {
                                throw runtime_error("pfmatrix(): can't solve equations of 2nd degree or higher");
                            }
                        });
                });
        }
    }
}

matrix &
pfmatrix::operator ()(const ex &p, int k)
{
    const auto key = make_pair(p, k);
    auto it = residues.find(key);
    if (it != residues.end()) {
        return it->second;
    }
    return residues.insert(make_pair(key, matrix(nrows, ncols))).first->second;
}

matrix
pfmatrix::to_matrix()
{
    matrix m(nrows, ncols);
    for (const auto &kv : residues) {
        const auto &pi = kv.first.first;
        const auto &ki = kv.first.second;
        const auto &Ci = kv.second;
        assert((ki < 0) || pi.is_zero());
        m = m.add(Ci.mul_scalar(pow(x - pi, ki)));
    }
    return m;
}

void
pfmatrix::normalize()
{
    for (auto it = residues.begin(); it != residues.end();) {
        auto &ci = it->second;
        ci = normal(ci);
        if (ci.is_zero_matrix()) {
            it = residues.erase(it);
        } else {
            it++;
        }
    }
}

void
pfmatrix::add(const matrix &C, const ex &p1, int k1)
{
    // M += C*(x-p1)^k1
    matrix &R = (*this)(p1, k1);
    R = R.add(C);
}

void
pfmatrix::add_div(const matrix &C, const ex &p1, int k1, const ex &p2)
{
    // M += C*(x-p1)/(x-p2)
    if (p1 == p2) {
        add(C, p1, k1 - 1);
    }
    else if (k1 >= 0) {
        assert(p1 == 0);
        for (int k = 0; k < k1; k++) {
            add(C.mul_scalar(pow(p2, k1 - 1 - k)), p1, k);
        }
        add(C.mul_scalar(pow(p2, k1)), p2, -1);
    }
    else {
        ex dx = p2 - p1;
        for (int k = k1; k < 0; k++) {
            add(C.mul_scalar(-pow(dx, k1 - 1 - k)), p1, k);
        }
        add(C.mul_scalar(pow(dx, k1)), p2, -1);
    }
}

void
pfmatrix::add_mul(const matrix &C, const ex &p1, int k1, const ex &p2)
{
    // M += C*(x-p1)^k1*(x-p2) = C*(x-p1)^k1*{(x-p1) + (p1-p2)}
    add(C.mul_scalar(p1 - p2), p1, k1);
    if (k1 == -1) {
        add(C, 0, 0);
    }
    else {
        add(C, p1, k1 + 1);
    }
}

void
pfmatrix::add_pow(const matrix &C, const ex &p1, int k1, const ex &p2, int k2)
{
    if (k1 >= 0) assert(p1 == 0);
    if (k2 >= 0) assert(p2 == 0);
    if (p1 == p2) {
        add(C, p1, k1 + k2);
    }
    else if (k1 == 0) {
        add(C, p2, k2);
    }
    else if (k2 == 0) {
        add(C, p1, k1);
    }
    else if ((k1 >= 0) && (k2 >= 0)) {
        add(C, p1, k1 + k2);
    }
    else if ((k1 < 0) && (k2 < 0)) {
        ex d = p1 - p2;
        add(C.mul_scalar(pow(d, k2)), p1, k1);
        for (int i = k2; i < 0; i++) {
            add_pow(C.mul_scalar(-pow(d, i)),
                    (k1 + 1 == 0) ? 0 : p1, k1 + 1,
                    p2, k2 - i - 1);
        }
    }
    else if ((k1 < 0) && (k2 > 0)) {
        long b = 1;
        for (int i = 0; i < min(-k1, k2 + 1); i++) {
            add(C.mul_scalar(b*pow(p1, k2 - i)), p1, k1 + i);
            b = b*(k2 - i)/(i + 1);
        }
        for (int i = min(-k1, k2 + 1); i < k2 + 1; i++) {
            long b2 = b;
            for (int j = 0; j < k1 + i + 1; j++) {
                add(C.mul_scalar(((k1 + i - j) % 2 == 0 ? b2 : -b2)*pow(p1, k1 + k2 - j)), 0, j);
                b2 = b2*(k1 + i - j)/(j + 1);
            }
            b = b*(k2 - i)/(i + 1);
        }
    }
    else if ((k1 > 0) && (k2 < 0)) {
        add_pow(C, p2, k2, p1, k1);
    }
    else {
        assert(false);
    }
}

pfmatrix
pfmatrix::with_constant_t(const matrix &T) const
{
    return with_constant_t(T.inverse(), T);
}

pfmatrix
pfmatrix::with_constant_t(const matrix &L, const matrix &R) const
{
    pfmatrix m(nrows, ncols, x);
    for (auto &kv : residues) {
        const auto &pi = kv.first.first;
        const auto &ki = kv.first.second;
        const auto &Ci = kv.second;
        m(pi, ki) = L.mul(Ci).mul(R);
    }
    return m;
}

pfmatrix
pfmatrix::with_balance_t(const matrix &P, const ex &x1, const ex &x2) const
{
    pfmatrix m(nrows, ncols, x);
    const matrix coP = ex_to_matrix(unit_matrix(P.rows()) - P);
    if (x1 == infinity) {
        const matrix neg_coP = coP.mul_scalar(-1);
        for (const auto &kv : residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &Ci = kv.second;
            // coP Ci coP (x-pi)^ki + P Ci P (x-pi)^ki
            m.add(coP.mul(Ci).mul(coP).add(P.mul(Ci).mul(P)), pi, ki);
            // coP Ci P -(x-x2) (x-pi)^ki
            m.add_mul(neg_coP.mul(Ci).mul(P), pi, ki, x2);
            // P Ci coP -1/(x-x2) (x-pi)^ki
            m.add_div(P.mul(Ci).mul(neg_coP), pi, ki, x2);
        }
        // -P/(x-x2)
        m.add(P.mul_scalar(-1), x2, -1);
    }
    else if (x2 == infinity) {
        const matrix neg_coP = coP.mul_scalar(-1);
        for (const auto &kv : residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &Ci = kv.second;
            // coP Ci coP (x-pi)^ki + P Ci P (x-pi)^ki
            m.add(coP.mul(Ci).mul(coP).add(P.mul(Ci).mul(P)), pi, ki);
            // P Ci coP -(x-x1) (x-pi)^ki
            m.add_mul(P.mul(Ci).mul(neg_coP), pi, ki, x1);
            // coP Ci P -1/(x-x1) (x-pi)^ki
            m.add_div(neg_coP.mul(Ci).mul(P), pi, ki, x1);
        }
        // P/(x-x1)
        m.add(P, x1, -1);
    }
    else {
        // Ci (x-pi)^ki
        m.residues = residues;
        const matrix neg_P = P.mul_scalar(-1);
        const matrix coP_x12 = coP.mul_scalar(x1-x2);
        for (const auto &kv : residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &Ci = kv.second;
            // coP Ci P (x1-x2)/(x-x1) (x-pi)^ki
            m.add_div(coP_x12.mul(Ci).mul(P), pi, ki, x1);
            // P Ci coP (x2-x1)/(x-x2) (x-pi)^ki
            m.add_div(neg_P.mul(Ci).mul(coP_x12), pi, ki, x2);
        }
        // P/(x-x1) - P/(x-x2)
        m.add(P, x1, -1);
        m.add(neg_P, x2, -1);
    }
    m.normalize();
    return m;
}

/* Apply a balance with projector P between x=x1 and x=x2 to a
 * given matrix m.
 */
matrix
with_balance_t(const matrix &m, const matrix &P, const ex &x1, const ex &x2, const ex &x)
{
    matrix coP = P.mul_scalar(-1);
    for (unsigned i = 0; i < coP.cols(); i++) {
        coP(i, i) += 1;
    }
    ex k, d;
    if (x1 == infinity) {
        k = -(x - x2);
        d = -1;
    }
    else if (x2 == infinity) {
        k = -1/(x - x1);
        d = 1/pow(x - x1, 2);
    }
    else {
        k = (x - x2)/(x - x1);
        d = (x2 - x1)/pow(x - x1, 2);
    }
    return ex_to_matrix((coP + 1/k*P)*m*(coP + k*P) - d/k*P);
}

/* Block-triangular permutation of a matrix.
 * This class computes it, and keeps the results.
 *
 * This one is a class instead of a function, because
 * recursive closures are awkward (impossible?) in C++.
 */
class block_triangular_permutation {
    public:
    block_triangular_permutation(const matrix &m);
    const matrix & t();
    const vector<int> &block_size();

    private:
    struct node {
        int index;
        int link;
        bool onstack;
        node() : index(-1) {}
    };
    const matrix &m;
    vector<node> nodes;
    stack<int> st;
    matrix result_t;
    int index_counter;
    int output_counter;
    vector<int> result_block_size;

    void visit(int i);
};

block_triangular_permutation::
block_triangular_permutation(const matrix &m)
    : m(m),
        nodes(m.rows()),
        result_t(m.rows(), m.cols()),
        index_counter(0),
        output_counter(0)
{
    for (int i = 0; i < (int)m.rows(); i++) {
        if (nodes[i].index == -1)
            visit(i);
    }
}

const matrix &
block_triangular_permutation::t()
{ return result_t; }

const vector<int> &
block_triangular_permutation::block_size()
{ return result_block_size; }

void
block_triangular_permutation::visit(int i)
{
    nodes[i].index = nodes[i].link = index_counter++;
    st.push(i);
    nodes[i].onstack = true;
    for (int j = 0; j < (int)m.rows(); j++) {
        if (i == j) continue;
        if (m(i, j).is_zero()) continue;
        if (nodes[j].index == -1) {
            visit(j);
            nodes[i].link = min(nodes[i].link, nodes[j].link);
        }
        else if (nodes[j].onstack) {
            nodes[i].link = min(nodes[i].link, nodes[j].index);
        }
    }
    if (nodes[i].link == nodes[i].index) {
        int size = 0;
        for (;;) {
            int j = st.top();
            st.pop();
            nodes[j].onstack = false;
            result_t(j, output_counter++) = 1;
            size++;
            if (j == i) break;
        }
        result_block_size.push_back(size);
    }
}

/* Compute and return characteristic polynomial of a matrix by
 * first permuting it into block-triangular shape.
 *
 * This is significantly faster than the naive approach for
 * larger (e.g. 10x10), sparse matrices.
 */
ex
charpoly_by_blocks(const matrix &m, const ex &lambda)
{
    block_triangular_permutation btp(m);
    matrix mm = btp.t().transpose().mul(m).mul(btp.t());
    ex cp = 1;
    int o = 0;
    for (int size : btp.block_size()) {
        matrix b = ex_to<matrix>(sub_matrix(mm, o, size, o, size));
        cp *= b.charpoly(lambda);
        o += size;
    }
    return cp;
}

template<typename F> void
factor_twice_and_iter(const ex &e, F yield)
{
    factor_iter(factor(e),
        [&](const ex &f1, int k1) {
            factor_iter(factor(f1),
                [&](const ex &f2, int k2) {
                    yield(f2, k1*k2);
                });
        });
}

/* Find and return eigenvalues of a matrix.
 * Throw an error, if unable.
 */
map<ex, unsigned, ex_is_less>
eigenvalues(const matrix &m)
{
    map<ex, unsigned, ex_is_less> eigenvalues;
    symbol lambda("λ");
    ex charpoly = charpoly_by_blocks(m, lambda);
    factor_and_iter(charpoly,
        [&](const ex &f, int k) {
            int deg = f.degree(lambda);
            if (deg == 0) {
                // No roots in λ here.
            }
            else if (deg == 1) {
                // f == c0 + c1*λ
                ex c0 = f.coeff(lambda, 0);
                ex c1 = f.coeff(lambda, 1);
                eigenvalues[ratcan(-c0/c1)] += k;
            }
            else {
                throw runtime_error("eigenvalues(): can't solve equations of 2nd degree or higher");
            }
        });
    return eigenvalues;
}

/* Since the internals of matrix class are protected, we need
 * to employ this hack to access them.
 */
class matrix_hack : public matrix {
    public:
    void append_rows(const matrix &src);
    exvector &mvec();
    void resize(unsigned nrows);
};

void
matrix_hack::append_rows(const matrix &src)
{
    assert(col == src.cols());
    row += src.rows();
    m.insert(m.end(), ((const matrix_hack*)&src)->m.begin(), ((const matrix_hack*)&src)->m.end());
}

exvector &
matrix_hack::mvec()
{
    ensure_if_modifiable();
    return m;
}

void
matrix_hack::resize(unsigned nrows)
{
    row = nrows;
    m.resize(row*col);
}

/* Multiply a submatrix by the LCM of all of its denominators.
 * Divide it by the GCD of all of its numerators.
 * Apply normal() on each element.
 */
void
rescale_submatrix(matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc)
{
    vector<ex> n(nr*nc), d(nr*nc);
    ex mul = 0;
    ex div = 0;
    auto it_n = n.begin();
    auto it_d = d.begin();
    for (unsigned i = r; i < r + nr; i++) {
        for (unsigned j = c; j < c + nc; j++) {
            ex nd = m(i, j).normal().numer_denom();
            ex numer = nd.op(0);
            ex denom = nd.op(1);
            *it_n++ = numer;
            *it_d++ = denom;
            div = div.is_zero() ? numer : gcd(div, numer);
            mul = mul.is_zero() ? denom : lcm(mul, denom);
        }
    }
    if (div.is_zero()) return;
    if (mul.is_zero()) return;
    // It would be tempting to exit here if div=mul=1, but we'd
    // then discard the normal() call results.
    it_n = n.begin();
    it_d = d.begin();
    for (unsigned i = r; i < r + nr; i++) {
        for (unsigned j = c; j < c + nc; j++) {
            ex nn, dd;
            bool ok1 = divide(*it_n++, div, nn);
            bool ok2 = divide(mul, *it_d++, dd);
            assert(ok1 && ok2);
            m(i, j) = nn*dd;
        }
    }
}

/* Transform a given matrix into upper-echelon form via Gauss
 * elimination.
 *
 * Requires O(N^3) GCD operations and about the same amount of
 * arithmetic ones for dense matrices, but about O(N*K) for
 * sparse matrices with K entries.
 */
void
echelon_form_gauss(matrix &m)
{
    unsigned nr = m.rows();
    unsigned nc = m.cols();
    exvector &mv = ((matrix_hack*)&m)->mvec();
    for (unsigned i = 0; i < nr*nc; i++) {
        mv[i] = mv[i].normal();
    }
    unsigned r0 = 0;
    unsigned c0 = 0;
    for (; (c0 < nc) && (r0 < nr - 1); c0++) {
        for (unsigned r = r0; r < nr; r++) {
            mv[r*nc + c0] = mv[r*nc + c0].normal();
        }
        // No normalization before is_zero() here, because
        // we maintain the matrix normalized throughout the
        // algorithm.
        unsigned pivot = r0;
        while ((pivot < nr) && mv[pivot*nc + c0].is_zero()) {
            pivot++;
        }
        if (pivot == nr) {
            // The whole column below r0:c0 is zero, let's skip
            // it.
            continue;
        }
        if (pivot > r0) {
            // Found a non-zero row somewhere below r0; let's
            // swap it in.
            for (unsigned c = c0; c < nc; c++) {
                mv[pivot*nc + c].swap(mv[r0*nc + c]);
            }
        }
        ex a = mv[r0*nc + c0];
        for (unsigned r = r0 + 1; r < nr; r++) {
            ex b = mv[r*nc + c0];
            if (!b.is_zero()) {
                ex k = b/a;
                mv[r*nc + c0] = 0;
                for (unsigned c = c0 + 1; c < nc; c++) {
                    mv[r*nc + c] = normal(mv[r*nc + c] - k*mv[r0*nc + c]);
                }
            }
        }
        r0++;
    }
    // Zero out the remaining rows (just in case).
    for (unsigned r = r0 + 1; r < nr; r++) {
        for (unsigned c = 0; c < nc; c++) {
            mv[r*nc + c] = 0;
        }
    }
}

/* Transform a given matrix into upper-echelon form via Bareiss
 * elimination.
 *
 * Requires O(N^2) GCD operations and O(N^3) arithmetic ones.
 * Generally faster than Gauss on dense matrices. Don't use this
 * on sparse matrices: it is slower than Gauss and you'll get
 * pretty big common factors on each row.
 */
void
echelon_form_bareiss_plain(matrix &m)
{
    unsigned nr = m.rows();
    unsigned nc = m.cols();
    for (unsigned r = 0; r < nr; r++) {
        rescale_submatrix(m, r, 1, 0, nc);
    }
    exvector &mv = ((matrix_hack*)&m)->mvec();
    ex z = 1;
    unsigned r0 = 0;
    for (unsigned c0 = 0; (c0 < nc) && (r0 < nr - 1); c0++) {
        unsigned pivot = r0;
        // No normalization before is_zero() here, because
        // we maintain the matrix normalized throughout the
        // algorithm.
        while ((pivot < nr) && mv[pivot*nc + c0].is_zero()) pivot++;
        if (pivot == nr) {
            // The whole column below r0:c0 is zero, let's skip
            // it.
            continue;
        }
        if (pivot > r0) {
            // Found a non-zero row somewhere below r0; let's
            // swap it in.
            for (unsigned c = c0; c < nc; c++) {
                mv[pivot*nc + c].swap(mv[r0*nc + c]);
            }
        }
        ex a = mv[r0*nc + c0];
        for (unsigned r = r0 + 1; r < nr; r++) {
            ex b = mv[r*nc + c0];
            for (unsigned c = c0 + 1; c < nc; c++) {
                ex v;
                bool divok = divide(a*mv[r*nc + c] - b*mv[r0*nc + c], z, v);
                assert(divok);
                // Since v is a polynomial, the call to normal()
                // here does not involve GCD, but rather just
                // expands sums while preserving common factors,
                // which is nomally preferable to full expand().
                m[r*nc + c] = normal(v);
            }
            mv[r*nc + c0] = 0;
        }
        z = a;
        r0++;
    }
    // Zero out the remaining rows (just in case).
    for (unsigned r = r0 + 1; r < nr; r++) {
        for (unsigned c = 0; c < nc; c++) {
            mv[r*nc + c] = 0;
        }
    }
}

/* Transform a given matrix into upper-echelon form via a version
 * of Bareiss elimination.
 *
 * Requires O(N^2) GCD operations and O(N^3) arithmetic ones.
 * Generally faster than Gauss on dense matrices, but slower on
 * sparse ones.
 */
void
echelon_form_bareiss(matrix &m)
{
    unsigned nr = m.rows();
    unsigned nc = m.cols();
    for (unsigned r = 0; r < nr; r++) {
        rescale_submatrix(m, r, 1, 0, nc);
    }
    exvector &mv = ((matrix_hack*)&m)->mvec();
    exvector extraf(nr);
    for (unsigned r = 0; r < nr; r++) { extraf[r] = 1; }
    ex z = 1;
    unsigned r0 = 0;
    for (unsigned c0 = 0; (c0 < nc) && (r0 < nr - 1); c0++) {
        unsigned pivot = r0;
        // No normalization before is_zero() here, because
        // we maintain the matrix normalized throughout the
        // algorithm.
        while ((pivot < nr) && mv[pivot*nc + c0].is_zero()) { pivot++; }
        if (pivot == nr) {
            // The whole column below r0:c0 is zero, we can skip
            // it.
            continue;
        }
        if (pivot > r0) {
            // Found a non-zero row somewhere below r0; let's
            // swap it in.
            for (unsigned c = c0; c < nc; c++) {
                mv[pivot*nc + c].swap(mv[r0*nc + c]);
            }
            extraf[pivot].swap(extraf[r0]);
        }
        ex a = mv[r0*nc + c0];
        for (unsigned r = r0 + 1; r < nr; r++) {
            ex b = mv[r*nc + c0];
            if (b.is_zero()) {
                // This branch is the same as the next one, but
                // is is simplified due to b being 0.
                ex f;
                assert(divide(a*extraf[r0]*extraf[r], z, f));
                extraf[r] = f;
            } else {
                ex g = gcd(a*extraf[r0], b*extraf[r]);
                ex f;
                gcd(extraf[r0]*g, z, &f);
                f = normal(f);
                ex d;
                assert(divide(z*f, extraf[r0]*extraf[r], d));
                for (unsigned c = c0 + 1; c < nc; c++) {
                    ex v;
                    assert(divide(a*mv[r*nc + c] - b*mv[r0*nc + c], d, v));
                    // Since v is a polynomial, the call to normal()
                    // here does not involve GCD, but rather just
                    // expands sums while preserving common factors,
                    // which is nomally preferable to full expand().
                    m[r*nc + c] = normal(v);
                }
                mv[r*nc + c0] = 0;
                extraf[r] = f;
            }
        }
        z = a*extraf[r0];
        r0++;
    }
    // Zero out the remaining rows (just in case).
    for (unsigned r = r0 + 1; r < nr; r++) {
        for (unsigned c = 0; c < nc; c++) {
            mv[r*nc + c] = 0;
        }
    }
    // At this point mv[r*nc + c]*extraf[r] is precisely what
    // plain Bareiss algorithm would have computed; we can now
    // discard this extra factor freely.
}

/* A vector (sub-)space represented by a set of basis vectors.
 */
struct vspace {
    matrix basis;
    vspace(unsigned n);
    vspace(const matrix &b);
    unsigned dim() const;
    matrix basis_col(unsigned i) const;
    matrix basis_row(unsigned i) const;
    const matrix basis_cols() const;
    const matrix &basis_rows() const;
    bool contains(const matrix &v) const;
    void add_rows(const matrix &v);
    void normalize();
};

vspace::vspace(unsigned n)
    : basis(0, n)
{
}

vspace::vspace(const matrix &b)
    : basis(b)
{
    for (unsigned i = 0; i < basis.rows(); i++) {
        rescale_submatrix(basis, i, 1, 0, basis.cols());
    }
    normalize();
}

void
vspace::add_rows(const matrix &v)
{
    ((matrix_hack*)&basis)->append_rows(v);
}

void
vspace::normalize()
{
    echelon_form_gauss(basis);
    //echelon_form_bareiss(basis);
    unsigned nrows = basis.rows();
    for (; nrows > 0; nrows--) {
        for (unsigned c = 0; c < basis.cols(); c++) {
            if (!basis(nrows - 1, c).normal().is_zero()) goto done;
        }
    }
done:;
    ((matrix_hack*)&basis)->resize(nrows);
    for (unsigned i = 0; i < basis.rows(); i++) {
        rescale_submatrix(basis, i, 1, 0, basis.cols());
    }
}

unsigned
vspace::dim() const
{
    return basis.rows();
}

matrix
vspace::basis_col(unsigned i) const
{
    matrix v(basis.cols(), 1);
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

matrix
vspace::basis_row(unsigned i) const
{
    matrix v(1, basis.cols());
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

const matrix &
vspace::basis_rows() const
{
    return basis;
}

const matrix
vspace::basis_cols() const
{
    return basis.transpose();
}

bool
vspace::contains(const matrix &v) const
{
    assert(v.nops() == basis.cols());
    matrix vv = v;
    rescale_submatrix(vv, 0, v.rows(), 0, v.cols());
    unsigned p = 0;
    // Division-free subtraction of basis vectors from v.
    for (unsigned i = 0; i < basis.rows(); i++, p++) {
        for (;;) {
            if (!basis(i, p).is_zero()) break;
            vv.let_op(p) = normal(vv.op(p));
            if (!vv.op(p).is_zero())
                return false;
            p++;
        }
        const ex &vv_p = vv.op(p);
        if (!vv_p.is_zero()) {
            const ex b_ip = basis(i, p);
            vv.let_op(p) = 0;
            for (unsigned j = p + 1; j < basis.cols(); j++) {
                vv.let_op(j) = normal(vv.op(j)*b_ip - basis(i, j)*vv_p);
            }
        }
    }
    for (unsigned i = p; i < basis.cols(); i++) {
        vv.let_op(i) = normal(vv.op(i));
        if (!vv.op(i).is_zero())
            return false;
    }
    return true;
}

/* Find the null space (a.k.a. the kernel) of a given matrix.
 */
vspace
nullspace(const matrix &m)
{
    unsigned nrows = m.rows();
    unsigned ncols = m.cols();
    matrix v(ncols, 1);
    std::vector<symbol> tmp;
    for (unsigned i = 0; i < ncols; i++) {
        symbol t;
        v(i, 0) = t;
        tmp.push_back(t);
    }
    // Solve M*V = 0
    matrix zero(nrows, 1);
    matrix s = normal(m).solve(v, zero, solve_algo::gauss);
    matrix coeff(ncols, ncols);
    for (unsigned k = 0; k < ncols; k++) {
        for (unsigned i = 0; i < ncols; i++) {
            coeff(k, i) = s(i, 0).coeff(tmp[k]);
        }
    }
    return vspace(coeff);
}

/* Find the (right) eigenspace corresponding to a given eigenvalue.
 */
vspace
eigenspace(const matrix &m, const ex &eval)
{
    matrix mm = m;
    for (unsigned i = 0; i < m.rows(); i++) {
        mm(i, i) -= eval;
    }
    return nullspace(mm);
}

/* Find the left eigenspace corresponding to a given eigenvalue.
 */
vspace
eigenspace_left(const matrix &m, const ex &eval)
{
    return eigenspace(m.transpose(), eval);
}

/* Find the set of right eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_right(const matrix &m, const ex &eval)
{
    vspace es = eigenspace(m, eval);
    vector<matrix> evectors(es.dim());
    for (unsigned i = 0; i < es.dim(); i++) {
        evectors[i] = es.basis_col(i);
    }
    return evectors;
}

/* Find the set of left eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_left(const matrix &m, const ex &eval)
{
    vspace es = eigenspace_left(m, eval);
    vector<matrix> evectors(es.dim());
    for (unsigned i = 0; i < es.dim(); i++) {
        evectors[i] = es.basis_row(i);
    }
    return evectors;
}

/* Compute the Jordan normal form J of a given matrix M.
 * Return a transformation matrix Q, such that Q^-1 M Q = J.
 * Also return a list of Jordan cell sizes.
 */
pair<matrix, vector<int>>
jordan(const matrix &m)
{
    assert(m.cols() == m.rows());
    unsigned n = m.rows();
    map<ex, unsigned, ex_is_less> eval2almul = eigenvalues(m);
    matrix q(n, n);
    vector<int> jcs;
    int idxq = 0;
    for (const auto &kv : eval2almul) {
        const auto &eval = kv.first;
        const auto &almul = kv.second;
        matrix mm = m;
        for (unsigned i = 0; i < n; i++) {
            mm(i, i) -= eval;
        }
        vector<vspace> nspace;
        vspace xspace(n);
        matrix mmk = mm;
        for (;;) {
            auto ns = nullspace(mmk);
            nspace.push_back(ns);
            if (ns.dim() == almul) break;
            mmk = mmk.mul(mm);
        }
        int kk = 0;
        for (int i = nspace.size() - 1; i >= 0; i--) {
            /* We're expecting kt - kk chains to start at rank i.
             */
            int kt = (i == 0) ?
                nspace[i].dim() : nspace[i].dim() - nspace[i-1].dim();
            for (int k = nspace[i].dim() - 1; kk < kt; k--) {
                /* Find a vector in nspace[i] that is not in
                 * nspace[i-1] (and thus, a chain head), and not
                 * in xspace (not excluded).
                 */
                assert(k >= 0);
                auto v = nspace[i].basis_col(k);
                if ((i != 0) && nspace[i-1].contains(v)) continue;
                if (xspace.contains(v)) continue;
                /* Exclude v from the future results (so no
                 * other vector is linearly dependent on it),
                 * and build a Jordan chain starting from v.
                 */
                kk++;
                xspace.add_rows(v.transpose());
                for (unsigned r = 0; r < n; r++) {
                    q(r, idxq + i) = v.op(r);
                }
                for (int j = i - 1; j >= 0; j--) {
                    v = mm.mul(v);
                    //assert(!xspace.contains(v));
                    xspace.add_rows(v.transpose());
                    for (unsigned r = 0; r < n; r++) {
                        q(r, idxq + j) = v.op(r);
                    }
                }
                /* Only rescale chains of length > 1; single
                 * eigenvectors are already well scaled by
                 * vspace().
                 */
                if (i != 0) {
                    rescale_submatrix(q, 0, n, idxq, i+1);
                }
                xspace.normalize();
                idxq += i + 1;
                jcs.push_back(i + 1);
            }
            assert(kk == kt);
        }
    }
    //assert(q.rank() == n);
    return make_pair(q, jcs);
}

/* Find matrix V, such that its rows span a left invariant
 * subspace of B and form a dual basis with columns of U
 * (that is, V*U=I).
 *
 * There may be zero, one, or many such matrices. Return a list
 * of them.
 */
vector<matrix>
dual_basis_spanning_left_invariant_subspace(const matrix &m, const matrix &u)
{
    LOGME;
    matrix mt = m.transpose();
    matrix lev(0, m.rows()); // left (row-)eigenvectors of m
    for (const auto &ev : eigenvalues(mt)) {
        const auto &eval = ev.first;
        vspace vs = eigenspace(mt, eval);
        ((matrix_hack*)&lev)->append_rows(vs.basis_rows());
        auto v = vs.basis_rows().mul(m);
    }
    matrix tmp(u.cols(), lev.rows());
    exmap tmpz;
    for (unsigned i = 0; i < tmp.nops(); i++) {
        symbol t;
        tmp.let_op(i) = t;
        tmpz[t] = 0;
    }
    matrix identity = identity_matrix(u.cols());
    vector<matrix> results;
    try {
        // Solve x*lev*u=1, and take x*lev as the result.
        matrix x = matrix_solve_left(lev.mul(u), tmp, identity);
        for (unsigned i = 0; i < x.nops(); i++) {
            x.let_op(i) = x.op(i).subs(tmpz);
        }
        results.push_back(x.mul(lev));
    } catch (const std::runtime_error &e) {
        if (e.what() != std::string("matrix::solve(): inconsistent linear system")) {
            throw;
        }
    }
    return results;
}

/* Fuchsification
 * ============================================================
 */

pair<vector<int>, matrix>
alg1(matrix l0, const vector<int> &jcs)
{
    LOGME;
    unsigned n = l0.rows();
    vector<int> s(n);
    matrix d(n, n);
    matrix ident = identity_matrix(n);
    unsigned nontrivialcells = 0;
    for (int cs : jcs) {
        if (cs > 1) nontrivialcells++;
    }
    for (;;) {
        matrix lx = l0;
        matrix tmp(n, 1);
        for (unsigned i = 0; i < n; i++) {
            tmp.let_op(i) = symbol();
            if (s[i]) {
                for (unsigned j = 0; j < n; j++) {
                    lx(i, j) = 0;
                    lx(j, i) = 0;
                }
            }
        }
        unsigned i = 0;
        matrix sol;
        for (;; i++) {
            assert(i < n);
            if (s[i]) continue;
            matrix lxb4i = ex_to_matrix(sub_matrix(lx, 0, n, 0, i));
            matrix tmpi = ex_to_matrix(sub_matrix(tmp, 0, i, 0, 1));
            matrix lxi = ex_to_matrix(sub_matrix(lx, 0, n, i, 1));
            try {
                sol = lxb4i.solve(tmpi, lxi);
                break;
            } catch (const std::runtime_error &e) {
                if (e.what() == std::string("matrix::solve(): inconsistent linear system")) {
                    continue;
                }
                throw;
            }
        }
        matrix d0(n, n);
        matrix invd0(n, n);
        for (unsigned j = 0; j < i; j++) {
            d0(j, i) = sol.op(j);
            invd0(j, i) = jcs[j] == jcs[i] ? -sol.op(j) : 0;
        }
        l0 = ident.sub(invd0).mul(l0).mul(ident.add(d0));
        d = d.add(d0).add(d.mul(d0));
        s[i] = true;
        if (i < nontrivialcells) break;
    }
    // Check Alg.1 promise
    for (unsigned j = 0; j < n; j++) {
        for (unsigned k = 0; k < n; k++) {
            if ((!s[j]) && (s[k])) {
                assert(l0(j, k) == 0);
            }
        }
    }
    return make_pair(s, d);
}

pair<matrix, matrix>
alg1x(const matrix &a0, const matrix &a1)
{
    LOGME;
    unsigned n = a0.rows();
    const auto &ucs = jordan(a0);
    const matrix &u = ucs.first;
    const vector<int> &jcs = ucs.second;
    matrix invu = u.inverse();
    unsigned ncells = jcs.size();
    vector<int> jce(ncells);
    vector<int> jcb(ncells);
    int nsimplecells = 0;
    for (unsigned i = 0; i < ncells; i++) {
        jce[i] = (i == 0) ? jcs[0] : jce[i - 1] + jcs[i];
        jcb[i] = (i == 0) ? 0 : jcb[i - 1] + jcs[i - 1];
        if (jcs[i] == 1) nsimplecells++;
    }
    matrix l0(ncells, ncells);
    matrix l1(ncells, ncells);
    for (unsigned k = 0; k < ncells; k++) {
        matrix v0t = ex_to_matrix(sub_matrix(invu, jce[k]-1, 1, 0, n));
        for (unsigned l = 0; l < ncells; l++) {
            matrix u0 = ex_to_matrix(sub_matrix(u, 0, n, jcb[l], 1));
            l0(k, l) = v0t.mul(a1).mul(u0).op(0);
            l1(k, l) = v0t.mul(u0).op(0);
        }
    }
    ex det = l0.sub(l1.mul_scalar(symbol("λ"))).determinant().normal();
    if (det != 0) {
        loge("det|l0 - λ*l1| = {}", det);
        throw fuchsia_error("matrix is Moser-irreducible");
    }
    const auto &sd = alg1(l0, jcs);
    const auto &s = sd.first;
    const auto &d = sd.second;
    matrix ie = identity_matrix(n);
    for (unsigned i = 0; i < ncells; i++) {
        for (unsigned j = 0; j < ncells; j++) {
            if (d(i, j) == 0) {
                int ni = jcb[i];
                int nj = jcb[j];
                for (int k = 0; k < min(jcs[i], jcs[j]); k++) {
                    ie(ni+k, nj+k) += d(i,j);
                }
            }
        }
    }
    matrix ut = u.mul(ie);
    matrix invut = ut.inverse();
    int ns = 0;
    for (unsigned i = 0; i < ncells; i++) {
        if (s[i]) ns++;
    }
    matrix ru(n, ns);
    matrix rv(ns, n);
    for (unsigned i = 0, ii = 0; i < ncells; i++) {
        if (s[i]) {
            for (unsigned j = 0; j < n; j++) {
                ru(j, ii) = ut(j, jcb[i]);
                rv(ii, j) = invut(jcb[i], j);
            }
            ii++;
        }
    }
    return make_pair(ru, rv);
}

matrix
c0_infinity(const pfmatrix &pfm)
{
    matrix m(pfm.nrows, pfm.ncols);
    for (const auto &kv : pfm.residues) {
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        if (ki == -1) m = m.sub(ci);
    }
    return m;
}

int
complexity(const ex &e)
{
    if (is_exactly_a<numeric>(e)) {
        const numeric n = ex_to<numeric>(e);
        if (n.is_zero()) return 0;
        return 8 + n.numer().int_length() + n.denom().int_length();
    } else {
        int c = 15;
        for (const auto &sube : e) {
            c += complexity(sube);
        }
        return c;
    }
}

int
complexity(const pfmatrix &pfm)
{
    int c = 0;
    for (const auto kv : pfm.residues) {
        const auto &pi = kv.first.first;
        const auto &ci = kv.second;
        c += complexity(pi);
        c += complexity(ci);
    }
    return c;
}

ex
balance_t(const matrix &p, const ex &x1, const ex &x2, const ex &x)
{
    matrix cop = identity_matrix(p.rows()).sub(p);
    if (x1 == infinity) {
        return cop - (x - x2)*p;
    }
    else if (x2 == infinity) {
        return cop - 1/(x - x1)*p;
    }
    else {
        return cop + (x - x2)/(x - x1)*p;
    }
}

pair<pfmatrix, ex>
fuchsify(const pfmatrix &m)
{
    LOGME;
    pfmatrix pfm = m;
    ex t = unit_matrix(m.nrows, m.ncols);
    // 1. Take a look at all {pi->pj} balances (where ki<>-1 and
    // kj<>-1), choose the least complex one.
    //
    // 2. Apply the found balance.
    //
    // 3. Repeat until Fuchsian.
    struct Reduction { ex pi; ex pj; matrix p; pfmatrix pfm; };
    for (;;) {
        logd("Current matrix complexity: {}", complexity(pfm));
        logd("Current matrix expansion:");
        for (const auto kv : pfm.residues) {
            const auto &xi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &c = kv.second;
            if (c.is_zero_matrix()) continue;
            logd("- pole of order {} at {}, complexity={}", ki, xi, complexity(c));
        }
        // Compute the highest powers of decomposition at each
        // singular point (including infinity).
        map<ex, int, ex_is_less> poincare_map;
        for (auto &kv : pfm.residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &ci = kv.second;
            if (ci.is_zero_matrix()) {
                continue;
            }
            if (ki <= -1) {
                if (poincare_map.find(pi) == poincare_map.end()) poincare_map[pi] = -1;
                poincare_map[pi] = min(poincare_map[pi], ki);
            }
            if (ki >= -1) {
                if (poincare_map.find(infinity) == poincare_map.end()) poincare_map[infinity] = -1;
                poincare_map[infinity] = max(poincare_map[infinity], ki);
            }
        }
        /* List all possible reductions between each pair of
         * singular points.
         */
        vector<Reduction> reductions;
        vector<Reduction> poor_reductions;
        bool done = true;
        for (auto &kvi : poincare_map) {
            const auto &pi = kvi.first;
            const auto &ki = kvi.second;
            if (ki == -1) continue;
            done = false;
            const auto a0 = (pi == infinity) ? pfm(0, ki) : pfm(pi, ki);
            if (a0.is_zero_matrix()) continue;
            const auto a1 = (pi == infinity) ? (ki == 0) ? c0_infinity(pfm) : pfm(0, ki - 1) : pfm(pi, ki + 1);
            const auto uv = alg1x(a0, a1);
            for (auto &kvj : poincare_map) {
                const auto &pj = kvj.first;
                const auto &kj = kvj.second;
                if (pi == pj) continue;
                /* To guarantee progress, we only reduce from
                 * points with higher Poincare rank to the ones
                 * with lower.
                 */
                if (abs(kj + 1) >= abs(ki + 1)) continue;
                //const auto b0 = (pj == infinity) ? pfm(0, kj) : pfm(pj, kj);
                const auto b0 = (pj == infinity) ? (kj == -1) ? c0_infinity(pfm) : pfm(0, kj) : pfm(pj, kj);
                if (b0.is_zero_matrix()) continue;
                for (auto dualb : dual_basis_spanning_left_invariant_subspace(b0, uv.first)) {
                    auto p = uv.first.mul(dualb);
                    reductions.push_back(Reduction {
                        pi,
                        pj,
                        p,
                        pfm.with_balance_t(p, pi, pj)
                    });
                }
            }
            ex randomp = rand()%20 - 10;
            poor_reductions.push_back(Reduction {
                pi,
                randomp,
                uv.first.mul(uv.second),
                pfm.with_balance_t(uv.first.mul(uv.second), pi, randomp)
            });
        }
        if (done) { break; }
        if (reductions.size() > 0) {
            // Find and apply the smallest reduction.
            size_t mini = 0;
            int minc = complexity(reductions[0].pfm);
            for (size_t i = 1; i < reductions.size(); i++) {
                int c = complexity(reductions[i].pfm);
                if (c < minc) {
                    mini = i;
                    minc = c;
                }
            }
            Reduction &r = reductions[mini];
            logi("use balance between {} and {} with projector:\n{}", r.pi, r.pj, r.p);
            logi("new complexity is: {}", minc);
            t = t * balance_t(r.p, r.pi, r.pj, pfm.x);
            pfm = r.pfm;
        } else {
            logd("* no good reductions found, looking at the bad ones");
            size_t mini = 0;
            int minc = complexity(poor_reductions[0].pfm);
            for (size_t i = 1; i < poor_reductions.size(); i++) {
                int c = complexity(poor_reductions[i].pfm);
                if (c < minc) {
                    mini = i;
                    minc = c;
                }
            }
            Reduction &r = poor_reductions[mini];
            logi("apply balance between {} and {} with projector:\n{}", r.pi, r.pj, r.p);
            t = t*balance_t(r.p, r.pi, r.pj, pfm.x);
            pfm = r.pfm;
        }
    }
    return make_pair(pfm, t);
}

/* Logging formatters
 * ============================================================
 */

template<> inline void
log_format(ostream &o, const matrix &m)
{
    save_matrix(o, m);
}

template<> inline void
log_format(ostream &o, const vector<int> &v)
{
    o << "[";
    for (size_t i = 0; i < v.size(); i++) {
        if (i != 0) o << ", ";
        o << v[i];
    }
    o << "]";
}

template<> inline void
log_format(ostream &o, const vspace &v)
{
    log_format(o, v.basis_rows());
}
