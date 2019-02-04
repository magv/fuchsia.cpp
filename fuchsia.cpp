#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <tuple>

using namespace GiNaC;
using namespace std;

static bool COLORS = !!isatty(STDOUT_FILENO);
static bool PARANOID = false;
static bool VERBOSE = true;

/* LOGGING
 * ============================================================
 *
 * As with everything in C++, you can "optimize" this piece of
 * code to e.g. use compile-only format string parsing, minimize
 * number of created functions, and so on, but at the expense
 * of kilolines of code, and your own sanity lost in the fight
 * versus byzantine template rules. Please don't.
 */

static auto _log_starttime = chrono::steady_clock::now();
static auto _log_lasttime = chrono::steady_clock::now();
static int _log_depth = 0;

const char*
log_adv(const char *fmt)
{
    for (int i = 0; ; i++) {
        if (fmt[i] == '{') {
            cerr.write(fmt, i);
            return fmt + i + 2;
        }
        if (fmt[i] == 0) {
            cerr.write(fmt, i);
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
log_print_start(const char *pre, const char *post)
{
    auto t = chrono::steady_clock::now();
    auto dt = chrono::duration_cast<chrono::duration<double>>(t - _log_starttime).count();
    cerr << pre << std::fixed << std::setprecision(4) << dt << "s +";
    cerr << chrono::duration_cast<chrono::duration<double>>(t - _log_lasttime).count() << "s";
    for (int i = 0; i < _log_depth; i++) {
        cerr << " *";
    }
    cerr << post;
    _log_lasttime = t;
}

template<typename T> const char *
log_print_one(const char *fmt, const T &value)
{
    fmt = log_adv(fmt);
    log_format(cerr, value);
    return fmt;
}

void
log_print_end(const char *fmt)
{
    cerr << fmt;
    if (COLORS) cerr << "\033[0m";
    cerr << endl;
}

struct _sequencehack {
    template<typename ...Args>
    _sequencehack(Args &&...) {}
};

template<typename ...Args> static void
log_fmt(const char *pre, const char *post, const char *fmt, const Args &...args)
{
    log_print_start(pre, post);
    (void) _sequencehack {
        (fmt = log_print_one(fmt, args), 0)
        ...
    };
    log_print_end(fmt);
}

/* Log an debug message. These can be suppressed by setting
 * VERBOSE to false.
 */
template<typename... Args> static inline void
logd(const char *fmt, const Args &...args)
{
    if (VERBOSE) {
        log_fmt(COLORS ? "\033[2;37m[dbg " : "[dbg ", "] ", fmt, args...);
    }
}

/* Log an information message.
 */
template<typename... Args> static inline void
logi(const char *fmt, const Args &...args)
{
    log_fmt(COLORS ? "\033[32m[inf " : "[inf ", COLORS ? "]\033[0m " : "] ", fmt, args...);
}

/* Log an error message.
 */
template<typename... Args> static inline void
loge(const char *fmt, const Args &...args)
{
    log_fmt(COLORS ? "\033[31m[err " : "[err ", "] ", fmt, args...);
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
        if (VERBOSE) { \
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

int
randint(int a, int b)
{
    return rand() % (b - a + 1) + a;
}

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
glue_matrix(const ex &a, const ex &b, unsigned offset)
{
    matrix aa = ex_to_matrix(a);
    matrix bb = ex_to_matrix(b);
    for (unsigned r = 0; r < bb.rows(); r++) {
        for (unsigned c = 0; c < bb.cols(); c++) {
            aa(offset + r, offset + c) = bb(r, c);
        }
    }
    return aa;
}

matrix
matrix_cut(const matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc)
{
    matrix res(nr, nc);
    for (unsigned i = 0; i < nr; i++) {
        for (unsigned j = 0; j < nc; j++) {
            res(i, j) = m(r + i, c + j);
        }
    }
    return res;
}

matrix
matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs)
{
    return m.transpose().solve(vars.transpose(), rhs.transpose()).transpose();
}

template <typename F> void
matrix_map_inplace(matrix &m, F f)
{
    for (unsigned i = 0; i < m.nops(); i++) {
        m.let_op(i) = f(m.op(i));
    }
}

template <typename F> matrix
matrix_map(const matrix &m, F f)
{
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.nops(); i++) {
        r.let_op(i) = f(m.op(i));
    }
    return r;
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
    for (int depth = 0;;) {
        int c = i.get();
        if (c == EOF)
            throw parse_error("got EOF before any data");
        if (isspace(c)) continue;
        if (c == '(') {
            int c = i.get();
            if (c == EOF)
                throw parse_error("got EOF before any data");
            if (c == '*') { depth++; continue; }
            throw parse_error("matrix file should start with a '{' or a (* comment *)");
        }
        if (depth == 0) {
            if (c == '{') { i.putback(c); break; }
            throw parse_error("matrix file should start with a '{' or a (* comment *)");
        } else {
            if (c == '*') {
                int c = i.get();
                if (c == EOF)
                    throw parse_error("EOF before any data");
                if (c == ')') { depth--; continue; }
            }
        }
    }
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
    if (!i) throw parse_error("the matrix file was not found");
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
    // Write to a temporary file, rename it to the correct location
    // at the end. This is done so that if the filesystem has
    // atomic renames, the resulting file would never contain
    // partial output, even if Fuchsia was Ctrl-C'ed.
    string tmpfilename = string(filename) + string(".tmp~");
    ofstream f(tmpfilename);
    save_matrix(f, m);
    f.close();
    int r = rename(tmpfilename.c_str(), filename);
    if (r != 0) {
        throw system_error(r, system_category()); 
    }
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
    return matrix_map(m, [&](auto &&e) { return ratcan(e); });
}

/* Normal form for matrices of rational expressions.
 *
 * This is similar to 'ratcan', but doesn't expand numerator and
 * denominator needlessly.
 */
matrix
normal(const matrix &m)
{
    return matrix_map(m, [&](auto &&e) { return normal(e); });
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
    ex sol = lsolve(eqlist, clist);
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
    matrix to_matrix() const;
    pfmatrix block(unsigned offset, unsigned size) const;
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
    // M' = B(x2, x1) M B(x1, x2) - P/(x - x2) + P/(x - x1)
    pfmatrix with_balance_t(const matrix &P, const ex &x1, const ex &x2) const;
    // M' = (1 - (x - p)^k D) M (1 + (x - p)^k D) - k (x - p)^(k-1) D
    pfmatrix with_off_diagonal_t(const matrix &D, const ex &p, int k) const;
};

bool
pfmatrix::key_is_less::operator ()(const key &k1, const key &k2) const
{
    if (k1.second < k2.second) return true;
    if (k1.second > k2.second) return false;
    return ex_is_less()(k1.first, k2.first);
}

pfmatrix::pfmatrix(unsigned nrows, unsigned ncols, const symbol &x)
    : nrows(nrows), ncols(ncols), x(x)
{ }

pfmatrix::pfmatrix(const matrix &m, const symbol &x)
    : nrows(m.rows()), ncols(m.cols()), x(x)
{
    LOGME;
    if (!m.has(x)) {
        logi("warning: matrix does not contain variable {}", x);
    }
    for (unsigned i = 0; i < nrows; i++) {
        logd("converting row {}", i);
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
    const auto it = residues.find(key);
    if (it != residues.end()) {
        return it->second;
    }
    return residues.insert(make_pair(key, matrix(nrows, ncols))).first->second;
}

matrix
pfmatrix::to_matrix() const
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

pfmatrix
pfmatrix::block(unsigned offset, unsigned size) const
{
    pfmatrix m(size, size, x);
    for (auto &kv : residues) {
        const auto &pi = kv.first.first;
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        matrix b = matrix_cut(ci, offset, size, offset, size);
        if (!b.is_zero_matrix())
            m(pi, ki) = b;
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
    m.normalize();
    return m;
}

pfmatrix
pfmatrix::with_balance_t(const matrix &P, const ex &x1, const ex &x2) const
{
    LOGME;
    pfmatrix m(nrows, ncols, x);
    const matrix coP = identity_matrix(P.rows()).sub(P);
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

pfmatrix
pfmatrix::with_off_diagonal_t(const matrix &D, const ex &p, int k) const
{
    // M' = (1 - (x - p)^k D) M (1 + (x - p)^k D) - k (x - p)^(k-1) D
    //    = M + (x-p)^k (M D - D M - 1/x k D)
    LOGME;
    assert(normal(D.mul(D)).is_zero_matrix());
    pfmatrix m = *this;
    for (auto &kv : residues) {
        const auto &pi = kv.first.first;
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        if (ci.is_zero_matrix()) continue;
        // (Ci D - D Ci) (x - pi)^ki (x - p)^k
        m.add_pow(ci.mul(D).sub(D.mul(ci)), pi, ki, p, k);
    }
    m(p, k - 1) = m(p, k - 1).add(D.mul_scalar(-k));
    // Only three blocks should be affected by this transformation:
    // * the block in D
    // * everything to the left of the D block
    // * everything to the bottom of the D block
    m.normalize();
    if (PARANOID) {
        auto t = identity_matrix(nrows).add(D.mul_scalar(pow(x - p, k)));
        auto tinv = identity_matrix(nrows).sub(D.mul_scalar(pow(x - p, k)));
        auto mm = tinv.mul(to_matrix().mul(t).sub(ex_to_matrix(t.diff(x))));
        assert(normal(m.to_matrix().sub(mm)).is_zero_matrix());
    }
    return m;
}

/* Transformations
 */

enum transformation_kind {
    tk_constant,
    tk_balance,
    tk_off_diagonal
};

struct transformation {
    typedef pair<transformation_kind, lst> component_t;
    typedef vector<component_t> list_t;
    list_t components;
    unsigned size;
    transformation(unsigned size) : components(), size(size) {}
    pfmatrix apply(const pfmatrix &m) const;
    matrix to_matrix() const;
    matrix to_inverse_matrix() const;
    transformation widen(unsigned size, unsigned offset) const;
    void add_constant_t(const matrix &l, const matrix &r);
    void add_balance_t(const matrix &p, const ex &x1, const ex &x2, const ex &x);
    void add_off_diagonal_t(const matrix &d, const ex &p, int k, const ex &x);
    void add(const transformation &t);
};

matrix
transform(const matrix &m, const matrix &t, const symbol &x)
{
    return t.inverse().mul(m.mul(t).sub(ex_to_matrix(t.diff(x))));
}

pfmatrix
transformation::apply(const pfmatrix &m) const
{
    LOGME;
    pfmatrix pfm = m;
            ex balance_t_matrix(const matrix &p, const ex &x1, const ex &x2, const ex &x);
    for (const auto &c : components) {
        switch (c.first) {
        case tk_constant:
            pfm = pfm.with_constant_t(
                    ex_to_matrix(c.second.op(0)),
                    ex_to_matrix(c.second.op(1)));
            break;
        case tk_balance:
            pfm = pfm.with_balance_t(
                    ex_to_matrix(c.second.op(0)),
                    c.second.op(1),
                    c.second.op(2));
            break;
        case tk_off_diagonal:
            pfm = pfm.with_off_diagonal_t(
                    ex_to_matrix(c.second.op(0)),
                    c.second.op(1),
                    ex_to<numeric>(c.second.op(2)).to_int());
            break;
        }
    }
    return pfm;
}

matrix
transformation::to_matrix() const
{
    ex balance_t_matrix(const matrix &p, const ex &x1, const ex &x2, const ex &x);
    ex m = unit_matrix(size);
    for (const auto &c : components) {
        switch (c.first) {
        case tk_constant:
            m *= c.second.op(1);
            break;
        case tk_balance:
            m *= balance_t_matrix(
                    ex_to_matrix(c.second.op(0)),
                    c.second.op(1),
                    c.second.op(2),
                    c.second.op(3));
            break;
        case tk_off_diagonal:
            m *= unit_matrix(size) +
                    pow(c.second.op(3) - c.second.op(1), c.second.op(2))*
                        c.second.op(0);
            break;
        }
    }
    return ex_to_matrix(m);
}

matrix
transformation::to_inverse_matrix() const
{
    ex balance_t_matrix(const matrix &p, const ex &x1, const ex &x2, const ex &x);
    ex m = unit_matrix(size);
    for (const auto &c : components) {
        switch (c.first) {
        case tk_constant:
            m = c.second.op(0)*m;
            break;
        case tk_balance:
            m = balance_t_matrix(
                    ex_to_matrix(c.second.op(0)),
                    c.second.op(2),
                    c.second.op(1),
                    c.second.op(3))*m;
            break;
        case tk_off_diagonal:
            m = (unit_matrix(size) -
                    pow(c.second.op(3) - c.second.op(1), c.second.op(2))*
                        c.second.op(0))*m;
            break;
        }
    }
    return ex_to_matrix(m);
}

transformation
transformation::widen(unsigned size, unsigned offset) const
{
    LOGME;
    transformation t(size);
    for (const auto &c : components) {
        switch (c.first) {
        case tk_constant:
            t.add_constant_t(
                glue_matrix(unit_matrix(size), c.second.op(0), offset),
                glue_matrix(unit_matrix(size), c.second.op(1), offset));
            break;
        case tk_balance:
            t.add_balance_t(
                glue_matrix(matrix(size, size), c.second.op(0), offset),
                c.second.op(1), c.second.op(2), c.second.op(3));
            break;
        case tk_off_diagonal:
            t.add_off_diagonal_t(
                glue_matrix(matrix(size, size), c.second.op(0), offset),
                c.second.op(1),
                ex_to<numeric>(c.second.op(2)).to_int(),
                c.second.op(3));
            break;
        }
    }
    return t;
}

void
transformation::add_constant_t(const matrix &l, const matrix &r)
{
    components.push_back(make_pair(tk_constant, lst{l, r}));
}

void
transformation::add_balance_t(const matrix &p, const ex &x1, const ex &x2, const ex &x)
{
    components.push_back(make_pair(tk_balance, lst{p, x1, x2, x}));
}

void
transformation::add_off_diagonal_t(const matrix &d, const ex &p, int k, const ex &x)
{
    components.push_back(make_pair(tk_off_diagonal, lst{d, p, k, x}));
}

void
transformation::add(const transformation &t)
{
    for (const auto &c : t.components) {
        components.push_back(c);
    }
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
    block_triangular_permutation(const pfmatrix &m);
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

matrix
boolean_matrix(const pfmatrix &m)
{
    matrix b(m.nrows, m.ncols);
    for (const auto &kv : m.residues) {
        const auto &ci = kv.second;
        for (unsigned k = 0; k < ci.nops(); k++) {
            if (!ci.op(k).is_zero()) b.let_op(k) = 1;
        }
    }
    return b;
}

block_triangular_permutation::
block_triangular_permutation(const pfmatrix &m)
    : block_triangular_permutation(boolean_matrix(m))
{ }

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
    LOGME;
    block_triangular_permutation btp(m);
    matrix mm = btp.t().transpose().mul(m).mul(btp.t());
    ex cp = 1;
    int o = 0;
    for (int size : btp.block_size()) {
        matrix b = matrix_cut(mm, o, size, o, size);
        cp *= b.charpoly(lambda);
        o += size;
    }
    return cp;
}

/* Compute and return determinant of a matrix by
 * first permuting it into block-triangular shape.
 */
ex
determinant_by_blocks(const matrix &m)
{
    LOGME;
    block_triangular_permutation btp(m);
    matrix mm = btp.t().transpose().mul(m).mul(btp.t());
    ex res = 1;
    int o = 0;
    for (int size : btp.block_size()) {
        matrix b = matrix_cut(mm, o, size, o, size);
        res *= b.determinant();
        o += size;
    }
    return res;
}

/* Find and return eigenvalues of a matrix.
 * Throw an error, if unable.
 */
map<ex, unsigned, ex_is_less>
eigenvalues(const matrix &m, bool skip_roots=false)
{
    LOGME;
    map<ex, unsigned, ex_is_less> eigenvalues;
    symbol lambda("L");
    ex charpoly = charpoly_by_blocks(m, lambda);
    factor_iter(factor_fixed(charpoly.numer()),
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
                if (skip_roots) {
                    logi("skipping eigenvalues with roots: RootOf[({})^{}, {}]", f, k, lambda);
                } else {
                    loge("could not factor this part of charpoly (in {}): {}", lambda, f);
                    throw fuchsia_error("eigenvalues(): eigenvalues contain roots of 2nd degree or higher");
                }
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

/* A vector (sub-)space represented by a set of basis vectors.
 *
 * The basis vectors are stored as row vectors, but can be viewed
 * as either row or column vectors; hence the *_row and *_col
 * set of functions.
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
    assert(i < basis.rows());
    matrix v(basis.cols(), 1);
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

matrix
vspace::basis_row(unsigned i) const
{
    assert(i < basis.rows());
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
        // Advance p to the first non-zero column of basis[i].
        for (;;) {
            // This assertion should only fail if the normalize()
            // was not called between add_rows() and contains().
            assert(p < basis.cols());
            if (!basis(i, p).is_zero()) break;
            vv.let_op(p) = normal(vv.op(p));
            // If vv has non-zero columns before p, it's not in
            // the basis.
            if (!vv.op(p).is_zero())
                return false;
            p++;
        }
        // Subtract basis[i] from vv, if vv[p] != 0.
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
    matrix s = normal(m).solve(v, zero);
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
        if (PARANOID) {
            assert(normal(m.mul(evectors[i]).sub(evectors[i].mul_scalar(eval))).is_zero_matrix());
        }
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
        if (PARANOID) {
            assert(normal(evectors[i].mul(m).sub(evectors[i].mul_scalar(eval))).is_zero_matrix());
        }
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
    LOGME;
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
                for (unsigned r = 0; r < n; r++) {
                    q(r, idxq + i) = v.op(r);
                }
                for (int j = i - 1; j >= 0; j--) {
                    v = mm.mul(v);
                    for (unsigned r = 0; r < n; r++) {
                        q(r, idxq + j) = v.op(r);
                    }
                }
                if (xspace.contains(v)) {
                    // The last element of a chain is in xspace.
                    // We'll skip this chain.
                    continue;
                }
                /* Only rescale chains of length > 1; single
                 * eigenvectors are already well scaled by
                 * vspace().
                 */
                if (i != 0) {
                    rescale_submatrix(q, 0, n, idxq, i+1);
                }
                xspace.add_rows(matrix_cut(q, 0, n, idxq, i+1).transpose());
                xspace.normalize();
                idxq += i + 1;
                jcs.push_back(i + 1);
                kk++;
            }
            assert(kk == kt);
        }
    }
    if (PARANOID) {
        assert(q.rank() == n);
    }
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
    for (const auto &ev : eigenvalues(mt, true)) {
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
        matrix x = matrix_solve_left(normal(lev.mul(u)), tmp, identity);
        matrix_map_inplace(x, [&](auto &&e) { return e.subs(tmpz); });
        matrix r = x.mul(lev);
        if (PARANOID) {
            matrix z = normal(r.mul(u)).sub(identity_matrix(u.cols()));
            assert(z.is_zero_matrix());
        }
        results.push_back(r);
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
    for (const auto &kv : pfm.residues) {
        const auto &pi = kv.first.first;
        const auto &ci = kv.second;
        c += complexity(pi);
        c += complexity(ci);
    }
    return c;
}

ex
balance_t_matrix(const matrix &p, const ex &x1, const ex &x2, const ex &x)
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

pair<pfmatrix, transformation>
fuchsify(const pfmatrix &m)
{
    LOGME;
    pfmatrix pfm = m;
    transformation t(m.nrows);
    // 1. Take a look at all {pi->pj} balances (where ki<>-1 and
    // kj<>-1), choose the least complex one.
    //
    // 2. Apply the found balance.
    //
    // 3. Repeat until Fuchsian.
    struct Reduction { ex pi; ex pj; matrix p; pfmatrix pfm; };
    for (;;) {
        logd("Current matrix complexity: {}", complexity(pfm));
        logd("Current expansion:");
        // Compute the highest (most distant from -1) powers
        // of decomposition at each singular point (including
        // infinity), place them in poincare_map.
        map<ex, int, ex_is_less> poincare_map;
        matrix c0inf = normal(c0_infinity(pfm));
        if (!c0inf.is_zero_matrix()) {
            poincare_map[infinity] = -1;
        }
        for (const auto &kv : pfm.residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &ci = kv.second;
            if (ci.is_zero_matrix()) continue;
            logd("- residue with power {} at {}={}, complexity={}",
                    ki, pfm.x, pi, complexity(ci));
            if (ki <= -1) {
                if (poincare_map.find(pi) == poincare_map.end())
                    poincare_map[pi] = -1;
                poincare_map[pi] = min(poincare_map[pi], ki);
            }
            if (ki > -1) {
                if (poincare_map.find(infinity) == poincare_map.end())
                    poincare_map[infinity] = -1;
                poincare_map[infinity] = max(poincare_map[infinity], ki);
            }
        }
        logd("Current poles:");
        for (auto &kvi : poincare_map) {
            const auto &xi = kvi.first;
            const auto &ki = kvi.second;
            const auto c = (xi != infinity) ? pfm(xi, ki) : (ki != -1) ? pfm(0, ki) : c0inf;
            if (c.is_zero_matrix()) continue;
            if (abs(ki + 1) > 0) {
                logd("- pole with power {} at {}={}, complexity={}, residue rank={}",
                        ki, pfm.x, xi, complexity(c), c.rank());
            } else {
                logd("- pole with power {} at {}={}, complexity={}",
                        ki, pfm.x, xi, complexity(c));
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
            const auto a0 = (pi != infinity) ? pfm(pi, ki) : pfm(0, ki);
            if (a0.is_zero_matrix()) continue;
            done = false;
            const auto a1 = (pi != infinity) ? pfm(pi, ki + 1) : (ki != 0) ? pfm(0, ki - 1) : c0inf;
            symbol lambda("L");
            // n = (a0 a1-l)
            //     (0  a0  )
            matrix n(2*a0.rows(), 2*a0.cols());
            for (unsigned i = 0; i < a0.rows(); i++) {
                for (unsigned j = 0; j < a0.cols(); j++) {
                    n(i, j) = a0(i, j);
                    n(i + a0.rows(), j + a0.cols()) = a0(i, j);
                }
            }
            for (unsigned i = 0; i < a1.rows(); i++) {
                for (unsigned j = 0; j < a1.cols(); j++) {
                    n(i, j + a0.cols()) = (i == j) ? a1(i, j) - lambda : a1(i, j);
                }
            }
            vspace ns = nullspace(n);
            vspace ws(a0.cols());
            // Find the span of coefficients of the second half of ns as a poly in lambda.
            for (unsigned i = 0; i < ns.dim(); i++) {
                int maxdeg = 0;
                for (unsigned j = 0; j < a0.cols(); j++) {
                    auto &&e = ns.basis_rows()(i, j + a0.cols());
                    maxdeg = max(maxdeg, e.degree(lambda));
                }
                matrix s(maxdeg + 1, a0.cols());
                for (unsigned j = 0; j < a0.cols(); j++) {
                    // It would be great to only expand by lambda here.
                    ex e = expand(ns.basis_rows()(i, j + a0.cols()));
                    for (int deg = 0; deg <= maxdeg; deg++) {
                        s(deg, j) = factor(normal(e.coeff(lambda, deg)));
                    }
                }
                ws.add_rows(s);
            }
            assert(!ws.basis_rows().has(lambda));
            ws.normalize();
            matrix wscols = ws.basis_cols();
            if (ws.dim() == 0) {
                throw fuchsia_error("matrix is Moser-irreducible");
            }
            if (PARANOID) {
                // ws is a subset of nullspace(a0)
                matrix z = normal(a0.mul(wscols));
                assert(z.is_zero_matrix());
            }
            for (auto &kvj : poincare_map) {
                const auto &pj = kvj.first;
                const auto &kj = kvj.second;
                if (pi == pj) continue;
                /* To guarantee progress, we only reduce from
                 * points with higher Poincare rank to the ones
                 * with lower.
                 */
                if (abs(kj + 1) >= abs(ki + 1)) continue;
                const auto b0 = (pj != infinity) ? pfm(pj, kj) : (kj != -1) ? pfm(0, kj) : c0inf;
                if (b0.is_zero_matrix()) continue;
                logi("Looking at reductions between {} and {}", pi, pj);
                bool nosinglevector = true;
                if (ws.dim() > 1) {
                    for (unsigned b = 0; b < ws.dim(); b++) {
                        logd("Looking at basis vector {}", b);
                        for (auto &&dualb : dual_basis_spanning_left_invariant_subspace(b0, ws.basis_col(b))) {
                            auto p = normal(ws.basis_col(b).mul(dualb));
                            auto pfm2 = pfm.with_balance_t(p, pi, pj);
                            logd("Complexity: {}", complexity(pfm2));
                            reductions.push_back(Reduction { pi, pj, p, pfm2 });
                            nosinglevector = false;
                        }
                    }
                }
                if (nosinglevector) {
                    for (auto &&dualb : dual_basis_spanning_left_invariant_subspace(b0, wscols)) {
                        auto p = normal(wscols.mul(dualb));
                        auto pfm2 = pfm.with_balance_t(p, pi, pj);
                        logd("Complexity: {}", complexity(pfm2));
                        reductions.push_back(Reduction { pi, pj, p, pfm2 });
                    }
                }
            }
            ex randomp = randint(-10, 10);
            matrix z(a0.rows(), a0.cols());
            for (auto &&dualb : dual_basis_spanning_left_invariant_subspace(z, wscols)) {
                auto p = wscols.mul(dualb);
                poor_reductions.push_back(Reduction {
                    pi,
                    randomp,
                    p,
                    pfm.with_balance_t(p, pi, randomp)
                });
                break;
            }
        }
        if (done) { break; }
        if (reductions.size() > 0) {
            // Find and apply the smallest reduction.
            size_t mini = 0;
            int minc = complexity(reductions[0].pfm);
            logi("Reduction between {} and {} can give complexity {}",
                    reductions[0].pi, reductions[0].pj, minc);
            logd("Projector:\n{}", reductions[0].p);
            for (size_t i = 1; i < reductions.size(); i++) {
                int c = complexity(reductions[i].pfm);
                if (c < minc) {
                    mini = i;
                    minc = c;
                }
                logi("Reduction between {} and {} can give complexity {}",
                        reductions[i].pi, reductions[i].pj, c);
                logd("Projector:\n{}", reductions[i].p);
            }
            Reduction &r = reductions[mini];
            logi("Use balance between {} and {} with projector:\n{}", r.pi, r.pj, r.p);
            t.add_balance_t(r.p, r.pi, r.pj, pfm.x);
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
            t.add_balance_t(r.p, r.pi, r.pj, pfm.x);
            pfm = r.pfm;
        }
    }
    return make_pair(pfm, t);
}

bool
is_fuchsian(const pfmatrix &pfm)
{
    for (const auto &kv : pfm.residues) {
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        if (ci.is_zero_matrix()) continue;
        if (ki != -1) return false;
    }
    return true;
}

/* Normalization
 * ============================================================
 */

ex
dot(const matrix &v1, const matrix &v2)
{
    assert(v1.nops() == v2.nops());
    ex res = 0;
    for (unsigned i = 0; i < v1.nops(); i++) {
        res += v1.op(i)*v2.op(i);
    }
    return res;
}

matrix
cross(const matrix &v1, const matrix &v2)
{
    matrix res(v1.nops(), v2.nops());
    for (unsigned i = 0; i < v1.nops(); i++) {
        const ex &a = v1.op(i);
        for (unsigned j = 0; j < v2.nops(); j++) {
            res(i, j) = a*v2.op(j);
        }
    }
    return res;
}

pair<pfmatrix, transformation>
normalize(const pfmatrix &m, const symbol &eps)
{
    LOGME;
    pfmatrix pfm = m;
    transformation t(m.nrows);
    for (;;) {
        logd("Current matrix complexity: {}", complexity(pfm));
        logd("Current expansion:");
        struct Residue { ex pi; matrix ci; ex eval; vector<matrix> evecs; };
        vector<Residue> residues1, residues2;
        for (const auto &kv : pfm.residues) {
            const auto &pi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &ci = kv.second;
            if (ci.is_zero_matrix()) continue;
            logd("- residue with power {} at {}={}, complexity={}",
                    ki, pfm.x, pi, complexity(ci));
            if (ki != -1) {
                throw fuchsia_error("normalize(): the matrix is not Fuchsian");
            }
            for (const auto kv : eigenvalues(ci, true)) {
                const auto &eval = kv.first;
                const auto &almul = kv.second;
                logd("  * eigenvalue^{}: {}", almul, eval);
                ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
                if (!is_a<numeric>(ev0)) {
                    logi("this eigenvalue can not be normalized: {}", eval);
                    throw fuchsia_error("normalize(): unsupported form of eigenvalue");
                }
                numeric ev0n = ex_to<numeric>(ev0);
                if (ev0n < numeric(-1, 2)) {
                    residues1.push_back(Residue { pi, ci, eval, eigenvectors_right(ci, eval) });
                }
                if (numeric(1, 2) <= ev0n) {
                    residues2.push_back(Residue { pi, ci, eval, eigenvectors_left(ci, eval) });
                }
            }
        }
        matrix c0inf = normal(c0_infinity(pfm));
        if (!c0inf.is_zero_matrix()) {
            logd("- effective residue with power {} at {}={}, complexity={}",
                    -1, pfm.x, infinity, complexity(c0inf));
            for (const auto kv : eigenvalues(c0inf, true)) {
                const auto &eval = kv.first;
                const auto &almul = kv.second;
                logd("  * eigenvalue^{}: {}", almul, eval);
                ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
                if (!is_a<numeric>(ev0)) {
                    logi("this eigenvalue can not be normalized: {}", eval);
                    throw fuchsia_error("normalize(): unsupported form of eigenvalue");
                }
                numeric ev0n = ex_to<numeric>(ev0);
                if (ev0n < numeric(-1, 2)) {
                    residues1.push_back(Residue { infinity, c0inf, eval, eigenvectors_right(c0inf, eval) });
                }
                if (numeric(1, 2) <= ev0n) {
                    residues2.push_back(Residue { infinity, c0inf, eval, eigenvectors_left(c0inf, eval) });
                }
            }
        }
        if (residues1.size() + residues2.size() == 0) break;
        logd("{} more eigevalues to normalize", residues1.size() + residues2.size());
        struct Reduction { ex pi; ex pj; matrix p; pfmatrix pfm; };
        vector<Reduction> reductions;
        for (const auto &r1 : residues1) {
            for (const auto evec1 : r1.evecs) {
                for (const auto &r2 : residues2) {
                    for (const auto evec2 : r2.evecs) {
                        ex scale = normal(dot(evec1, evec2));
                        if (scale.is_zero()) continue;
                        matrix p = normal(cross(evec1, evec2).mul_scalar(1/scale));
                        reductions.push_back(Reduction {
                            r1.pi,
                            r2.pi,
                            p,
                            pfm.with_balance_t(p, r1.pi, r2.pi)
                        });
                    }
                }
            }
        }
        if (reductions.size() > 0) {
            size_t mini = 0;
            logd("{} suitable reductions found", reductions.size());
            int minc = complexity(reductions[0].pfm);
            logi("Reduction between {} and {} can give complexity {}",
                    reductions[0].pi, reductions[0].pj, minc);
            for (size_t i = 1; i < reductions.size(); i++) {
                int c = complexity(reductions[i].pfm);
                if (c < minc) {
                    mini = i;
                    minc = c;
                }
                logi("Reduction between {} and {} can give complexity {}",
                        reductions[i].pi, reductions[i].pj, c);
            }
            Reduction &r = reductions[mini];
            logi("Use balance between {} and {} with projector:\n{}", r.pi, r.pj, r.p);
            t.add_balance_t(r.p, r.pi, r.pj, pfm.x);
            pfm = r.pfm;
        } else {
            loge("no suitable reductions found");
            throw fuchsia_error("normalize(): no suitable reductions found");
        }
    }
    return make_pair(pfm, t);
}

bool
is_normal(const pfmatrix &pfm, const symbol &eps)
{
    for (const auto &kv : pfm.residues) {
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        if (ci.is_zero_matrix()) continue;
        if (ki != -1) return false;
        for (const auto kv : eigenvalues(ci)) {
            const auto &eval = kv.first;
            ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
            if (!is_a<numeric>(ev0)) return false;
            numeric ev0n = ex_to<numeric>(ev0);
            if (ev0n < numeric(-1, 2)) return false;
            if (numeric(1, 2) <= ev0n) return false;
        }
    }
    matrix c0inf = normal(c0_infinity(pfm));
    if (!c0inf.is_zero_matrix()) {
        for (const auto kv : eigenvalues(c0inf)) {
            const auto &eval = kv.first;
            ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
            if (!is_a<numeric>(ev0)) return false;
            numeric ev0n = ex_to<numeric>(ev0);
            if (ev0n < numeric(-1, 2)) return false;
            if (numeric(1, 2) <= ev0n) return false;
        }
    }
    return true;
}


/* Factorization
 * ============================================================
 */

pair<pfmatrix, transformation>
factorize(const pfmatrix &m, const symbol &eps)
{
    LOGME;
    lst tmp;
    for (unsigned i = 0; i < m.nrows*m.nrows; i++) {
        tmp.append(symbol());
    }
    matrix t = matrix(m.nrows, m.nrows, tmp);
    symbol mu("MU");
    lst eqs;
    int nresidues = 0;
    for (const auto &kv : m.residues) {
        const auto &ki = kv.first.second;
        const auto &ci = kv.second;
        if (ci.is_zero_matrix()) continue;
        assert(ki == -1);
        nresidues++;
        matrix ci_eps = ci.mul_scalar(1/eps);
        matrix ci_mu(ci_eps.rows(), ci_eps.cols());
        for (unsigned i = 0; i < ci_eps.nops(); i++) {
            ci_mu.let_op(i) = ci_eps.op(i).subs(exmap{{eps, mu}});
        }
        matrix eq = ci_eps.mul(t).sub(t.mul(ci_mu));
        for (unsigned i = 0; i < eq.nops(); i++) {
            const ex &e = eq.op(i);
            if (!e.is_zero()) eqs.append(e == 0);
        }
    }
    if (nresidues == 1) {
        logd("Only one residue, let's try using Jordan form");
        for (const auto &kv : m.residues) {
            const auto &ci = kv.second;
            if (ci.is_zero_matrix()) continue;
            matrix j = jordan(ci.mul_scalar(1/eps)).first;
            matrix ij = j.inverse();
            logi("Use constant transformation:\n{}", j);
            transformation t_(m.nrows);
            t_.add_constant_t(ij, j);
            return make_pair(m.with_constant_t(ij, j), t_);
        }
    }
    logd("solving {} linear equations in {} variables", eqs.nops(), tmp.nops());
    ex sol = lsolve(eqs, tmp);
    logd("found a solution");
    matrix_map_inplace(t, [&](auto &&e) { return e.subs(sol); });
    tmp.append(mu);
    for (int range = 0;; range += 1 + range/16) {
        logd("substituting free variables, range={}", range);
        try {
            exmap map;
            for (unsigned i = 0; i < tmp.nops(); i++) {
                map[tmp.op(i)] = randint(-range, range);
            }
            matrix st = matrix_map(t, [&](auto &&e) {
                return normal(e.subs(map));
            });
            matrix invst = normal(st.inverse());
            transformation t_(m.nrows);
            t_.add_constant_t(invst, st);
            logi("Use constant transformation:\n{}", st);
            return make_pair(m.with_constant_t(invst, st), t_);
        } catch (const pole_error & e) {
            logd("got error: {}; retrying", e.what());
            continue;
        } catch (const std::runtime_error & e) {
            if (e.what() == std::string("matrix::inverse(): singular matrix")) {
                logd("got singular matrix; retrying");
                continue;
            }
            throw;
        }
    }
    assert(false);
}

bool
is_factorized(const pfmatrix &pfm, const symbol &eps)
{
    for (const auto &kv : pfm.residues) {
        const auto &ci = kv.second;
        if (ci.is_zero_matrix()) continue;
        exmap map = {{eps, 2*eps}};
        for (unsigned i = 0; i < ci.nops(); i++) {
            const ex &e = ci.op(i);
            if (!normal(e.subs(map) - e*2).is_zero()) return false;
        }
    }
    return true;
}

/* Block reduction
 * ============================================================
 */

pair<pfmatrix, transformation>
reduce_diagonal_blocks(const pfmatrix &m, const symbol &eps)
{
    LOGME;
    block_triangular_permutation btp(m);
    transformation t(m.nrows);
    t.add_constant_t(btp.t().transpose(), btp.t());
    // TODO: avoid normalize() call in with_constant_t()
    pfmatrix pfm = m.with_constant_t(btp.t().transpose(), btp.t());
    int offs = 0;
    for (int size : btp.block_size()) {
        logi("Reducing {}x{} diagonal block at offset {}", size, size, offs);
        pfmatrix b = pfm.block(offs, size);
        auto mt1 = fuchsify(b);
        assert(is_fuchsian(mt1.first));
        auto pfmx1 = t.apply(m);
        t.add(mt1.second.widen(pfm.nrows, offs));
        auto mt2 = normalize(mt1.first, eps);
        assert(is_normal(mt2.first, eps));
        t.add(mt2.second.widen(pfm.nrows, offs));
        auto mt3 = factorize(mt2.first, eps);
        assert(is_factorized(mt3.first, eps));
        t.add(mt3.second.widen(pfm.nrows, offs));
        offs += size;
    }
    logd("Applying the combined reduction to the whole system");
    pfm = t.apply(m);
    pfm.normalize();
    offs = 0;
    for (int size : btp.block_size()) {
        logd("Double-checking {}x{} block at offset {}", size, size, offs);
        pfmatrix b = pfm.block(offs, size);
        assert(is_fuchsian(b));
        assert(is_normal(b, eps));
        assert(is_factorized(b, eps));
        offs += size;
    }
    return make_pair(pfm, t);
}

pair<pfmatrix, transformation>
fuchsify_off_diagonal_blocks(const pfmatrix &m)
{
    LOGME;
    block_triangular_permutation btp(m);
    transformation t(m.nrows);
    t.add_constant_t(btp.t().transpose(), btp.t());
    // TODO: avoid normalize() call in with_constant_t()
    pfmatrix pfm = m.with_constant_t(btp.t().transpose(), btp.t());
    auto bs = btp.block_size();
    unsigned offs1 = 0;
    for (int i1 = 0; i1 < (int)bs.size(); i1++) {
        unsigned size1 = bs[i1];
        unsigned offs2 = offs1;
        for (int i2 = i1 - 1; i2 >= 0; i2--) {
            unsigned size2 = bs[i2];
            offs2 -= size2;
loop:;
            logd("Looking at {}x{} block at {}:{}", size1, size2, offs1, offs2);
            assert(offs1 > offs2);
            // Sort the residues.
            std::vector<pfmatrix::key> keys;
            keys.reserve(pfm.residues.size());
            for(auto &&kv: pfm.residues) {
                if (kv.first.second != -1) keys.push_back(kv.first);
            }
            sort(keys.begin(), keys.end(), [](auto &&k1, auto &&k2) {
                // Reduction of the residues with higher absolute
                // Poincare rank spoils the residues with lower
                // ones, so we should start with the largest.
                int p1 = k1.second + 1;
                int p2 = k2.second + 1;
                if (abs(p1) > abs(p2)) return true;
                if (abs(p1) < abs(p2)) return false;
                // The order beyond this point doesn't matter.
                if (p1 < p2) return true;
                if (p1 > p2) return false;
                return ex_is_less()(k1.first, k2.first);
            });
            for (auto &&kk : keys) {
                // It's useful that pi and ki do not point to
                // pfm, because we'll overwrite that, but we'll
                // still need pi and ki afterwards.
                const auto &pi = kk.first;
                const auto &ki = kk.second;
                const auto &ci = pfm(pi, ki);
                if (ci.is_zero_matrix()) continue;
                auto b = matrix_cut(ci, offs1, size1, offs2, size2);
                if (b.is_zero_matrix()) continue;
                logi("Reducing {}x{} block at {}:{}, at {}={}, k={}", size1, size2, offs1, offs2, pfm.x, pi, ki);
                lst d_vars;
                matrix d(size1, size2);
                for (unsigned i = 0; i < size1*size2; i++) {
                    symbol t;
                    d_vars.append(t);
                    d.let_op(i) = t;
                }
                matrix eqmx;
                if (ki < 0) {
                    // M = {A 0} = {a/x   0  }
                    //     {B C}   {b x^k c/x}
                    // D = {0         0}
                    //     {d x^(k+1) 0}
                    // M' = off_diagonal_t(M, D, p, k+1)
                    //    = M + x^(k+1) (M D - D M - 1/x (k+1) D)
                    // B' = B + x^k (c d - d a - (k+1) d)
                    auto a = matrix_cut(pfm(pi, -1), offs2, size2, offs2, size2);
                    auto c = matrix_cut(pfm(pi, -1), offs1, size1, offs1, size1);
                    eqmx = b.add(c.mul(d).sub(d.mul(a)).sub(d.mul_scalar(ki + 1)));
                } else {
                    assert(pi == 0);
                    // The case of ki >= 0 is special, because all pfm(*, -1)
                    // residues contribute to the transformed off-diagonal
                    // cell. This is because (x-p1)^0 = (x-0)^0.
                    //
                    // We could reproduce the correct calculations of
                    // the transformed cell here, but just reusing
                    // with_off_diagonal_t() code is equally valid, even if
                    // slower.
                    matrix D(m.nrows, m.ncols);
                    for (unsigned i = 0; i < size1; i++) {
                        for (unsigned j = 0; j < size2; j++) {
                            D(offs1 + i, offs2 + j) = d(i, j);
                        }
                    }
                    pfmatrix pfm2 = pfm.with_off_diagonal_t(D, pi, ki + 1);
                    eqmx = matrix_cut(pfm2(pi, ki), offs1, size1, offs2, size2);
                }
                lst eq;
                for (unsigned i = 0; i < eqmx.nops(); i++) {
                    eq.append(numer(eqmx.op(i)) == 0);
                }
                ex sol = lsolve(eq, d_vars);
                assert(sol.nops() == d_vars.nops());
                matrix D(m.nrows, m.ncols);
                for (unsigned i = 0; i < size1; i++) {
                    for (unsigned j = 0; j < size2; j++) {
                        D(offs1 + i, offs2 + j) = d(i, j).subs(sol);
                    }
                }
                logi("Use off-diagonal transformation, p={}, k={}, D=\n{}", pi, ki + 1, D);
                pfm = pfm.with_off_diagonal_t(D, pi, ki + 1);
                pfm.normalize();
                assert(matrix_cut(pfm(pi, ki), offs1, size1, offs2, size2).is_zero_matrix());
                t.add_off_diagonal_t(D, pi, ki + 1, pfm.x);
                goto loop;
            }
        }
        offs1 += size1;
    }
    return make_pair(pfm, t);
}

pair<pfmatrix, transformation>
reduce(const pfmatrix &m, const symbol &eps)
{
    LOGME;
    auto mt1 = reduce_diagonal_blocks(m, eps);
    auto mt2 = fuchsify_off_diagonal_blocks(mt1.first);
    auto mt3 = factorize(mt2.first, eps);
    transformation t(m.nrows);
    t.add(mt1.second);
    t.add(mt2.second);
    t.add(mt3.second);
    assert(is_fuchsian(mt3.first));
    assert(is_normal(mt3.first, eps));
    assert(is_factorized(mt3.first, eps));
    return make_pair(mt3.first, t);
}

/* Simplification
 * ============================================================
 */

int
nonzero_count(const pfmatrix &pfm)
{
    int nnz = 0;
    for (const auto &kv : pfm.residues) {
        const auto &ci = kv.second;
        for (unsigned i = 0; i < ci.nops(); i++) {
            nnz += !ci.op(i).is_zero();
        }
    }
    return nnz;
}

/* This simplification routine tries to increase the number of
 * zero terms in off-diagonal parts of the matrix with constant
 * transformations.
 *
 * Most useful after factorization, because factors of the form
 * eps/(1+eps) are skipped, so it's best if epsilon cancels
 * out, which is guaranteed in the epsilon form.
 */
pair<pfmatrix, transformation>
simplify_off_diagonal_blocks(const pfmatrix &m)
{
    LOGME;
    block_triangular_permutation btp(m);
    pfmatrix pfm = m.with_constant_t(btp.t().transpose(), btp.t());
    logi("Shuffle into block-diagonal form with:\n{}", btp.t());
    logd("Block sizes: {}", btp.block_size());
    transformation tr(m.nrows);
    tr.add_constant_t(btp.t().transpose(), btp.t());
    auto bs = btp.block_size();
    int nterms1 = nonzero_count(pfm);
    logd("Initial number of non-zero terms: {}", nterms1);
    for (;;) {
        bool done = true;
        // Traversing rows from the bottom to the top. This way
        // the bottom row collects the least amount of garbage
        // coefficients.
        for (int row = pfm.nrows, ir = bs.size() - 1; ir >= 0; ir--) {
            int rowend = row;
            row -= bs[ir];
            for (int r = rowend - 1; r >= row; r--) {
                for (int c = 0; c < row; c++) {
                    logd("Looking at {}:{}", r, c);
                    // Off-diagonal transformation 1+K*1(r,c), with r>c, has this effect:
                    // 1) m[rr, c] += K*m[rr, r]
                    // 2) m[r, cc] -= K*m[c, cc]
                    //
                    // Depending on K some terms become zero, others
                    // become (or remain) non-zero. Our strategy is to
                    // find K that results in the most zeros.
                    map<ex, int, ex_is_less> kzerocnt;
                    for (const auto &kv : pfm.residues) {
                        const auto &ci = kv.second;
                        // ci'(r, c) = ci(r, c) + K*(ci(r, r) - ci(c, c))
                        ex dc = normal(ci(r, r) - ci(c, c));
                        if (!dc.is_zero()) {
                            ex k = ratcan(-ci(r, c)/dc);
                            if (k.info(info_flags::rational)) kzerocnt[k] += 1;
                        }
                        for (int rr = row; rr < (int)pfm.nrows; rr++) {
                            if (rr == r) continue;
                            // ci'(rr, c) = ci(rr, c) + K*ci(rr, r);
                            if (ci(rr, r).is_zero()) continue;
                            ex k = ratcan(-ci(rr, c)/ci(rr, r));
                            if (k.info(info_flags::rational)) kzerocnt[k] += 1;
                        }
                        for (int cc = 0; cc < row; cc++) {
                            if (cc == c) continue;
                            // ci'(r, cc) = ci(r, cc) - K*ci(c, cc);
                            if (ci(c, cc).is_zero()) continue;
                            ex k = ratcan(ci(r, cc)/ci(c, cc));
                            if (k.info(info_flags::rational)) kzerocnt[k] += 1;
                        }
                    }
                    for (auto &&kv : kzerocnt) {
                        logd("* k={} removes {} terms", kv.first, kv.second);
                    }
                    if (kzerocnt.size() == 0) continue;
                    ex k = max_element(begin(kzerocnt), end(kzerocnt), [] (auto &&kv1, auto &&kv2) {
                        if (kv1.second < kv2.second) return true;
                        if (kv1.second > kv2.second) return false;
                        int c1 = complexity(kv1.first);
                        int c2 = complexity(kv2.first);
                        if (c1 > c2) return true;
                        if (c1 < c2) return false;
                        return ex_is_less()(kv1.first, kv2.first);
                    })->first;
                    if (kzerocnt[k] > kzerocnt[0]) {
                        logd("* best k here is: {} (should eliminate {} terms)", k, kzerocnt[k] - kzerocnt[0]);
                        matrix t = identity_matrix(pfm.nrows);
                        matrix invt = identity_matrix(pfm.nrows);
                        t(r, c) = k;
                        invt(r, c) = -k;
                        logi("Use constant transformation:\n{}", t);
                        int n1 = nonzero_count(pfm);
                        pfm = pfm.with_constant_t(invt, t);
                        tr.add_constant_t(invt, t);
                        int n2 = nonzero_count(pfm);
                        logd("Term reduction: {} ({} vs {})", n1 - n2, n2, n1);
                        done = false;
                    }
                }
            }
        }
        if (done) break;
        logd("We'll now make an additional pass");
    }
    int nterms2 = nonzero_count(pfm);
    logd("Total terms eliminated: {} ({} vs {})", nterms1 - nterms2, nterms2, nterms1);
    return make_pair(pfm, tr);
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
