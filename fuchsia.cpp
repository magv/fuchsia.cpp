#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <fstream>
#include <tuple>
#include <assert.h>

using namespace GiNaC;
using namespace std;

matrix
ex_to_matrix(const ex &m)
{
    return ex_to<matrix>(m.evalm());
}

/* Infinity object.
 */
struct infinity_s { infinity_s() {} };
typedef structure<infinity_s> infinity_t;

template<> void
infinity_t::print(const print_context & c, unsigned level) const
{
    if (is_a<print_tree>(c))
        inherited::print(c, level);
    c.s << "Infinity";
}

const ex infinity = infinity_t();

/* Load matrix from a file in Mathematica format.
 */
pair<matrix, symtab>
load_matrix(const char *filename, const symtab &table)
{
    parser reader(table);
    ifstream i(filename);
    ex x = reader(i);
    i.close();
    // TODO: check that the input is indeed a matrix of rational
    // expressions, with no imaginary or irrational numbers; signal
    // errors if this is not the case.
    matrix m = ex_to<matrix>(lst_to_matrix(ex_to<lst>(x)));
    return make_pair(m, reader.get_syms());
}

/* Save matrix to a file in Mathematica format.
 */
void
save_matrix(const char *filename, const matrix &m)
{
    ofstream f(filename);
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
    f.close();
}

/* Find the set of right eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_right(matrix m, ex eigenvalue)
{
    unsigned n = m.cols();
    // Construct eigenvector of symbol temporaries.
    std::vector<symbol> tmp;
    matrix ev(n, 1);
    for (unsigned i = 0; i < n; i++) {
        symbol t;
        tmp.push_back(t);
        ev(i, 0) = t;
    }
    // Construct mm = m - eigenvalue*I
    matrix mm(n, n);
    for (unsigned i = 0; i < n; i++) {
        mm(i, i) = eigenvalue;
    }
    mm = m.sub(mm);
    // Solve mm*v = 0
    matrix rhs(n, 1);
    matrix s = mm.solve(ev, rhs);
    // Find eigenvectors
    vector<matrix> evectors;
    for (unsigned k = 0; k < n; k++) {
        matrix ev(n, 1);
        for (unsigned i = 0; i < n; i++) {
            ev(i, 0) = normal(s(i, 0).coeff(tmp[k]));
        }
        // is_zero_matrix() is only sufficient here because
        // of the normal() call above.
        if (not ev.is_zero_matrix()) {
            evectors.push_back(ev);
        }
    }
    return evectors;
}

/* Find the set of left eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_left(matrix m, ex eigenvalue)
{
    auto v = eigenvectors_right(m.transpose(), eigenvalue);
    for (auto &e : v) e = e.transpose();
    return v;
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
    ex nd = numer_denom(e);
    return nd.op(0).expand()/nd.op(1).expand();
}

/* Canonical form for matrices of rational expressions.
 */
matrix
ratcan(const matrix &m)
{
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.rows(); i++) {
        for (unsigned j = 0; j < m.cols(); j++) {
            r(i, j) = ratcan(m(i, j));
        }
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

/* Iterate through terms of a partial fraction decomposition
 * of a rational expression in x. Call yield(p, q, k) for
 * each p*q^k term.
 */
template<typename F> void
partial_fraction_iter(const ex &e, const symbol &x, F yield)
{
    ex nd = numer_denom(e);
    ex numer = nd.op(0);
    ex denom = nd.op(1);
    auto qr = poly_divmod(numer, denom, x);
    ex q = qr.first;
    ex r = qr.second;
    denom = factor(denom);
    struct Term { ex f; int n; };
    vector<Term> factors;
    factor_iter(denom,
        [&](const ex &f, int i) {
            factors.push_back(Term{f, i});
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
            eq += expand(pfnum*(denom/pow(f.f, k)));
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
    int idx = 0;
    for (const auto &f : factors) {
        int deg = f.f.degree(x);
        // For each power of factor...
        for (int k = 1; k <= f.n; k++) {
            ex pfnum = 0;
            // For each power of the parfrac term numerator...
            for (int l = 0; l < deg; l++) {
                pfnum += sol.op(idx).rhs()*pow(x, l);
                idx++;
            }
            if (!pfnum.is_zero()) yield(pfnum, f.f, -k);
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

    pfmatrix(unsigned nrows, unsigned ncols);
    pfmatrix(const matrix &m, const symbol &x);
    matrix &operator ()(const ex &p, int k);
    matrix to_matrix(const ex &x);
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
    pfmatrix with_balance_t(const matrix &P, const ex &x1, const ex &p2) const;
};

bool
pfmatrix::key_is_less::operator ()(const key &k1, const key &k2) const
{
    return (k1.second < k2.second) || ex_is_less()(k1.first, k2.first);
}

pfmatrix::pfmatrix(unsigned nrows, unsigned ncols)
    : nrows(nrows), ncols(ncols)
{
}

pfmatrix::pfmatrix(const matrix &m, const symbol &x)
    : nrows(m.rows()), ncols(m.cols())
{
    for (unsigned i = 0; i < nrows; i++) {
        for (unsigned j = 0; j < ncols; j++) {
            partial_fraction_iter(m(i, j), x,
                [&](const auto &p, const auto &q, int k) {
                    int deg = q.degree(x);
                    if (k >= 0) {
                        // q == x
                        (*this)(0, k)(i, j) = p;
                    }
                    else if (deg == 1) {
                        // q == c0 + c1*x
                        ex c0 = q.coeff(x, 0);
                        ex c1 = q.coeff(x, 1);
                        (*this)(normal(-c0/c1), k)(i, j) = p*pow(c1, k);
                    }
                    else {
                        throw runtime_error("pfmatrix(): can't solve equations of 2nd degree or higher");
                    }
                });
        }
    }
}

matrix &
pfmatrix::operator ()(const ex &p, int k)
{
    auto key = make_pair(p, k);
    decltype(residues)::iterator it;
    tie(it, ignore) = residues.emplace(
            piecewise_construct,
            make_tuple(key),
            std::make_tuple(nrows, ncols));
    return (*it).second;
}

matrix
pfmatrix::to_matrix(const ex &x)
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
pfmatrix::add(const matrix &C, const ex &p1, int k1)
{
    // M += C*(x-p1)^k1
    matrix &R = (*this)(p1, k1);
    R = R.add(C);
}

void
pfmatrix::add_div(const matrix &C, const ex &p1, int k1, const ex &p2)
{
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
    pfmatrix m(nrows, ncols);
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
    pfmatrix m(nrows, ncols);
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
    // TODO: normalize and clean residues here maybe?
    return m;
}

matrix
balance_t(const matrix &m, const matrix &P, const ex &x1, const ex &x2, const ex &x)
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
vector<ex>
eigenvalues(const matrix &m)
{
    vector<ex> eigenvalues;
    symbol lambda("λ");
    ex charpoly = charpoly_by_blocks(m, lambda);
    factor_iter(factor(charpoly),
        [&](const ex &f, int k) {
            int deg = f.degree(lambda);
            if (deg == 0) {
                // No roots in λ here.
            }
            else if (deg == 1) {
                // f == c0 + c1*λ
                ex c0 = f.coeff(lambda, 0);
                ex c1 = f.coeff(lambda, 1);
                for (int i = 0; i < k; i++)
                    eigenvalues.push_back(normal(-c0/c1));
            }
            else {
                throw runtime_error("eigenvalues(): can't solve equations of 2nd degree or higher");
            }
        });
    return eigenvalues;
}
