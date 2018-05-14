#include <unistd.h>
#include "fuchsia.cpp"

static const char usagetext[] = R"(
Ss{NAME}
    Nm{fuchsia} -- transform linear differential equations into epsilon form.

Ss{SYNOPSYS}
    Nm{fuchsia} [options] Cm{command} Ar{args} ...

Ss{COMMANDS}
    Cm{show} [Fl{-x} Ar{name}] Ar{matrix}
        Show a description of a given matrix.

    Cm{reduce} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] Ar{matrix}
        Find an epsilon form of the given matrix.

    Cm{fuchsify} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] Ar{matrix}
        Find a transformation that will transform a given matrix into Fuchsian
        form.

    Cm{normalize} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] Ar{matrix}
        Find a transformation that will transform a given Fuchsian matrix into
        normalized form.

    Cm{factorize} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] Ar{matrix}
        Find a transformation that will make a given normalized matrix
        proportional to the infinitesimal parameter.

    Cm{sort} [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] Ar{matrix}
        Find a block-triangular form of the given matrix by shuffling.

    Cm{transform} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] Ar{matrix} Ar{transform} ...
        Transform a given matrix using a given transformation.

    Cm{changevar} [Fl{-x} Ar{name}] [Fl{-y} Ar{name}] [Fl{-m} Ar{path}] Ar{matrix} Ar{expr}
        Perform a change of variable from x to y, such that x=expr(y).

Ss{OPTIONS}
    Fl{-h}         Show this help message.
    Fl{-v}         Print a more verbose log.
    Fl{-x} Ar{name}    Use this name for the free variable (default: x).
    Fl{-y} Ar{name}    Use this name for the new free variable (default: y).
    Fl{-e} Ar{name}    Use this name for the infinitesimal parameter (default: eps).
    Fl{-m} Ar{path}    Save the resulting matrix into this file.
    Fl{-t} Ar{path}    Save the resulting transformation into this file.

Ss{ARGUMENTS}
    Ar{matrix}     Read the input matrix from this file.
    Ar{transform}  Read the transformation matrix from this file.
    Ar{expr}       Arbitrary expression.

Ss{AUTHORS}
    Vitaly Magerya <vitaly.magerya@tx97.net>
)";

static bool COLORS = !!isatty(STDOUT_FILENO);

void
usage()
{
    const char *p = strchr(usagetext, '\n') + 1;
    for (;;) {
        const char *l = strchr(p + 2, '{');
        if (l == NULL) break;
        const char *r = strchr(l, '}');
        if (r == NULL) break;
        const char *a = "", *b = "\033[0m";
        if (l[-2] == 'S' && l[-1] == 's') { a = "\033[1m"; }
        if (l[-2] == 'N' && l[-1] == 'm') { a = "\033[1;35m"; }
        if (l[-2] == 'F' && l[-1] == 'l') { a = "\033[33m"; }
        if (l[-2] == 'C' && l[-1] == 'm') { a = "\033[1m"; }
        if (l[-2] == 'A' && l[-1] == 'r') { a = "\033[32m"; }
        cout.write(p, l - p - 2);
        if (COLORS) cout << a;
        cout.write(l + 1, r - l - 1);
        if (COLORS) cout << b;
        p = r + 1;
    }
    cout << p;
}

int
main(int argc, char *argv[])
{
    const char *var_x_name = "x";
    const char *var_y_name = "y";
    const char *var_eps_name = "eps";
    const char *matrix_m_path = NULL;
    const char *matrix_t_path = NULL;
    const char *matrix_i_path = NULL;
    for (int opt; (opt = getopt(argc, argv, "hvx:e:y:m:t:i:s:")) != -1;) {
        switch (opt) {
        case 'h': usage(); return 0;
        case 'v': log_verbose = true; break;
        case 'x': var_x_name = optarg; break;
        case 'y': var_y_name = optarg; break;
        case 'e': var_eps_name = optarg; break;
        case 'm': matrix_m_path = optarg; break;
        case 't': matrix_t_path = optarg; break;
        case 'i': matrix_i_path = optarg; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    matrix matrix_m(0, 0);
    matrix matrix_t(0, 0);
    matrix matrix_i(0, 0);
    symbol x(var_x_name);
    symbol y(var_y_name);
    symbol eps(var_eps_name);
    symtab vars = {{var_x_name, x}, {var_y_name, y}, {var_eps_name, eps}};
    if ((argc == 1) && !strcmp(argv[0], "help")) {
        usage();
        return 0;
    }
    else if ((argc == 2) && !strcmp(argv[0], "show")) {
        auto ms = load_matrix(argv[1], vars);
        cout << "Matrix size: " << ms.first.rows() << "x" << ms.first.cols() << endl;
        pfmatrix pfm(ms.first, x);
        cout << "Matrix complexity: " << complexity(pfm) << endl;
        cout << "Matrix expansion:" << endl;
        for (const auto kv : pfm.residues) {
            const auto &xi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &c = kv.second;
            if (c.is_zero_matrix()) continue;
            cout << "  pole of power " << ki << " at " << x << "=" << (ki >= 0 ? infinity : xi) << endl;
            cout << "    complexity: " << complexity(c) << endl;
            for (auto ev : eigenvalues(c, true)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
        }
        matrix c = normal(c0_infinity(pfm));
        if (!c.is_zero_matrix()) {
            cout << "  effective pole of power " << -1 << " at " << x << "=" << infinity << endl;
            cout << "    complexity: " << complexity(c) << endl;
            for (auto ev : eigenvalues(c, true)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
        }
    }
    else if ((argc == 2) && !strcmp(argv[0], "sort")) {
        auto ms = load_matrix(argv[1], vars);
        block_triangular_permutation btp(ms.first);
        matrix_m = btp.t().transpose().mul(ms.first).mul(btp.t());
        matrix_t = btp.t();
        matrix_i = btp.t().transpose();
    }
    else if ((argc == 2) && !strcmp(argv[0], "fuchsify")) {
        auto ms = load_matrix(argv[1], vars);
        pfmatrix pfm(ms.first, x);
        auto r = fuchsify(pfm);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else if ((argc == 2) && !strcmp(argv[0], "normalize")) {
        auto ms = load_matrix(argv[1], vars);
        pfmatrix pfm(ms.first, x);
        auto r = normalize(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else if ((argc == 2) && !strcmp(argv[0], "factorize")) {
        auto ms = load_matrix(argv[1], vars);
        pfmatrix pfm(ms.first, x);
        auto r = factorize(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else if ((argc == 2) && !strcmp(argv[0], "reduce")) {
        auto ms = load_matrix(argv[1], vars);
        pfmatrix pfm(ms.first, x);
        auto r = reduce(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else if ((argc >= 3) && !strcmp(argv[0], "transform")) {
        matrix m;
        tie(m, vars) = load_matrix(argv[1], vars);
        for (int i = 2; i < argc; i++) {
            matrix t;
            tie(t, vars) = load_matrix(argv[i], vars);
            m = matrix_inverse(t).mul(m.mul(t).sub(ex_to_matrix(t.diff(x))));
        }
        pfmatrix pfm(m, x);
        matrix_m = pfm.to_matrix();
    }
    else if ((argc == 3) && !strcmp(argv[0], "changevar")) {
        auto ms = load_matrix(argv[1], vars);
        parser reader(ms.second);
        auto xsubs = reader(argv[2]);
        matrix m2 = ex_to_matrix(ms.first.subs(exmap{{x, xsubs}})).mul_scalar(xsubs.diff(y));
        matrix_m = pfmatrix(m2, y).to_matrix();
    }
    else if (argc == 0) {
        cerr << "fuchsia: no command provided (use -h to see usage)" << endl;
        return 1;
    }
    else {
        cerr << "fuchsia: unrecognized command '"
             << argv[0]
             << "' (use -h to see usage)" << endl;
        return 1;
    }
    if (matrix_m.nops() > 0) {
        if (matrix_m_path != NULL) {
            logi("Saving the matrix to {}", matrix_m_path);
            save_matrix(matrix_m_path, matrix_m);
        } else {
            save_matrix(cout, matrix_m);
            cout << endl;
        }
    }
    if (matrix_t.nops() > 0) {
        if (matrix_t_path != NULL) {
            logi("Saving the transformation to {}", matrix_t_path);
            save_matrix(matrix_t_path, matrix_t);
        } else {
            logi("Not saving the transformation matrix (no -t argument)");
        }
    }
    if (matrix_i.nops() > 0) {
        if (matrix_i_path != NULL) {
            logi("Saving the inverse transformation to {}", matrix_i_path);
            save_matrix(matrix_i_path, matrix_i);
        } else {
            logi("Not saving the inverse transformation matrix (no -i argument)");
        }
    }
    return 0;
}
