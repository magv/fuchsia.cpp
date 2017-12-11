#include <unistd.h>
#include "fuchsia.cpp"

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

void
usage()
{
    cerr <<
        "Usage:\n"
        "    fuchsia [options] <command> <args>...\n"
        "\n"
        "Commands:\n"
        "    show [-x <name>] <matrix>\n"
        "        show a description of a given matrix\n"
        "\n"
        "    sort [-m <path>] [-t <path>] <matrix>\n"
        "        find a block-triangular form of the given matrix\n"
        "\n"
        "    transform [-x <name>] [-m <path>] <matrix> <transform>\n"
        "        transform a given matrix using a given transformation\n"
        "\n"
        "Options:\n"
        "    -h          show this help message\n"
        "    -v          produce a more verbose log\n"
        "    -x <name>   use this name for the free variable (default: x)\n"
        "    -e <name>   use this name for the infinitesimal parameter (default: eps)\n"
        "    -m <path>   save the resulting matrix into this file\n"
        "    -t <path>   save the resulting transformation into this file\n"
        "\n"
        "Arguments:\n"
        "    <matrix>    read the input matrix from this file\n"
        "    <transform> read the transformation matrix from this file\n"
        "    <expr>      arbitrary expression\n";
}

int
main(int argc, char *argv[])
{
    const char *var_x_name = "x";
    const char *var_eps_name = "eps";
    const char *matrix_m_path = NULL;
    const char *matrix_t_path = NULL;
    for (int opt; (opt = getopt(argc, argv, "hvx:e:m:t:s:")) != -1;) {
        switch (opt) {
        case 'h': usage(); return 0;
        case 'v': log_verbose = true; break;
        case 'x': var_x_name = optarg; break;
        case 'e': var_eps_name = optarg; break;
        case 'm': matrix_m_path = optarg; break;
        case 't': matrix_t_path = optarg; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    matrix matrix_m(0, 0);
    matrix matrix_t(0, 0);
    if ((argc == 1) && !strcmp(argv[0], "help")) {
        usage();
        return 0;
    }
    else if ((argc == 2) && !strcmp(argv[0], "show")) {
        symtab s;
        matrix m;
        tie(m, s) = load_matrix(argv[1], s);
        symbol x = ex_to<symbol>(s[var_x_name]);
        cout << "Matrix size: " << m.rows() << "x" << m.cols() << endl;
        cout << "Matrix expansion:" << endl;
        pfmatrix pfm(m, x);
        for (const auto kv : pfm.residues) {
            const auto &xi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &c = kv.second;
            if (c.is_zero_matrix()) continue;
            cout << "  pole of power " << ki << " at " << x << "=" << (ki >= 0 ? infinity : xi) << endl;
            for (auto ev : eigenvalues(c)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
        }
        matrix c = normal(c0_infinity(pfm));
        if (!c.is_zero_matrix()) {
            cout << "  effective pole of power " << -1 << " at " << x << "=" << infinity << endl;
            for (auto ev : eigenvalues(c)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
        }
    }
    else if ((argc == 2) && !strcmp(argv[0], "sort")) {
        symtab s;
        matrix m;
        tie(m, s) = load_matrix(argv[1], s);
        block_triangular_permutation btp(m);
        matrix_m = btp.t().transpose().mul(m).mul(btp.t());
        matrix_t = btp.t();
    }
    else if ((argc == 3) && !strcmp(argv[0], "transform")) {
        symtab s;
        matrix m, t;
        tie(m, s) = load_matrix(argv[1], s);
        tie(t, s) = load_matrix(argv[2], s);
        symbol x = ex_to<symbol>(s[var_x_name]);
        pfmatrix pfm(t.inverse().mul(m.mul(t).sub(ex_to_matrix(t.diff(x)))), x);
        matrix_m = pfm.to_matrix();
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
            save_matrix(matrix_m_path, matrix_m);
        } else {
            save_matrix(cout, matrix_m);
            cout << endl;
        }
    }
    if (matrix_t.nops() > 0) {
        if (matrix_t_path != NULL) {
            save_matrix(matrix_t_path, matrix_t);
        } else {
            logi("not saving the transformation matrix (no -t argument)");
        }
    }
    return 0;
}
