#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <chrono>
#include <numeric>
#include <algorithm>

using namespace std;

// Hamming(15, 11) 基础参数
const int BASE_N = 15;
const int BASE_K = 11;
const int BASE_M = 4; // 校验位 N - K

/* =========================================================
   1. 基础数学工具 (GF(2) 矩阵运算)
   ========================================================= */

using Matrix = vector<vector<int>>;

// 矩阵乘法 A * B
Matrix mat_mul(const Matrix& A, const Matrix& B) {
    int r = A.size();
    int c = B[0].size();
    int k_dim = B.size();
    Matrix C(r, vector<int>(c, 0));
    for (int i = 0; i < r; i++) {
        for (int k = 0; k < k_dim; k++) {
            if (A[i][k] == 1) {
                for (int j = 0; j < c; j++) {
                    C[i][j] ^= B[k][j];
                }
            }
        }
    }
    return C;
}

// 向量 * 矩阵
vector<int> vec_mul(const vector<int>& v, const Matrix& M) {
    int n = v.size();
    int m = M[0].size();
    vector<int> res(m, 0);
    for (int i = 0; i < n; i++) {
        if (v[i] == 1) {
            for (int j = 0; j < m; j++) {
                res[j] ^= M[i][j];
            }
        }
    }
    return res;
}

// 矩阵转置
Matrix mat_transpose(const Matrix& M) {
    int r = M.size();
    int c = M[0].size();
    Matrix T(c, vector<int>(r));
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            T[j][i] = M[i][j];
    return T;
}

// 生成单位矩阵
Matrix identity(int n) {
    Matrix I(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) I[i][i] = 1;
    return I;
}

// 高斯消元求逆矩阵
Matrix mat_inv(Matrix A) {
    int n = A.size();
    if (A.empty() || A[0].size() != n) return {};
    Matrix I = identity(n);
    for (int i = 0; i < n; i++) {
        int pivot = -1;
        for (int j = i; j < n; j++) {
            if (A[j][i] == 1) { pivot = j; break; }
        }
        if (pivot == -1) return {};
        if (pivot != i) {
            swap(A[i], A[pivot]);
            swap(I[i], I[pivot]);
        }
        for (int j = 0; j < n; j++) {
            if (j != i && A[j][i] == 1) {
                for (int k = 0; k < n; k++) {
                    A[j][k] ^= A[i][k];
                    I[j][k] ^= I[i][k];
                }
            }
        }
    }
    return I;
}

/* =========================================================
   2. Hamming(15,11) 核心组件
   ========================================================= */

void generate_base_hamming_matrices(Matrix& G_sys, Matrix& H_sys) {
    Matrix H_raw(BASE_M, vector<int>(BASE_N));
    for (int j = 0; j < BASE_N; j++) {
        int val = j + 1;
        for (int i = 0; i < BASE_M; i++) {
            H_raw[i][j] = (val >> i) & 1;
        }
    }

    Matrix tempH = H_raw;
    int n = BASE_N;
    int m = BASE_M;
    int k = BASE_K;

    for (int i = 0; i < m; i++) {
        int pivot_col = k + i;
        int pivot_row = -1;
        for (int r = i; r < m; r++) {
            if (tempH[r][pivot_col] == 1) { pivot_row = r; break; }
        }
        if (pivot_row == -1) {
            for (int c = 0; c < k; c++) {
                int r2 = -1;
                for (int rr = i; rr < m; rr++) if (tempH[rr][c] == 1) { r2 = rr; break; }
                if (r2 != -1) {
                    for (int row = 0; row < m; row++) swap(tempH[row][c], tempH[row][pivot_col]);
                    pivot_row = r2;
                    break;
                }
            }
        }
        if (pivot_row != i) swap(tempH[i], tempH[pivot_row]);
        for (int r = 0; r < m; r++) {
            if (r != i && tempH[r][pivot_col] == 1) {
                for (int c = 0; c < n; c++) tempH[r][c] ^= tempH[i][c];
            }
        }
    }
    H_sys = tempH;

    Matrix PT(m, vector<int>(k));
    for (int r = 0; r < m; r++)
        for (int c = 0; c < k; c++)
            PT[r][c] = H_sys[r][c];

    Matrix P_mat = mat_transpose(PT);
    G_sys.assign(k, vector<int>(n, 0));
    for (int i = 0; i < k; i++) G_sys[i][i] = 1;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < m; j++)
            G_sys[i][k + j] = P_mat[i][j];
}

Matrix create_block_diagonal(const Matrix& small_M, int L) {
    int r = small_M.size();
    int c = small_M[0].size();
    Matrix BigM(r * L, vector<int>(c * L, 0));
    for (int k = 0; k < L; k++) {
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                BigM[k * r + i][k * c + j] = small_M[i][j];
            }
        }
    }
    return BigM;
}

Matrix generate_invertible_S(int size) {
    Matrix S, S_inv;
    do {
        S = identity(size);
        for (int i = 0; i < size * 2; i++) {
            int r1 = rand() % size;
            int r2 = rand() % size;
            if (r1 != r2) for (int j = 0; j < size; j++) S[r1][j] ^= S[r2][j];
        }
        S_inv = mat_inv(S);
    } while (S_inv.empty());
    return S;
}

Matrix generate_permutation_P(int size) {
    Matrix P(size, vector<int>(size, 0));
    vector<int> perm(size);
    for (int i = 0; i < size; i++) perm[i] = i;
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap(perm[i], perm[j]);
    }
    for (int i = 0; i < size; i++) P[i][perm[i]] = 1;
    return P;
}

vector<int> decode_block(vector<int> block, const Matrix& H) {
    Matrix HT = mat_transpose(H);
    vector<int> s = vec_mul(block, HT);
    bool error = false;
    for (int bit : s) if (bit) error = true;
    if (!error) return block;

    int n = H[0].size();
    int m = H.size();
    int error_pos = -1;
    for (int j = 0; j < n; j++) {
        bool match = true;
        for (int i = 0; i < m; i++) {
            if (H[i][j] != s[i]) { match = false; break; }
        }
        if (match) { error_pos = j; break; }
    }
    if (error_pos != -1) block[error_pos] ^= 1;
    return block;
}

/* =========================================================
   3. 综合测试逻辑 (遵循“解密出错的核心”推导)
   ========================================================= */

void run_test(int L) {
    cout << "\n========================================================" << endl;
    cout << "  运行测试: 分块数量 L = " << L << endl;
    cout << "  方案: 每分块 50% 概率产生随机错误 (全局混淆模式)" << endl;
    cout << "========================================================" << endl;

    int total_K = L * BASE_K;
    int total_N = L * BASE_N;

    Matrix G_base, H_base;
    generate_base_hamming_matrices(G_base, H_base);

    Matrix G_total = create_block_diagonal(G_base, L);
    Matrix S = generate_invertible_S(total_K);
    Matrix P = generate_permutation_P(total_N);

    Matrix G_pub = mat_mul(mat_mul(S, G_total), P);
    Matrix S_inv = mat_inv(S);
    Matrix P_inv = mat_transpose(P);

    int rounds = 1000, success_count = 0;
    double encrypt_time_total = 0, decrypt_time_total = 0;

    for (int r = 0; r < rounds; r++) {
        vector<int> m(total_K);
        for (int i = 0; i < total_K; i++) m[i] = rand() % 2;

        auto start = chrono::high_resolution_clock::now();
        vector<int> c = vec_mul(m, G_pub);

        // --- 核心修改：遵循您的方案注入错误 ---
        // 1. 遍历 L 个分块，每个分块 50% 概率产生 1 个错
        vector<int> e(total_N, 0);
        for (int i = 0; i < L; i++) {
            if ((rand() % 100) < 50) {
                int bit_pos = i * BASE_N + (rand() % BASE_N);
                e[bit_pos] = 1;
            }
        }
        // 此时 e 是直接在加密后空间叠加的随机噪声，不与 P 耦合
        for (int i = 0; i < total_N; i++) c[i] ^= e[i];

        auto end = chrono::high_resolution_clock::now();
        encrypt_time_total += (chrono::duration<double, milli>(end - start)).count();

        // --- 解密过程 ---
        start = chrono::high_resolution_clock::now();

        // 逆置换：c' = c * P^-1
        // 这里会发生 e * P^-1，导致错误分布被打乱，可能在某些块聚集
        vector<int> c_prime = vec_mul(c, P_inv);

        vector<int> m_prime_total;
        for (int k = 0; k < L; k++) {
            vector<int> block(BASE_N);
            for (int i = 0; i < BASE_N; i++) block[i] = c_prime[k * BASE_N + i];
            vector<int> corrected_block = decode_block(block, H_base);
            for (int i = 0; i < BASE_K; i++) m_prime_total.push_back(corrected_block[i]);
        }

        vector<int> m_decrypted = vec_mul(m_prime_total, S_inv);
        end = chrono::high_resolution_clock::now();
        decrypt_time_total += (chrono::duration<double, milli>(end - start)).count();

        if (m == m_decrypted) success_count++;
    }

    cout << "[3] 测试结果:" << endl;
    cout << "    解密成功率: " << (double)success_count / rounds * 100.0 << "%" << endl;
    cout << "    (注：成功率 < 100% 证明了 P 的全局随机性导致 e*P^-1 错误超限)" << endl;
    cout << "    平均加密时间: " << encrypt_time_total / rounds << " ms" << endl;
    cout << "    平均解密时间: " << decrypt_time_total / rounds << " ms" << endl;
    cout << "    密文扩张率: " << fixed << setprecision(2) << (double)total_N / total_K << "x" << endl;
}

int main() {
    srand(time(NULL));
    cout << "=== Problem A: McEliece (Hamming 15,11) ===\n";
    run_test(10);
    run_test(20);
    return 0;
}