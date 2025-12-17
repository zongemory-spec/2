#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <numeric>

using namespace std;

// ==========================================
// 1. 全局参数配置 (BCH 15, 7, t=2)
// ==========================================
const int BLOCK_N = 15;  // 码字长度 n
const int BLOCK_K = 7;   // 信息位长度 k
const int BLOCK_R = 8;   // 校验位长度 n-k
const int BLOCK_T = 2;   // 纠错能力 t

// GF(2^4) 表
vector<int> gf_log(16);
vector<int> gf_exp(32);
const int PRIMITIVE_POLY = 19; // x^4 + x + 1 = 10011(bin)

// ==========================================
// 2. 有限域 GF(2^4) 运算
// ==========================================

void init_gf_tables() {
    int x = 1;
    for (int i = 0; i < 15; i++) {
        gf_exp[i] = x;
        gf_log[x] = i;
        x <<= 1;
        if (x & 16) {
            x ^= PRIMITIVE_POLY;
        }
    }
    for (int i = 15; i < 32; i++) {
        gf_exp[i] = gf_exp[i % 15];
    }
    gf_log[0] = -1;
}

int gf_mult(int a, int b) {
    if (a == 0 || b == 0) return 0;
    return gf_exp[(gf_log[a] + gf_log[b]) % 15];
}

int gf_div(int a, int b) {
    if (a == 0) return 0;
    if (b == 0) return 0;
    int res = gf_log[a] - gf_log[b];
    if (res < 0) res += 15;
    return gf_exp[res];
}

int gf_pow(int a, int n) {
    if (a == 0) return 0;
    if (n == 0) return 1;
    int log_val = (gf_log[a] * n) % 15;
    if (log_val < 0) log_val += 15;
    return gf_exp[log_val];
}

// ==========================================
// 3. 矩阵运算工具 (GF(2) 域)
// ==========================================

typedef vector<vector<int>> Matrix;

// 向量 * 矩阵
vector<int> vec_mat_mul(const vector<int>& v, const Matrix& M) {
    int rows = M.size();
    int cols = M[0].size();
    vector<int> res(cols, 0);
    for (int i = 0; i < rows; i++) {
        if (v[i] == 1) {
            for (int j = 0; j < cols; j++) res[j] ^= M[i][j];
        }
    }
    return res;
}

// 矩阵 * 矩阵
Matrix mat_mul(const Matrix& A, const Matrix& B) {
    int rA = A.size();
    int cA = A[0].size();
    int cB = B[0].size();
    Matrix C(rA, vector<int>(cB, 0));
    for (int i = 0; i < rA; i++) {
        for (int k = 0; k < cA; k++) {
            if (A[i][k] == 1) {
                for (int j = 0; j < cB; j++) C[i][j] ^= B[k][j];
            }
        }
    }
    return C;
}

// 矩阵转置
Matrix mat_transpose(const Matrix& M) {
    int r = M.size();
    int c = M[0].size();
    Matrix T(c, vector<int>(r));
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) T[j][i] = M[i][j];
    }
    return T;
}

// 矩阵求逆 (高斯消元)
Matrix mat_inverse(Matrix A) {
    int n = A.size();
    Matrix I(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) I[i][i] = 1;

    for (int i = 0; i < n; i++) {
        int pivot = i;
        while (pivot < n && A[pivot][i] == 0) pivot++;
        if (pivot == n) return {};
        swap(A[i], A[pivot]);
        swap(I[i], I[pivot]);

        for (int k = 0; k < n; k++) {
            if (k != i && A[k][i] == 1) {
                for (int j = 0; j < n; j++) {
                    A[k][j] ^= A[i][j];
                    I[k][j] ^= I[i][j];
                }
            }
        }
    }
    return I;
}

// 生成随机可逆矩阵 S
Matrix gen_S(int size) {
    while (true) {
        Matrix S(size, vector<int>(size));
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) S[i][j] = rand() % 2;
        if (!mat_inverse(S).empty()) return S;
    }
}

// 修改后的全局置换矩阵 P 生成函数
Matrix gen_P(int total_N, int block_size) {
    Matrix P(total_N, vector<int>(total_N, 0));
    vector<int> p(total_N);
    iota(p.begin(), p.end(), 0); // 生成 0 到 N-1 的序列
    random_shuffle(p.begin(), p.end()); // 全局随机打乱
    for (int i = 0; i < total_N; i++) {
        P[i][p[i]] = 1;
    }
    return P;
}

// ==========================================
// 4. BCH 核心
// ==========================================

vector<int> encode_bch_block(const vector<int>& msg) {
    vector<int> temp(15, 0);
    for (int i = 0; i < 7; i++) temp[i] = msg[i];
    int g[] = { 1, 1, 1, 0, 1, 0, 0, 0, 1 };
    for (int i = 0; i < 7; i++) {
        if (temp[i] == 1) {
            for (int j = 0; j < 9; j++) temp[i + j] ^= g[j];
        }
    }
    vector<int> codeword = msg;
    for (int i = 7; i < 15; i++) codeword.push_back(temp[i]);
    return codeword;
}

Matrix gen_H_sys_block() {
    Matrix P_part(7, vector<int>(8));
    for (int i = 0; i < 7; i++) {
        vector<int> u(7, 0); u[i] = 1;
        vector<int> c = encode_bch_block(u);
        for (int j = 0; j < 8; j++) P_part[i][j] = c[7 + j];
    }
    Matrix PT = mat_transpose(P_part);
    Matrix H(8, vector<int>(15));
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 7; j++) H[i][j] = PT[i][j];
        H[i][7 + i] = 1;
    }
    return H;
}

Matrix GLOBAL_H_BLOCK;

vector<int> decode_bch_block(const vector<int>& r) {
    int s1 = 0, s3 = 0;
    for (int i = 0; i < 15; i++) {
        if (r[i] == 1) {
            s1 ^= gf_pow(2, 14 - i);
            s3 ^= gf_pow(2, (14 - i) * 3);
        }
    }
    if (s1 == 0 && s3 == 0) {
        vector<int> msg(7);
        for (int i = 0; i < 7; i++) msg[i] = r[i];
        return msg;
    }

    vector<int> err_locs;
    if (s1 != 0) {
        int delta = s3 ^ gf_pow(s1, 3);
        if (delta == 0) {
            err_locs.push_back(14 - gf_log[s1]);
        }
        else {
            int sigma1 = s1;
            int sigma2 = gf_div(delta, s1);
            for (int j = 0; j < 15; j++) {
                int inv_X = gf_exp[(15 - j) % 15];
                if ((1 ^ gf_mult(sigma1, inv_X) ^ gf_mult(sigma2, gf_pow(inv_X, 2))) == 0) {
                    err_locs.push_back(14 - j);
                }
            }
        }
    }

    vector<int> corrected = r;
    for (int idx : err_locs) {
        if (idx >= 0 && idx < 15) corrected[idx] ^= 1;
    }
    vector<int> msg(7);
    for (int i = 0; i < 7; i++) msg[i] = corrected[i];
    return msg;
}

vector<int> batch_decode(const vector<int>& c, int L) {
    vector<int> m;
    for (int i = 0; i < L; i++) {
        vector<int> chunk(c.begin() + i * 15, c.begin() + (i + 1) * 15);
        vector<int> decoded_chunk = decode_bch_block(chunk);
        m.insert(m.end(), decoded_chunk.begin(), decoded_chunk.end());
    }
    return m;
}

// ==========================================
// 5. McEliece 系统
// ==========================================

Matrix gen_G_sys(int L) {
    Matrix G(7 * L, vector<int>(15 * L, 0));
    Matrix G_sub(7, vector<int>(15));
    for (int i = 0; i < 7; i++) {
        vector<int> u(7, 0); u[i] = 1;
        G_sub[i] = encode_bch_block(u);
    }
    for (int k = 0; k < L; k++) {
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 15; j++) {
                G[k * 7 + i][k * 15 + j] = G_sub[i][j];
            }
        }
    }
    return G;
}

// 加密过程：按照要求在每个分块以 50% 概率注入 1 个错误
vector<int> encrypt(const vector<int>& m, const Matrix& G_pub, int L) {
    vector<int> c = vec_mat_mul(m, G_pub);
    for (int k = 0; k < L; k++) {
        if ((rand() % 100) < 50) {
            int err_pos = k * 15 + (rand() % 15);
            c[err_pos] ^= 1;
        }
    }
    return c;
}

vector<int> decrypt(const vector<int>& c, const Matrix& S, const Matrix& P, int L) {
    Matrix P_inv = mat_transpose(P);
    vector<int> c_prime = vec_mat_mul(c, P_inv);
    vector<int> mS = batch_decode(c_prime, L);
    Matrix S_inv = mat_inverse(S);
    return vec_mat_mul(mS, S_inv);
}

// ==========================================
// 6. 测试流程
// ==========================================

void run_test_case(int L) {
    int K = 7 * L;
    int N = 15 * L;

    cout << "\n==========================================" << endl;
    cout << "  运行测试方案: L=" << L << " (N=" << N << ", K=" << K << ")" << endl;
    cout << "==========================================" << endl;

    Matrix S = gen_S(K);
    Matrix G_sys = gen_G_sys(L);
    Matrix P = gen_P(N, 15);
    Matrix G_pub = mat_mul(mat_mul(S, G_sys), P);

    int rounds = 1000;
    int success = 0;
    clock_t start_t = clock();

    for (int r = 0; r < rounds; r++) {
        vector<int> m(K);
        for (int i = 0; i < K; i++) m[i] = rand() % 2;

        vector<int> c = encrypt(m, G_pub, L);
        vector<int> dec = decrypt(c, S, P, L);

        if (dec == m) success++;

        if (r == 0 && L == 10) {
            cout << "[调试] 方案 L=10 第一行 G_sys 验证:" << endl;
            cout << "Expected: 1 0 0 0 0 0 0 1 1 1 0 1 0 0 0" << endl;
            cout << "Actual:   ";
            for (int j = 0; j < 15; j++) cout << G_sys[0][j] << " ";
            cout << endl;
        }
    }

    clock_t end_t = clock();
    double total_ms = (double)(end_t - start_t) / CLOCKS_PER_SEC * 1000.0;

    cout << "测试次数: " << rounds << endl;
    cout << "成功次数: " << success << endl;
    cout << "成功率:   " << fixed << setprecision(2) << (double)success / rounds * 100.0 << "%" << endl;
    cout << "总时间:   " << total_ms << " ms" << endl;
    cout << "平均时间: " << total_ms / rounds << " ms/次" << endl;
    cout << "扩张率:   " << fixed << setprecision(2) << (double)N / K << endl;
}

int main() {
    srand((unsigned)time(0));
    init_gf_tables();
    GLOBAL_H_BLOCK = gen_H_sys_block();

    run_test_case(10);
    run_test_case(20);

    return 0;
}