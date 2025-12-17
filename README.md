# McEliece Cryptosystem Implementations with BCH and Hamming Codes

[![GitHub license](https://img.shields.io/github/license/zongemory-spec/2?style=flat-square)](https://github.com/zongemory-spec/2/blob/main/LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/zongemory-spec/2?style=flat-square)](https://github.com/zongemory-spec/2/stargazers)
[![GitHub repo size](https://img.shields.io/github/repo-size/zongemory-spec/2?style=flat-square)](https://github.com/zongemory-spec/2)

This project provides two distinct implementations of the McEliece public-key cryptosystem, leveraging different linear block codes: BCH(15,7, t=2) and Hamming(15,11) codes. Both implementations adhere to the classic McEliece framework, utilizing the robust error-correcting capabilities of these codes to establish secure encryption and decryption. The BCH-based implementation is designed to correct up to 2 errors per block, while the Hamming-based variant corrects up to 1 error per block. The repository encompasses a comprehensive suite of functionalities, including code generation, encoding/decoding, matrix operations over GF(2) (and GF(2⁴) for BCH), key generation, encryption, decryption, and integrated performance testing (evaluating success rate, time consumption, and ciphertext expansion rate).

## Table of Contents
- [🚀 Features](#-features)
- [🛠 Dependencies & Installation](#-dependencies--installation)
- [💻 Usage & Running Experiments](#-usage--running-experiments)
- [⚙️ Configuration & Parameters](#️-configuration--parameters)
- [🔬 Reproducibility](#-reproducibility)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)

## 🚀 Features

*   **Dual Error-Correcting Code Support**: Implements the McEliece cryptosystem with both BCH(15,7, t=2) and Hamming(15,11) codes, allowing adaptation to scenarios with varying error correction requirements.
*   **Complete Cryptosystem Workflow**: Covers the entire cryptographic process, from key generation (including invertible matrix S, permutation matrix P, and public key) to encryption (message multiplied by public key with error injection) and decryption (inverse permutation, decoding, and inverse matrix operation).
*   **GF(2) & GF(2⁴) Operations**: Provides essential matrix operations (multiplication, transposition, inversion) over GF(2), and specialized Galois field operations (multiplication, division, exponentiation) over GF(2⁴) specifically for the BCH code implementation.
*   **Performance Testing**: Includes built-in test functions to rigorously verify decryption success rate, measure average encryption/decryption time, and calculate the ciphertext expansion rate under different block counts (L).

## 🛠 Dependencies & Installation

This project is self-contained and relies solely on the C++ Standard Library, requiring no external dependencies.

### Dependencies

| Category          | Dependency          | Version / Notes                               |
| :---------------- | :------------------ | :-------------------------------------------- |
| **Programming Language** | C++                 | C++11 or newer standard compliant compiler |
| **Standard Libraries**   | C++ Standard Library | `iostream`, `vector`, `ctime`, `cstdlib`, `algorithm`, etc. |

### Installation

To get started, clone the repository and compile the source files using a C++ compiler (e.g., `g++`).

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/zongemory-spec/2.git
    cd 2
    ```

2.  **Compile the source files**:
    *   **For the BCH-based implementation**:
        ```bash
        g++ B\(1\).cpp -o bch_mceliece -std=c++11
        ```
    *   **For the Hamming-based implementation**:
        ```bash
        g++ A\(2\).cpp -o hamming_mceliece -std=c++11
        ```

## 💻 Usage & Running Experiments

After successful compilation, you can run the executables directly. The test cases for each cryptosystem will execute automatically and print the performance results to the console.

### Running the Tests

```bash
# Run BCH-based McEliece cryptosystem test
./bch_mceliece

# Run Hamming-based McEliece cryptosystem test
./hamming_mceliece
```

### Expected Output

Upon execution, the programs will output detailed test results, including:
*   **Decryption success rate**: This should ideally be close to 100% when the injected errors are within the code's correction capability.
*   **Average encryption/decryption time per round**: Provides insight into the computational efficiency.
*   **Ciphertext expansion rate**: Calculated as `total_N / total_K`, indicating how much the message length expands after encryption.

### Modifying Test Conditions

To adjust the experimental conditions, you can modify the following parameters directly within the source code:
*   **`L` (Number of blocks)**: This parameter determines the total message length (`total_K = L × k`) and ciphertext length (`total_N = L × n`). You can find and adjust `L` in the `run_test_case` function for the BCH implementation (`B(1).cpp`) or the `run_test` function for the Hamming implementation (`A(2).cpp`). The provided tests use `L=10` and `L=20`.
*   **Error Injection Probability**: The probability of injecting an error into each block can be adjusted within the test functions. The current default is 50% per block.

## ⚙️ Configuration & Parameters

The behavior of the cryptosystems and their performance tests are governed by several key parameters, primarily defined within the source code.

### Key Parameters

| Parameter Name     | Code File(s) | Type   | Description                                                                                                                                                                                                                                                                          |
| :----------------- | :----------- | :----- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `L`                | `A(2).cpp`, `B(1).cpp` | `int`  | **Number of blocks**. Determines the total information bits (`total_K = L × k`) and total ciphertext bits (`total_N = L × n`) processed in a test run. Higher `L` values lead to longer messages and more extensive tests.                                                              |
| `BLOCK_N`          | `B(1).cpp`   | `int`  | **Code length** for the BCH(15,7, t=2) code. The total number of bits in a codeword block. (Fixed at 15 for BCH).                                                                                                                                                                     |
| `BLOCK_K`          | `B(1).cpp`   | `int`  | **Information bits** for the BCH(15,7, t=2) code. The number of message bits encoded into each block. (Fixed at 7 for BCH).                                                                                                                                                           |
| `BLOCK_T`          | `B(1).cpp`   | `int`  | **Error-correcting capability** for the BCH(15,7, t=2) code. The maximum number of errors per block that the decoder can correct. (Fixed at 2 for BCH).                                                                                                                               |
| `BASE_N`           | `A(2).cpp`   | `int`  | **Code length** for the Hamming(15,11) code. The total number of bits in a codeword block. (Fixed at 15 for Hamming).                                                                                                                                                               |
| `BASE_K`           | `A(2).cpp`   | `int`  | **Information bits** for the Hamming(15,11) code. The number of message bits encoded into each block. (Fixed at 11 for Hamming).                                                                                                                                                         |
| `BASE_M`           | `A(2).cpp`   | `int`  | **Parity bits** for the Hamming(15,11) code. The number of parity bits added to each block (`BASE_N - BASE_K`). (Fixed at 4 for Hamming).                                                                                                                                              |
| `Error Probability` | `A(2).cpp`, `B(1).cpp` | `double` | The probability (e.g., 0.5 for 50%) that an error will be injected into a bit during the encryption process. This parameter directly influences the success rate of decryption and is typically set within the test functions. (Currently 50% per block in the provided code). |

## 🔬 Reproducibility

To ensure the reproducibility of experiments, especially those involving random error injection or key generation, it is crucial to manage the random number generator's seed.

The C++ standard library's `rand()` function is typically seeded using `srand()`. By default, many applications use `srand(time(NULL))` to initialize the random number generator with a time-dependent seed, leading to different results on each run.

To achieve reproducible results:

1.  **Locate `srand()` calls**: Examine the source files (`A(2).cpp` and `B(1).cpp`) for calls to `srand()`.
2.  **Set a fixed seed**: Change `srand(time(NULL))` to `srand(fixed_seed)` where `fixed_seed` is a constant integer (e.g., `srand(42)`). This will ensure that the sequence of "random" numbers generated by `rand()` is identical across multiple runs.

    Example modification:
    ```cpp
    // Original (non-reproducible)
    // srand(time(NULL));

    // Modified (reproducible)
    srand(42); // Use any fixed integer value
    ```
    By setting a fixed seed, the error patterns injected and other random aspects of the cryptosystem will be consistent, allowing for exact replication of experimental outcomes.

## 🤝 Contributing

This repository is primarily for an academic assignment. While direct contributions in the form of pull requests are not actively sought, we welcome and encourage you to open issues for any clarifications regarding the experiments, code logic, or potential improvements. Your feedback is valuable for understanding and enhancing the project.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/zongemory-spec/2/blob/main/LICENSE) file for details.

---
**Author**: zongemory-spec
