# Wiener Filter in MIPS Assembly

A Wiener filter written in pure **MIPS32 assembly** for the MARS simulator. It reads a noisy signal and a desired signal from text files, finds the optimal filter coefficients, filters the signal, and reports the MMSE. All of the math is done in hand-written floating-point code on Coprocessor 1.

A **C++ reference model** follows the same steps and rounds at the same points, so the two outputs can be compared value by value.

> Team assignment for *Computer Architecture Lab (CO2008)*, HCMUT, Semester 252 (Spring 2026). **Graded 10/10.**
> My part: file I/O and the ASCII → double parser.

---

## What it does

```
input.txt, desired.txt
   │  parse ASCII → double
   ▼
correlation  γxx, γdx
   │
   ▼
Toeplitz R  →  solve R·h = γdx   (Gaussian elimination)
   │
   ▼
y = h * x   →   MMSE   →   console + output.txt
```

| Step | MIPS procedure | What happens |
|---|---|---|
| 1 | `loadInputFile` | Reads the file with syscalls 13/14/16 and parses ASCII byte by byte into IEEE-754 doubles (sign, integer and fraction parts) |
| 2 | `build_Correlation` | Finite-sample estimates of the autocorrelation γxx(m) and the cross-correlation γdx(m) |
| 3 | `build_R_matrix` | Builds the 10×10 Toeplitz matrix R[i][j] = γxx(\|i−j\|), stored row-major in a flat 800-byte block |
| 4 | `calculate_H_matrix` | Gaussian elimination, then back substitution → optimal coefficients h |
| 5 | `calculate_output` | Convolution y[n] = Σ h[k]·x[n−k] |
| 6 | `calculate_MMSE` | Mean squared error between y and d |
| 7 | `write_double_1digit_file` / `print_double_1digit_console` | Formats doubles to 1 decimal by hand (MARS has no syscall that writes floats to a file) |

### Low-level details worth noting
- **Floating point on Coprocessor 1:** `l.d` / `s.d`, `mul.d` / `div.d`, and conversions between integer and float registers with `mtc1` / `cvt.d.w` / `round.w.d` / `mfc1`.
- **Manual memory layout:** 2-D matrix indexing computed as `(i*10 + j) << 3`, with `.align 3` so doubles are 8-byte aligned.
- **Procedure calls:** `$ra` and caller-saved temporaries are saved on the stack around nested `jal` calls.
- **Output formatting:** numbers are written character by character with syscall 15, because MARS can only write raw bytes to a file.

---

## Run it

**MIPS (MARS 4.5):**
1. Copy `input.txt` and `desired.txt` from a test case (e.g. `tests/tc6/`) into the folder you launch MARS from.
2. Open `src/wiener_filter.asm`, assemble it, and run it.
3. The result is printed to the console and written to `output.txt`.

Example (`tests/tc6`):
```
Filtered output: 1.1 2.1 3.6 1.5 -0.6 -3.8 -2.3 3.4 6.4 5.3
MMSE: 1.2
```

**C++ reference + test runner:**
```bash
cmake -B build && cmake --build build
tests/run_reference_tests.sh          # runs all test cases, prints PASS/FAIL
```

---

## Validation

8 test cases (`tests/tc1` to `tests/tc8`); TC3 checks the error path, because its input has 9 values instead of 10.

| Test | MMSE | Note |
|---|---|---|
| TC1 | 0.1 | |
| TC2 | 0.1 | |
| TC3 | — | `Error : size not match` |
| TC4 | 0.0 | all-negative signal |
| TC5 | 0.1 | |
| TC6 | 1.2 | TC6–TC8 share one desired signal |
| TC7 | 2.5 | |
| TC8 | 0.5 | |

**How it was checked:** `tests/run_reference_tests.sh` runs the C++ model automatically on all 8 cases and checks both the output and the exit code (0 for valid input, 1 for the TC3 error case). The MIPS outputs were checked by hand in MARS. For TC2 and TC4–TC8, the saved MARS output matches the C++ output value for value. TC1's saved MARS output was produced by an older build, from before the rounding fix, and still needs to be re-run.

**Error path difference:** on a size mismatch, the MIPS program prints `Error : size not match` to the console and exits before it opens `output.txt`. The C++ model writes the message to both the console and the file.

---

## Design notes & limitations
- **Problem → fix: rounding.** The first version rounded only the final results, and some outputs came out 0.1 off from the expected values. Rounding the running sum to 1 decimal after every accumulation, in both the convolution and the MMSE, made the results match. The C++ model rounds at exactly the same points, which is what makes a value-by-value comparison possible.
- **No pivoting:** the Gaussian elimination does not pivot, and it does not check for a zero pivot. Here R is built from the *biased* autocorrelation estimate, so R = AᵀA / N, where A is the zero-padded convolution matrix of x. That makes R symmetric positive definite whenever x is not all zeros, and elimination without pivoting is numerically stable for such matrices. Its pivots are positive in exact arithmetic. However, R can still be ill-conditioned, an all-zero input would make it singular, and the solver has no guard for either case. A general-purpose solver would add a pivot check or partial pivoting, or use Levinson–Durbin to take advantage of the Toeplitz structure in O(M²).
- **Tie-breaking in rounding:** MIPS `round.w.d` rounds ties to even, while C++ `std::round` rounds ties away from zero. The two can disagree on exact .x5 values, for example 0.25.
- **Fixed size:** the filter length and signal length are both 10, and the parser reads at most 1024 bytes. It does not check bounds, so an input file with more than 10 numbers would write past the end of its array before the size check runs.

## Repository layout
```
src/wiener_filter.asm       MIPS implementation (MARS)
reference/wiener_filter.cpp C++ reference model
tests/tcN/                  input.txt, desired.txt, expected_output.txt
tests/run_reference_tests.sh  runs the C++ model on all test cases
```
