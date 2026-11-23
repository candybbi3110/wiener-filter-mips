# Wiener Filter in MIPS Assembly

A Wiener filter written in pure **MIPS32 assembly** for the MARS simulator. It reads a noisy signal and a desired signal from text files, finds the optimal filter coefficients, filters the signal, and reports the MMSE. All of the math is done in hand-written floating-point code on Coprocessor 1.

A **C++ reference model** follows the same steps and rounds at the same points, so the two outputs can be compared value by value.

> Team assignment for *Computer Architecture Lab (CO2008)*, HCMUT, Semester 252 (Spring 2026). **Graded 10/10.**
> My part: file I/O and the ASCII → double parser.

---

## What it does

```
input.txt   (x[n], 10 values) ┐
desired.txt (d[n], 10 values) ┴─► parse ─► correlate ─► Toeplitz R ─► solve R·h = γ ─► y = h * x ─► MMSE
                                                                                         │
                                                                    console + output.txt ◄┘
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
1. Put `input.txt` and `desired.txt` in the folder you launch MARS from, e.g. copy them from `tests/sample/`.
2. Open `src/wiener_filter.asm`, assemble it, and run it.
3. The result is printed to the console and written to `output.txt`.

**C++ reference:**
```bash
cmake -B build && cmake --build build
./build/wiener_ref tests/sample/input.txt tests/sample/desired.txt output.txt
```

Example (`tests/sample`):
```
Filtered output: -0.4 -2.6 2.9 2.7 -0.2 -1.4 -0.8 3.1 4.8 5.0
MMSE: 4.7
```
If the two files do not hold the same number of values, the output is `Error : size not match` (see `tests/size-mismatch`).

---

## Design notes & limitations
- **Rounding to 1 decimal at every step:** the assignment's expected results are given to 1 decimal place, so intermediate values in the convolution and the MMSE are rounded at the same points in both versions.
- **No pivoting:** the Gaussian elimination does not pivot. That is fine here because R is a symmetric Toeplitz autocorrelation matrix, which is positive semi-definite. A general-purpose solver would need partial pivoting, or Levinson–Durbin to take advantage of the Toeplitz structure in O(M²).
- **Fixed size:** the filter length and signal length are both 10, and the parser reads at most 1024 bytes. It does not check bounds, so an input file with more than 10 numbers would write past the end of its array before the size check runs.

## Repository layout
```
src/wiener_filter.asm       MIPS implementation (MARS)
reference/wiener_filter.cpp C++ reference model
tests/sample/               input.txt, desired.txt, expected_output.txt
tests/size-mismatch/        error-path test
```
