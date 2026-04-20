# Ratio_attack_against_G-G

This repository contains an implementation of the **G+G** digital signature scheme and the **ratio attack** against the G+G scheme.

## Programs

- **`GG_main`**: Runs a single end-to-end demo of **G+G = {Keygen, Sign, Verify}**.
- **`GG_correlation`**: Runs the **ratio attack** experiment (correlation-based analysis).
  - The experimental results can be found in **`correlation_data.csv`**.

## Build & Run (macOS / Homebrew)

> These commands assume your dependencies are installed via Homebrew and available under `/opt/homebrew`.

### 1) `GG_main`

**Compile:**

```bash
rm -f GG_main && g++ -Wall -I/opt/homebrew/include/ -I/opt/homebrew/include/eigen3 -std=gnu++14 -Ofast \
  $(ls *.cc | grep -v -E '^(IBE|GG_correlation|GG_klein_check)\.cc$') \
  -o GG_main -L/opt/homebrew/lib/ -lntl -lgmp -lssl -lcrypto
```

**Run:**

```bash
./GG_main
```

### 2) `GG_correlation`

**Compile:**

```bash
rm -f GG_correlation
g++ -Wall -I/opt/homebrew/include/ -I/opt/homebrew/include/eigen3 -std=gnu++14 -Ofast \
  $(ls *.cc | grep -v -E '^(IBE|GG_main|GG_klein_check)\.cc$') \
  -o GG_correlation \
  -L/opt/homebrew/lib/ -lntl -lgmp -lssl -lcrypto
```

**Run:**

```bash
./GG_correlation
```
