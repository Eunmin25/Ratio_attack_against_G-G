#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

echo "[1/4] Rebuild core objects (C++14 for Eigen)"
g++ -Wall -I/opt/homebrew/include/ -std=c++14 -Ofast -c GG_scheme.cc -o GG_scheme.o
g++ -Wall -I/opt/homebrew/include/ -std=c++14 -Ofast -c GG_sampler.cc -o GG_sampler.o
g++ -Wall -I/opt/homebrew/include/ -I/opt/homebrew/include/eigen3 -std=c++14 -Ofast -c GG_klein.cc -o GG_klein.o
g++ -Wall -I/opt/homebrew/include/ -std=c++14 -Ofast -c GG_klein_check.cc -o GG_klein_check.o
g++ -Wall -I/opt/homebrew/include/ -std=gnu++0x -Ofast -c GG_correlation.cc -o GG_correlation.o

echo "[2/4] Link binaries (without GG_decomp_stub)"
g++ -o GG_klein_check GG_klein_check.o GG_scheme.o GG_klein.o GG_sampler.o GG_poly.o Sampling.o Random.o Scheme.o io.o Algebra.o FFT.o -L/opt/homebrew/lib/ -lntl -lgmp -lssl -lcrypto
g++ -o GG_correlation GG_correlation.o Algebra.o FFT.o GG_klein.o GG_poly.o GG_sampler.o GG_scheme.o Sampling.o Random.o Scheme.o io.o -L/opt/homebrew/lib/ -lntl -lgmp -lssl -lcrypto

echo "[3/4] y-covariance check (mode off)"
GG_PI_PARAM_MODE=0 GG_STRICT_KLEIN=1 GG_CHECK_SAMPLES="${GG_CHECK_SAMPLES:-200}" GG_CHECK_COEFF="${GG_CHECK_COEFF:-0}" \
  ./GG_klein_check > /tmp/gg_klein_check_mode0.out 2>/tmp/gg_klein_check_mode0.err
echo "  output: /tmp/gg_klein_check_mode0.out"
echo "  stderr: /tmp/gg_klein_check_mode0.err"

echo "[4/4] y-covariance check (mode on: /2pi)"
GG_PI_PARAM_MODE=1 GG_STRICT_KLEIN=1 GG_CHECK_SAMPLES="${GG_CHECK_SAMPLES:-200}" GG_CHECK_COEFF="${GG_CHECK_COEFF:-0}" \
  ./GG_klein_check > /tmp/gg_klein_check_mode1.out 2>/tmp/gg_klein_check_mode1.err
echo "  output: /tmp/gg_klein_check_mode1.out"
echo "  stderr: /tmp/gg_klein_check_mode1.err"

echo
echo "Done. Compare mode0 vs mode1 covariance fit from:"
echo "  /tmp/gg_klein_check_mode0.out"
echo "  /tmp/gg_klein_check_mode1.out"

