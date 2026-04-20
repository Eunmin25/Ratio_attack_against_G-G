#include <iostream>
#include <ctime>
#include "GG_scheme.h"

using namespace std;
using namespace GG;

int main() {
    // Initialize random seed for IBE's Sampling functions
    srand((unsigned)time(nullptr));
    
    cout << "==============================================\n";
    cout << "G+G Signature Scheme (Toy Implementation)\n";
    cout << "==============================================\n\n";
    
    cout << "Parameters:\n";
    cout << "  Ring degree N = " << GG_N << "\n";
    cout << "  Module dimensions: m = " << GG_m << ", k = " << GG_k << "\n";
    cout << "  Modulus q = " << GG_q << "\n";
    cout << "  Secret key bound η = " << GG_ETA << "\n";
    cout << "  Secret key distribution: χ_η = U({y ∈ R | ||y||_∞ ≤ η})\n";
    cout << "  Covariance base σ = " << GG_SIGMA_BASE << "\n";
    cout << "  Using Klein sampling: y ~ D_{Z^m, Σ(S)}\n";
    cout << "  where Σ(S) = σ² I_m - s² SS^T, S = rot(ζs)\n\n";
    
    // Key Generation
    cout << "Generating keys...\n";
    PublicKey pk;
    SecretKey sk;
    keygen(pk, sk);
    cout << "Keys generated successfully!\n";
    
    // Print public key
    print_public_key(pk);
    
    // Print secret key
    print_secret_key(sk);
    
    // Sign a message
    string message = "Hello, G+G signature over module lattices!";
    cout << "\nSigning message: \"" << message << "\"\n";
    Signature sig = sign(pk, sk, message);
    
    // Print signature
    print_signature(sig);
    
    // Verify signature
    cout << "Verifying signature...\n";
    bool valid = verify(pk, message, sig);
    cout << "  Verification: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    cout << "==============================================\n";
    cout << "Note: This is a TOY implementation for\n";
    cout << "educational purposes only. Not secure!\n";
    cout << "==============================================\n";
    
    return 0;
}
