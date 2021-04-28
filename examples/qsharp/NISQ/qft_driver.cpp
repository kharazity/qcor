#include <iostream> 
#include <vector>
#include "import_kernel_utils.hpp"

using QPEOracleSignature = KernelSignature<qubit>;
qcor_import_qsharp_kernel(QCOR__InverseQFT);

__qpu__ void qpe(qreg q, QPEOracleSignature oracle) {
  // Extract the counting qubits and the state qubit
  auto counting_qubits = q.extract_range({0,3});
  auto state_qubit = q[3];
  // Put it in |1> eigenstate
  X(state_qubit);
  // Create uniform superposition on all 3 qubits
  H(counting_qubits);

  // run ctr-oracle operations
  for (auto i : range(counting_qubits.size())) {
    const int nbCalls = 1 << i;
    for (auto j : range(nbCalls)) {
      oracle.ctrl(counting_qubits[i], state_qubit);
    }
  }

  // Run Inverse QFT on counting qubits
  // Using the Q# Kernel (wrapped as a QCOR kernel)
  QCOR__InverseQFT(parent_kernel, counting_qubits);

  // Measure the counting qubits
  Measure(counting_qubits);
}

// Oracle to consider
__qpu__ void oracle(qubit q) { T(q); }

int main(int argc, char **argv) {
  auto q = qalloc(4);
  qpe::print_kernel(std::cout, q, oracle);
}
