#pragma once

#include "qcor_observable.hpp"
#include "qcor_utils.hpp"
#include "qrt.hpp"

namespace qcor {
enum class QrtType { NISQ, FTQC };

// The QuantumKernel represents the super-class of all qcor
// quantum kernel functors. Subclasses of this are auto-generated
// via the Clang Syntax Handler capability. Derived classes
// provide a destructor implementation that builds up and
// submits quantum instructions to the specified backend. This enables
// functor-like capability whereby programmers instantiate temporary
// instances of this via the constructor call, and the destructor is
// immediately called. More advanced usage is of course possible for
// qcor developers.
//
// This class works by taking the Derived type (CRTP) and the kernel function
// arguments as template parameters. The Derived type is therefore available for
// instantiation within provided static methods on QuantumKernel. The Args...
// are stored in a member tuple, and are available for use when evaluating the
// kernel. Importantly, QuantumKernel provides static adjoint and ctrl methods
// for auto-generating those circuits.
//
// The Syntax Handler will take kernels like this
// __qpu__ void foo(qreg q) { H(q[0]); }
// and create a derived type of QuantumKernel like this
// class foo : public qcor::QuantumKernel<class foo, qreg> {...};
// with an appropriate implementation of constructors and destructors.
// Users can then call for adjoint/ctrl methods like this
// foo::adjoint(q); foo::ctrl(1, q);
template <typename Derived, typename... Args>
class QuantumKernel {
 protected:
  // Tuple holder for variadic kernel arguments
  std::tuple<Args...> args_tuple;

  // Parent kernel - null if this is the top-level kernel
  // not null if this is a nested kernel call
  std::shared_ptr<qcor::CompositeInstruction> parent_kernel;

  // Default, submit this kernel, if parent is given
  // turn this to false
  bool is_callable = true;

  // Turn off destructor execution, useful for
  // qcor developers, not to be used by clients / programmers
  bool disable_destructor = false;

 public:
  // Flag to indicate we only want to
  // run the pass manager and not execute
  bool optimize_only = false;
  QrtType runtime_env = QrtType::NISQ;
  // Default constructor, takes quantum kernel function arguments
  QuantumKernel(Args... args) : args_tuple(std::forward_as_tuple(args...)) {
    runtime_env = (__qrt_env == "ftqc") ? QrtType::FTQC : QrtType::NISQ;
  }

  // Internal constructor, provide parent kernel, this
  // kernel now represents a nested kernel call and
  // appends to the parent kernel
  QuantumKernel(std::shared_ptr<qcor::CompositeInstruction> _parent_kernel,
                Args... args)
      : args_tuple(std::forward_as_tuple(args...)),
        parent_kernel(_parent_kernel),
        is_callable(false) {
    runtime_env = (__qrt_env == "ftqc") ? QrtType::FTQC : QrtType::NISQ;
  }

  // Static method for printing this kernel as a flat qasm string
  static void print_kernel(std::ostream &os, Args... args) {
    Derived derived(args...);
    derived.disable_destructor = true;
    derived(args...);
    xacc::internal_compiler::execute_pass_manager();
    os << derived.parent_kernel->toString() << "\n";
  }

  static void print_kernel(Args... args) {
    Derived derived(args...);
    derived.disable_destructor = true;
    derived(args...);
    xacc::internal_compiler::execute_pass_manager();
    std::cout << derived.parent_kernel->toString() << "\n";
  }

  // Static method to query how many instructions are in this kernel
  static std::size_t n_instructions(Args... args) {
    Derived derived(args...);
    derived.disable_destructor = true;
    derived(args...);
    return derived.parent_kernel->nInstructions();
  }

  // Create the Adjoint of this quantum kernel
  static void adjoint(std::shared_ptr<CompositeInstruction> parent_kernel,
                      Args... args) {
    // instantiate and don't let it call the destructor
    Derived derived(args...);
    derived.disable_destructor = true;

    // run the operator()(args...) call to get the parent_kernel
    derived(args...);

    // get the instructions
    auto instructions = derived.parent_kernel->getInstructions();
    std::shared_ptr<CompositeInstruction> program = derived.parent_kernel;

    // Assert that we don't have measurement
    if (!std::all_of(
            instructions.cbegin(), instructions.cend(),
            [](const auto &inst) { return inst->name() != "Measure"; })) {
      error(
          "Unable to create Adjoint for kernels that have Measure operations.");
    }

    auto provider = qcor::__internal__::get_provider();
    for (int i = 0; i < instructions.size(); i++) {
      auto inst = derived.parent_kernel->getInstruction(i);
      // Parametric gates:
      if (inst->name() == "Rx" || inst->name() == "Ry" ||
          inst->name() == "Rz" || inst->name() == "CPHASE" ||
          inst->name() == "U1" || inst->name() == "CRZ") {
        inst->setParameter(0, -inst->getParameter(0).template as<double>());
      }
      // Handles T and S gates, etc... => T -> Tdg
      else if (inst->name() == "T") {
        auto tdg = provider->createInstruction("Tdg", inst->bits());
        program->replaceInstruction(i, tdg);
      } else if (inst->name() == "S") {
        auto sdg = provider->createInstruction("Sdg", inst->bits());
        program->replaceInstruction(i, sdg);
      }
    }

    // We update/replace instructions in the derived.parent_kernel composite,
    // hence collecting these new instructions and reversing the sequence.
    auto new_instructions = derived.parent_kernel->getInstructions();
    std::reverse(new_instructions.begin(), new_instructions.end());

    // add the instructions to the current parent kernel
    parent_kernel->addInstructions(new_instructions);

    // no measures, so no execute
  }

  // Create the controlled version of this quantum kernel
  static void ctrl(std::shared_ptr<CompositeInstruction> parent_kernel,
                   int ctrlIdx, Args... args) {
    // instantiate and don't let it call the destructor
    Derived derived(args...);
    derived.disable_destructor = true;

    // run the operator()(args...) call to get the the functor
    // as a CompositeInstruction (derived.parent_kernel)
    derived(args...);

    // Use the controlled gate module of XACC to transform
    auto tempKernel = qcor::__internal__::create_composite("temp_control");
    tempKernel->addInstruction(derived.parent_kernel);

    auto ctrlKernel = qcor::__internal__::create_ctrl_u();
    ctrlKernel->expand({
        std::make_pair("U", tempKernel),
        std::make_pair("control-idx", ctrlIdx),
    });

    // std::cout << "HELLO\n" << ctrlKernel->toString() << "\n";
    for (int instId = 0; instId < ctrlKernel->nInstructions(); ++instId) {
      parent_kernel->addInstruction(
          ctrlKernel->getInstruction(instId)->clone());
    }
    // Need to reset and point current program to the parent
    quantum::set_current_program(parent_kernel);
  }

  // Create the controlled version of this quantum kernel
  static void ctrl(std::shared_ptr<CompositeInstruction> parent_kernel,
                   qubit ctrl_qbit, Args... args) {
    int ctrl_bit = (int)ctrl_qbit.second;

    // instantiate and don't let it call the destructor
    Derived derived(args...);
    derived.disable_destructor = true;

    // run the operator()(args...) call to get the the functor
    // as a CompositeInstruction (derived.parent_kernel)
    derived(args...);

    // Use the controlled gate module of XACC to transform
    auto tempKernel = qcor::__internal__::create_composite("temp_control");
    tempKernel->addInstruction(derived.parent_kernel);

    auto ctrlKernel = qcor::__internal__::create_ctrl_u();
    ctrlKernel->expand({{"U", tempKernel},
                        {"control-idx", ctrl_bit},
                        {"control-buffer", ctrl_qbit.first}});

    for (int instId = 0; instId < ctrlKernel->nInstructions(); ++instId) {
      parent_kernel->addInstruction(
          ctrlKernel->getInstruction(instId)->clone());
    }
    // Need to reset and point current program to the parent
    quantum::set_current_program(parent_kernel);
  }

  static Eigen::MatrixXcd as_unitary_matrix(Args... args) {
    Derived derived(args...);
    derived.disable_destructor = true;
    derived(args...);
    qcor::KernelToUnitaryVisitor visitor(derived.parent_kernel->nLogicalBits());
    xacc::InstructionIterator iter(derived.parent_kernel);
    while (iter.hasNext()) {
      auto inst = iter.next();
      if (!inst->isComposite() && inst->isEnabled()) {
        inst->accept(&visitor);
      }
    }
    return visitor.getMat();
  }

  static double observe(Observable &obs, Args... args) {
    // instantiate and don't let it call the destructor
    Derived derived(args...);
    derived.disable_destructor = true;

    // run the operator()(args...) call to get the the functor
    // as a CompositeInstruction (derived.parent_kernel)
    derived(args...);

    auto instructions = derived.parent_kernel->getInstructions();
    // Assert that we don't have measurement
    if (!std::all_of(
            instructions.cbegin(), instructions.cend(),
            [](const auto &inst) { return inst->name() != "Measure"; })) {
      error("Unable to observe kernels that already have Measure operations.");
    }

    xacc::internal_compiler::execute_pass_manager();

    // Will fail to compile if more than one qreg is passed.
    std::tuple<Args...> tmp(std::forward_as_tuple(args...));
    auto q = std::get<qreg>(tmp);
    return qcor::observe(derived.parent_kernel, obs, q);
  }

  static double observe(std::shared_ptr<Observable> obs, Args... args) {
    // instantiate and don't let it call the destructor
    Derived derived(args...);
    derived.disable_destructor = true;

    // run the operator()(args...) call to get the the functor
    // as a CompositeInstruction (derived.parent_kernel)
    derived(args...);

    auto instructions = derived.parent_kernel->getInstructions();
    // Assert that we don't have measurement
    if (!std::all_of(
            instructions.cbegin(), instructions.cend(),
            [](const auto &inst) { return inst->name() != "Measure"; })) {
      error("Unable to observe kernels that already have Measure operations.");
    }

    xacc::internal_compiler::execute_pass_manager();

    // Will fail to compile if more than one qreg is passed.
    std::tuple<Args...> tmp(std::forward_as_tuple(args...));
    auto q = std::get<qreg>(tmp);
    return qcor::observe(derived.parent_kernel, obs, q);
  }

  static std::string openqasm(Args... args) {
    Derived derived(args...);
    derived.disable_destructor = true;
    derived(args...);
    xacc::internal_compiler::execute_pass_manager();
    return __internal__::translate("staq", derived.parent_kernel);
  }

  virtual ~QuantumKernel() {}
};

template <typename... Args>
using callable_function_ptr =
    void (*)(std::shared_ptr<xacc::CompositeInstruction>, Args...);

template <typename... Args>
class CallableKernel {
 protected:
  callable_function_ptr<Args...> &function_pointer;

 public:
  CallableKernel(callable_function_ptr<Args...> &&f) : function_pointer(f) {}
  void operator()(std::shared_ptr<xacc::CompositeInstruction> ir,
                  Args... args) {
    function_pointer(ir, args...);
  }
  void ctrl(std::shared_ptr<xacc::CompositeInstruction> ir, int ctrl_qbit,
            Args... args) {
    auto tempKernel = qcor::__internal__::create_composite("temp_control");
    function_pointer(tempKernel, args...);

    auto ctrlKernel = qcor::__internal__::create_ctrl_u();
    ctrlKernel->expand({
        std::make_pair("U", tempKernel),
        std::make_pair("control-idx", ctrl_qbit),
    });

    for (int instId = 0; instId < ctrlKernel->nInstructions(); ++instId) {
      ir->addInstruction(ctrlKernel->getInstruction(instId)->clone());
    }
  }

  void ctrl(std::shared_ptr<xacc::CompositeInstruction> ir, qubit ctrl_qbit,
            Args... args) {
    int ctrl_bit = (int)ctrl_qbit.second;
    ctrl(ir, ctrl_bit, args...);
  }

  void adjoint(std::shared_ptr<CompositeInstruction> ir, Args... args) {
    auto tempKernel = qcor::__internal__::create_composite("temp_adjoint");
    function_pointer(tempKernel, args...);
 
    // get the instructions
    auto instructions = tempKernel->getInstructions();
    std::shared_ptr<CompositeInstruction> program = tempKernel;

    // Assert that we don't have measurement
    if (!std::all_of(
            instructions.cbegin(), instructions.cend(),
            [](const auto &inst) { return inst->name() != "Measure"; })) {
      error(
          "Unable to create Adjoint for kernels that have Measure operations.");
    }

    auto provider = qcor::__internal__::get_provider();
    for (int i = 0; i < instructions.size(); i++) {
      auto inst = tempKernel->getInstruction(i);
      // Parametric gates:
      if (inst->name() == "Rx" || inst->name() == "Ry" ||
          inst->name() == "Rz" || inst->name() == "CPHASE" ||
          inst->name() == "U1" || inst->name() == "CRZ") {
        inst->setParameter(0, -inst->getParameter(0).template as<double>());
      }
      // Handles T and S gates, etc... => T -> Tdg
      else if (inst->name() == "T") {
        auto tdg = provider->createInstruction("Tdg", inst->bits());
        program->replaceInstruction(i, tdg);
      } else if (inst->name() == "S") {
        auto sdg = provider->createInstruction("Sdg", inst->bits());
        program->replaceInstruction(i, sdg);
      }
    }

    // We update/replace instructions in the derived.parent_kernel composite,
    // hence collecting these new instructions and reversing the sequence.
    auto new_instructions = tempKernel->getInstructions();
    std::reverse(new_instructions.begin(), new_instructions.end());

    // add the instructions to the current parent kernel
    ir->addInstructions(new_instructions);
  }
};

}  // namespace qcor