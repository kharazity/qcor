#pragma once

#include "FermionOperator.hpp"
#include "Observable.hpp"
#include "qcor_utils.hpp"

namespace qcor{
    using Observable = xacc::Observable;
class ObservableProvider : public xacc::Identifiable{
    public:
        virtual std::shared_ptr<Observable> createOperator(HeterogeneousMap &options) = 0;        
};//ObservableProvider
}//namespace qcor
