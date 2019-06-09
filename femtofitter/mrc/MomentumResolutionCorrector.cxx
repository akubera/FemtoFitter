
#include "MomentumResolutionCorrector.hpp"

/// \cond classimp
ClassImp(MomentumResolutionCorrector);
/// \endcond

#include <iostream>


MomentumResolutionCorrector::MomentumResolutionCorrector()
 : TNamed("mrc", "MomentumResolutionCorrection")
 , data()
{
    // std::map<int, float> m;
    // m = {{4, 5.4}, {7, 1.1}};
    // std::cout << m[4] << "\n";
    // m.emplace(std::make_pair<int, float>(4, 5.5));
    // m.emplace(4, 5.5);
    // std::array<int, 3> q = {1,1,1};
    // a.emplace(q, m);

    I3 key = {{1,1,1}};
    std::map<u8, u8> val {{5, 5}};

    data.emplace(key, val); // {{3, {{3, {{3, 5}}}}}});
}
