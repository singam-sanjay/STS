#include <cstdlib>
#include <iostream>

#ifndef _CHECK_ASSERTS_
#define _CHECK_ASSERTS_

#define CHECK(COND, ERR_MSG)                                          \
  {                                                                   \
    if (!(COND)) {                                                    \
      std::cerr << __FILE__ << "::" << __LINE__ << ":: " << (ERR_MSG) \
                << std::endl;                                         \
      exit(1);                                                        \
    }                                                                 \
  }

#define CHECK_EQ(VAL1, VAL2)                                                  \
  {                                                                           \
    auto _val1 = VAL1;                                                        \
    auto _val2 = VAL2;                                                        \
    CHECK(_val1 == _val2, (#VAL1 + std::string("(") + std::to_string(_val1) + \
                           std::string(") != ") + #VAL2 + std::string("(") +  \
                           std::to_string(_val2) + std::string(")")));        \
  }

#define CHECK_LT(VAL1, VAL2)                                                 \
  {                                                                          \
    auto _val1 = VAL1;                                                       \
    auto _val2 = VAL2;                                                       \
    CHECK(_val1 < _val2, (#VAL1 + std::string("(") + std::to_string(_val1) + \
                          std::string(") >= ") + #VAL2 + std::string("(") +  \
                          std::to_string(_val2) + std::string(")")));        \
  }

#define REPORT(x) \
  { std::cout << #x << " = " << x << '\n'; }

#endif
