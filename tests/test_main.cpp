#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include "test.hpp"

void print_unit_msg(const std::string& msg)
{
    BOOST_LOG_TRIVIAL(info) << "\n" <<
    std::string(msg.length() + 15, '=') <<
    "\n==> UNIT TEST: " << msg << "\n" <<
    std::string(msg.length() + 15, '=') << "\n";
}