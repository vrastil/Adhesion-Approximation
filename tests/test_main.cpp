#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include "test.hpp"

void print_unit_msg(const std::string& msg)
{
    std::cout << std::endl;
    std::cout << std::string(msg.length() + 15, '=') << std::endl;
    std::cout << "==> UNIT TEST: " << msg << std::endl;
    std::cout << std::string(msg.length() + 15, '=') << std::endl;
}