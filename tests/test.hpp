#pragma once

inline void print_unit_msg(const std::string& msg)
{
    std::cout << std::endl;
    std::cout << std::string(msg.length() + 15, '=') << std::endl;
    std::cout << "==> UNIT TEST: " << msg << std::endl;
    std::cout << std::string(msg.length() + 15, '=') << std::endl;
}