#pragma once

#include <cstddef>


namespace rapidgzip::blockfinder
{
class Interface
{
public:
    virtual
    ~Interface() = default;

    [[nodiscard]] virtual size_t
    find() = 0;
};
}  // rapidgzip::blockfinder
