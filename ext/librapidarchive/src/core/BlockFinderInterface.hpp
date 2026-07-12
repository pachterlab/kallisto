#pragma once

#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>

namespace rapidgzip
{
class BlockFinderInterface
{
public:
    enum class GetReturnCode
    {
        SUCCESS,
        TIMEOUT,
        FAILURE,
    };

public:
    virtual
    ~BlockFinderInterface() = default;

    [[nodiscard]] virtual size_t
    size() const = 0;

    [[nodiscard]] virtual bool
    finalized() const = 0;

    [[nodiscard]] virtual std::pair<std::optional<size_t>, GetReturnCode>
    get( size_t blockIndex,
         double timeoutInSeconds ) = 0;

    [[nodiscard]] virtual std::optional<size_t>
    get( size_t blockIndex )
    {
        const auto [result, returnCode] = get( blockIndex, std::numeric_limits<double>::infinity() );
        assert( returnCode != GetReturnCode::TIMEOUT );
        return result;
    }

    [[nodiscard]] virtual size_t
    find( size_t encodedBlockOffsetInBits ) const = 0;
};
}  // namespace rapidgzip
