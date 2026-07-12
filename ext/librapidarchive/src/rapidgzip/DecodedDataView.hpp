#pragma once

#include <array>
#include <cstdint>

#include <core/VectorView.hpp>


namespace rapidgzip::deflate
{
/**
 * Only one of the two will contain non-empty VectorViews depending on whether marker bytes might appear.
 * @ref dataWithMarkers will be empty when @ref setInitialWindow has been called.
 */
struct DecodedDataView
{
public:
    [[nodiscard]] constexpr size_t
    size() const noexcept
    {
        return dataWithMarkers[0].size() + dataWithMarkers[1].size() + data[0].size() + data[1].size();
    }

    [[nodiscard]] size_t
    dataSize() const noexcept
    {
        return data[0].size() + data[1].size();
    }

    [[nodiscard]] size_t
    dataWithMarkersSize() const noexcept
    {
        return dataWithMarkers[0].size() + dataWithMarkers[1].size();
    }

    [[nodiscard]] constexpr bool
    containsMarkers() const noexcept
    {
        return !dataWithMarkers[0].empty() || !dataWithMarkers[1].empty();
    }

public:
    std::array<VectorView<uint16_t>, 2> dataWithMarkers{};
    std::array<VectorView<uint8_t>, 2> data{};
};
}  // namespace rapidgzip::deflate
