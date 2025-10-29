#include <type_traits>
#include <limits>
#include <stdexcept>
#include <string>

// written by chatgpt
template <typename ToType, typename FromType>
constexpr ToType safe_num_cast(FromType value) {
    static_assert(std::is_integral_v<ToType> && std::is_integral_v<FromType>,
                  "shrink_int requires integral types");

    using ToTypeLimits = std::numeric_limits<ToType>;
    // using FromTypeLimits = std::numeric_limits<FromType>;

    // Handle signedness properly
    if constexpr (std::is_signed_v<FromType> && std::is_signed_v<ToType>) {
        if (value < ToTypeLimits::min() || value > ToTypeLimits::max()) {
            throw std::runtime_error("shrink_int: signed value out of range");
        }
    } else if constexpr (!std::is_signed_v<FromType> && std::is_signed_v<ToType>) {
        if (value > static_cast<std::make_unsigned_t<ToType>>(ToTypeLimits::max())) {
            throw std::runtime_error("shrink_int: unsigned value too FromType for signed target");
        }
    } else if constexpr (std::is_signed_v<FromType> && !std::is_signed_v<ToType>) {
        if (value < 0 || static_cast<std::make_unsigned_t<FromType>>(value) > ToTypeLimits::max()) {
            throw std::runtime_error("shrink_int: negative or out-of-range value for unsigned target");
        }
    } else { // both unsigned
        if (value > ToTypeLimits::max()) {
            throw std::runtime_error("shrink_int: unsigned value out of range");
        }
    }

    return static_cast<ToType>(value);
}