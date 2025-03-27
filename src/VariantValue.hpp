#pragma once
#include <variant>
#include <string>
#include <complex>
#include <optional>
#include <iostream>
#include <type_traits>

class VariantValue {
public:
    using ValueType = std::variant<
        std::monostate,
        bool,
        int,
        float,
        double,
        std::string,
        std::complex<double>
    >;

    // Construct/input ====================================
    VariantValue() = default;
    
    // General constructor（C++17）
    template <typename T, typename = std::enable_if_t<
        std::is_constructible_v<ValueType, T>
    >>
    VariantValue(T&& value) : data_(std::forward<T>(value)) {}

    // General setter ==================================
    template <typename T, typename = std::enable_if_t<
        std::is_constructible_v<ValueType, T>
    >>
    void set(T&& value) {
        data_ = std::forward<T>(value);
    }

    // Type-safe access =================================
    template <typename T>
    std::optional<T> try_get() const noexcept {
        const T* ptr = std::get_if<T>(&data_);
        return ptr ? std::optional(*ptr) : std::nullopt;
    }

    template <typename T>
    T get() const {
        return std::get<T>(data_);
    }

    // Type conversion =====================================
    template <typename To>
    std::optional<To> get_as() const noexcept {  
        return std::visit([](const auto& v) -> std::optional<To> {
            using From = std::decay_t<decltype(v)>;
            if constexpr (std::is_convertible_v<From, To>) {
                return static_cast<To>(v);
            }
            // Convert string to number
            else if constexpr (std::is_same_v<From, std::string>) {
                try {
                    if constexpr (std::is_same_v<To, int>) {
                        return std::stoi(v);
                    } else if constexpr (std::is_same_v<To, float>) {
                        return std::stof(v);
                    } else if constexpr (std::is_same_v<To, double>) {
                        return std::stod(v);
                    } else {
                        return std::nullopt;
                    }
                } catch (...) {
                    return std::nullopt;
                }
            }
            else {
                return std::nullopt;
            }
        }, data_);  // Call std::visit directly and pass data_
    }

    // Visitor Mode ===================================
    template <typename Visitor>
    decltype(auto) visit(Visitor&& vis) const {
        return std::visit(std::forward<Visitor>(vis), data_);
    }

    // Comparison operator ===================================
    bool operator==(const VariantValue& rhs) const noexcept {
        return data_ == rhs.data_;
    }

    bool operator!=(const VariantValue& rhs) const noexcept {
        return !(*this == rhs);
    }

    // Stream output =======================================
    friend std::ostream& operator<<(std::ostream& os, const VariantValue& var) {
        var.visit([&os](const auto& v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, std::monostate>) {
                os << "null";
            } else if constexpr (std::is_same_v<T, std::string>) {
                os << '"' << v << '"';
            } else {
                os << v;
            }
        });
        return os;
    }

private:
    ValueType data_;
};