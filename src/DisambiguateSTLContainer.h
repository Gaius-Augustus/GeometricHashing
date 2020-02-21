#ifndef DISAMBIGUATESTLCONTAINER_H
#define DISAMBIGUATESTLCONTAINER_H

// https://stackoverflow.com/questions/21660325/disambiguate-template-specialization-between-map-like-and-vector-like-containers
inline constexpr auto is_container_impl(...) {
    return std::false_type{};
}

template <typename C>
constexpr auto is_container_impl(C const & c)
-> decltype(begin(c), end(c), std::true_type{}) {   // https://stackoverflow.com/questions/16044514/what-is-decltype-with-two-arguments
    return std::true_type{};
}

template <typename C>
constexpr auto is_container(C const & c) {
    return is_container_impl(c);
}

inline constexpr auto is_map_container_impl(...) {
    return std::false_type{};
}

template <typename C, typename = typename C::key_type, typename = typename C::value_type::first_type>
constexpr auto is_map_container_impl(C const &) {
    return std::true_type{};
}

template <typename C>
constexpr auto is_map_container(C const & c) {
    return is_map_container_impl(c);
}

inline constexpr auto is_associative_container_impl(...) {
    return std::false_type{};
}

template <typename C, typename = typename C::key_type>
constexpr auto is_associative_container_impl(C const &) {
    return std::true_type{};
}

template <typename C>
constexpr auto is_associative_container(C const & c) {
    return is_associative_container_impl(c);
}
#endif // DISAMBIGUATESTLCONTAINER_H
