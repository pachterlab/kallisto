/*
 * Copyright (c) 2023-2025, NVIDIA CORPORATION.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

template <typename T1, typename T2>
struct tuple_size<cuco::pair<T1, T2>> : integral_constant<size_t, 2> {};

template <typename T1, typename T2>
struct tuple_size<const cuco::pair<T1, T2>> : tuple_size<cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_size<volatile cuco::pair<T1, T2>> : tuple_size<cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_size<const volatile cuco::pair<T1, T2>> : tuple_size<cuco::pair<T1, T2>> {};

template <std::size_t Index, typename T1, typename T2>
struct tuple_element<Index, cuco::pair<T1, T2>> {
  using type = void;
};

template <typename T1, typename T2>
struct tuple_element<0, cuco::pair<T1, T2>> {
  using type = T1;
};

template <typename T1, typename T2>
struct tuple_element<1, cuco::pair<T1, T2>> {
  using type = T2;
};

template <typename T1, typename T2>
struct tuple_element<0, const cuco::pair<T1, T2>> : tuple_element<0, cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_element<1, const cuco::pair<T1, T2>> : tuple_element<1, cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_element<0, volatile cuco::pair<T1, T2>> : tuple_element<0, cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_element<1, volatile cuco::pair<T1, T2>> : tuple_element<1, cuco::pair<T1, T2>> {};

template <typename T1, typename T2>
struct tuple_element<0, const volatile cuco::pair<T1, T2>> : tuple_element<0, cuco::pair<T1, T2>> {
};

template <typename T1, typename T2>
struct tuple_element<1, const volatile cuco::pair<T1, T2>> : tuple_element<1, cuco::pair<T1, T2>> {
};

template <std::size_t Index, typename T1, typename T2>
__host__ __device__ constexpr auto get(cuco::pair<T1, T2>& p) ->
  typename tuple_element<Index, cuco::pair<T1, T2>>::type&
{
  static_assert(Index < 2);
  if constexpr (Index == 0) {
    return p.first;
  } else {
    return p.second;
  }
}

template <std::size_t Index, typename T1, typename T2>
__host__ __device__ constexpr auto get(cuco::pair<T1, T2>&& p) ->
  typename tuple_element<Index, cuco::pair<T1, T2>>::type&&
{
  static_assert(Index < 2);
  if constexpr (Index == 0) {
    return cuda::std::move(p.first);
  } else {
    return cuda::std::move(p.second);
  }
}

template <std::size_t Index, typename T1, typename T2>
__host__ __device__ constexpr auto get(cuco::pair<T1, T2> const& p) ->
  typename tuple_element<Index, cuco::pair<T1, T2>>::type const&
{
  static_assert(Index < 2);
  if constexpr (Index == 0) {
    return p.first;
  } else {
    return p.second;
  }
}

template <std::size_t Index, typename T1, typename T2>
__host__ __device__ constexpr auto get(cuco::pair<T1, T2> const&& p) ->
  typename tuple_element<Index, cuco::pair<T1, T2>>::type const&&
{
  static_assert(Index < 2);
  if constexpr (Index == 0) {
    return cuda::std::move(p.first);
  } else {
    return cuda::std::move(p.second);
  }
}
