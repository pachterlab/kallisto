/*
 * Copyright (c) 2020-2023, NVIDIA CORPORATION.
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

#pragma once

#include <cuco/utility/error.hpp>

#include <cuda_runtime_api.h>

#define STRINGIFY_DETAIL(x) #x
#define CUCO_STRINGIFY(x)   STRINGIFY_DETAIL(x)

/**
 * @brief Error checking macro for CUDA runtime API functions.
 *
 * Invokes a CUDA runtime API function call. If the call does not return
 * `cudaSuccess`, invokes cudaGetLastError() to clear the error and throws an
 * exception detailing the CUDA error that occurred
 *
 * Defaults to throwing `cuco::cuda_error`, but a custom exception may also be
 * specified.
 *
 * Example:
 * ```c++
 *
 * // Throws `cuco::cuda_error` if `cudaMalloc` fails
 * CUCO_CUDA_TRY(cudaMalloc(&p, 100));
 *
 * // Throws `std::runtime_error` if `cudaMalloc` fails
 * CUCO_CUDA_TRY(cudaMalloc(&p, 100), std::runtime_error);
 * ```
 *
 */
#define CUCO_CUDA_TRY(...)                                               \
  GET_CUCO_CUDA_TRY_MACRO(__VA_ARGS__, CUCO_CUDA_TRY_2, CUCO_CUDA_TRY_1) \
  (__VA_ARGS__)
#define GET_CUCO_CUDA_TRY_MACRO(_1, _2, NAME, ...) NAME
#define CUCO_CUDA_TRY_2(_call, _exception_type)                                                    \
  do {                                                                                             \
    cudaError_t const error = (_call);                                                             \
    if (cudaSuccess != error) {                                                                    \
      cudaGetLastError();                                                                          \
      throw _exception_type{std::string{"CUDA error at: "} + __FILE__ + CUCO_STRINGIFY(__LINE__) + \
                            ": " + cudaGetErrorName(error) + " " + cudaGetErrorString(error)};     \
    }                                                                                              \
  } while (0);
#define CUCO_CUDA_TRY_1(_call) CUCO_CUDA_TRY_2(_call, cuco::cuda_error)

/**
 * @brief Error checking macro for CUDA runtime API that asserts the result is
 * equal to `cudaSuccess`.
 *
 */
#define CUCO_ASSERT_CUDA_SUCCESS(expr) \
  do {                                 \
    cudaError_t const status = (expr); \
    assert(cudaSuccess == status);     \
  } while (0)

/**
 * @brief Macro for checking (pre-)conditions that throws an exception when
 * a condition is violated.
 *
 * Defaults to throwing `cuco::logic_error`, but a custom exception may also be
 * specified.
 *
 * Example usage:
 * ```
 * // throws cuco::logic_error
 * CUCO_EXPECTS(p != nullptr, "Unexpected null pointer");
 *
 * // throws std::runtime_error
 * CUCO_EXPECTS(p != nullptr, "Unexpected nullptr", std::runtime_error);
 * ```
 * @param ... This macro accepts either two or three arguments:
 *   - The first argument must be an expression that evaluates to true or
 *     false, and is the condition being checked.
 *   - The second argument is a string literal used to construct the `what` of
 *     the exception.
 *   - When given, the third argument is the exception to be thrown. When not
 *     specified, defaults to `cuco::logic_error`.
 * @throw `_exception_type` if the condition evaluates to 0 (false).
 */
#define CUCO_EXPECTS(...)                                             \
  GET_CUCO_EXPECTS_MACRO(__VA_ARGS__, CUCO_EXPECTS_3, CUCO_EXPECTS_2) \
  (__VA_ARGS__)

#define GET_CUCO_EXPECTS_MACRO(_1, _2, _3, NAME, ...) NAME

#define CUCO_EXPECTS_3(_condition, _reason, _exception_type)                    \
  do {                                                                          \
    static_assert(std::is_base_of_v<std::exception, _exception_type>);          \
    (_condition) ? static_cast<void>(0)                                         \
                 : throw _exception_type /*NOLINT(bugprone-macro-parentheses)*/ \
      {"CUCO failure at: " __FILE__ ":" CUCO_STRINGIFY(__LINE__) ": " _reason}; \
  } while (0)

#define CUCO_EXPECTS_2(_condition, _reason) CUCO_EXPECTS_3(_condition, _reason, cuco::logic_error)

/**
 * @brief Indicates that an erroneous code path has been taken.
 *
 * Example usage:
 * ```c++
 * // Throws `cuco::logic_error`
 * CUCO_FAIL("Unsupported code path");
 *
 * // Throws `std::runtime_error`
 * CUCO_FAIL("Unsupported code path", std::runtime_error);
 * ```
 *
 * @param ... This macro accepts either one or two arguments:
 *   - The first argument is a string literal used to construct the `what` of
 *     the exception.
 *   - When given, the second argument is the exception to be thrown. When not
 *     specified, defaults to `cuco::logic_error`.
 * @throw `_exception_type` if the condition evaluates to 0 (false).
 */
#define CUCO_FAIL(...)                                       \
  GET_CUCO_FAIL_MACRO(__VA_ARGS__, CUCO_FAIL_2, CUCO_FAIL_1) \
  (__VA_ARGS__)

#define GET_CUCO_FAIL_MACRO(_1, _2, NAME, ...) NAME

#define CUCO_FAIL_2(_what, _exception_type)      \
  /*NOLINTNEXTLINE(bugprone-macro-parentheses)*/ \
  throw _exception_type { "CUCO failure at:" __FILE__ ":" CUCO_STRINGIFY(__LINE__) ": " _what }

#define CUCO_FAIL_1(_what) CUCO_FAIL_2(_what, cuco::logic_error)
