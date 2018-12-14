//  Copyright: (c) lesya89
#ifndef MODULES_TASK_3_ZOLOTAREVA_TASK3_LSD_RADIX_SORT_INCLUDE_RADIX_SORT_H_
#define MODULES_TASK_3_ZOLOTAREVA_TASK3_LSD_RADIX_SORT_INCLUDE_RADIX_SORT_H_

#include <algorithm>  // For swap
#include <cstdint>    // For system independent base types
#include <cstring>    // For memset
#include <vector>
#include <utility>

template <typename ValueType>
void RadixSortMSD(ValueType* array, size_t N);

template <size_t Base, typename ValueType>
void RadixSortLSD(ValueType* array, size_t N);

namespace details {

using Radixable8 = uint8_t;
using Radixable16 = uint16_t;
using Radixable32 = uint32_t;
using Radixable64 = uint64_t;

template <size_t Bytes>
struct Radixable {};

template <>
struct Radixable<1> {
  using type = Radixable8;
};
template <>
struct Radixable<2> {
  using type = Radixable16;
};
template <>
struct Radixable<4> {
  using type = Radixable32;
};
template <>
struct Radixable<8> {
  using type = Radixable64;
};

template <class ValueType, class RadixableType>
struct RadixTransform {
  RadixableType TransformFrom(ValueType);
  ValueType TransformTo(RadixableType);
};

/*
 * The first group of transformation rules are by far the easiest and most
 * straight forward. Unsigned numbers are already in radixable form.
 */
template <>
struct RadixTransform<uint8_t, Radixable8> {
  Radixable8 TransformFrom(uint8_t value) { return value; }
  uint8_t TransformTo(Radixable8 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint16_t, Radixable16> {
  Radixable16 TransformFrom(uint16_t value) { return value; }
  uint16_t TransformTo(Radixable16 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint32_t, Radixable32> {
  Radixable32 TransformFrom(uint32_t value) { return value; }
  uint32_t TransformTo(Radixable32 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint64_t, Radixable64> {
  Radixable64 TransformFrom(uint64_t value) { return value; }
  uint64_t TransformTo(Radixable64 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<int8_t, Radixable8> {
  Radixable8 TransformFrom(int8_t value) {
    Radixable8 rdxble = *reinterpret_cast<Radixable8 *>(&value);
    rdxble ^= (1 << 7);
    return rdxble;
  }

  int8_t TransformTo(Radixable8 rdxble) {
    int8_t value = *reinterpret_cast<int8_t*>(&rdxble);
    value ^= (1 << 7);
    return value;
  }
};

template <>
struct RadixTransform<int16_t, Radixable16> {
  Radixable16 TransformFrom(int16_t value) {
    Radixable16 rdxble = *reinterpret_cast<Radixable16*>(&value);
    rdxble ^= (1 << 15);
    return rdxble;
  }

  int16_t TransformTo(Radixable16 rdxble) {
    int16_t value = *reinterpret_cast<int16_t *>(&rdxble);
    value ^= (1 << 15);
    return value;
  }
};

template <>
struct RadixTransform<int32_t, Radixable32> {
  Radixable32 TransformFrom(int32_t value) {
    Radixable32 rdxble = *reinterpret_cast<Radixable32 *>(&value);
    rdxble ^= (1 << 31);
    return rdxble;
  }

  int32_t TransformTo(Radixable32 rdxble) {
    int32_t value = *reinterpret_cast<int32_t *>(&rdxble);
    value ^= (1 << 31);
    return value;
  }
};

template <>
struct RadixTransform<int64_t, Radixable64> {
  Radixable64 TransformFrom(int64_t value) {
    Radixable64 rdxble = *reinterpret_cast<Radixable64 *>(&value);
    rdxble ^= (((Radixable64)1) << 63);
    return rdxble;
  }

  int64_t TransformTo(Radixable64 rdxble) {
    int64_t value = *reinterpret_cast<int64_t *>(&rdxble);
    value ^= (((uint64_t)1) << 63);
    return value;
  }
};

template <>
struct RadixTransform<float, Radixable32> {
  Radixable32 TransformFrom(float value) {
    union _ { float value; Radixable32 cg;};
    Radixable32 rdxble = ((union _ *)&value)->cg;
    // Two's complement negative numbers, but keep them negative.
    if ((rdxble >> 31) == 1) {
      rdxble *= -1;
      rdxble ^= (1 << 31);
    }
    // Flip the high bit on all numbers to swap the positive and negative
    // regions.
    rdxble ^= (1 << 31);
    return rdxble;
  }

  float TransformTo(Radixable32 rdxble) {
    rdxble ^= (1 << 31);
    if ((rdxble >> 31) == 1) {
      rdxble ^= (1 << 31);
      rdxble *= -1;
    }
    union __ { float fg; Radixable32 rdxble;};
    float rdxdbl = ((union __ *)&rdxble)->fg;
    return rdxdbl;
  }
};

template <>
struct RadixTransform<double, Radixable64> {
  Radixable64 TransformFrom(double value) {
    union _ {double value; Radixable64 cg;};
    Radixable64 rdxble = ((union _ *)&value)->cg;
    // Two's complement negative numbers, but keep them negative.
    if ((rdxble >> 63) == 1) {
      rdxble *= -1;
      rdxble ^= (((Radixable64)1) << 63);
    }
    // Flip the high bit on all numbers to swap the positive and negative
    // regions.
    rdxble ^= (((Radixable64)1) << 63);
    return rdxble;
  }

  double TransformTo(Radixable64 rdxble) {
    rdxble ^= (((Radixable64)1) << 63);
    if ((rdxble >> 63) == 1) {
      rdxble ^= (((Radixable64)1) << 63);
      rdxble *= -1;
    }
    union __ { double fg; Radixable64 rdxble;};
    double rdxdbl = ((union __ *)&rdxble)->fg;
    return rdxdbl;
  }
};


template <typename RadixableType>
void RadixSortMSD(RadixableType* data, size_t N, size_t bit) {
  if (N <= 1 || std::is_sorted(data, data + N)) {
    return;
  }

  int left = 0;
  int right = N;
  while (left < right) {
    size_t left_bit = ((data[left] >> bit) & 1);

    if (left_bit == 0) {
      ++left;
      continue;
    }

    size_t right_bit = ((data[right - 1] >> bit) & 1);
    if (right_bit == 1) {
      --right;
      continue;
    }

    std::swap(data[left], data[right - 1]);
    ++left;
    --right;
  }

  /*
   * As long as more bits exist, call both sides of the pivot recursively.
   */
  if (bit > 0) {
    RadixSortMSD(data, left, bit - 1);
    RadixSortMSD(data + left, N - left, bit - 1);
  }
}

template <size_t Base, typename RadixableType>
void RadixSortLSD(RadixableType* data, size_t N, size_t) {
  std::vector<RadixableType[Base]> buckets(N);
  size_t pointers[Base];

  for (size_t div = 1; !std::is_sorted(data, data + N); div *= Base) {
    std::memset(pointers, 0, sizeof(pointers));

    for (size_t i = 0; i < N; ++i) {
      size_t bkt = (data[i] / div) % Base;
      buckets[pointers[bkt]++][bkt] = data[i];
    }

    int index = 0;
    for (size_t b = 0; b < Base; ++b) {
      for (size_t p = 0; p < pointers[b]; ++p) {
        data[index++] = buckets[p][b];
      }
    }
  }
}

template <typename RadixableType, typename ValueType>
void TransformAndInvoke(ValueType* array, size_t N,
                        void (*RadixFunc)(RadixableType*, size_t, size_t)) {
  RadixTransform<ValueType, RadixableType> transformer;

  RadixableType* radixable_array = reinterpret_cast<RadixableType*>(array);

  for (size_t i = 0; i < N; ++i) {
    radixable_array[i] = transformer.TransformFrom(array[i]);
  }

  RadixFunc(radixable_array, N, (8 * sizeof(RadixableType)) - 1);

  for (size_t i = 0; i < N; ++i) {
    array[i] = transformer.TransformTo(radixable_array[i]);
  }
}

}  // namespace details

/* Actual implementation of MSD Radix Sort. */
template <typename ValueType>
void RadixSortMSD(ValueType* array, size_t N) {
  using RadixableType = typename details::Radixable<sizeof(ValueType)>::type;
  details::TransformAndInvoke<RadixableType>(array, N, details::RadixSortMSD);
}

/* Actual implementation of LSD Radix Sort. */
template <size_t Base, typename ValueType>
void RadixSortLSD(ValueType* array, size_t N) {
  using RadixableType = typename details::Radixable<sizeof(ValueType)>::type;
  details::TransformAndInvoke<RadixableType>(array, N,
                                             details::RadixSortLSD<Base>);
}

#endif  // MODULES_TASK_3_ZOLOTAREVA_TASK3_LSD_RADIX_SORT_INCLUDE_RADIX_SORT_H_
