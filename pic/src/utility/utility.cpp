// src/utility/utility.cpp

#include <cmath>

namespace maxwell {

//==============================================================================
//  Binary logarithm
//==============================================================================
uint Binary_logarithm(const uint num)
{
    uint res = 0;
    
    for (uint n = num; !(n & 1); n >>= 1) { ++res; }

    return res;
}

/// useful /// uint32_t Ceil_binary_power(uint32_t num)
/// useful /// {
/// useful ///     --num;
/// useful /// 
/// useful ///     num |= num >> 1;
/// useful ///     num |= num >> 2;
/// useful ///     num |= num >> 4;
/// useful ///     num |= num >> 8;
/// useful ///     num |= num >> 16;
/// useful /// 
/// useful ///     return ++num;
/// useful /// }

//==============================================================================
//  Power
//==============================================================================
uint Power(const uint base, const uint exponent)
{
    uint res = 1;

    for (uint bas = base, exp = exponent; exp; exp >>= 1, bas *= bas)
    {
        if (exp & 1) { res *= bas; }
    }

    return res;
}

//==============================================================================
//  Sub hypercube count
//==============================================================================
// recursive
uint Sub_hypercube_count(const int sub_dim, const int dim)
{
    if (sub_dim < 0 || dim < 0 || sub_dim >= dim) { return 0; }
    else if (sub_dim == dim == 0) { return 1; }
    else
    {
        return (Sub_hypercube_count(sub_dim, dim - 1) << 1)
            + Sub_hypercube_count(sub_dim - 1, dim - 1);
    }
}


//==============================================================================
//
//  WIKIPEDIA MERGE SORT
//
//==============================================================================
#include <algorithm>
#include <iterator>
#include <omp.h>
#include <memory>

template <typename Iterator>
void mergesort(Iterator from, Iterator to)
{
#pragma omp parallel
    {
#pragma omp single nowait
        static_assert(!std::is_same<typename std::iterator_traits<Iterator>::value_type, void>::value);

        auto n = std::distance(from, to);

        if (1 < n)
        {
#pragma omp task firstprivate (from, to, n)
            {
                Iterator l_from = from;
                Iterator l_to = l_from;
                std::advance(l_to, n/2);
                mergesort(l_from, l_to);
            }
#pragma omp task firstprivate (from, to, n)
            {
                Iterator r_from = from;
                std::advance(r_from, n/2);
                Iterator r_to = r_from;
                std::advance(r_to, n-(n/2));
                mergesort(r_from, r_to);
            }
#pragma omp taskwait

            auto tmp_array = std::make_unique<typename Iterator::value_type[]>(n);
            Iterator l_iter = from;
            Iterator l_end = l_iter;
            std::advance(l_end, n/2);
            Iterator r_iter = l_end;
            Iterator& r_end = to;

            auto tmp_iter = tmp_array.get();

            while (l_iter != l_end || r_iter != r_end)
            {
                if (*l_iter < *r_iter)
                {
                    *tmp_iter = std::move(*l_iter);
                    ++l_iter;
                    ++tmp_iter;
                }
                else
                {
                    *tmp_iter = std::move(*r_iter);
                    ++r_iter;
                    ++tmp_iter;
                }

                if (l_iter == l_end)
                {
                    std::copy(
                                std::make_move_iterator(r_iter),
                                std::make_move_iterator(r_end),
                                tmp_iter
                    );

                    break;
                }

                if (r_iter == r_end)
                {
                    std::copy(
                                std::make_move_iterator(l_iter),
                                std::make_move_iterator(l_end),
                                tmp_iter
                    );

                    break;
                }
            }

            std::copy(
                        std::make_move_iterator(tmp_array.get()),
                        std::make_move_iterator(&tmp_array[n]),
                        from
            );
        }
    }
}

Итеративная реализация на языке C++:

template<typename T>
void MergeSort(T a[], size_t l)
{
    size_t BlockSizeIterator;
    size_t BlockIterator;
    size_t LeftBlockIterator;
    size_t RightBlockIterator;
    size_t MergeIterator;

    size_t LeftBorder;
    size_t MidBorder;
    size_t RightBorder;
    for (BlockSizeIterator = 1; BlockSizeIterator < l; BlockSizeIterator *= 2)
    {
        for (BlockIterator = 0; BlockIterator < l - BlockSizeIterator; BlockIterator += 2 * BlockSizeIterator)
        {
            //Производим слияние с сортировкой пары блоков начинающуюся с элемента BlockIterator
            //левый размером BlockSizeIterator, правый размером BlockSizeIterator или меньше
            LeftBlockIterator = 0;
            RightBlockIterator = 0;
            LeftBorder = BlockIterator;
            MidBorder = BlockIterator + BlockSizeIterator;
            RightBorder = BlockIterator + 2 * BlockSizeIterator;
            RightBorder = (RightBorder < l) ? RightBorder : l;
            int* SortedBlock = new int[RightBorder - LeftBorder];

            //Пока в обоих массивах есть элементы выбираем меньший из них и заносим в отсортированный блок
            while (LeftBorder + LeftBlockIterator < MidBorder && MidBorder + RightBlockIterator < RightBorder)
            {
                if (a[LeftBorder + LeftBlockIterator] < a[MidBorder + RightBlockIterator])
                {
                    SortedBlock[LeftBlockIterator + RightBlockIterator] = a[LeftBorder + LeftBlockIterator];
                    LeftBlockIterator += 1;
                }
                else
                {
                    SortedBlock[LeftBlockIterator + RightBlockIterator] = a[MidBorder + RightBlockIterator];
                    RightBlockIterator += 1;
                }
            }
            //После этого заносим оставшиеся элементы из левого или правого блока
            while (LeftBorder + LeftBlockIterator < MidBorder)
            {
                SortedBlock[LeftBlockIterator + RightBlockIterator] = a[LeftBorder + LeftBlockIterator];
                LeftBlockIterator += 1;
            }
            while (MidBorder + RightBlockIterator < RightBorder)
            {
                SortedBlock[LeftBlockIterator + RightBlockIterator] = a[MidBorder + RightBlockIterator];
                RightBlockIterator += 1;
            }

            for (MergeIterator = 0; MergeIterator < LeftBlockIterator + RightBlockIterator; MergeIterator++)
            {
                a[LeftBorder + MergeIterator] = SortedBlock[MergeIterator];
            }
            delete SortedBlock;
        }
    }
}

} // namespace maxwell

// src/utility/utility.cpp
