

#ifndef POWERGRID_PERMUTE_HPP
#define POWERGRID_PERMUTE_HPP
#include <armadillo>
#include <tuple>
#include <algorithm>
#include <vector>

using namespace arma;
typedef std::tuple<std::size_t,std::size_t,std::size_t> D3tuple;

template <typename T>
inline Col<uword> shape (const Cube<T>& x)
{
    return { x.n_rows, x.n_cols, x.n_slices };
}
template <typename T>
static Cube<T> permute (Cube<T>& cube, const D3tuple& order)
{
    uword idx1 = std::get<0>(order);
    uword idx2 = std::get<1>(order);
    uword idx3 = std::get<2>(order);

    Col<uword> dimension = shape(cube);

    uword rows = dimension(idx1 - 1);
    uword cols = dimension(idx2 - 1);
    uword slis = dimension(idx3 - 1);

    Cube<T> output;
    output.zeros(rows, cols, slis);

    uword perm = idx1*100 + idx2*10 + idx3;

    switch (perm)
    {
        case 123:
        {
            output = cube; // identity
        }
        break;
        case 132:
        {
            for (int c = 0; c < cube.n_cols; ++c)
                for (int r = 0; r < cube.n_rows; ++r)
                    for (int s = 0; s < cube.n_slices; ++s)
                        output(r, s, c) = cube(r, c, s);
        }
        break;
        case 213:
        {
            for (int c = 0; c < cube.n_cols; ++c)
                for (int r = 0; r < cube.n_rows; ++r)
                    for (int s = 0; s < cube.n_slices; ++s)
                        output(c, r, s) = cube(r, c, s);
        }
        break;
        case 231:
        {
            for (int c = 0; c < cube.n_cols; ++c)
                for (int r = 0; r < cube.n_rows; ++r)
                    for (int s = 0; s < cube.n_slices; ++s)
                        output(c, s, r) = cube(r, c, s);
        }
        break;
        case 312:
        {
            for (int c = 0; c < cube.n_cols; ++c)
                for (int r = 0; r < cube.n_rows; ++r)
                    for (int s = 0; s < cube.n_slices; ++s)
                        output(s, r, c) = cube(r, c, s);
        }
        break;
        case 321:
        {
            for (int c = 0; c < cube.n_cols; ++c)
                for (int r = 0; r < cube.n_rows; ++r)
                    for (int s = 0; s < cube.n_slices; ++s)
                        output(s, c, r) = cube(r, c, s);
        }
        break;
    }

    return output;
}

#endif //POWERGRID_PERMUTE_HPP
