#include <mpi.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <algebra/finite_fields/finite_fields.hpp>
#include <algebra/matrices/matrices.hpp>

namespace {

using Field = algebra::FiniteField<11, 3>;
using Element = Field::Element;
constexpr std::size_t dimension = ALGEBRA_ASSIGNMENT_SIZE;
constexpr std::uint32_t seed = 123456u;
using Matrix = algebra::Matrix<Element, dimension, dimension>;
using Clock = std::chrono::steady_clock;

static_assert(algebra::FieldElement<Element>);
static_assert(algebra::RingElement<Matrix>);

double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::uint32_t draw_coefficient(std::mt19937& generator, std::uint32_t prime) {
    if (prime < 2) {
        throw std::invalid_argument("The prime must be at least two");
    }
    constexpr std::uint64_t range = std::uint64_t{1} << 32;
    const std::uint64_t limit = range - range % prime;
    for (;;) {
        const std::uint64_t value = generator();
        if (value < limit) {
            return static_cast<std::uint32_t>(value % prime);
        }
    }
}

Element next_element(std::mt19937& generator) {
    Element::polynomial_type polynomial;
    for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        const auto coefficient = draw_coefficient(generator, Field::characteristic);
        polynomial.set_coefficient(degree, Element::coefficient_type::from_integer(coefficient));
    }
    return Element::from_polynomial(polynomial);
}

void populate(Matrix& a, Matrix& b) {
    std::mt19937 generator(seed);
    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            a(i, j) = next_element(generator);
            b(i, j) = next_element(generator);
        }
    }
}

Element determinantOLD(Matrix matrix) {
    Element result = Element::one();
    for (std::size_t k = 0; k < dimension; ++k) {
        std::size_t pivot_row = k;
        while (pivot_row < dimension && matrix(pivot_row, k).is_zero()) {
            ++pivot_row;
        }
        if (pivot_row == dimension) {
            return Element::zero();
        }
        if (pivot_row != k) {
            matrix.swap_rows(pivot_row, k);
            result = result.negate();
        }
        const Element pivot = matrix(k, k);
        result = result.multiply(pivot);
        Element inverse = Element::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        for (std::size_t i = k + 1; i < dimension; ++i) {
            if (matrix(i, k).is_zero()) {
                continue;
            }
            const Element multiplier = matrix(i, k).multiply(inverse);
            matrix(i, k) = Element::zero();
            for (std::size_t j = k + 1; j < dimension; ++j) {
                const Element correction = multiplier.multiply(matrix(k, j));
                matrix(i, j) = matrix(i, j).add(correction.negate());
            }
        }
    }
    return result;
}

uint64_t encode(const Element& element) {
    return element.index();
}

Element decode(uint64_t index) {
    return Element::from_index(index);
}

Element determinant_sequential(Matrix matrix) {
    Element result = Element::one();
    for (std::size_t k = 0; k < dimension; ++k) {
        std::size_t pivot_row = k;
        while (pivot_row < dimension && matrix(pivot_row, k).is_zero()) {
            ++pivot_row;
        }
        if (pivot_row == dimension) {
            return Element::zero();
        }
        if (pivot_row != k) {
            matrix.swap_rows(pivot_row, k);
            result = result.negate();
        }
        const Element pivot = matrix(k, k);
        result = result.multiply(pivot);
        Element inverse = Element::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        for (std::size_t i = k + 1; i < dimension; ++i) {
            if (matrix(i, k).is_zero()) {
                continue;
            }
            const Element multiplier = matrix(i, k).multiply(inverse);
            matrix(i, k) = Element::zero();
            for (std::size_t j = k + 1; j < dimension; ++j) {
                const Element correction = multiplier.multiply(matrix(k, j));
                matrix(i, j) = matrix(i, j).add(correction.negate());
            }
        }
    }
    return result;
}

Element determinant_mpi(Matrix matrix, std::size_t block_size) {
    int rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    Element result = Element::one();

    std::cout << "rank " << rank << " computing determinant of matrix\n";

    std::vector<size_t> P(dimension); 
    Matrix L;
    Matrix U;

    if (rank == 0) {
        std::cout << "rank 0 initializing...\n";
        for (size_t i = 0; i < dimension; ++i) {P[i] = i;}
        L = Matrix::one();
        U = matrix;
    }

    // Iterate over columns in blocks.
    if (block_size == 0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (size_t k = 0; k < dimension; k += block_size) {

        const size_t end = std::min(k + block_size, dimension);
        int singular = 0;

        // PART 1: Sequential panel factorization
        if (rank == 0) {
            for (size_t l = k; l < end; ++l) {
                // Find first nonzero pivot in column l
                size_t pivot_row = l;

                while (pivot_row < dimension &&
                    U(pivot_row, l).is_zero()) {
                    ++pivot_row;
                }

                if (pivot_row == dimension) {
                    singular = 1;
                    break;
                }

                // Swap rows
                if (pivot_row != l) {
                    U.swap_rows(pivot_row, l);

                    std::swap(P[pivot_row], P[l]);

                    // Swap previously computed L multipliers.
                    for (size_t j = 0; j < l; ++j) {
                        std::swap(L(pivot_row, j), L(l, j));
                    }

                    result = result.negate();
                }

                // Calculate inverse once per pivot.
                Element inverse = U(l, l).reciprocal();

                // Eliminate within the current panel.
                for (size_t i = l + 1; i < dimension; ++i) {

                    Element factor = U(i, l) * inverse;

                    L(i, l) = factor;
                    U(i, l) = Element::zero();

                    for (size_t j = l + 1; j < end; ++j) {
                        U(i, j) -= factor * U(l, j);
                    }
                }
            }

            if (!singular) {
                // Complete U12 using forward substitution.
                for (size_t i = k; i < end; ++i) {
                    for (size_t j = end; j < dimension; ++j) {
                        for (size_t t = k; t < i; ++t) {
                            U(i, j) -= L(i, t) * U(t, j);
                        }
                    }
                }
            }
        }

        // Communicate singularity before continuing
        MPI_Bcast(&singular, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (singular) {
            return Element::zero();
        }

        // PART 2: Parallel trailing-matrix update

        size_t b = end - k;
        size_t remaining = dimension - end;

        // No trailing matrix after the last panel.
        if (remaining == 0) {
            continue;
        }

        // Divide rows between all ranks.
        size_t base = remaining / num_ranks;
        size_t extra = remaining % num_ranks;

        size_t rows = base + (rank < extra ? 1 : 0);

        // Shared buffer: everyone needs all of U12.
        std::vector<uint64_t> packed_U12(b * remaining);

        // Rank 0 owns these complete packed regions.
        std::vector<uint64_t> packed_L21;
        std::vector<uint64_t> packed_U22;

        // Every rank owns its local pieces.
        std::vector<uint64_t> local_L21(rows * b);
        std::vector<uint64_t> local_U22(rows * remaining);

        // Scatter/gather metadata, used by rank 0.
        std::vector<int> counts_L(num_ranks);
        std::vector<int> displs_L(num_ranks);
        std::vector<int> counts_U(num_ranks);
        std::vector<int> displs_U(num_ranks);

        if (rank == 0) {
            packed_L21.resize(remaining * b);
            packed_U22.resize(remaining * remaining);

            // Pack U12.
            for (size_t i = 0; i < b; ++i) {
                for (size_t j = 0; j < remaining; ++j) {
                    packed_U12[i * remaining + j] =
                        encode(U(k + i, end + j));
                }
            }

            // Pack L21.
            for (size_t i = 0; i < remaining; ++i) {
                for (size_t j = 0; j < b; ++j) {
                    packed_L21[i * b + j] =
                        encode(L(end + i, k + j));
                }
            }

            // Pack U22.
            for (size_t i = 0; i < remaining; ++i) {
                for (size_t j = 0; j < remaining; ++j) {
                    packed_U22[i * remaining + j] =
                        encode(U(end + i, end + j));
                }
            }

            // Prepare counts and displacements.
            for (int r = 0; r < num_ranks; ++r) {
                size_t r_rows =
                    base + (static_cast<size_t>(r) < extra ? 1 : 0);

                size_t r_offset =
                    static_cast<size_t>(r) * base +
                    std::min<size_t>(r, extra);

                counts_L[r] = static_cast<int>(r_rows * b);
                displs_L[r] = static_cast<int>(r_offset * b);

                counts_U[r] = static_cast<int>(r_rows * remaining);
                displs_U[r] = static_cast<int>(r_offset * remaining);
            }
        }

        // 1. Broadcast U12 to every rank.
        MPI_Bcast(
            packed_U12.data(),
            static_cast<int>(packed_U12.size()),
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );

        // 2. Scatter L21 rows.
        MPI_Scatterv(
            packed_L21.data(),
            counts_L.data(),
            displs_L.data(),
            MPI_UINT64_T,
            local_L21.data(),
            static_cast<int>(local_L21.size()),
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );

        // 3. Scatter U22 rows.
        MPI_Scatterv(
            packed_U22.data(),
            counts_U.data(),
            displs_U.data(),
            MPI_UINT64_T,
            local_U22.data(),
            static_cast<int>(local_U22.size()),
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );

        // 4. EVERY rank performs its own trailing update.
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < remaining; ++j) {

                Element value = decode(local_U22[i * remaining + j]);

                for (size_t t = 0; t < b; ++t) {
                    Element l = decode(local_L21[i * b + t]);
                    Element u = decode(packed_U12[t * remaining + j]);

                    value -= l * u;
                }

                local_U22[i * remaining + j] = encode(value);
            }
        }

        // 5. Gather updated U22 rows back to rank 0.
        MPI_Gatherv(
            local_U22.data(),
            static_cast<int>(local_U22.size()),
            MPI_UINT64_T,
            packed_U22.data(),
            counts_U.data(),
            displs_U.data(),
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );

        // 6. Rank 0 writes the results into global U.
        if (rank == 0) {
            for (size_t i = 0; i < remaining; ++i) {
                for (size_t j = 0; j < remaining; ++j) {
                    U(end + i, end + j) =
                        decode(packed_U22[i * remaining + j]);
                }
            }
        }

    }

    if (rank == 0) {
        for (size_t i = 0; i < dimension; ++i) {
            result = result.multiply(U(i, i));
        }
    }

    return result;
}

void print_determinant(const char* name, const Element& value, double elapsed) {
    std::cout << name << " = " << value << '\n';
    std::cout << "Coefficients [";
    for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        if (degree != 0) {
            std::cout << ", ";
        }
        std::cout << value.coefficient(degree);
    }
    std::cout << "]\nSeconds " << elapsed << '\n' << std::flush;
}

}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string mode = argc > 1 ? argv[1] : "mpi";
    size_t block_size = argc > 2 ? std::stoull(argv[2]) : 32;

    if (rank == 0) {
        std::cout << "Block size: " << block_size << '\n';
        std::cout << "Mode: " << mode << '\n';
    }

    try {
        const auto total_start = Clock::now();

        Matrix a;
        Matrix b;
        Matrix c;

        Field::initialize();

        std::cout << "I am rank " << rank << " of " << size << std::endl;

        Element det_a, det_b, det_c;
        double a_seconds, b_seconds, c_seconds, pop_seconds, mult_seconds;

        if (rank == 0) {
            auto pop_start = Clock::now();
            populate(a, b);
            pop_seconds = seconds(pop_start);

            auto mult_start = Clock::now();
            c = a.multiply(b);
            mult_seconds = seconds(mult_start);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (mode == "mpi") {
            auto a_start = Clock::now();
            det_a = determinant_mpi(a, block_size);
            a_seconds = seconds(a_start);

            auto b_start = Clock::now();
            det_b = determinant_mpi(b, block_size);
            b_seconds = seconds(b_start);

            auto c_start = Clock::now();
            det_c = determinant_mpi(c, block_size);
            c_seconds = seconds(c_start);
        }
        else if (mode == "sequential") {
            if (rank == 0) {
                auto a_start = Clock::now();
                det_a = determinant_sequential(a);
                a_seconds = seconds(a_start);

                auto b_start = Clock::now();
                det_b = determinant_sequential(b);
                b_seconds = seconds(b_start);

                auto c_start = Clock::now();
                det_c = determinant_sequential(c);
                c_seconds = seconds(c_start);
            }
        }
        else {
            throw std::invalid_argument("Invalid mode: " + mode);
        }
        
        if (rank == 0) {
            print_determinant("det(A)", det_a, a_seconds);
            print_determinant("det(B)", det_b, b_seconds);
            print_determinant("det(C)", det_c, c_seconds);

            const auto expected = det_a.multiply(det_b);
            const bool verified = det_c == expected;

            std::cout << std::fixed << std::setprecision(2)
                      << "det(A) * det(B) = " << expected << '\n'
                      << "Identity verified " << std::boolalpha
                      << verified << "\n\n"
                      << "Total populate seconds: " << pop_seconds << " \n"
                      << "Total multiply seconds: " << mult_seconds << " \n"
                      << "Total determinant seconds: " << a_seconds + b_seconds + c_seconds << " \n"
                      << "Total seconds " << seconds(total_start)
                      << '\n';
        }
    }
    catch (const std::exception& error) {
        std::cerr << "Rank " << rank
                  << ": " << error.what() << '\n';

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}

