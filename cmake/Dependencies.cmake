include(FetchContent)

if(NOT TARGET algebra::finite_fields)
    FetchContent_Declare(
        algebra_finite_fields
        GIT_REPOSITORY git@github.com:MohamedAKonsowa/algebra-finite-fields-group2.git
        GIT_TAG 56840c9ead8e781fd18a4662dfc8dca73e612c86
    )
    FetchContent_MakeAvailable(algebra_finite_fields)
endif()

#if(NOT TARGET algebra::finite_fields)
#    FetchContent_Declare(
#        algebra_finite_fields
#        GIT_REPOSITORY git@github.com:nathanjrussell/algebra-finite-fields.git
#        GIT_TAG 74b3d569fb48818f4c10f3c32889f2128485f1a7
#    )
#    FetchContent_MakeAvailable(algebra_finite_fields)
#endif()

if(NOT TARGET algebra::matrices)
    FetchContent_Declare(
        algebra_matrices
        GIT_REPOSITORY git@github.com:nathanjrussell/algebra-matrices.git
        GIT_TAG 84abddac26d2805072efcb1c44bbe7faac45fa6b
    )
    FetchContent_MakeAvailable(algebra_matrices)
endif()
